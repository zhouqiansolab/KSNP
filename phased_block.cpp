//
// Created by ixiaohu on 2022/1/30.
//

#include <cassert>
#include <iostream>
#include <zlib.h>
#include <algorithm>
#include <map>
#include <set>
#include <cstring>
#include <fstream>

#include "phased_block.h"

PostProcess::PostProcess(std::vector<SNP> &snp, const std::vector<Read_Allele> &read,
                         std::vector<Phased_Block> &b, std::vector<int8_t> &h) :
                         snp_column(snp), read_row(read),
                         blocks(b), hint(h) {
	fprintf(stderr, "Post-process haplotype blocks from DBG\n");
}

void PostProcess::de_overlap() {
	// There must be no gap within SNP block from DBG
	for (const auto &b : blocks) assert(b.snp_idx.back() - b.snp_idx.front() + 1 == b.size());
	std::sort(blocks.begin(), blocks.end());

	// Remove short blocks fully covered by large block
	std::vector<Phased_Block> big_blocks;
	int shadowed_n = 0, length_sum = 0;
	for (int i = 0; i < blocks.size(); ) {
		big_blocks.push_back(blocks[i]);

		const auto &bi = blocks[i];
		int begin_i = bi.snp_idx.front(), end_i = bi.snp_idx.back();
		int next = blocks.size();
		for (int j = i+1; j < blocks.size(); j++) {
			const auto &bj = blocks[j];
			int begin_j = bj.snp_idx.front(), end_j = bj.snp_idx.back();
			if (end_j > end_i) {
				next = j;
				break;
			} else {
				assert(begin_j >= begin_i);
				shadowed_n++;
				length_sum += bj.size();
			}
		}
		i = next;
	}

	// Remove overlaps between blocks
	int over_sum = 0, snp_sum = 0;
	blocks.clear();
	for (int i = 0; i < big_blocks.size(); ) {
		auto &bi = big_blocks[i];

		int beg_i = bi.snp_idx.front(), end_i = bi.snp_idx.back() + 1;
		int next = big_blocks.size();
		auto *tree_array = new int[bi.size() + 1];
		memset(tree_array, 0, (bi.size() + 1) * sizeof(int));
		for (int j = i + 1; j < big_blocks.size(); j++) {
			const auto &bj = big_blocks[j];
			int beg_j = bj.snp_idx.front(), end_j = bj.snp_idx.back() + 1;
			if (beg_j >= end_i) {
				next = j;
				break;
			}
//			fprintf(stderr, "block %d [%d,%d) overlaps with block %d [%d,%d)\n", j, beg_j, end_j, i, beg_i, end_i);
			assert(end_j > end_i); // Blocks must be tiled
			// So far, I just append the non-overlap rear part of next block to the previous block.
			// There might be some better solutions to merge. Please check out the necessary of optimization.
			for (int k = end_i-beg_j; k < bj.size(); k++) {
				bi.push_back(bj.snp_idx[k], bj.phase[k]);
			}
			tree_array[beg_j-bi.snp_idx.front()]++; tree_array[end_i-bi.snp_idx.front()]--;
		}
		blocks.push_back(bi);
		snp_sum += bi.size();
		for (int j = 0, t = 0; j < bi.size(); j++) {
			t += tree_array[j];
			if (t > 0) over_sum++;
		}
		delete [] tree_array;
		i = next;
	}
	fprintf(stderr, "    Remove %.2f %% overlaps between blocks\n", 100.0 * over_sum / snp_sum);
}

int PostProcess::mec_of_template(Phased_Block &b, int p, int TLEN, uint tmp) {
	std::set<int> read_set; // Reads on b[p, p+TLEN)
	int max_r = -1;
	for (int i = 0; i < TLEN; i++) {
		const auto &s = snp_column[b.ps + p + i];
		for (int j = 0; j < s.size(); j++) {
			int rid = s.rid[j];
			int l = read_row[rid].front().snp_idx;
			int r = read_row[rid].back().snp_idx;
			if (b.snp_idx.front() <= l and r <= b.snp_idx.back()) {
				// Select reads within this block
				max_r = std::max(max_r, r - b.ps);
				read_set.insert(rid);
			}
		}
	}
	if (read_set.empty()) return INT32_MAX;

	// Set bits in template indicates SNPs should be reversed. For example
	// block 0101010101, template
	//         00100    is used for correcting error, then the block becomes
	//       0101101010. When set number is odd, all SNPs after template is reversed.
	uint8_t set_odd = 0U; // Whether set number in template is odd.
	for (uint i = 0; i < TLEN; i++) {
		if ((tmp & (1U << i)) != 0) set_odd ^= 1U;
		if (set_odd) b.phase[p + i] ^= 1;
	}
	if (set_odd) for (int i = p + TLEN; i <= max_r; i++) b.phase[i] ^= 1;

	// Calculating MEC
	int mec = 0;
	for (const auto rid : read_set) {
		int h0 = 0, h1 = 0;
		for (const auto &v : read_row[rid]) {
			if (v.allele == -1) continue;
			if (v.allele != b.phase[v.snp_idx - b.ps]) h0++;
			else h1++;
		}
		mec += std::min(h0, h1);
	}

	set_odd = 0U; // Rollback
	for (uint i = 0; i < TLEN; i++) {
		if ((tmp & (1U << i)) != 0) set_odd ^= 1U;
		if (set_odd) b.phase[p + i] ^= 1;
	}
	if (set_odd) for (int i = p + TLEN; i <= max_r; i++) b.phase[i] ^= 1;
	return mec;
}

bool PostProcess::try_to_switch(const Phased_Block &b, int p) {
	const auto &pivot_snp = snp_column[b.ps + p]; // [0,p] kept, [p+1,n) reversed
	int old_mec = 0, new_mec = 0;
	for (int i = 0; i < pivot_snp.size(); i++) {
		const auto &read = read_row[pivot_snp.rid[i]];
		int l = read.front().snp_idx, r = read.back().snp_idx;
		int old_c0 = 0, old_c1 = 0; // For original haplotype
		int new_c0 = 0, new_c1 = 0; // For reversed haplotype
		if (b.snp_idx.front() <= l and r <= b.snp_idx.back()) {
			for (const auto &g : read) {
				if (g.allele == -1) continue;
				if (g.allele != b.phase[g.snp_idx - b.ps]) {
					old_c0++;
					if (g.snp_idx - b.ps > p) new_c1++;
					else new_c0++;
				} else {
					old_c1++;
					if (g.snp_idx - b.ps > p) new_c0++;
					else new_c1++;
				}
			}
		}
		old_mec += std::min(old_c0, old_c1);
		new_mec += std::min(new_c0, new_c1);
	}
	return new_mec < old_mec;
}

bool PostProcess::try_to_flip(const Phased_Block &b, int p) {
	const auto &pivot_snp = snp_column[b.ps + p]; // [0,p] kept, [p+1,n) reversed
	int old_mec = 0, new_mec = 0;
	for (int i = 0; i < pivot_snp.size(); i++) {
		const auto &read = read_row[pivot_snp.rid[i]];
		int l = read.front().snp_idx, r = read.back().snp_idx;
		int old_c0 = 0, old_c1 = 0; // For original haplotype
		int new_c0 = 0, new_c1 = 0; // For reversed haplotype
		if (b.snp_idx.front() <= l and r <= b.snp_idx.back()) {
			for (const auto &g : read) {
				if (g.allele == -1) continue;
				if (g.allele != b.phase[g.snp_idx - b.ps]) {
					old_c0++;
					if (g.snp_idx - b.ps == p) new_c1++;
					else new_c0++;
				} else {
					old_c1++;
					if (g.snp_idx - b.ps == p) new_c0++;
					else new_c1++;
				}
			}
		}
		old_mec += std::min(old_c0, old_c1);
		new_mec += std::min(new_c0, new_c1);
	}
	return new_mec < old_mec;
}

void PostProcess::error_correction() {
	int need_n = 0; for (int i = 0; i < snp_column.size(); i++) need_n += (hint[i] == 1);
	fprintf(stderr, "    Local correctness for %.2f %% phased SNPs from graph\n", 100.0 * need_n / snp_column.size());

	const int LOCAL_SIZE = 5; // Template size
	const int HT_N = 14; // High-quality CUT templates
	const uint HIGH_TEMP[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 15, 20, 24};
	const int LT_N = 17;
	const uint LOW_TEMP[] = {11, 13, 14, 16, 17, 18, 19, 21, 22, 23, 25, 26, 27, 28, 29, 30, 31};

	for (auto &b : blocks) {
		// s(i)=true represents the phased result from i+1 to end is reversed
		auto *set_bit = (bool*) calloc(b.size(), sizeof(bool));
		for (uint t : HIGH_TEMP) {
			memset(set_bit, false, b.size() * sizeof(bool));
			for (int i = 0; i + LOCAL_SIZE <= b.size(); i++) {
				if (!hint[b.ps + i]) continue;
				int mec0 = mec_of_template(b, i, LOCAL_SIZE, 0);
				int mec1 = mec_of_template(b, i, LOCAL_SIZE, t);
				if (mec0 <= mec1) { // If do not update, skip 1
					continue;
				}
				for (int j = 0; j < LOCAL_SIZE; j++) {
					if ((t & (1<<j)) != 0) {
						set_bit[i + j] = true;
					}
				}
				i += LOCAL_SIZE - 1; // If update, skip 5.
			}

			// An cumulative variable to achieve this:
			// Old block={0,0,1,0,0,1}
			// Flag=     {F,T,F,F,T,F}
			// Cum=      {0,1,1,1,0,0}
			// New block={0,1,0,1,0,1}
			// Old block is reversed at Cum=1, and kept at Cum=0.
			uint8_t cum = 0U;
			for (int i = 0; i < b.size(); i++) {
				if (set_bit[i]) cum ^= 1U;
				if (cum == 1U) b.phase[i] ^= 1;
			}
		}
		free(set_bit);
	}

	// Try to flip each SNV to lower switch error
	for (auto &b : blocks) {
		for (int i = 0; i < b.size(); i++) {
			if (try_to_flip(b, i)) {
				b.phase[i] ^= 1;
				i += 64;
			}
		}
	}

	// Try to switch each SNV to lower hamming error
	for (auto &b : blocks) {
		auto *set_bit = (bool*) calloc(b.size(), sizeof(bool));
		for (int i = 1; i < b.size(); i++) {
			set_bit[i] = try_to_switch(b, i - 1);
			if (set_bit[i]) i += 64;
		}
		uint8_t cum = 0U;
		for (int i = 0; i < b.size(); i++) {
			if (set_bit[i]) cum ^= 1U;
			if (cum == 1U) b.phase[i] ^= 1;
		}
		free(set_bit);
	}
}

int8_t PostProcess::connect_block(const Phased_Block &prev, const Phased_Block &curr) {
	const int8_t BROKE = -1, KEPT = 0, REVERSE = 1;

	std::set<int> read_set;
	int p_beg = prev.snp_idx.front(), p_end = prev.snp_idx.back(); // [beg, end]
	int c_beg = curr.snp_idx.front(), c_end = curr.snp_idx.back();
	for (int i = 0; i < snp_column[p_end].size(); i++) {
		int rid = snp_column[p_end].rid[i];
		const auto &v = read_row[rid];
		// Skip reads overflow two blocks, or within the left side
		if (v.front().snp_idx < p_beg or v.back().snp_idx > c_end or v.back().snp_idx < c_beg) continue;
		read_set.insert(rid);
	}
	for (int i = 0; i < snp_column[c_beg].size(); i++) {
		int rid = snp_column[c_beg].rid[i];
		const auto &v = read_row[rid];
		// Skip reads overflow two blocks, or within the right side
		if (v.front().snp_idx < p_beg or v.back().snp_idx > c_end or v.front().snp_idx > p_end) continue;
		read_set.insert(rid);
	}
	if (read_set.empty()) return BROKE; // No reads crossing two consecutive blocks

	// Calculate MEC for direct and reversed connecting two blocks
	int kept_mec = 0, reverse_mec = 0;
	bool active = false; // An active read has non-gap alleles on two blocks simultaneously.
	for (const auto rid : read_set) {
		int kept_h0 = 0, kept_h1 = 0;
		int reverse_h0 = 0, reverse_h1 = 0;
		bool prev_active = false, curr_active = false;
		for (const auto &v : read_row[rid]) {
			if (v.allele == -1) continue;
			if (v.snp_idx <= p_end) {
				prev_active = true;
				if (v.allele == prev.phase[v.snp_idx - prev.ps]) {
					kept_h0++; reverse_h0++;
				} else {
					kept_h1++; reverse_h1++;
				}
			} else if (v.snp_idx >= c_beg) {
				curr_active = true;
				if (v.allele == curr.phase[v.snp_idx - curr.ps]) {
					kept_h0++; reverse_h1++;
				} else {
					kept_h1++; reverse_h0++;
				}
			} // else in the gap between two blocks
		}
		if (!prev_active or !curr_active) continue; // not an active read
		active = true;
		kept_mec += std::min(kept_h0, kept_h1);
		reverse_mec += std::min(reverse_h0, reverse_h1);
	}
	if (!active) return BROKE;
	if (kept_mec > reverse_mec) return REVERSE;
	if (kept_mec < reverse_mec) return KEPT;
	if (kept_mec == reverse_mec) return BROKE; // An conservative strategy
}

void PostProcess::merge() {
	// Absolute state against the original block
	auto *absolute = new int8_t[blocks.size()];
	const int8_t BROKE = -1, KEPT = 0, REVERSE = 1;

	absolute[0] = KEPT;
	for (int i = 1; i < blocks.size(); i++) {
		const auto &prev = blocks[i-1], &curr = blocks[i];
		// Relative state between the previous and current block
		int8_t relative = connect_block(prev, curr);
		absolute[i] = relative;
		if (relative == BROKE or absolute[i-1] == BROKE) continue; // No connected, do nothing
		// If previous block is absolutely reversed, KEPT becomes to REVERSE, and REVERSE becomes gto KEPT.
		if (absolute[i-1] == REVERSE) absolute[i] = (relative ^ 1);
	}

	std::vector<Phased_Block> merged_blocks;
	merged_blocks.push_back(blocks[0]);
	int n = 0;
	for (int i = 1; i < blocks.size(); i++) {
		const auto &b = blocks[i];
		if (absolute[i] == BROKE) {
			merged_blocks.push_back(b); // Add a new block to the merged set
		} else if (absolute[i] == KEPT) {
			// Append the block to the last one in the merged set.
			for (int j = 0; j < b.size(); j++)
				merged_blocks.back().push_back(b.snp_idx[j], b.phase[j]);
			n++;
		} else {
			// Append in reversed way
			for (int j = 0; j < b.size(); j++)
				merged_blocks.back().push_back(b.snp_idx[j], b.phase[j] ^ 1);
			n++;
		}
	}
	delete [] absolute;

	fprintf(stderr, "    Finally, %ld merged haplotype blocks\n", merged_blocks.size());
	blocks = merged_blocks;
}

void PostProcess::resolve() {
	de_overlap(); // TODO: remove contained blocks or overlaps in DBG part
	error_correction();
	merge();

	// Set phased information for original input variants
	for (const auto &b : blocks) {
		for (int i = 0; i < b.size(); i++) {
			auto &snp = snp_column[b.snp_idx[i]];
			snp.ps = snp_column[b.ps].pos;
			snp.gt = b.phase[i];
		}
	}
	blocks.clear();
}

void output_vcf_header(VCF_Header &header, const char *fn) {
	FILE *f_out = fn ?fopen(fn, "w") :stdout;
	if (f_out == nullptr) {
		fprintf(stderr, "ERR: can not open file %s\n", fn);
		std::abort();
	}
	header.annotate_for_phasing();
	for (const auto &line : header.to_str()) {
		fprintf(f_out, "%s\n", line.c_str());
	}
	if (fn) fclose(f_out);
}

void output_haplotype_block(const std::vector<SNP> &vcf_snp, const char *dst_fn) {
	FILE *f_out = dst_fn ?fopen(dst_fn, "a") :stdout;
	for (const auto &snp : vcf_snp) {
		auto fields = split_str(snp.line, '\t');
		auto fmt = split_str(fields[VCF_FORMAT].c_str(), ':');
		auto smp = split_str(fields[VCF_SAMPLE].c_str(), ':');
		assert(fmt.size() == smp.size());
		int pos = atoi(fields[VCF_POS].c_str()); assert(snp.pos == pos);

		// The arrangement of GT:GQ:PS is recommended.
		bool has_gt = false, has_ps = false;
		for (const auto &k : fmt) { has_gt |= (k == "GT"); has_ps |= (k == "PS"); }
		if (!has_gt) { fmt.emplace_back("GT"); smp.emplace_back("."); }
		if (!has_ps) { fmt.emplace_back("PS"); smp.emplace_back("."); }
		for (int i = 0; i < fmt.size(); i++) {
			if (fmt[i] == "GT") { std::swap(fmt[i], fmt[0]); std::swap(smp[i], smp[0]); }
			if (fmt[i] == "GQ") { std::swap(fmt[i], fmt[1]); std::swap(smp[i], smp[1]); }
			if (fmt[i] == "PS") { std::swap(fmt[i], fmt[2]); std::swap(smp[i], smp[2]); }
		}

		// See whether this SNP from original VCF file is phased by KSNP.
		if (snp.gt != -1) {
			smp[0] = snp.gt == 0 ?"0|1" :"1|0";
			smp[2] = std::to_string(snp.ps); assert(snp.ps != -1);
		}

		fields[VCF_FORMAT] = join_str(fmt, ':');
		fields[VCF_SAMPLE] = join_str(smp, ':');
		auto o = join_str(fields, '\t');
		fprintf(f_out, "%s\n", o.c_str());
		free(snp.line);
	}
	fflush(f_out);
	if (dst_fn) fclose(f_out);
	fprintf(stderr, "Phased results are written into the output file\n");
	fprintf(stderr, "\n");
}
