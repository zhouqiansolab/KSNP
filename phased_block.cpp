//
// Created by ixiaohu on 2022/1/30.
//

#include <cassert>
#include <iostream>
#include <fstream>
#include <zlib.h>

#include "phased_block.h"

std::vector<SNP_Block> remove_overlap(std::vector<SNP_Block> &blocks) {
	bool *shadowed = (bool*) calloc(blocks.size(), sizeof(bool));
	for (int i = 1; i < blocks.size(); i++) {
		assert(blocks[i].anchor >= blocks[i-1].anchor);
	}
	int cut_off_n = 0;
	for (int i = 0; i < blocks.size(); i++) {
		const auto &bi = blocks[i];
		int begin_i = bi.anchor, end_i = bi.anchor + bi.n, len_i = end_i - begin_i;
		for (int j = i+1; j < blocks.size(); j++) {
			const auto &bj = blocks[j];
			int begin_j = bj.anchor, end_j = bj.anchor + bj.n, len_j = end_j - begin_j;
			if (begin_j >= end_i) break;
			int overlap = std::min(end_i, end_j) - std::max(begin_i, begin_j);
			if (overlap == len_j) shadowed[j] = true;
			else if (overlap == len_i) {
				if (len_j > len_i) {
					shadowed[i] = true;
					break; // Don't worry, the blocks shadowed by $i can also be shadowed by $j
				}
			} else {
				// Cut off the rear part of block $i
				blocks[i].n -= overlap;
				cut_off_n++;
				break;
			}
		}
	}
	int cnt = 0;
	for (int i = 0; i < blocks.size(); i++) {
		if (shadowed[i]) continue;
		cnt++;
		const auto &bi = blocks[i]; assert(bi.n > 0);
		int end_i = bi.anchor + bi.n;
		for (int j = i+1; j < blocks.size(); j++) {
			if (shadowed[j]) continue;
			const auto &bj = blocks[j];
			int begin_j = bj.anchor;
			assert(begin_j >= end_i);
		}
	}
	std::cerr << "Remove " << blocks.size() - cnt << " shadowed block pairs, ";
	std::cerr << "and cut off " << cut_off_n << " block pairs" << std::endl;

	std::vector<SNP_Block> ret;
	for (int i = 0; i < blocks.size(); i++) {
		const auto &b = blocks[i];
		if (shadowed[i]) free(blocks[i].bit);
		else ret.push_back(b);
	}
	blocks.clear();
	free(shadowed);
	return ret;
}

extern std::string global_chrom;
extern std::vector<SNP> global_snps;

std::vector<SNP_Block> input_blocks_from_file(const char *fn) {
	freopen(fn, "r", stdin);
	std::vector<SNP_Block> ret;
	while (true) {
		std::string block_id; std::cin >> block_id;
		if (std::cin.eof()) break;
		std::string chrom; std::cin >> chrom; global_chrom = chrom;
		std::string pos_list; std::cin >> pos_list;
		std::string base_list; std::cin >> base_list;
		std::vector<int> pos_v;
		int num = 0;
		for (auto c : pos_list) {
			if (c == '_') {
				pos_v.push_back(num);
				num = 0;
			} else {
				num *= 10;
				num += c - '0';
			}
		}
		pos_v.push_back(num);
		assert(pos_v.size() == base_list.size());
		int p = binary_search(global_snps, pos_v[0]);
		assert(p != -1 && pos_v[0] == global_snps[p].pos);
		int n = pos_v.size();
		SNP_Block b(p);
		for (int i = 0; i < n; i++) {
			if (global_snps[p + i].alt == base_list[i]) b.add(1);
			else b.add(0);
		}
		ret.push_back(b);

		std::cin >> block_id >> chrom >> pos_list >> base_list; // Skip the mate line
	}
	std::cerr << "Input " << ret.size() << " blocks from stdin" << std::endl;
	fclose(stdin);
	return ret;
}

void view_merged_blocks() {
	int blocks_n = 0, merge_n = 0;
	while (true) {
		blocks_n++;
		std::string block_id; std::cin >> block_id;
		if (std::cin.eof()) break;
		std::string chrom; std::cin >> chrom;
		std::string pos_list; std::cin >> pos_list;
		std::string base_list; std::cin >> base_list;
		std::vector<int> pos_v;
		int num = 0;
		for (auto c : pos_list) {
			if (c == '_') {
				pos_v.push_back(num);
				num = 0;
			} else {
				num *= 10;
				num += c - '0';
			}
		}
		pos_v.push_back(num);
		assert(pos_v.size() == base_list.size());
		for (int i = 1; i < pos_v.size(); i++) {
			int p1 = binary_search(global_snps, pos_v[i-1]);
			int p2 = binary_search(global_snps, pos_v[i]);
			if (p2 != p1 + 1) {
				merge_n++;
			}
		}
		std::cin >> block_id >> chrom >> pos_list >> base_list; // Skip the mate line
	}
	std::cerr << "Collect " << blocks_n << " block pairs" << std::endl;
	std::cerr << "Found " << merge_n << " merged blocks" << std::endl;
}


void Phased_Result::add(const SNP_Block &b) {
	for (int i = 0; i < b.n; i++) {
		index[n] = b.anchor + i;
		bit[n] = b.bit[i];
		head[n] = b.anchor;
		n++;
	}
}

void Phased_Result::merge(const SNP_Block &b, bool reversed) {
	for (int i = 0; i < b.n; i++) {
		index[n] = b.anchor + i;
		bit[n] = reversed ?(b.bit[i] ^ 1) :b.bit[i];
		head[n] = head[n-1];
		n++;
	}
}

void Phased_Result::write_vcf(const char *src, const char *dst) {
	gzFile f_in = gzopen(src, "r");
	FILE *f_out = fopen(dst, "w");
	int line_n = 0, pid = 0;
	const int SNP_BUF_SIZE = 4 * 1024 * 1024;
	char *buf = (char*) malloc(SNP_BUF_SIZE * sizeof(char));
	std::vector<int> pos;
	bool format_GT = false, format_GQ = false, format_PS = false;
	bool format_fields = false;
	while (gzgets(f_in, buf, SNP_BUF_SIZE)) {
		// VCF headers
		if (buf[0] == '#') {
			std::string header(buf);
			if (header.find("##FORMAT") != std::string::npos) format_fields = true;
			if (header.find("##FORMAT=<ID=GT") != std::string::npos) format_GT = true;
			if (header.find("##FORMAT=<ID=GQ") != std::string::npos) format_GQ = true;
			if (header.find("##FORMAT=<ID=PS") != std::string::npos) format_PS = true;
			if (header.find("##FORMAT") == std::string::npos && format_fields && !format_GT)
				fprintf(f_out, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
			if (header.find("##FORMAT") == std::string::npos && format_fields && !format_GQ)
				fprintf(f_out, "##FORMAT=<ID=GQ,Number=1,Type=Float,Description=\"Genotype Quality\">\n");
			if (header.find("##FORMAT") == std::string::npos && format_fields && !format_PS)
				fprintf(f_out, "##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase Set\">\n");
			if (header.find("##FORMAT") == std::string::npos && format_fields) format_fields = false;
			fprintf(f_out, "%s", buf);
			continue;
		}

		// VCF columns
		const char *CHROM, *POS, *ID;
		const char *REF, *ALT, *QUAL, *FILTER;
		const char *INFO, *FORMAT, *SAMPLE;
		int tab_n = 0;
		for (int i = 0; buf[i]; i++) {
			if (i == 0) {
				CHROM = buf + i;
			} else if (buf[i] == '\t') {
				tab_n++;
				if (tab_n == 1) POS = buf + i + 1;
				else if (tab_n == 2) ID = buf + i + 1;
				else if (tab_n == 3) REF = buf + i + 1;
				else if (tab_n == 4) ALT = buf + i + 1;
				else if (tab_n == 5) QUAL = buf + i + 1;
				else if (tab_n == 6) FILTER = buf + i + 1;
				else if (tab_n == 7) INFO = buf + i + 1;
				else if (tab_n == 8) FORMAT = buf + i + 1;
				else if (tab_n == 9) SAMPLE = buf + i + 1;
				buf[i] = '\0';
			} else if (buf[i] == '\n') buf[i] = '\0';
		}
		pos.push_back(atoi(POS));

		// (format, sample) key-value pairs
		std::vector<std::string> key, value;
		std::string temp;
		for (int i = 0; true; i++) {
			if (FORMAT[i] == '\0') {
				key.push_back(temp);
				break;
			} else if (FORMAT[i] == ':') {
				key.push_back(temp);
				temp = "";
			} else {
				temp += FORMAT[i];
			}
		}
		temp = "";
		for (int i = 0; true; i++) {
			if (SAMPLE[i] == '\0') {
				value.push_back(temp);
				break;
			} else if (SAMPLE[i] == ':') {
				value.push_back(temp);
				temp = "";
			} else {
				temp += SAMPLE[i];
			}
		}
		assert(key.size() == value.size());
		bool has_gt = false, has_gq = false, has_ps = false;
		for (const auto &k : key) {
			if (k == "GT") has_gt = true;
			else if (k == "GQ") has_gq = true;
			else if (k == "PS") has_ps = true;
		}
		if (!has_gt) { key.emplace_back("GT"); value.emplace_back("."); }
		if (!has_gq) { key.emplace_back("GQ"); value.emplace_back("."); }
		if (!has_ps) { key.emplace_back("PS"); value.emplace_back("."); }
		// It is recommended that re-arrange GT:GQ:PS ahead.
		for (int i = 0; i < key.size(); i++) {
			if (key[i] == "GT") { std::swap(key[i], key[0]); std::swap(value[i], value[0]); }
			if (key[i] == "GQ") { std::swap(key[i], key[1]); std::swap(value[i], value[1]); }
			if (key[i] == "PS") { std::swap(key[i], key[2]); std::swap(value[i], value[2]); }
		}

		while (pid < n && index[pid] < line_n) pid++;
		if (pid < n && index[pid] == line_n) {
			if (bit[pid] == 1) value[0] = "1|0"; else value[0] = "0|1";
			assert(head[pid] < pos.size());
			value[2] = std::to_string(pos[head[pid]]);
		}

		// Output phased VCF
		fprintf(f_out, "%s\t%s\t%s\t", CHROM, POS, ID);
		fprintf(f_out, "%s\t%s\t%s\t%s\t%s\t", REF, ALT, QUAL, FILTER, INFO);
		for (int i = 0; i < key.size(); i++)
			if (i == 0) fprintf(f_out, "%s", key[i].c_str());
			else fprintf(f_out, ":%s", key[i].c_str());
		fprintf(f_out, "\t");
		for (int i = 0; i < value.size(); i++)
			if (i == 0) fprintf(f_out, "%s", value[i].c_str());
			else fprintf(f_out, ":%s", value[i].c_str());
		fprintf(f_out, "\n");
		line_n++;
	}
	free(buf);
	pos.clear();
	gzclose(f_in); fclose(f_out);
	std::cerr << "The phased results are written into VCF file" << std::endl;
}

static uint8_t locate_first_snp(const SNP_Block &block, const std::vector<SNP> &ros, int rid) {
	for (int i = 0; i < block.n; i++) {
		const auto &reads_on_snp = ros[block.anchor + i];
		int l = 0, r = reads_on_snp.n - 1, ans = -1;
		while (l <= r) {
			int mid = (l + r) / 2;
			if (reads_on_snp.rid[mid] < rid) {
				l = mid + 1;
			} else if (reads_on_snp.rid[mid] > rid) {
				r = mid - 1;
			} else {
				ans = mid;
				break;
			}
		}
		if (ans == -1) break;
		if (reads_on_snp.bit[ans] == 2) continue;
		if (reads_on_snp.bit[ans] == block.bit[i]) return 1;
		else return 0;
	}
	return 2;
}

static uint8_t locate_last_snp(const Phased_Result &block, const std::vector<SNP> &ros, int rid) {
	for (int i = block.n-1; i >= 0; i--) {
		const auto &reads_on_snp = ros[block.index[i]];
		int l = 0, r = reads_on_snp.n - 1, ans = -1;
		while (l <= r) {
			int mid = (l + r) / 2;
			if (reads_on_snp.rid[mid] < rid) {
				l = mid + 1;
			} else if (reads_on_snp.rid[mid] > rid) {
				r = mid - 1;
			} else {
				ans = mid;
				break;
			}
		}
		if (ans == -1) break;
		if (reads_on_snp.bit[ans] == 2) continue;
		if (reads_on_snp.bit[ans] == block.bit[i]) return 1;
		else return 0;
	}
	return 2;
}

Phased_Result merge_blocks(const std::vector<SNP_Block> &blocks, const std::vector<SNP> &ros) {
	for (int i = 1; i < blocks.size(); i++) {
		assert(blocks[i].anchor >= blocks[i-1].anchor + blocks[i-1].n);
	}
	int memory_size = 0;
	for (const auto &b : blocks) memory_size += b.n;
	Phased_Result ret(memory_size);
	ret.add(blocks[0]);
	int merged_n = 1;
	for (int i = 1; i < blocks.size(); i++) {
		const auto &R1 = ros[ret.index[ret.n - 1]]; // The reads sit on the last SNP of the block i-1
		const auto &R2 = ros[blocks[i].anchor]; // The reads sit on the first SNP of the block i
		// Both the two arrays are sorted. It is easy to find the overlapped elements in a linear way
		int cnt_forward = 0, cnt_reverse = 0;
		int j = 0, k = 0;
		while (j < R1.n && k < R2.n) {
			if (R1.rid[j] < R2.rid[k]) j++;
			else if (R1.rid[j] > R2.rid[k]) k++;
			else {
				// In this loop, we process every read cross two blocks
				uint8_t judge1 = locate_last_snp(ret, ros, R1.rid[j]);
				uint8_t judge2 = locate_first_snp(blocks[i], ros, R2.rid[k]);
				if ((judge1 == 1 && judge2 == 1) || (judge1 == 0 && judge2 == 0)) cnt_forward++;
				if ((judge1 == 0 && judge2 == 1) || (judge1 == 1 && judge2 == 0)) cnt_reverse++;
				j++; k++;
			}
		}
		if (cnt_forward == 0 && cnt_reverse == 0) ret.add(blocks[i]);
		else if (cnt_forward >= cnt_reverse) ret.merge(blocks[i], false);
		else ret.merge(blocks[i], true);

		if (cnt_reverse == 0 && cnt_forward == 0) merged_n++;
	}
	std::cerr << "In the end, there are " << merged_n << " merged blocks" << std::endl;

	return ret;
}

void Phased_Result::output() {
	int blocks_n = 1;
	for (int i = 1; i < n; i++) if (head[i] != head[i-1]) blocks_n++;

	std::vector<int> pos[blocks_n];
	std::string base[blocks_n], rev_base[blocks_n];

	blocks_n = 0;
	for (int i = 0; i < n; i++) {
		if (i > 0 && head[i] != head[i-1]) blocks_n++;

		pos[blocks_n].push_back(global_snps[index[i]].pos);
		if (bit[i] == 1) {
			base[blocks_n] += global_snps[index[i]].alt;
			rev_base[blocks_n] += global_snps[index[i]].ref;
		} else {
			base[blocks_n] += global_snps[index[i]].ref;
			rev_base[blocks_n] += global_snps[index[i]].alt;
		}
	}

	for (int i = 0; i <= blocks_n; i++) {
		std::cout << "block" << 2*i+1 << "\t" << global_chrom.c_str() << "\t";
		for (int j = 0; j < pos[i].size(); j++) {
			if (j > 0) std::cout << "_";
			std::cout << pos[i][j];
		}
		std::cout << "\t" << base[i] << std::endl;

		std::cout << "block" << 2*i+2 << "\t" << global_chrom.c_str() << "\t";
		for (int j = 0; j < pos[i].size(); j++) {
			if (j > 0) std::cout << "_";
			std::cout << pos[i][j];
		}
		std::cout << "\t" << rev_base[i] << std::endl;
	}
}

void output_blocks(const std::vector<SNP_Block> &blocks) {
	for (int i = 0; i < blocks.size(); i++) {
		fprintf(stdout, "block%d\t", 2*i+1);
		fprintf(stdout, "%s\t", global_chrom.c_str());
		const auto &b = blocks[i];
		for (int j = 0; j < b.n; j++) {
			if (j > 0) fprintf(stdout, "_");
			assert(b.anchor + j < global_snps.size());
			fprintf(stdout, "%d", global_snps[b.anchor + j].pos);
		}
		fprintf(stdout, "\t");
		for (int j = 0; j < b.n; j++) {
			assert(b.anchor + j < global_snps.size());
			if (b.bit[j] == 0) fprintf(stdout, "%c", global_snps[b.anchor + j].ref);
			else fprintf(stdout, "%c", global_snps[b.anchor + j].alt);
		}
		fprintf(stdout, "\n");

		fprintf(stdout, "block%d\t", 2*i+2);
		fprintf(stdout, "%s\t", global_chrom.c_str());
		for (int j = 0; j < b.n; j++) {
			if (j > 0) fprintf(stdout, "_");
			assert(b.anchor + j < global_snps.size());
			fprintf(stdout, "%d", global_snps[b.anchor + j].pos);
		}
		fprintf(stdout, "\t");
		for (int j = 0; j < b.n; j++) {
			assert(b.anchor + j < global_snps.size());
			if ((b.bit[j] ^ 1) == 0) fprintf(stdout, "%c", global_snps[b.anchor + j].ref);
			else fprintf(stdout, "%c", global_snps[b.anchor + j].alt);
		}
		fprintf(stdout, "\n");
	}
}
