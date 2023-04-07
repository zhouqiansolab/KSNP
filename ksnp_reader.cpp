//
// Created by ixiaohu on 2022/11/16.
//

#include <sstream>
#include <algorithm>
#include <memory.h>
#include <cassert>
#include <fstream>

#include "ksnp_reader.h"
#include "htslib/sam.h"
#include "realignment.h"
#include "time_stamp.h"

void VCF_Header::add_line(const char *buf) {
	if (buf[0] == '#' and buf[1] == '#') {
		std::string key, value;
		int i = 2; for (; buf[i] != '='; i++) key += buf[i];
		i++; for (; buf[i] !='\0'; i++) value += buf[i];
		keys.push_back(key);
		values.push_back(value);
		line_n++;
	} else if (buf[0] == '#') {
		auto s = split_str(buf, '\t');
		if (s.size() - VCF_SAMPLE != 1) {
			fprintf(stderr, "ERR: %ld samples are found, but KSNP is designed for a single individual\n",
		        s.size() - VCF_SAMPLE);
		}
		sample_line = buf;
	}
}

void VCF_Header::annotate_for_phasing() {
	bool has_gt = false, has_ps = false;
	for (int i = 0; i < line_n; i++) {
		if (keys[i] == "FORMAT") {
			if (values[i].find("ID=GT") != std::string::npos) has_gt = true;
			if (values[i].find("ID=PS") != std::string::npos) has_ps = true;
		}
	}
	if (not has_gt) add_line("<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
	if (not has_ps) add_line("<ID=PS,Number=1,Type=Integer,Description=\"Phase Set\">");
}

std::vector<std::string> VCF_Header::to_str() {
	std::vector<std::string> ret;
	for (int i = 0; i < line_n; i++) {
		ret.push_back("##" + keys[i] + "=" + values[i]);
	}
	ret.push_back(sample_line);
	return ret;
}

Variant_Table input_vcf(const char *fn, const char *req_chrom) {
	gzFile in = gzopen(fn, "r");
	if (in == nullptr) {
		fprintf(stderr, "ERR: Can not open VCF file %s\n", fn);
		std::abort();
	}

	Variant_Table table; table.size = 0;
	std::map<std::string,int> dict;
	const int SNP_BUF_SIZE = 4 * 1024 * 1024;
	char *buf = new char[SNP_BUF_SIZE];
	while (gzgets(in, buf, SNP_BUF_SIZE) != nullptr) {
		int len = strlen(buf); if (buf[len-1] == '\n') buf[len-1] = '\0';
		if (buf[0] == '#') { table.header.add_line(buf); continue; }
		auto fields = split_str(buf, '\t');
		const auto &chrom = fields[VCF_CHROM];
		if (req_chrom and chrom != std::string(req_chrom)) continue;
		int pos = atoi(fields[VCF_POS].c_str());
		if (fields[VCF_REF].size() != 1 or fields[VCF_ALT].size() != 1) continue; // Not a SNV
		char ref = fields[VCF_REF][0];
		char alt = fields[VCF_ALT][0];
		if (dict.find(chrom) == dict.end()) {
			table.chromosome.push_back(chrom);
			table.variants.emplace_back(std::vector<SNP>());
			table.size++;
			dict[chrom] = table.size - 1;
		}
		table.variants[dict[chrom]].emplace_back(SNP(pos, ref, alt, buf));
	}
	delete [] buf;
	if (table.size == 0) {
		fprintf(stderr, "ERR: input no variants to phase\n");
		std::abort();
	}

	int snp_n = 0;
	for (auto &variant: table.variants) {
		snp_n += variant.size();
		std::sort(variant.begin(), variant.end());
	}
	return table;
}

std::vector<SNP> load_phased_snp(const char *fn) {
	std::vector<SNP> ret;
	gzFile in = gzopen(fn, "r");
	if (in == nullptr) return ret;
	const int SNP_BUF_SIZE = 4 * 1024 * 1024;
	char *buf = new char[SNP_BUF_SIZE];
	while (gzgets(in, buf, SNP_BUF_SIZE) != nullptr) {
		if (buf[0] == '#') continue;
		int len = strlen(buf);
		if (buf[len-1] == '\n') buf[len-1] = '\0';
		auto fields = split_str(buf, '\t');
		int pos = atoi(fields[VCF_POS].c_str());
		char ref = fields[VCF_REF][0];
		char alt = fields[VCF_ALT][0];
		SNP s(pos, ref, alt, nullptr);
		auto fmt = split_str(fields[VCF_FORMAT].c_str(), ':');
		auto smp = split_str(fields[VCF_SAMPLE].c_str(), ':');
		for (int i = 0; i < fmt.size(); i++) {
			if (fmt[i] == "GT") {
				if (smp[i] == "0|1") s.gt = 0;
				else if (smp[i] == "1|0") s.gt = 1;
			} else if (fmt[i] == "PS") {
				if (smp[i] != ".") s.ps = atoi(smp[i].c_str());
			}
		}
		if (s.ps == -1 or s.gt == -1) continue;
		ret.push_back(s);
	}
	delete [] buf;
	return ret;
}

std::vector<std::string> split_str(const char *s, char sep) {
	std::vector<std::string> ret;
	std::string temp;
	for (int i = 0; true; i++) {
		if (s[i] == '\0') {
			if (!temp.empty()) ret.push_back(temp);
			break;
		} else if (s[i] == sep) {
			if (!temp.empty()) ret.push_back(temp);
			temp = "";
		} else temp += s[i];
	}
	return ret;
}

std::string join_str(const std::vector<std::string> &a, char sep) {
	std::string ret;
	for (int i = 0; i < a.size(); i++) {
		if (i > 0) ret += sep;
		ret += a[i];
	}
	return ret;
}

FASTA_Reader::FASTA_Reader(const std::string &fn) {
	ref_fp = gzopen(fn.c_str(), "r");
	if (ref_fp == nullptr) {
		fprintf(stderr, "ERR: open reference file %s failed\n", fn.c_str());
		std::abort();
	}
	std::ifstream in(fn + ".fai");
	if (!in.is_open()) {
		fprintf(stderr, "ERR: reference sequence is not indexed. Use `samtools faidx`\n");
		std::abort();
	}
	auto *line_buffer = new char[10 * 1024];
	while (in.getline(line_buffer, 10 * 1024)) {
		std::stringstream ss(line_buffer);
		std::string chr_name;
		FAIDX_Contig t;
		ss >> chr_name >> t.length >> t.offset >> t.line_bases >> t.line_width;
		dict[chr_name] = t;
	}
	delete [] line_buffer;
	in.close();
}

char* FASTA_Reader::get_contig(const std::string &chr_name) {
	const auto &c = dict[chr_name];
	gzseek(ref_fp, c.offset, SEEK_SET);
	auto *sequence = new char[c.length + 5]; int length = 0;
	auto *buf = new char[c.line_width + 5];
	while (length < c.length) {
		 gzgets(ref_fp, buf, c.line_width + 5);
		 int n = std::min(c.line_bases, c.length - length);
		 memcpy(sequence + length, buf, n);
		 length += n;
	}
	delete [] buf;
	assert(length == c.length);
	return sequence;
}

FASTQ_Reader::FASTQ_Reader(const char *fn) {
	fp = gzopen(fn, "r");
	if (fp == nullptr) {
		fprintf(stderr, "ERR: open file %s failed\n", fn);
		std::abort();
	}
	begin = end = is_eof = 0;
	stream_buf = (char*) malloc(CHUNK_SIZE);
	last_header = 0;
}

int FASTQ_Reader::getc() {
	if (eof()) return -1;
	if (begin >= end) { // If queue buffer is empty, read a chunk of data again
		begin = 0;
		end = gzread(fp, stream_buf, CHUNK_SIZE);
		if (end == 0) {
			is_eof = 1;
			return -1;
		}
	}
	return stream_buf[begin++];
}

static inline int roundup(uint32_t x) { // Rounded to the next closest 2^k
	x--;
	x |= (x>>1u); x |= (x>>2u);
	x |= (x>>4u); x |= (x>>8u);
	x |= (x>>16u);
	x++;
	return x;
}

int FASTQ_Reader::gets(String_Buffer &str, int delimiter, int *found_del, int mode) {
	bool got_any = false;
	if (mode == MODE_WRITE) str.l = 0; // Reset string buffer
	while (true) {
		if (begin >= end) {
			if (!is_eof) { // Reset the queue
				begin = 0;
				end = gzread(fp, stream_buf, CHUNK_SIZE);
				if (end == 0) { is_eof = 1; break; }
			} else break;
		}
		got_any = true;
		int i;
		if (delimiter == SEP_LINE) {
			for (i = begin; i < end; i++) if (stream_buf[i] == '\n') break;
		}  else if (delimiter == SEP_SPACE) {
			for (i = begin; i < end; i++) if (isspace(stream_buf[i])) break;
		} else if (delimiter == SEP_TAB) {
			for (i = begin; i < end; i++)
				if (isspace(stream_buf[i] and stream_buf[i] != ' ')) break;
		} else { // the character delimiter itself
			for (i = begin; i < end; i++) if (stream_buf[i] == delimiter) break;
		}
		int input_n = (i - begin); // Do not include the delimiter
		if (str.m - str.l < input_n + 1) { // Expand the string buffer
			str.m = str.l + input_n + 1;
			str.m = roundup(str.m);
			str.s = (char*) realloc(str.s, str.m * sizeof(char));
		}
		memcpy(str.s + str.l, stream_buf + begin, input_n);
		str.l = str.l + input_n;
		begin = i + 1; // Move the queue front forward
		if (i < end) { // Meet the delimiter
			// Return the found delimiter among several possible characters
			if (found_del) *found_del = ((int)stream_buf[i]);
			break;
		}
	}
	if (!got_any and eof()) return -1; // Reach the end of file
	if (str.s == nullptr) {
		// Input file does not contain any other characters except the delimiter
		// Return an empty string
		str.m = 1;
		str.s = (char*) calloc(str.m, sizeof(char));
	} else if (delimiter == SEP_LINE and str.l > 1 and str.s[str.l-1] == '\r') {
		// This is for Windows OS '\r\n'
		str.l--;
	}
	str.s[str.l] = '\0'; // The end of string
	return str.l;
}

Sequence FASTQ_Reader::dump() const {
	int m = name.l + comment.l + seq.l + qual.l + 4, p = 0; // Total memory space

	Sequence ret;
	ret.data = (char*) calloc(m, sizeof(char)); // a chunk of memory allocated
	ret.name = ret.data + p; memcpy(ret.name, name.s, name.l); p += name.l + 1;
	ret.comment = ret.data + p; memcpy(ret.comment, comment.s, comment.l); p += comment.l + 1;
	ret.seq = ret.data + p; memcpy(ret.seq, seq.s, seq.l); p += seq.l + 1;
	ret.qual = ret.data + p; memcpy(ret.qual, qual.s, qual.l);
	ret.len = seq.l;
	if (comment.l == 0) ret.comment = nullptr; // comment and qual coule be NULL
	if (qual.l == 0) ret.qual = nullptr;
	return ret;
}

Sequence FASTQ_Reader::read1() {
	Sequence null;
	int c;
	if (last_header == 0) { // jump to the next QNAME line
		// There might be empty lines between FASTA/Q records
		while ((c = getc()) != -1 and c != '>' and c != '@');
		if (c == -1) return null;
		last_header = c;
	}
	name.l = comment.l = seq.l = qual.l = 0; // reset all members
	/** QNAME */
	int found_del;
	if (gets(name, SEP_SPACE, &found_del, MODE_WRITE) == -1) return null; // normal exit: EOF

	/** COMMENT */
	if (found_del != '\n') gets(comment, SEP_LINE, nullptr, MODE_WRITE); // FASTA/Q comment

	/** BASE */
	// Read until the header of next FASTA/Q record or quality string, which are >@+
	// Note that there are possibly many '\n' within the base string, especially in reference FASTA file
	if (seq.s == nullptr) { // initially allocate a chunk of memory
		seq.m = 256;
		seq.s = (char*) malloc(seq.m * sizeof(char));
	}
	while ((c = getc()) != -1 and c != '>' and c != '@' and c != '+') {
		if (c == '\n') continue; // skip empty lines
		seq.s[seq.l++] = (char)c; // this is safe because the '\0' space allows directly putting this character
		gets(seq, SEP_LINE, nullptr, MODE_APPEND); // read the rest of the line
	}
	if (seq.l + 1 >= seq.m) { // seq->seq.s[seq->seq.l] below may be out of boundary
		seq.m = seq.l + 2;
		seq.m = roundup(seq.m);
		seq.s = (char*) realloc(seq.s, seq.m);
	}
	seq.s[seq.l] = '\0';

	// Character + indicates FASTQ file
	if (c == '>' or c == '@') last_header = c; // meet the header of next FASTA record
	if (c != '+') return dump(); // FASTA header or EOF

	/** QUAL */
	if (qual.m < seq.m) {
		qual.m = seq.m;
		qual.s = (char*) realloc(qual.s, qual.m * sizeof(char));
	}
	while ((c = getc()) != -1 and c != '\n'); // skip until '\n'
	if (c == -1) {
		fprintf(stderr, "[%s] Error: absence of quality string\n", __func__);
		std::abort();
	}
	while (gets(qual, SEP_LINE, nullptr, MODE_APPEND) != -1 and qual.l < seq.l);
	if (qual.l != seq.l) {
		fprintf(stderr, "[%s] Error: %d bases but only %d quality scores found\n", __func__, seq.l, qual.l);
		std::abort();
	}

	last_header = 0; // have not meet the header of next FASTQ record
	while ((c = getc()) != -1 and c != '>' and c != '@'); // this sets eof()=true before returning an empty read
	if (c != -1) last_header = c;
	return dump();
}

void FASTQ_Reader::close() {
	free(name.s); free(comment.s); free(seq.s); free(qual.s);
	free(stream_buf);
	gzclose(fp);
}

static int binary_search_snp(const std::vector<SNP> &snps, int pos) {
	int l = 0, r = snps.size() - 1, ret = -1;
	while (l <= r) {
		int mid = (l + r) / 2;
		if (snps[mid].pos < pos) l = mid + 1;
		else {
			ret = mid;
			r = mid - 1;
		}
	}
	return ret;
}

std::vector<Read_Allele> detect_allele(const char *bam_fn, const std::string &chr_name,
									   std::vector<SNP> &snps, int len, const char *seq) {
	Time1 bam_time, dal_time, tick_time;
	std::vector<Read_Allele> ret; int total_allele = 0;
	tick_time = TimeStamp::get_time();
	samFile *bam_fp = sam_open(bam_fn, "r");
	if (bam_fn == nullptr) {
		fprintf(stderr, "ERR: can not open BAM file %s\n", bam_fn);
		std::abort();
	}
	bam_hdr_t *bam_header = sam_hdr_read(bam_fp);
	if (bam_header == nullptr) {
		fprintf(stderr, "ERR: can not read header in BAM file\n");
		std::abort();
	}
	hts_idx_t *bam_idx = bam_index_load(bam_fn);
	if (bam_idx == nullptr) {
		fprintf(stderr, "ERR: can not load BAM index; use `samtools index`\n");
		std::abort();
	}
	bam1_t *aln = bam_init1();
	hts_itr_t *iter = sam_itr_querys(bam_idx, bam_header, chr_name.c_str());
	if (iter == nullptr) {
		fprintf(stderr,"ERR: invalid region for chromosome %s\n", chr_name.c_str());
		std::abort();
	}
	bam_time += TimeStamp::get_time() - tick_time;
	Realignment realign(len, seq);
	while (true) {
		tick_time = TimeStamp::get_time();
		int got_any = sam_itr_next(bam_fp, iter, aln);
		bam_time += TimeStamp::get_time() - tick_time;
		if (got_any < 0) break;

		tick_time = TimeStamp::get_time();
		int ref_start = aln->core.pos + 1;
		int bs = binary_search_snp(snps, ref_start);
		if (bs == -1) continue;

		Read_Allele real;
		int read_len = 0, ref_len = 0;
		const uint32_t *cigar_array = bam_get_cigar(aln);
		for (int i = 0; i < aln->core.n_cigar; i++) {
			char op_chr = bam_cigar_opchr(cigar_array[i]);
			int op_len = bam_cigar_oplen(cigar_array[i]);
			if (op_chr != 'D' && op_chr != 'H') read_len += op_len;
			if (op_chr != 'I' && op_chr != 'S' && op_chr != 'H') ref_len += op_len;
		}

		int cid = 0; // Cigar iterator
		int que_pointer = 0, ref_pointer = ref_start; // Pointers sliding the aligned window
		for (int i = bs; i < snps.size(); i++) {
			auto &snp = snps[i];
			if (snp.pos >= ref_start + ref_len) break;

			// Find the cigar interval overlapping the SNP.
			// The cigar intervals are consecutive.
			// For example, 12D2X19=3I has cigar intervals [0,12), [12, 14), [14, 33), [33, 33).
			while (cid < aln->core.n_cigar) {
				char op_chr = bam_cigar_opchr(cigar_array[cid]);
				int op_len = bam_cigar_oplen(cigar_array[cid]);
				if (op_chr == 'H') { cid++; continue; }
				if (op_chr == 'I' || op_chr == 'S') {
					que_pointer += op_len;
					cid++;
					continue;
				}
				if (ref_pointer + op_len > snp.pos) break; // Cigar interval overlapping the SNP found
				ref_pointer += op_len;
				if (op_chr != 'D') que_pointer += op_len;
				cid++;
			}

			char op_chr = bam_cigar_opchr(cigar_array[cid]);
			int que_pos = op_chr == 'D' ?que_pointer - 1 : que_pointer + snp.pos - ref_pointer;
			auto pair = realign.bit_vector_dp(aln, que_pos, snp.pos - 1, snp.alt);
			int allele;
			if (pair.first < pair.second) allele = 0;
			else if (pair.first > pair.second) allele = 1;
			else allele = -1;
			real.emplace_back(Allele_Call(que_pos, i, allele));
		}

		// Remove marginal gaps
		int l_active = -1, r_active = -2;
		for (int i = 0; i < real.size(); i++) {
			const auto &v = real[i];
			if (v.allele == -1) continue;
			if (l_active == -1) l_active = i;
			r_active = i;
		}
		if (l_active == -1) continue;
		Read_Allele no_gap_real;
		for (int i = l_active; i <= r_active; i++) no_gap_real.push_back(real[i]);
		real = no_gap_real;

		// Only push back informative reads
		if (real.size() >= 2) {
			total_allele += real.size();
			for (const auto &v : real) snps[v.snp_idx].add_read(ret.size(), v.allele);
			ret.push_back(real);
		}
		dal_time += TimeStamp::get_time() - tick_time;
	}

	bam_destroy1(aln);
	bam_hdr_destroy(bam_header);
	sam_close(bam_fp);
	fprintf(stderr, "Detected %d alleles on %ld informative reads\n", total_allele, ret.size());
	fprintf(stderr, "    Decompress BAM costs %.3f CPU and %.3f real seconds\n", bam_time.cpu, bam_time.real);
	fprintf(stderr, "    Realignment costs %.3f CPU and %.3f real seconds\n", dal_time.cpu, dal_time.real);
	return ret;
}

std::vector<Read_Allele> input_allele(const char *fn, std::vector<SNP> &snps) {
	std::vector<Read_Allele> ret; int total_allele = 0;
	std::ifstream in(fn); assert(in.is_open());
	int n = 0; in >> n;
	for (int i = 0; i < n; i++) {
		Read_Allele read, temp;
		std::string qname; int m;
		in >> qname >> m;
		for (int j = 0; j < m; j++) {
			int pos, allele;
			in >> pos >> allele;
			pos++;
			int snp_idx = binary_search_snp(snps, pos); assert(snp_idx != -1);
			temp.emplace_back(Allele_Call(-1, snp_idx, allele));
		}

		int l_active = -1, r_active = -2; // Remove marginal gaps
		for (int j = 0; j < temp.size(); j++) {
			const auto &v = temp[j];
			if (v.allele == -1) continue;
			if (l_active == -1) l_active = j;
			r_active = j;
		}
		for (int j = l_active; j <= r_active; j++) read.push_back(temp[j]);

		if (read.size() >= 2) {
			assert(read.front().allele != -1 and read.back().allele != -1);
			total_allele += read.size();
			for (auto x : read) snps[x.snp_idx].add_read(ret.size(), x.allele);
			ret.push_back(read);
		}
	}
	in.close();
	fprintf(stderr, "Detected %d alleles on %ld informative reads\n", total_allele, ret.size());
	return ret;
}