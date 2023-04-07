//
// Created by ixiaohu on 2022/11/16.
//

#ifndef KSNP_DEV_KSNP_READER_H
#define KSNP_DEV_KSNP_READER_H

#include <vector>
#include <iostream>
#include <map>
#include <zlib.h>
#include <cstring>

/**************
 *    VCF     *
 **************/
const int VCF_CHROM  = 0;
const int VCF_POS    = 1;
const int VCF_ID     = 2;
const int VCF_REF    = 3;
const int VCF_ALT    = 4;
const int VCF_QUAL   = 5;
const int VCF_FILTER = 6;
const int VCF_INFO   = 7;
const int VCF_FORMAT = 8;
const int VCF_SAMPLE = 9;

class VCF_Header {
private:
	int line_n;
	std::vector<std::string> keys; /** INFO, FILTER, FOMRAT, contig ... */
	std::vector<std::string> values;
	std::string sample_line; /** The last line includes sample information */

public:
	VCF_Header(): line_n(0) {}

	/** Parse a header line in VCF file */
	void add_line(const char *buf);

	/** Append FORMAT=<ID=GT> and FORMAT=<ID=PS> if absent in original file */
	void annotate_for_phasing();

	/** Return parsed content in text format */
	std::vector<std::string> to_str();
};

struct SNP {
	int pos;
	char ref;
	char alt;
	char *line;

	/** Detected alleles from aligned reads */
	std::vector<int> rid; /** Read indices on this SNP */
	std::vector<int> allele; /** Alleles of reads. 0:REF; 1:ALT; -1:gap */

	/** Phased result */
	int ps; /** Phased set: the first SNP position of the haplotype block */
	int gt; /** GT=0 for 0|1; GT=1 for 1|0; GT=-1 for unknown */

	SNP() { pos = -1; ref = alt = 'N'; ps = gt = 0; } // FIXME: this should not exist

	SNP(int p, char r, char a, const char *b):pos(p), ref(r), alt(a) {
		ps = -1;
		gt = -1;
		if (b) line = strdup(b);
		else line = nullptr;
	}

	inline void add_read(int r, int a) {
		rid.push_back(r);
		allele.push_back(a);
	}

	inline int size() const { return (int)rid.size(); }

	bool operator < (const SNP &o) const { return this->pos < o.pos; }
};

/** Table of chromosomes and variants */
struct Variant_Table {
	int size;
	VCF_Header header;
	std::vector<std::string> chromosome;
	std::vector<std::vector<SNP>> variants;
};
Variant_Table input_vcf(const char *fn, const char *req_chrom);

std::vector<SNP> load_phased_snp(const char *fn); /** Load phased SNPs for debug KSNP */

std::vector<std::string> split_str(const char *s, char sep);
std::string join_str(const std::vector<std::string> &a, char sep);

/**************
 *   FASTA    *
 **************/
struct FAIDX_Contig {
	int length; /** Length of the contig (bp) */
	long offset; /** Offset of the first base (0-based) */
	int line_bases; /** Bases in each line (except to last line) */
	int line_width; /** Width of each line (except to last line) */
	FAIDX_Contig():length(0), offset(0), line_bases(0), line_width() {}
};

class FASTA_Reader {
private:
	gzFile_s *ref_fp;
	std::map<std::string, FAIDX_Contig> dict; /** Dictionary from chromosome to offset */

public:
	explicit FASTA_Reader(const std::string &fn);

	inline int get_length(const std::string &chr_name) {
		return dict.find(chr_name) != dict.end() ?dict[chr_name].length :-1;
	}

	char* get_contig(const std::string &chr_name);

	inline void close() { gzclose(ref_fp); }
};

struct Sequence {
	char *data;
	char *name, *comment, *seq, *qual;
	int len; // length of sequencing fragment
	Sequence() { data = name = comment = seq = qual = nullptr; len = 0; }
	void destroy() const { if (data) free(data); }
};

class FASTQ_Reader {
private:
	const int SEP_SPACE = 0; // isspace(): \t, \n, \v, \f, \r
	const int SEP_TAB   = 1; // isspace() && !' '
	const int SEP_LINE  = 2; // line separator: "\n" (Unix) or "\r\n" (Windows)
	const int MODE_WRITE = 0;
	const int MODE_APPEND = 1;
	const int CHUNK_SIZE = 16 * 1024; // 16KB buffer size
	char *stream_buf; // fixed input queue stream buffer
	int begin, end, is_eof; // [begin, end) is the stream interval pointer
	gzFile_s *fp; // the input file; is_eof marks its EOF

	struct String_Buffer {
		int l, m;
		char *s;
		String_Buffer(): l(0), m(0), s(nullptr){}
	};
	String_Buffer name, comment, seq, qual; // recycled buffers for reading FASTA/Q
	int last_header; // The header character has been in the previous call, @ for FASTQ, > for FASTA

	// fetch one character and return it, return -1 if EOF
	inline int getc();

	// fetch a string until meeting specified delimiter; if EOF return -1, or return string length.
	// @found_del is the found delimiter, as space could be any character among \t \n \r ...
	int gets(String_Buffer &str, int delimiter, int *found_del, int mode);

	Sequence dump() const;

public:
	explicit FASTQ_Reader(const char *fn);

	// end of input file and stream buffer is also empty
	bool eof() const { return (is_eof and begin >= end); }

	// Return data=NULL if eof()=true
	Sequence read1();

	void close();
};

/**************
 *    BAM     *
 **************/
struct Allele_Call {
	int que_pos; /** Allele position on query read */
	int snp_idx; /** Variant(reference and alternative allele) index in VCF */
	int allele;  /** 0: read supports REF allele, 1: read support ALT allele, -1: gap or unknown. */
	Allele_Call (int q, int s, int a): que_pos(q), snp_idx(s), allele(a) {}
};
typedef std::vector<Allele_Call> Read_Allele; /** All alleles on a read */

/**
 * Detect alleles by realignment
 * @param bam_fn aligned reads
 * @param chr_name chromosome name
 * @param snps variants to phase
 * @param len reference sequence length
 * @param seq reference sequence
 * @return all alleles on informative reads
 */
std::vector<Read_Allele> detect_allele(const char *bam_fn, const std::string &chr_name,
									   std::vector<SNP> &snps, int len, const char *seq);

std::vector<Read_Allele> input_allele(const char *fn, std::vector<SNP> &snps);

#endif //KSNP_DEV_KSNP_READER_H
