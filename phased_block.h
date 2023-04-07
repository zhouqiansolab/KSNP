//
// Created by ixiaohu on 2022/1/30.
//

#ifndef KSNP_PHASED_BLOCK_H
#define KSNP_PHASED_BLOCK_H

#include <cstdlib>
#include <cstdint>
#include <vector>

#include "snp_dbg.h"

class PostProcess {
private:
	std::vector<SNP> &snp_column;
	const std::vector<Read_Allele> &read_row;
	std::vector<Phased_Block> &blocks;
	std::vector<int8_t> &hint;

	void de_overlap();

	/** Calculate MEC score for a seed template TENC on block [p, p+TLEN) */
	int mec_of_template(Phased_Block &b, int p, int TLEN, uint tmp);
	void error_correction();

	/** Connect two neighbouring blocks in which way.
	 * @return 0  Connect directly
	 * @return 1  Connect with the second block reversed
	 * @return -1 Do NOT connect them */
	int8_t connect_block(const Phased_Block &prev, const Phased_Block &curr);
	void merge();

	bool try_to_switch(const Phased_Block &b, int p);
	bool try_to_flip(const Phased_Block &b, int p);

public:
	PostProcess(std::vector<SNP> &snp, const std::vector<Read_Allele> &read,
        std::vector<Phased_Block> &b, std::vector<int8_t> &h);

	void resolve();
};

void output_vcf_header(VCF_Header &header, const char *fn);

void output_haplotype_block(const std::vector<SNP> &vcf_snp, const char *dst_fn);

#endif //KSNP_PHASED_BLOCK_H
