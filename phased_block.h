//
// Created by ixiaohu on 2022/1/30.
//

#ifndef KSNP_PHASED_BLOCK_H
#define KSNP_PHASED_BLOCK_H

#include <cstdlib>
#include <cstdint>
#include <vector>

#include "vcf_reader.h"

struct SNP_Block {
	int anchor;
	int n, m;
	uint8_t *bit;

	SNP_Block(int pos_0) {
		anchor = pos_0;
		n = 0;
		m = 64;
		bit = (uint8_t*) malloc(m * sizeof(uint8_t));
	}

	void add(uint8_t b) {
		if (n >= m) {
			m <<= 1;
			bit = (uint8_t*) realloc(bit, m * sizeof(uint8_t));
		}
		bit[n++] = b;
	}
};

struct Phased_Result {
	int n;
	int *index; // The index of SNP in VCF file
	uint8_t *bit;
	int *head;

	Phased_Result(int m) {
		n = 0;
		index = (int*) malloc(m * sizeof(int));
		bit = (uint8_t*) malloc(m * sizeof(uint8_t));
		head = (int*) malloc(m * sizeof(int));
	}

	void add(const SNP_Block &b);

	void merge(const SNP_Block &b, bool reversed);

	void write_vcf(const char *src, const char *dst);

	void output();

	void destroy() const {
		free(index); free(bit); free(head);
	}
};

std::vector<SNP_Block> remove_overlap(std::vector<SNP_Block> &blocks);

std::vector<SNP_Block> load_blocks_from_stdin();

void view_merged_blocks();

Phased_Result merge_blocks(const std::vector<SNP_Block> &blocks, const std::vector<SNP> &ros);

void output_blocks(const std::vector<SNP_Block> &blocks);

#endif //KSNP_PHASED_BLOCK_H
