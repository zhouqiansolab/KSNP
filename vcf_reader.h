//
// Created by ixiaohu on 2022/1/28.
//

#ifndef KSNP_VCF_READER_H
#define KSNP_VCF_READER_H

#include <vector>
#include <string>
#include <iostream>

struct SNP {
	int pos;
	char ref, alt;

	int n, m; // The reads on this SNP
	int *rid;
	uint8_t *bit;

	SNP(int p, char r, char a) {
		this->pos = p;
		this->ref = r;
		this->alt = a;

		this->n = 0;
		this->m = 16;
		this->rid = (int*) malloc(m * sizeof(int));
		this->bit = (uint8_t*) malloc(m * sizeof(uint8_t));
	}

	bool operator < (const SNP &o) const {
		return this->pos < o.pos;
	}

	void add(int id, uint8_t b) {
		if (n >= m) {
			m <<= 1;
			rid = (int*) realloc(rid, m * sizeof(int));
			bit = (uint8_t*) realloc(bit, m * sizeof(uint8_t));
		}
		rid[n] = id;
		bit[n] = b;
		n++;
	}

	void destroy() {
		free(rid);
		free(bit);
	}
};

std::vector<SNP> input_snp(const char *fn);

int binary_search(const std::vector<SNP> &snps, int pos);

#endif //KSNP_VCF_READER_H
