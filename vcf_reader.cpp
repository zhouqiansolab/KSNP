//
// Created by ixiaohu on 2022/1/28.
//

#include <fstream>
#include <sstream>
#include <zlib.h>

#include "vcf_reader.h"

extern std::string global_chrom;

std::vector<SNP> input_snp(const char *fn) {
	std::vector<SNP> ret;
	gzFile in = gzopen(fn, "r");
	if (in == nullptr) {
		std::cerr << "ERR: Open " << fn << " failed" << std::endl;
		std::abort();
	}
	const int SNP_BUF_SIZE = 4 * 1024 * 1024;
	char *buf = new char[SNP_BUF_SIZE];
	while (gzgets(in, buf, SNP_BUF_SIZE)) {
		if (buf[0] == '#') continue;
		std::istringstream iss(buf);
		std::string chrom; iss >> chrom;
		int pos; iss >> pos;
		std::string dummy; iss >> dummy;
		char ref; iss >> ref;
		char alt; iss >> alt;
		if (global_chrom.empty()) global_chrom = chrom;
		else if (global_chrom != chrom) {
			std::cerr << "ERR: Multiple chromosomes found" << std::endl;
			std::abort();
		}
		ret.emplace_back(SNP(pos, ref, alt));
	}
	delete [] buf;
	return ret;
}

int binary_search(const std::vector<SNP> &snps, int pos) {
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