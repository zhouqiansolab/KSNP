#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>
#include <sys/time.h>
#include <getopt.h>

#include "snp_dbg.h"

std::string global_chrom; // The program only processes one chromosome each time

static int usage() {
	std::cerr << "Usage: samtools view in.bam | ksnp -k <k-mer size> -v <VCF file> -o <output file>" << std::endl;
	std::cerr << "K-mer size supports INT value from 2 to 10 (default value is 2)" << std::endl;
	return 1;
}

double realtime(){
	struct timeval tp;
	struct timezone tzp;
	gettimeofday(&tp, &tzp);
	return tp.tv_sec + tp.tv_usec * 1e-6;
}

int main(int argc, char *argv[]) {
	if (argc == 1) return usage();
	double time_s = realtime();
	int c, K = 2;
	const char *vcf_fn = nullptr, *output_fn = nullptr;
	while ((c = getopt(argc, argv, "k:v:o:")) >= 0) {
		if (c == 'k') {
			K = atoi(optarg);
		} else if (c == 'v') {
			vcf_fn = optarg;
		} else if (c == 'o') {
			output_fn = optarg;
		} else return usage();
	}

	if (K < 2 || K > 10) return usage();
	if (vcf_fn == nullptr || output_fn == nullptr) return usage();

	auto snps = input_snp(vcf_fn);
	std::sort(snps.begin(), snps.end());
	std::cerr << "Input " << snps.size() << " SNPs" << std::endl;

	SNP_DBG dbg(K, snps.size());
	dbg.construct_graph(snps);
	dbg.trim_graph();
	dbg.remove_tip();
	dbg.remove_bubble();
	dbg.remove_tip_again();
	auto blocks = dbg.snp_haplotype();
	dbg.destroy();

	auto ro_blocks = remove_overlap(blocks);
	auto results = merge_blocks(ro_blocks, snps);
	results.write_vcf(vcf_fn, output_fn);

	// Free memory for blocks and resolution reads
	for (auto &b : ro_blocks) free(b.bit); ro_blocks.clear();
	for (auto &s : snps) s.destroy(); snps.clear();
	results.destroy();

	double time_e = realtime();
	std::cerr << "In total, " << time_e-time_s << " seconds elapsed" << std::endl;
	return 0;
}
