#include <getopt.h>
#include <cassert>

#include "ksnp_reader.h"
#include "snp_dbg.h"
#include "phased_block.h"
#include "time_stamp.h"

#define VERSION  "1.3"
#define DEV_MODE "RELEASE"

static int usage() {
	fprintf(stderr, "Usage: ksnp -b <BAM> -r <FASTA> -v <VCF> -o <Output>\n");
	fprintf(stderr, "  -b aligned reads in BAM format\n");
	fprintf(stderr, "  -r reference sequence for allele realignment in FASTA format\n");
	fprintf(stderr, "  -v heterozygous variants to phase in VCF format\n");
	fprintf(stderr, "  -o output file that phased results are written to (stdout)\n");
	fprintf(stderr, "  -k k-mer size to construct DBG, currently supporting 2,3,4,5 (2)\n");
	fprintf(stderr, "Version: %s (%s)\n", VERSION, DEV_MODE);
	return 1;
}

int main(int argc, char *argv[]) {
	if (argc == 1) return usage();
	TimeStamp stamp; stamp.push();

	int c, K = 2; // default value K = 2
	const char *bam_fn = nullptr, *vcf_fn = nullptr,  *ref_fn = nullptr, *output_fn = nullptr;
	const char *request_chromosome = nullptr;
	const char *resolution_fn = nullptr;
	const char *og_prefix = nullptr;
	while ((c = getopt(argc, argv, "b:v:o:c:k:r:l:R:")) >= 0) {
		if (c == 'b') {
			bam_fn = optarg;
		} else if (c == 'v') {
			vcf_fn = optarg;
		} else if (c == 'o') {
			output_fn = optarg;
		} else if (c == 'c') {
			request_chromosome = optarg;
		} else if (c == 'k') {
			K = atoi(optarg);
		} else if (c == 'r') {
			ref_fn = optarg;
		}
		// Development options
		else if (c == 'l') {
			og_prefix = optarg; // output graph of each step
		} else if (c == 'R') {
			resolution_fn = optarg; // input allele instead of BAM (for faster input)
		} else return usage();
	}

	if (bam_fn == nullptr) { fprintf(stderr, "ERR: please assign a BAM file of aligned reads\n"); return 1; }
	if (vcf_fn == nullptr) { fprintf(stderr, "ERR: please choose a VCF file to phase\n"); return 1; }
	if (ref_fn == nullptr) {
		fprintf(stderr, "ERR: please provide a reference sequence for detecting alleles by realignment\n" );
		return 1;
	}
	if (K < 2 or K > 5) {
		fprintf(stderr, "MSG: invalid K value is given. K is set to 2.\n");
		K = 2;
	}

	// Input all variants to phase
	stamp.push();
	auto variant_table = input_vcf(vcf_fn, request_chromosome);
	stamp.record("1 Input_VCF");

	// Phase SNPs of each chromosome
	FASTA_Reader ref_reader(ref_fn);
	output_vcf_header(variant_table.header, output_fn);
	for (int i = 0; i < variant_table.size; i++) {
		const auto &chr_name = variant_table.chromosome[i];
		auto &snp_column = variant_table.variants[i];
		fprintf(stderr, "Phase %ld SNPs on chromosome %s\n", snp_column.size(), chr_name.c_str());

		// Detecting alleles
		stamp.push();
		int length = ref_reader.get_length(chr_name); assert(length > 0);
		char *sequence = ref_reader.get_contig(chr_name);
		std::vector<Read_Allele> read_row;
		if (!resolution_fn) read_row = detect_allele(bam_fn, chr_name, snp_column, length, sequence);
		else read_row = input_allele(resolution_fn, snp_column);
		delete [] sequence;
		stamp.record("2 Detect_Allele");

		// DBG
		stamp.push();
		SNP_DBG dbg(K, snp_column, read_row);
		if (og_prefix) dbg.output_graph_to(og_prefix);
		dbg.resolve();
		auto hint_post = dbg.get_hint_post();
		auto blocks = dbg.haplotype_block();
		dbg.destroy();
		// Post Process
		PostProcess post(snp_column, read_row, blocks, hint_post);
		post.resolve();
		stamp.record("3 DBG");

		// Output phased result
		stamp.push();
		output_haplotype_block(snp_column, output_fn);
		stamp.record("4 Output_VCF");
	}
	ref_reader.close();

	stamp.record("5 Total");
//	stamp.output();
	return 0;
}
