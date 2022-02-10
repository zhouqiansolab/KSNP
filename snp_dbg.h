//
// Created by ixiaohu on 2022/1/28.
//

#ifndef KSNP_SNP_DBG_H
#define KSNP_SNP_DBG_H

#include <vector>

#include "vcf_reader.h"
#include "phased_block.h"

struct SNP_DBG {
	int snp_n;
	int K;
	int *node_cnt;
	int *edge_cnt;
	int *weight;

	SNP_DBG (int k, int n) {
		snp_n = n;
		K = k;
		node_cnt = (int*) calloc((1<<k) * n, sizeof(int));
		edge_cnt = (int*) calloc((1<<(k+1)) * n, sizeof(int));
		weight = (int*) calloc((1<<(k+1)) * n, sizeof(int));
	}

	void add_node(int p, uint32_t e) const {
		(node_cnt + (1<<K) * p)[e]++;
	}

	void add_edge(int p, uint32_t e) const {
		(edge_cnt + (1<<(K+1)) * p)[e]++;
	}

	void construct_graph(std::vector<SNP> &snps) const;

	void load_graph_from_stdin(const std::vector<SNP> &snps) const;

	inline bool filter_edge(int heaviest, int edge_weight) ;

	void trim_graph();

	void remove_tip();

	void remove_bubble();

	void remove_tip_again();

	void check_graph_symmetry();

	inline int out_degree(int p, int node) const {
		if (p == snp_n-1) return 0;
		int ret = 0;
		if (weight[(p+1) * (1<<(K+1)) + ((node<<1) | 0)] != 0) ret++;
		if (weight[(p+1) * (1<<(K+1)) + ((node<<1) | 1)] != 0) ret++;
		return ret;
	}

	inline int in_degree(int p, int node) const {
		if (p == 0) return 0;
		int ret = 0;
		if (weight[p * (1<<(K+1)) + ((0<<K) | node)] != 0) ret++;
		if (weight[p * (1<<(K+1)) + ((1<<K) | node)] != 0) ret++;
		return ret;
	}

	SNP_Block traversal_block(int p, int u);

	std::vector<SNP_Block> snp_haplotype();

	void output(const std::vector<SNP> &snps);

	void destroy() const {
		free(node_cnt);
		free(edge_cnt);
		free(weight);
	}
};

// The SAM structure only keeps necessary fields
struct SAM_Record {
	char *que_name;
	int flag;
	char *ref_name;
	int ref_start;
	char *cigar;
	char *bases;

	SAM_Record(char *s);
};

struct Parsed_CIGAR {
	int n, m;
	int *num;
	char *op;

	Parsed_CIGAR() {
		n = 0;
		m = 512;
		num = (int*) malloc(m * sizeof(int));
		op = (char*) malloc(m * sizeof(char));
	}

	void parse(const char *cigar);

	void destroy() const { free(num); free(op); }
};

#endif //KSNP_SNP_DBG_H
