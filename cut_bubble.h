//
// Created by ixiaohu on 2023/2/28.
//

#ifndef KSNP_DEV_CUT_BUBBLE_H
#define KSNP_DEV_CUT_BUBBLE_H

#include <string>
#include <vector>
#include "ksnp_reader.h"

struct Fragment {
	std::vector<char> allele_call; /** Alleles calls, 0: reference allele, 1: alternative allele */
	std::vector<int> variant_idx; /** Corresponding variant indexes in VCF file */
	std::vector<int> phred_quality; /** Phred-scaled quality score of each allele call */
	std::vector<double> probability; /** p = 10 ^ (phred / -10) */
	std::vector<double> prob_miscall; /** log10(p), log is low, precalculating for faster */
	std::vector<double> prob_correct; /** log10(1-p) */

	double prob_RH0, prob_RH1; /** P(R | H), P(R | complement(H)) */
	double prob_RC0, prob_RC1; /** P(R | H_cut), P(R | complement(H_cut)) */

	int size() const { return (int)allele_call.size(); }
};

/** Edge from variant(node) u to v */
struct Edge {
	int v; /** Node v */
	int frag_idx; /** Index to fragment providing this edge */
	char allele_u, allele_v; /** Allele u and allele v */
	double prob_u, prob_v; /** Quality of two allele calls */
	Edge(int nv, int f, char au, char av, double qu, double qv):
			v(nv), frag_idx(f), allele_u(au), allele_v(av), prob_u(qu), prob_v(qv) {}
	bool operator < (const Edge &e) const {
		return this->v < e.v; /** Cluster edges pointing to a same node together */
	}
};

/** Edges used for calculate MAX-CUT */
struct Cut_Edge {
	int u, v; /** Node u -> v */
	double w; /** Weight for a given haplotype */
	Cut_Edge(int u, int v, double w): u(u), v(v), w(w) {}
};

/** A variant also a node in graph */
struct Variant {
	std::vector<int> read_idx; /** Indexes to fragments covering the variant */
	std::vector<int> allele_idx; /** Indexes to alleles on fragments */
	std::vector<Edge> edges; /** Edges connected to the variant */
	std::vector<Cut_Edge> cut_edges; /** Deduplicated edges which sum up all edge weight for given haplotype */
	double score; /** Node score for a given haplotype */
	int parent; /** Which side of cut the node belongs to, p=S when score>0 and T when score<0 */
	int heap_location; /** The node position in heap has to be record as modifications are required */
	bool flipped; /** Whether this node belongs to S side of a cut (it should be flipped) */

	Variant():score(0) {}

	int size() const { return (int)read_idx.size(); }
};

/** Build max heap for variants/nodes */
class Max_Heap {
private:
	std::vector<int> e;
	std::vector<Variant> &v;
	void heapify(int f);

public:
	Max_Heap(std::vector<Variant> &v, std::vector<int> &e);

	/** Return top element */
	int top();
	/** Remove top element */
	void pop();
	/** set the ith element in heap to val */
	void update(int i, double val);
};

class CUT_Bubble {
private:
	std::vector<Fragment> fragments;
	std::vector<Variant> variants;

	/** Log sum approximation */
	const int LOG_SUM10_SIZE = 16000;
	const double LOG_SUM10_PCS = 1000.0;
	double log_sum10_lookup[16000];
	double add_logs(double a, double b);

	/** Build read-variant graph */
	void build_read_variant_graph();

	const int MIN_ALLELE_QUAL = 6; /** Minimum allele quality used to construct graph */
	const int CUT_CONVERGE = 5; /** Stop CUT if exceed this many iterations since improvement. */
	const int MAX_CUT_ROUND = 100; /** Maximum CUT rounds */

	/**
	 * Calculate the score of fragment driven from hap, i.e., likelihoods P(read|haplotype)
	 * @return log10(P(read | haplotype))
	 */
	double calc_fragment_score(const Fragment &f, const char *hap);

	/**
	 * The likelihood-based cut calculation is performed on fragments instead of edges.
	 *
	 * For each fragment calculate four values P(R|H), P(R|complement(H)), P(R|H_cut), P(R|complement(H_cut)),
	 * where likelihood is only for nodes added to the cut till now.
	 * H_cut is new haplotype formed by flipping the phase of vertices in shore T relative to shore S.
	 * Score of a node is log[P(R|H) + P(R|complement(H)] - log[P(R|H_cut) + P(R|complement(H_cut))].
	 */
	void init_node_score(const char *hap, int S, int T);

	/**
	 * Update fragment scores and variant scores as a result of adding a new node to the growing cut
	 */
	void update_node_score(const char *hap, int S, int T, int a, Max_Heap &max_heap);

	/** Optimize MAX-likelihood-CUT */
	double max_likelihood_cut(const char *hap);

public:
	CUT_Bubble();

	/** Construct a variant X fragment matrix for [l, r) interval */
	void haplotyping_bubble(const std::vector<SNP> &snp_column, const std::vector<Read_Allele> &read_row, int L, int R, char *hap1);
};

#endif //KSNP_DEV_CUT_BUBBLE_H

