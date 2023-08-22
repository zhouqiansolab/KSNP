//
// Created by ixiaohu on 2022/1/28.
//

#ifndef KSNP_SNP_DBG_H
#define KSNP_SNP_DBG_H

#include <vector>
#include <set>
#include <stdint.h>
#include "ksnp_reader.h"
#include "cut_bubble.h"

struct Phased_Block {
	int ps;
	std::vector<int> snp_idx;
	std::vector<int8_t> phase;

	explicit Phased_Block(int p): ps(p) {}

	size_t size() const { return snp_idx.size(); }

	void push_back(int id, int8_t pa) {
		snp_idx.push_back(id);
		phase.push_back(pa);
	}

	/** Sort by phased set of haplotype block */
	bool operator < (const Phased_Block &b) const {
		return this->ps < b.ps || (this->ps == b.ps && this->size() > b.size());
	}
};


class SNP_DBG {
private:
	const std::vector<SNP> &snp_column; /** SNPs to phase. */
	const std::vector<Read_Allele> &read_row; /** Reads with detected alleles. */
	int32_t K; /** K-mer size. */
	int32_t NODE_SHIFT; /** 2^K; the number of k-mer nodes in same rank on DBG. */
	int32_t NODE_MASK; /** 2^K-1; a k-mer node is encoded by a K-bit word. */
	int32_t EDGE_SHIFT; /** 2^(K+1); for diploid, each node has two out-going edge at most. */
	int32_t EDGE_MASK; /** 2^(K+1)-1 */
	int32_t KMER_N; /** #SNPs - K + 1; the number of horizontal k-mers: (i, i+1, i+K-1) ... */
	int *node_cnt; /** The number of a k-mer node occurrences in all reads. */
	int *edge_cnt; /** The number of reads supporting an edge. */
	int *weight; /** The assigned weight of an edge; The three arrays above are accessed
                   * by an index of shift times and an encoded word. */
	std::vector<int8_t> hint_post; /** Hint for post-processing (Caller is the placeholder).
                                    * 0: Phase of high confidence.
                                    * 1: Phase of low quality.
                                    * 3: SNPs from tangle bubbles. */

	/** Index in array of node(w, b) */
	inline int32_t VID(int32_t s, int32_t w) const {
		return s * NODE_SHIFT + w;
	}

	/** Hash value of edge u(s, w[0:K-1]) -> v(s+1, w[1:K]) */
	inline int32_t EID(int32_t s, int32_t e) const {
		return s * EDGE_SHIFT + e;
	}

	/** An occurrence of a node corresponding to SNPs(s, s+1, ..., s+K-1) and K-Mer(wK, w(k-1), ..., w1).
	 *  w is the K-bit word encoding the phase case, with w(i) denote its i-th bit. */
	inline void add_node(int32_t s, int32_t w) {
		node_cnt[VID(s,  w)]++;
		node_cnt[VID(s, w ^ NODE_MASK)]++;
	}

	/** An occurrence of an edge, which is encoded by (K+1)-bit word.
	 *  Same as adding node, this function is called only when building graph. */
	inline void add_edge(int32_t s, int32_t e) {
		edge_cnt[EID(s, e)]++;
		edge_cnt[EID(s, e ^ EDGE_MASK)]++;
	}

	/** Set value for edge: u(s, ..., s+K-1)[w(K), ..., w(1)] -> v(s+1, ... s+K)[w(K-1), ..., b] */
	inline void set_out_weight(int32_t s, int32_t w, int32_t b, int val) {
		int32_t e = (w << 1 | b);
		weight[EID(s, e)] = val;
		weight[EID(s, e ^ EDGE_MASK)] = val;
	}

	/** Get the weight of outgoing edge b(either 0 or 1) of a node. */
	inline int get_out_weight(int32_t s, int32_t w, int32_t b) {
		return weight[EID(s, w << 1 | b)];
	}

	/** Out-degree of node(s, w). */
	inline int out_degree(int32_t s, int32_t w) {
		return (s == KMER_N-1) ?0 :(get_out_weight(s, w, 0) != 0) + (get_out_weight(s, w, 1) != 0);
	}

	/** Set value for edge: u(s-1, ..., s+K-2)[b, ..., w(2)] <- v(s, ... s+K-1)[w(K-1), ..., w(1)] */
	inline void set_in_weight(int32_t s, int32_t w, int32_t b, int val) {
		int32_t e = (b << K | w);
		weight[EID(s - 1, e)] = val;
		weight[EID(s - 1, e ^ EDGE_MASK)] = val;
	}

	/** Get the weight of ingoing edge b(either 0 or 1) of a node. */
	inline int get_in_weight(int32_t s, int32_t w, int32_t b) {
		return weight[EID(s - 1, b << K | w)];
	}

	/** In-degree of node(s, w). */
	inline int in_degree(int32_t s, int32_t w) {
		return (s == 0) ?0 :(get_in_weight(s, w, 0) != 0) + (get_in_weight(s, w, 1) != 0);
	}

	/** @return K-bit word: [b, w(K), ..., w(2)] */
	inline int prev_node(int32_t w, int32_t b) const {
		return ((b << K | w) >> 1) & NODE_MASK;
	}

	/** @return K-bit word: [w(K-1), ..., w(1), b] */
	inline int next_node(int32_t w, int32_t b) const {
		return (w << 1 | b) & NODE_MASK;
	}

	/** Display the read-variant matrix between the range from l to r */
	void display_matrix(int L, int R);

	/** This function is deprecated. Because different K results in varied advantages in metrics
	 * of precision, recall and N50. It is hard to choose a K to balance them. */
	int automatic_k();

	/** Construct DBG by given K, SNPs and linked alleles (reads). */
	void construct_graph();

	/** Trim poor edges which are strongly likely to be wrong in graph. */
	void trim_graph();

	/** Remove short tips away from long and accurate haplotype path. */
	void remove_small_tips();

	/** Bubble path definition. */
	struct Bubble_Path {
		int32_t mid; // one path in bubble area
		int32_t elen; // forward linear path length
		int64_t forw; // bitwise encoded linear path
		Bubble_Path(int32_t m, int32_t e, int64_t f):mid(m), elen(e), forw(f) {}
	};

	/** Collect reads covering bubble area [s, s+blen) */
	std::vector<int> read_on_bubble(int l, int32_t s, int blen, const std::vector<Bubble_Path> &paths);

	/** Display reads on bubble area */
	void display_reads(int l, int64_t back, int32_t s, int blen,
	                   const std::vector<Bubble_Path> &paths, const std::vector<int> &read_set);

	/** Calculate MEC score for a path in bubble area. */
	int mec_bubble_path(const std::vector<int> &read_set, int base, const std::vector<int8_t> &path);

	/** Break small bubble by optimizing MEC score. */
	void topology_bubble();

	struct Bubble {
		int path_n; /** Complexity: Path number in the bubble */
		std::set<int> read_set; /** Reads in bubble */
		int l_bound, r_bound; /** Boundary of reads */
		int connect_n; /** Connected edge number */
		int best_mec; /** Minimum MEC score of the best path */

		Bubble() {
			path_n = 0;
			l_bound = INT32_MAX;
			r_bound = INT32_MIN;
			connect_n = 0;
			best_mec = INT32_MAX;
		}
	};

	/** Within a complexity bound, a bubble can be processed by traceback. */
	void solve_small_bubble(Bubble &bubble);

	/** Solve tangle area by dynamic programming. */
	void dp_tangle_area(const Bubble &bubble);

	/** Select a path with the maximum weight in tangle area. */
	void greedy_tangle_area(const Bubble &bubble);

	/** CUT for tangle bubble */
	CUT_Bubble max_cut;
	void cut_tangle_are(const Bubble &bubble, char *hap1);

	void solve_bubbles();

	/** Keep sure DBG is symmetrical so that meets requirement of processing diploid. */
	void check_graph_symmetry();

	const char *og_prefix = nullptr; /** Prefix of filename which graphs are outputted into. */
	void output_graph(const char *suffix); /** Output graph at each step. */

public:
	SNP_DBG (int k, const std::vector<SNP> &snp, const std::vector<Read_Allele> &read);

	/** Resolve DBG. */
	void resolve();

	/** Find out ambiguous phased SNPs.
	 * @return b(i) represents if the phase between SNP i and i-1 is undetermined. */
	std::vector<int8_t> get_hint_post() { return hint_post; }

	std::vector<Phased_Block> haplotype_block();

	void destroy();

	void output_graph_to(const char *prefix) { og_prefix = prefix; }
};

#endif //KSNP_SNP_DBG_H
