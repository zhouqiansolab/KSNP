//
// Created by ixiaohu on 2022/1/28.
//

#include <cassert>
#include <cmath>
#include <algorithm>
#include <set>
#include <cstring>
#include <queue>

#include "snp_dbg.h"
#include "time_stamp.h"

int SNP_DBG::automatic_k() {
	int ret = 2; double max_score = 0;
	for (int k = 2; k <= 5; k++) {
		EDGE_SHIFT = (1<<(k+1)); EDGE_MASK = EDGE_SHIFT - 1;
		edge_cnt = (int*) calloc(EDGE_SHIFT * snp_column.size(), sizeof(int));
		// K-mer utilization ratio is improtant
		int kmer_n = 0, effective_n = 0;
		for (const auto &read_variant : read_row) {
			int gap_n = 0; uint32_t temp = 0;
			for (int i = 0; i < k; i++) {
				gap_n += (read_variant[i].allele == -1);
				temp = temp << 1U | (read_variant[i].allele == 1);
			}
			kmer_n += (read_variant.size() >= k);
			effective_n += (gap_n == 0);
			for (int i = k; i < read_variant.size(); i++) {
				gap_n -= (read_variant[i-k].allele == -1);
				gap_n += (read_variant[i].allele == -1);
				temp = (temp << 1U | (read_variant[i].allele == 1)) & (uint32_t)EDGE_MASK;
				kmer_n++;
				effective_n += (gap_n == 0);
				if (gap_n == 0) edge_cnt[EID(read_variant[i-k].snp_idx, temp)]++;
			}
		}
		// Smaller complexity and higher connectivity, the graph is better.
		int graph_n = EDGE_SHIFT * ((int)snp_column.size() - k), edge_n = 0;
		int connect_n = 0;
		for (int s = 0; s + k < snp_column.size(); s++) {
			bool v = false;
			for (int e = 0; e < EDGE_SHIFT; e++) {
				edge_n += (edge_cnt[EID(s, e)] > 0);
				v |= (edge_cnt[EID(s, e)] > 0);
			}
			connect_n += v;
		}
		free(edge_cnt);

		double eff = 1.0 * effective_n / kmer_n;
		double com = 1.0 * edge_n / graph_n;
		double con = 1.0 * (snp_column.size() - k - connect_n) / (snp_column.size() - k);
		double sco = 100 * (eff - com) / (5 * con);
		fprintf(stderr, "    When K=%d, effective k-mer: %.3f, graph complexity: %.3f, ",
		        k, eff, com);
		fprintf(stderr, "graph disconnection: %.3f, estimated score: %.0f\n",
		        con, sco);
		if (sco > max_score) {
			max_score = sco;
			ret = k;
		}
	}
	fprintf(stderr, "    KSNP automatically chooses K=%d to construct DBG\n", ret);
	return ret;
}

SNP_DBG::SNP_DBG (int k, const std::vector<SNP> &snp, const std::vector<Read_Allele> &read):
	K(k), snp_column(snp), read_row(read) {
	fprintf(stderr, "Resolve haplotypes on DBG\n");
	int n = snp.size();
	KMER_N = n - K + 1;
	NODE_SHIFT = (1<<K); NODE_MASK = NODE_SHIFT - 1;
	EDGE_SHIFT = (1<<(K+1)); EDGE_MASK = EDGE_SHIFT - 1;
	node_cnt = (int*) calloc(NODE_SHIFT * n, sizeof(int));
	edge_cnt = (int*) calloc(EDGE_SHIFT * n, sizeof(int));
	weight   = (int*) calloc(EDGE_SHIFT * n, sizeof(int));
	for (int i = 0; i < n; i++) hint_post.push_back(0);
}

void SNP_DBG::display_matrix(int L, int R) {
	std::set<int> read_set;
	for (int i = L; i <= R; i++) {
		const auto &s = snp_column[i];
		for (int j = 0; j < s.size(); j++) {
			read_set.insert(s.rid[j]);
		}
	}
	for (auto rid : read_set) {
		std::string o;
		for (int i = L; i < read_row[rid].front().snp_idx; i++)
			o += '-';
		for (const auto &v: read_row[rid]) {
			if (L <= v.snp_idx and v.snp_idx <= R) {
				if (v.allele == -1) o += '-';
				else if (v.allele == 0) o += '0';
				else o += '1';
			}
		}
		for (int i = read_row[rid].back().snp_idx + 1; i <= R; i++)
			o += '-';
		int no_gap = 0;
		for (auto c : o) no_gap += (c != '-');
		if (no_gap > 1) fprintf(stderr, "%s\n", o.c_str());
	}
}

void SNP_DBG::construct_graph() {
	// Adding k-mers of reads to graph (always consider nodes of k=2)
	int32_t k2_shift = 1U << 2U, k2_mask = k2_shift - 1;
	auto *k2_node = (int32_t*) calloc(k2_shift * snp_column.size(), sizeof(int32_t));
	for (auto &read : read_row) {
		int32_t bit_count = 0, node_word = 0, edge_word = 0;
		for (auto &v : read) {
			int32_t s = v.snp_idx; int8_t allele = v.allele;
			if (allele == -1) {
				// Gap breaks k-mer
				bit_count = node_word = edge_word = 0;
			} else {
				bit_count++;
				node_word = ((node_word << 1) | allele) & NODE_MASK;
				edge_word = ((edge_word << 1) | allele) & EDGE_MASK;
				if (bit_count >= K) {
					// (s-K+1, s-K+2, ..., s)
					add_node(s - K + 1, node_word);
				}
				if (bit_count >= K + 1) {
					// u = (s-K, s-K+1, ..., s-1)
					// v = (s-K+1, s-K+2, ..., s)
					add_edge(s - K, edge_word);
				}
				if (bit_count >= 2) {
					int32_t k2_node_word = node_word & k2_mask;
					k2_node[(s - 2 + 1) * k2_shift + k2_node_word]++;
					k2_node[(s - 2 + 1) * k2_shift + (k2_node_word ^ k2_mask)]++;
				}
			}
		}
	}

	// Patch unconnected part to extend the graph connectivity
	for (int32_t s = 0; s < KMER_N - 1; s++) {
		bool connected = false;
		for (int32_t e = 0; e < EDGE_SHIFT; e++) {
			if (edge_cnt[EID(s, e)] > 0) {
				connected = true;
				break;
			}
		}
		if (connected) continue;
		int max_weight = -1, max_e = -1;
		for (int32_t e = 0; e < EDGE_SHIFT; e++) {
			int wu = node_cnt[VID(s, (e >> 1) & NODE_MASK)];
			int wv = node_cnt[VID(s + 1, e & NODE_MASK)];
			if (wu == 0 || wv == 0) continue;
			int we = int(0.5 * (wu + wv) + 0.499);
			if (we > max_weight) {
				max_weight = we;
				max_e = e;
			}
		}
		if (max_weight == -1) continue;
		// Set an determined edge so that does not introduce errors for bubbles.
		edge_cnt[EID(s, max_e)] = max_weight;
		edge_cnt[EID(s, max_e ^ EDGE_MASK)] = max_weight;
	}

	// Initialize weight of edges
	for (int32_t s = 0; s < KMER_N - 1; s++) {
		// u = (s, s+1, ..., s+K-1)
		// v =    (s+1, ..., s+K-1, s+K)
		for (int32_t e = 0; e < EDGE_SHIFT; e++) {
			weight[EID(s, e)] = edge_cnt[EID(s, e)];
		}
	}

	int nodes_n = 0, edges_n = 0;
	for (int32_t s = 0; s < KMER_N; s++) {
		for (int32_t w = 0; w < NODE_SHIFT; w++) {
			if (node_cnt[VID(s, w)] > 0) {
				nodes_n++;
				edges_n += out_degree(s, w);
			}
		}
	}
	fprintf(stderr, "    Construct DBG(K=%d) of %d nodes and %d edges\n", K, nodes_n, edges_n);

	// Evaluate quality of SNP phase
	const int DEPTH_THRESHOLD = 5; // A node with occurrence < 5 is poor quality
	const double MAX_INCLUDE = 0.05; // Maximum SNPs could be thought as poor quality
	const int MAX_DIFF = 20; // Difference is the number of |Phase(00)-Phase(01)|
	int cnt_ld = 0;
	int cnt_diff[MAX_DIFF + 5]; memset(cnt_diff, 0, sizeof(cnt_diff));
	for (int s = 0; s < snp_column.size() - 1; s++) {
		int total = k2_node[s * k2_shift + 0] + k2_node[s * k2_shift + 1];
		int diff = std::abs(k2_node[s * k2_shift + 0] - k2_node[s * k2_shift + 1]);
		if (total < DEPTH_THRESHOLD) { hint_post[s] = hint_post[s + 1] = 1; cnt_ld++; }
		if (diff < MAX_DIFF) cnt_diff[diff]++;
	}
	for (int i = 1; i < MAX_DIFF; i++) cnt_diff[i] += cnt_diff[i-1];
	int DIFF_THRESHOLD = 5;
	for (int d = 6; d < MAX_DIFF; d++) {
		if (cnt_ld + cnt_diff[d] < MAX_INCLUDE * snp_column.size()) {
			DIFF_THRESHOLD = d;
		}
	}
	for (int s = 0; s < snp_column.size() - 1; s++) {
		int diff = std::abs(k2_node[s * k2_shift + 0] - k2_node[s * k2_shift + 1]);
		if (diff < DIFF_THRESHOLD) { hint_post[s] = hint_post[s + 1] = true; }
	}
	free(k2_node);

	check_graph_symmetry(); // The initial DBG must have two symmetrical subgraph
	output_graph("1_init");
}

void SNP_DBG::trim_graph() {
	// Remove poor-quality competing edges compared with the most high-confident edge
	int trim_n = 0;
	const int CONFIDENT = 15, ABOVE_FILTER = 2, BELOW_FILTER = 5;
	for (int32_t s = 0; s < KMER_N - 1; s++) {
		int32_t max_weight = 0;
		for (int32_t w = 0; w < EDGE_SHIFT; w++) {
			max_weight = std::max(max_weight, weight[s * EDGE_SHIFT + w]);
		}
		if (max_weight >= CONFIDENT) {
			for (int32_t w = 0; w < EDGE_SHIFT; w++) {
				if (weight[s * EDGE_SHIFT + w] * ABOVE_FILTER < max_weight) {
					trim_n++; weight[s * EDGE_SHIFT + w] = 0;
				}
			}
		} else {
			for (int32_t w = 0; w < EDGE_SHIFT; w++) {
				if (weight[s * EDGE_SHIFT + w] * BELOW_FILTER < max_weight) {
					trim_n++; weight[s * EDGE_SHIFT + w] = 0;
				}
			}
		}
	}
	for (int32_t s = 0; s < KMER_N; s++) {
		for (int32_t w = 0; w < NODE_SHIFT; w++) {
			if (in_degree(s, w) == 0 and out_degree(s, w) == 0) {
				node_cnt[VID(s, w)] = 0;
			}
		}
	}
	fprintf(stderr, "    Trim out %d poor-quality edges from the graph\n", trim_n);
	output_graph("2_trim");
}

void SNP_DBG::remove_small_tips() {
	// Backward pass
	const int SHORT_LENGTH = 5, LONG_LENGTH = 15;
	const double WEIGHT_COEF = 1.5;
	int backward_remove = 0;
	for (int32_t s = KMER_N-1; s >= 0; s--) {
		for (int32_t w = 0; w < NODE_SHIFT/2; w++) {
			if (out_degree(s, w) != 2) continue;
			int32_t ns, nw;
			// Whether path leaving from edge 1 is a simple tip: nodes on this path
			// has one incoming or outgoing edge at most.
			ns = s + 1; nw = next_node(w, 1);
			bool p1_is_tip = true; int len1 = 0;
			while (ns < KMER_N - 1) {
				if (in_degree(ns, nw) == 2) p1_is_tip = false; // has two incoming edges
				int w0 = get_out_weight(ns, nw, 0);
				int w1 = get_out_weight(ns, nw, 1);
				if (w0 == 0 and w1 == 0) break; // reach end of the tip
				if (w0 > 0 and w1 > 0) p1_is_tip = false; // has two outgoing edges
				ns++; nw = next_node(nw, w0 > w1 ?0 :1); len1++;
				if (len1 > LONG_LENGTH) break;
			}

			// Whether path leaving from edge 0 is a simple tip
			ns = s + 1; nw = next_node(w, 0);
			bool p0_is_tip = true; int len0 = 0;
			while (ns < KMER_N - 1) {
				if (in_degree(ns, nw) == 2) p0_is_tip = false;
				int w0 = get_out_weight(ns, nw, 0);
				int w1 = get_out_weight(ns, nw, 1);
				if (w0 == 0 and w1 == 0) break;
				if (w0 > 0 and w1 > 0) p0_is_tip = false;
				ns++; nw = next_node(nw, w0 > w1 ?0 :1); len0++;
				if (len0 > LONG_LENGTH) break;
			}

			if (!p0_is_tip and !p1_is_tip) continue; // No tips found

			int tip_to_remove = -1;
			if (p0_is_tip and len0 < SHORT_LENGTH and len1 > LONG_LENGTH) tip_to_remove = 0;
			if (p1_is_tip and len1 < SHORT_LENGTH and len0 > LONG_LENGTH) tip_to_remove = 1;

			if (tip_to_remove == -1) continue;
			int min_len = std::min(len0, len1);
			int total_weight0 = get_out_weight(s, w, 0);
			ns = s + 1; nw = next_node(w, 0);
			for (int i = 0; i < min_len; i++) {
				int w0 = get_out_weight(ns, nw, 0);
				int w1 = get_out_weight(ns, nw, 1);
				total_weight0 += w0 > 0 ?w0 :w1;
				ns++; nw = next_node(nw, w0 > w1 ?0 :1);
			}
			int total_weight1 = get_out_weight(s, w, 1);
			ns = s + 1; nw = next_node(w, 1);
			for (int i = 0; i < min_len; i++) {
				int w0 = get_out_weight(ns, nw, 0);
				int w1 = get_out_weight(ns, nw, 1);
				total_weight1 += w0 > 0 ?w0 :w1;
				ns++; nw = next_node(nw, w0 > w1 ?0 :1);
			}
			bool low_quality_tip = false;
			low_quality_tip |= (p0_is_tip and total_weight1 < WEIGHT_COEF * total_weight0);
			low_quality_tip |= (p1_is_tip and total_weight0 < WEIGHT_COEF * total_weight1);

			backward_remove++;
			if(low_quality_tip) hint_post[s] = 2; // FIXME: better way can solve mistake of tips
			set_out_weight(s, w, tip_to_remove, 0);
			ns = s + 1; nw = next_node(w, tip_to_remove);
			while (ns < KMER_N - 1) {
				if(low_quality_tip) hint_post[ns] = 2;
				int w0 = get_out_weight(ns, nw, 0);
				int w1 = get_out_weight(ns, nw, 1);
				if (w0 == 0 and w1 == 0) break;
				if (w0 > 0) { set_out_weight(ns, nw, 0, 0); ns++; nw = next_node(nw, 0); }
				if (w1 > 0) { set_out_weight(ns, nw, 1, 0); ns++; nw = next_node(nw, 1); }
			}
		}
	}

	// Forward pass
	int forward_remove = 0;
	for (int32_t s = 1; s < KMER_N; s++) {
		for (int32_t w = 0; w < NODE_SHIFT/2; w++) {
			if (in_degree(s, w) != 2) continue;
			int32_t ns, nw;
			ns = s - 1; nw = prev_node(w, 1);
			bool p1_is_tip = true; int len1 = 0;
			while (ns > 0) {
				if (out_degree(ns, nw) == 2) p1_is_tip = false;
				int w0 = get_in_weight(ns, nw, 0);
				int w1 = get_in_weight(ns, nw, 1);
				if (w0 == 0 and w1 == 0) break;
				if (w0 > 0 and w1 > 0) p1_is_tip = false;
				ns--; nw = prev_node(nw, w0 > w1 ?0 :1); len1++;
				if (len1 > LONG_LENGTH) break;
			}

			ns = s - 1; nw = prev_node(w, 0);
			bool p0_is_tip = true; int len0 = 0;
			while (ns > 0) {
				if (out_degree(ns, nw) == 2) p0_is_tip = false;
				int w0 = get_in_weight(ns, nw, 0);
				int w1 = get_in_weight(ns, nw, 1);
				if (w0 == 0 and w1 == 0) break;
				if (w0 > 0 and w1 > 0) p0_is_tip = false;
				ns--; nw = prev_node(nw, w0 > w1 ?0 :1); len0++;
				if (len0 > LONG_LENGTH) break;
			}

			if (!p0_is_tip and !p1_is_tip) continue;

			int tip_to_remove = -1;
			if (p0_is_tip and len0 < SHORT_LENGTH and len1 > LONG_LENGTH) tip_to_remove = 0;
			if (p1_is_tip and len1 < SHORT_LENGTH and len0 > LONG_LENGTH) tip_to_remove = 1;

			if (tip_to_remove == -1) continue;
			int min_len = std::min(len0, len1);
			int total_weight0 = get_in_weight(s, w, 0);
			ns = s - 1; nw = prev_node(w, 0);
			for (int i = 0; i < min_len; i++) {
				int w0 = get_in_weight(ns, nw, 0);
				int w1 = get_in_weight(ns, nw, 1);
				total_weight0 += w0 > 0 ?w0 :w1;
				ns--; nw = prev_node(nw, w0 > w1 ?0 :1);
			}
			int total_weight1 = get_in_weight(s, w, 1);
			ns = s - 1; nw = prev_node(w, 1);
			for (int i = 0; i < min_len; i++) {
				int w0 = get_in_weight(ns, nw, 0);
				int w1 = get_in_weight(ns, nw, 1);
				total_weight1 += w0 > 0 ?w0 :w1;
				ns--; nw = prev_node(nw, w0 > w1 ?0 :1);
			}
			bool low_quality_tip = false;
			low_quality_tip |= (p0_is_tip and total_weight1 < WEIGHT_COEF * total_weight0);
			low_quality_tip |= (p1_is_tip and total_weight0 < WEIGHT_COEF * total_weight1);

			forward_remove++;
			if(low_quality_tip) hint_post[s] = 2;
			set_in_weight(s, w, tip_to_remove, 0);
			ns = s - 1; nw = prev_node(w, tip_to_remove);
			while (ns > 0) {
				if (low_quality_tip) hint_post[ns] = 2;
				int w0 = get_in_weight(ns, nw, 0);
				int w1 = get_in_weight(ns, nw, 1);
				if (w0 == 0 and w1 == 0) break;
				if (w0 > 0) { set_in_weight(ns, nw, 0, 0); ns--; nw = prev_node(nw, 0); }
				if (w1 > 0) { set_in_weight(ns, nw, 1, 0); ns--; nw = prev_node(nw, 1); }
			}
		}
	}
	for (int32_t s = 0; s < KMER_N; s++) {
		for (int32_t w = 0; w < NODE_SHIFT; w++) {
			if (in_degree(s, w) == 0 and out_degree(s, w) == 0) {
				node_cnt[VID(s, w)] = 0;
			}
		}
	}
	fprintf(stderr, "    Remove %d backward and %d forward short tips\n", backward_remove, forward_remove);
	output_graph("3_tip");
}

std::vector<int> SNP_DBG::read_on_bubble(int l, int32_t s, int blen, const std::vector<Bubble_Path> &paths) {
	uint8_t diff[blen]; memset(diff, 0, sizeof(diff));
	for (const auto &b: paths) {
		int8_t prev = b.mid >> (blen-1) & 1;
		for (int i = 1; i < blen; i++) {
			int8_t curr = b.mid >> (blen-1-i) & 1;
			if ((prev ^ curr) == 1) {
				diff[i] |= 1U; diff[i-1] |= 1U;
			} else {
				diff[i] |= 2U; diff[i-1] |= 2U;
			}
			prev = curr;
		}
	}
	std::set<int> read_v;
	for (int i = 0; i < blen; i++) {
		if (diff[i] == 3) {
			const auto &snp = snp_column[s + i];
			for (int j = 0; j < snp.rid.size(); j++) {
				if (snp.allele[j] != -1) {
					read_v.insert(snp.rid[j]);
				}
			}
		}
	}
	// Calculate boundary
	int L = s - l, R = s + blen + 128;
	for (const auto path: paths) R = std::min(R, s + blen + path.elen);
	std::vector<int> read_set;
	for (auto rid : read_v) {
		const auto &read = read_row[rid];
		// The bubble overlaps with another one.
		if (read.front().snp_idx < L or read.back().snp_idx >= R) {
			read_set.clear();
			return read_set;
		}
		read_set.push_back(rid);
	}
	return read_set;
}


int SNP_DBG::mec_bubble_path(const std::vector<int> &read_set, int base, const std::vector<int8_t> &path) {
	int ret = 0;
	for (auto rid : read_set) {
		int h0 = 0, h1 = 0; bool beyond = false;
		for (const auto &x : read_row[rid]) {
			if (x.allele == -1) continue;
			int i = x.snp_idx - base;
			if (i < 0 or i >= path.size()) {
				// This allele is beyond the boundary
				beyond = true;
				break;
			}
			if (path[i] == x.allele) h0++;
			else h1++;
		}
		if (beyond) return -1;
		if (!beyond) ret += std::min(h0, h1);
	}
	return ret;
}

void SNP_DBG::display_reads(int l, int64_t back, int32_t s, int blen, const std::vector<Bubble_Path> &paths,
                            const std::vector<int> &read_set) {
	fprintf(stderr, "Key_SNP:%d, Backward:%d, Bubble:%d, Path:%ld, Reads:%ld\n",
		s, l, blen, paths.size(), read_set.size());
	for (const auto &p: paths) {
		std::vector<int8_t> decode;
		fprintf(stderr, " ");
		for (int i = l-1; i >= 0; i--) {
			decode.push_back(back >> i & 1);
			fprintf(stderr, "%d", back >> i & 1);
		}
		fprintf(stderr, " ");
		for (int i = blen-1; i >= 0; i--) {
			decode.push_back(p.mid >> i & 1);
			fprintf(stderr, "%d", p.mid >> i & 1);
		}
		fprintf(stderr, " ");
		for (int i = 0; i < p.elen; i++) {
			decode.push_back(p.forw >> i & 1);
			fprintf(stderr, "%d", p.forw >> i & 1);
		}
		fprintf(stderr, "    MEC=%d\n", mec_bubble_path(read_set, s-l, decode));
	}
	for (const auto &rid : read_set) {
		const auto &read = read_row[rid];
		bool beyond = (read.front().snp_idx < s - l);
		fprintf(stderr, "%c", beyond ?'-' :' ');
		int output_n = std::max(read.front().snp_idx - (s-l), 0);
		if (output_n == l or output_n == l + blen) fprintf(stderr, " ");
		for (int i = s - l; i < read.front().snp_idx; i++) fprintf(stderr, " ");
		for (const auto &x : read) {
			if (x.snp_idx < s - l) continue;
			if (x.allele == 0) fprintf(stderr, "0");
			else if (x.allele == 1) fprintf(stderr, "1");
			else fprintf(stderr, "-");
			output_n++;
			if (output_n == l or output_n == l + blen) fprintf(stderr, " ");
		}

		fprintf(stderr, "\n");
	}
	fprintf(stderr, "\n");
}

void SNP_DBG::topology_bubble() {
	// Calculating maximum backward linear length. This implementation seeks for O(E) time complexity.
	// Dynamic programming in forward enumerating loop: f(v) = f(u) + 1 when ONLY u is connected to v.
	auto *back_linear = new int[KMER_N * NODE_SHIFT]; // Maximum length of backward linear path from node v
	memset(back_linear, 0, KMER_N * NODE_SHIFT * sizeof(int));
	auto *back_path = new int64_t[KMER_N * NODE_SHIFT]; // Bitwise encoded linear path
	memset(back_path, 0, KMER_N * NODE_SHIFT * sizeof(int64_t));
	for (int32_t s = 0; s < KMER_N; s++) {
		for (int32_t w = 0; w < NODE_SHIFT; w++) {
			int v = s * NODE_SHIFT + w;
			if (back_linear[v] != 0) continue;
			int32_t cs = s, cw = w; // Current node u
			while (cs < KMER_N) {
				if (out_degree(cs, cw) != 1) break;
				int w0 = get_out_weight(cs, cw, 0);
				int w1 = get_out_weight(cs, cw, 1);
				int ns = cs + 1, nw = next_node(cw, w0>0 ?0 :1); // Next node v
				if (in_degree(ns, nw) != 1) break; // v has another incoming edge
				back_linear[ns * NODE_SHIFT + nw] = back_linear[cs * NODE_SHIFT + cw] + 1;
				back_path[ns * NODE_SHIFT + nw] = (back_path[cs * NODE_SHIFT + cw] << 1 | ((cw>>(K-1))&1));
				cs = ns; cw = nw;
			}
		}
	}

	// Calculating maximum forward linear length.
	auto *forw_linear = new int[KMER_N * NODE_SHIFT]; // Maximum length of forward linear path
	memset(forw_linear, 0, KMER_N * NODE_SHIFT * sizeof(int));
	auto *forw_path = new int64_t[KMER_N * NODE_SHIFT];
	memset(forw_path, 0, KMER_N * NODE_SHIFT * sizeof(int64_t));
	for (int32_t s = KMER_N-1; s >= 0; s--) {
		for (int32_t w = 0; w < NODE_SHIFT; w++) {
			int v = s * NODE_SHIFT + w;
			if (forw_linear[v] != 0) continue;
			int32_t cs = s, cw = w;
			while (cs >= 0) {
				if (in_degree(cs, cw) != 1) break;
				int w0 = get_in_weight(cs, cw, 0);
				int w1 = get_in_weight(cs, cw, 1);
				int ns = cs - 1, nw = prev_node(cw, w0>0 ?0 :1);
				if (out_degree(ns, nw) != 1) break;
				forw_linear[ns * NODE_SHIFT + nw] = forw_linear[cs * NODE_SHIFT + cw] + 1;
				forw_path[ns * NODE_SHIFT + nw] = (forw_path[cs * NODE_SHIFT + cw] << 1 | (cw&1));
				cs = ns; cw = nw;
			}
		}
	}

	// Resolve bubbles
	const int MIN_LINEAR_LEN = 10; // Minimum linear path at both sides of a bubble,
								  // or minimum distance between two neighboring bubbles.
	const int MAX_LINEAR_LEN = 64; // Maximum linear path to extend, to avoid massive unrelated
	                               // nodes included. So, a path can be fit into a 64-bit number.
    const int MAX_BUBBLE_PATH = 256; // Maximum paths in bubble area; constraint the bubble complexity.
    for (int32_t s = 0; s < KMER_N; s++) {
    	for (int32_t w = 0; w < NODE_SHIFT; w++) {
    		int32_t u = s * NODE_SHIFT + w;
    		forw_linear[u] = std::min(forw_linear[u], MAX_LINEAR_LEN);
		    back_linear[u] = std::min(back_linear[u], MAX_LINEAR_LEN);
    	}
    }
	const int MAX_BUBBLE_LEN = 30; // Maximum bubble size; only consider small bubbles here.
	int resolve_n = 0;
	auto *vis = new bool[KMER_N * NODE_SHIFT]; // Whether a node is visited
	memset(vis, false, sizeof(bool) * KMER_N * NODE_SHIFT);
	int max_bubble_n = 0;
	for (int32_t s = 0; s < KMER_N; s++) {
		for (int32_t w = 0; w < NODE_SHIFT/2; w++) {
			if (vis[VID(s, w)] or out_degree(s, w) != 2 or back_linear[VID(s, w)] < MIN_LINEAR_LEN)
				continue;

			/* Go through the bubble area */
			int bubble_len = 0; // length of the bubble area (=SNP numbers in this area)
			std::vector<int32_t> prev, curr; // bitwise encoded bubble paths, for example,
			// bubble {x->y0, y0->z, x->y1, y1->z} has two paths x,y0,z and x,y1,z.
			std::vector<Bubble_Path> not_bubble; // Tips or incomplete bubble
			std::vector<Bubble_Path> bubble_paths; // Including linear path in two sides beside bubble paths.
			prev.push_back(w); // start node (s, w)
			vis[VID(s, w)] = true;
			for (int l = 1; l <= MAX_BUBBLE_LEN and s + l < KMER_N; l++) {
				curr.clear();
				int8_t merge_n[NODE_SHIFT]; memset(merge_n, 0, sizeof(merge_n));
				for (auto p : prev) {
					int32_t cs = s + l - 1, cw = (p & NODE_MASK);
					int w0 = get_out_weight(cs, cw, 0);
					int w1 = get_out_weight(cs, cw, 1);
					if (w0 > 0) { curr.push_back(p << 1 | 0); }
					if (w1 > 0) { curr.push_back(p << 1 | 1); }
					vis[VID(s+l, next_node(cw, w0>0 ?0 :1))] = true;
					// If one path reach its end (has no outgoing edge), ignore it, though, maybe it is
					// okay to delete it.
					if (w0 == 0 and w1 == 0) not_bubble.emplace_back(Bubble_Path(p, K+l-1, 0));
					else merge_n[next_node(cw, w0>0 ?0 :1)]++;
				}
				if (curr.empty()) break; // Only find branches but bubble
				if (curr.size() > MAX_BUBBLE_PATH) break; // The bubble area is too complex
				std::swap(prev, curr);

				// Require that at least one closed bubble is found
//				bool one_close = false;
//				for (int32_t m = 0; m < NODE_SHIFT; m++) {
//					// Two paths joins a same node
//					one_close |= (merge_n[m] == 2);
//				}
//				if (!one_close) continue;

				// There might be overlapping bubbles. Keep sure the biggest bubble is closed.
				bool all_close = true;
				for (auto p : prev) {
					int32_t cs = s + l, cw = (p & NODE_MASK);
					all_close &= (forw_linear[VID(cs, cw)] >= MIN_LINEAR_LEN);
				}
				if (all_close) { // End of the bubble area.
					bubble_len = K + l; // K is the first node and l is the extending length
					for (auto p : prev) {
						int32_t cs = s + l, cw = (p & NODE_MASK);
						bubble_paths.emplace_back(
							Bubble_Path(p, forw_linear[VID(cs, cw)], forw_path[VID(cs, cw)])
						);
					}
					break;
				}
			}
			if (bubble_len == 0) continue;

			/* Resolve this bubble */
			assert(bubble_len >= K + K); // From branching node, two paths join each other after at least K nodes.
			assert(!bubble_paths.empty());
			max_bubble_n = std::max(max_bubble_n, (int)bubble_paths.size());
			std::vector<int8_t> path; // Decoding bitwise bubble
			for (const auto &b : bubble_paths) {
				path.clear();
				int l = back_linear[VID(s, w)], m = bubble_len, r = b.elen;
				for (int j = l-1; j >= 0; j--) path.push_back(back_path[VID(s, w)] >> j & 1);
				for (int j = bubble_len-1; j >= 0; j--) path.push_back(b.mid >> j & 1);
				for (int j = 0; j < b.elen; j++) path.push_back(b.forw >> j & 1);

				for (int j = 0; j + K < path.size(); j++) {
					int32_t e = 0;
					for (int k = 0; k <= K; k++) e = e << 1 | path[j + k];
					assert(weight[EID(s-l+j, e)] > 0);
				}
			}

			for (const auto &p : not_bubble) { // First remove non-bubble path
				path.clear();
				for (int i = p.elen-1; i >= 0; i--) path.push_back(p.mid >> i & 1);
				for (int i = 0; i + K < path.size(); i++) {
					int32_t e = 0;
					for (int j = 0; j <= K; j++) e = e << 1 | path[i+j];
					weight[EID(s + i, e)] = 0;
					weight[EID(s + i, e ^ EDGE_MASK)] = 0;
				}
			}

			auto read_set = read_on_bubble(back_linear[VID(s, w)], s, bubble_len, bubble_paths);
			int left_boundary = snp_column.size()-1, right_boundary = -1;
			for (auto rid : read_set) {
				left_boundary = std::min(left_boundary, read_row[rid].front().snp_idx);
				right_boundary = std::max(right_boundary, read_row[rid].back().snp_idx);
			}
			if (read_set.empty()) continue; // This bubble overlaps with another one
//			display_reads(back_linear[VID(s, w)], back_path[VID(s, w)], s, bubble_len, bubble_paths, read_set);

			int min_mec = INT32_MAX, min_path = -1;
			for (int i = 0; i < bubble_paths.size(); i++) { // Choose the best path with optimized MEC score
				const auto &b = bubble_paths[i];

				path.clear();
				int l = back_linear[VID(s, w)], m = bubble_len, r = b.elen;
				// Backward linear path
				for (int j = l-1; j >= 0; j--) path.push_back(back_path[VID(s, w)] >> j & 1);
				// Middle bubble path (only consider reads on bubble part)
				for (int j = bubble_len-1; j >= 0; j--) path.push_back(b.mid >> j & 1);
				// Forward linear path (is of equal length and symmetrical)
				for (int j = 0; j < b.elen; j++) path.push_back(b.forw >> j & 1);

				int mec = mec_bubble_path(read_set, s-l, path);
				if (mec < min_mec) {
					min_mec = mec;
					min_path = i;
				}
				// Delete all paths for the moment
				for (int j = 0; j + K < path.size(); j++) {
					int32_t e = 0;
					for (int k = 0; k <= K; k++) e = e << 1 | path[j + k];
					weight[EID(s-l+j, e)] = 0;
					weight[EID(s-l+j, e^EDGE_MASK)] = 0;
				}
			}
			if (min_path != -1) {
				// Only keep the chosen path with minimum MEC score.
				const auto &b = bubble_paths[min_path];

				path.clear();
				int l = back_linear[VID(s, w)], m = bubble_len, r = b.elen;
				for (int j = l-1; j >= 0; j--) path.push_back(back_path[VID(s, w)] >> j & 1);
				for (int j = bubble_len-1; j >= 0; j--) path.push_back(b.mid >> j & 1);
				for (int j = 0; j < b.elen; j++) path.push_back(b.forw >> j & 1);

				for (int j = 0; j + K < path.size(); j++) {
					int32_t e = 0;
					for (int k = 0; k <= K; k++) e = e << 1 | path[j + k];
					weight[EID(s-l+j, e)] = 1;
					weight[EID(s-l+j, e^EDGE_MASK)] = 1;
				}
				resolve_n++;
			}
		}
	}
	delete [] vis;
	fprintf(stderr, "    Resolve %d small bubbles, max-complexity is %d\n", resolve_n, max_bubble_n);

	delete [] back_linear;
	delete [] forw_linear;
}

void SNP_DBG::solve_small_bubble(Bubble &b) {
	struct QueElement {
		int s, w;
		std::vector<int8_t> path;
		QueElement(int ss, int ww): s(ss), w(ww) {}
	};

	std::vector<std::vector<int8_t>> all_paths;
	std::queue<QueElement> que;
	for (int w = 0; w < NODE_SHIFT/2; w++) {
		if (node_cnt[VID(b.l_bound, w)] > 0) {
			QueElement u(b.l_bound, w);
			for (int i = 0; i < K; i++) u.path.push_back(w >> (K-1-i) & 1);
			que.push(u);
		}
	}
	while (!que.empty()) {
		auto curr = que.front(); que.pop();
		if (curr.s + K - 1 == b.r_bound) {
			assert(curr.path.size() == b.r_bound - b.l_bound + 1);
			all_paths.push_back(curr.path);
			continue;
		}
		for (int i = 0; i < 2; i++) {
			if (get_out_weight(curr.s, curr.w, i) > 0) {
				auto next = curr;
				next.s++; next.w = next_node(curr.w, i);
				next.path.push_back(i);
				que.push(next);
			}
		}
	}
	if (all_paths.empty()) {
		b.best_mec = -1;
		return;
	}
	assert(not all_paths.empty());

	int min_mec = INT32_MAX, min_id = -1;
	for (int i = 0; i < all_paths.size(); i++) {
		const auto &path = all_paths[i];
		int mec = 0;
		for (auto rid : b.read_set) {
			int h0 = 0, h1 = 0;
			for (const auto &v : read_row[rid]) {
				if (v.allele == -1) continue;
				if (v.allele != path[v.snp_idx - b.l_bound]) h0++;
				else h1++;
			}
			mec += std::min(h0, h1);
		}
		if (mec < min_mec) {
			min_mec = mec;
			min_id = i;
		}
	}

	for (int i = b.l_bound; i+K <= b.r_bound; i++) {
		for (int e = 0; e < EDGE_SHIFT; e++) {
			weight[EID(i, e)] = 0;
		}
	}
	for (int i = b.l_bound; i+K <= b.r_bound; i++) {
		int e = 0;
		for (int k = 0; k <= K; k++) {
			e = e << 1 | all_paths[min_id][i+k-b.l_bound];
		}
		weight[EID(i, e)] = edge_cnt[EID(i, e)];
		weight[EID(i, e^EDGE_MASK)] = edge_cnt[EID(i, e^EDGE_MASK)];
	}
	b.best_mec = min_mec;
}

void SNP_DBG::greedy_tangle_area(const Bubble &bubble) {
	int n = bubble.r_bound - bubble.l_bound; // Number of K-mer columns
	auto *f = new int[n * NODE_SHIFT * sizeof(int)]; // f(u): maximum weight until node u
	memset(f, 0, n * NODE_SHIFT * sizeof(int));
	auto *t = new int[n * NODE_SHIFT * sizeof(int)]; // t(u): back trace of node u
	memset(t, -1, n * NODE_SHIFT * sizeof(int));
	for (int s = bubble.l_bound + 1; s+K-1 <= bubble.r_bound; s++) {
		for (int w = 0; w < NODE_SHIFT; w++) {
			for (int b = 0; b < 2; b++) {
				int we = get_in_weight(s, w, b);
				if (we == 0) continue;
				// F(u) = max{F(u), F(v) + W(u, v)}
				int pw = prev_node(w, b);
				int temp = f[VID(s-1-bubble.l_bound, pw)] + we;
				if (temp > f[VID(s-bubble.l_bound, w)]) {
					f[VID(s-bubble.l_bound, w)] = temp;
					t[VID(s-bubble.l_bound, w)] = VID(s-1-bubble.l_bound, pw);
				}
			}
		}
	}
	for (int s = bubble.l_bound; s + K <= bubble.r_bound; s++) {
		for (int e = 0; e < EDGE_SHIFT; e++) {
			weight[EID(s, e)] = 0;
		}
	}
	int es = n-1;
	while (es > 0) {
		int ew = 0, max_w = 0;
		for (int w = 0; w < NODE_SHIFT; w++) {
			if (f[VID(es, w)] > max_w) {
				max_w = f[VID(es, w)];
				ew = w;
			}
		}

		int e = VID(es, ew);
		while (t[e] != -1) {
			int us = t[e] / NODE_SHIFT + bubble.l_bound, uw = t[e] % NODE_SHIFT; // node u
			int vs = e / NODE_SHIFT + bubble.l_bound, vw = e % NODE_SHIFT; // node v
			assert(us == vs - 1);
			for (int i = 0; i < K-1; i++) {
				assert((uw >> (K - 1 - (i + 1)) & 1) == (vw >> (K - 1 - i) & 1));
			}
			int x = uw << 1 | (vw & 1); // edge u->v
			weight[EID(us, x)] = edge_cnt[EID(us, x)];
			weight[EID(us, x^EDGE_MASK)] = edge_cnt[EID(us, x^EDGE_MASK)];
			e = t[e];
		}
		es = e / NODE_SHIFT; es -= K; // remove overlap
	}
	delete [] t;
	delete [] f;
}

void SNP_DBG::dp_tangle_area(const Bubble &bubble) {
	// Greedy selecting reads
	struct Read_Interval {
		int rid; int l, r;
		Read_Interval (int r, int ll, int rr): rid(r), l(ll), r(rr) {}
		bool operator < (const Read_Interval &x) const {
			if (this->l != x.l) return this->l < x.l;
			else return this->r > x.r;
		}
	};
	std::vector<Read_Interval> intvs;
	for (auto rid : bubble.read_set) {
		const auto &r = read_row[rid];
		assert(r.front().snp_idx >= bubble.l_bound and r.back().snp_idx <= bubble.r_bound);
		intvs.emplace_back(Read_Interval(rid, r.front().snp_idx, r.back().snp_idx));
	}
	std::sort(intvs.begin(), intvs.end()); // Sort by intervals

	const int MAX_COVERAGE = 15;
	auto *coverage = new int[bubble.r_bound - bubble.l_bound + 1];
	memset(coverage, 0, (bubble.r_bound-bubble.l_bound+1) * sizeof(int));
	auto *used = new bool[intvs.size()];
	std::vector<int> selected;
	while (!intvs.empty()) { // Until no read left or exceeding maximum coverage
		// Each iteration selects the longest slice linked by reads
		bool choose_any = false;
		memset(used, false, intvs.size() * sizeof(bool));
		int right_upto = -1;
		for (int i = 0, next; i < intvs.size(); i = next) {
			next = i + 1;
			const auto &p = intvs[i];
			if (p.r <= right_upto) continue; // Covered by the previous read
			bool exceed = false;
			for (const auto &v : read_row[p.rid]) {
				if (coverage[v.snp_idx-bubble.l_bound] + 1 > MAX_COVERAGE) {
					exceed = true;
					break;
				}
			}
			if (exceed) continue; // Exceeding the coverage limitation

			choose_any = true;
			right_upto = p.r;
			used[i] = true; selected.push_back(p.rid); // Set coverage
			for (const auto &v : read_row[p.rid]) coverage[v.snp_idx-bubble.l_bound]++;
			int max_r = right_upto; int jump = intvs.size();
			for (int j = i + 1; j < intvs.size(); j++) {
				const auto &c = intvs[j];
				if (c.l <= right_upto) { // Overlapped with previous read
					if (c.r > max_r) { // Select one with rightmost interval
						max_r = c.r;
						next = j;
					}
				} else {
					jump = j; // If no overlapping read, jump to it.
					break;
				}
			}
			if (max_r == p.r) next = jump; // No linked read, jump to a new slice.
		}

		if (!choose_any) break; // No read can be chosen

		// Remove used reads
		std::vector<Read_Interval> temp;
		for (int i = 0; i < intvs.size(); i++) {
			if (!used[i]) temp.push_back(intvs[i]);
		}
		intvs = temp; // Next round
	}
	delete [] used;

	int m = bubble.r_bound - bubble.l_bound + 1;
	int max_cov = 0;
	for (int i = bubble.l_bound; i <= bubble.r_bound; i++)
		max_cov = std::max(max_cov, coverage[i - bubble.l_bound]);

	auto *idx_column = new std::vector<int>[m];
	auto *allele_column = new std::vector<int8_t>[m];
	for (int i = 0; i < selected.size(); i++) {
		int rid = selected[i];
		for (const auto &v : read_row[rid]) {
			idx_column[v.snp_idx - bubble.l_bound].push_back(i);
			allele_column[v.snp_idx - bubble.l_bound].push_back(v.allele);
		}
	}
	for (int i = 0; i < m; i++) {
		assert(idx_column[i].size() == coverage[i]);
		assert(allele_column[i].size() == coverage[i]);
	}

	const int SHIFT_COLUMN = 1 << max_cov;
	auto *dp_table = new int[m * SHIFT_COLUMN];
	auto *min_extend = new int[m * SHIFT_COLUMN];
	memset(min_extend, 0x3f, m * SHIFT_COLUMN * sizeof(int));
	auto *back_trace = new int[m * SHIFT_COLUMN];
	fprintf(stderr, "Bubble=[%d,%d]; Select=%ld\n", bubble.l_bound, bubble.r_bound, selected.size());
	for (int i = 0; i < m; i++) {
		int intersection_n = 0;
		if (i > 0) {
			for (int j = 0, k = 0; j < idx_column[i-1].size() and k < idx_column[i].size(); ) {
				if (idx_column[i-1][j] < idx_column[i][k]) j++;
				else if (idx_column[i-1][j] < idx_column[i][k]) k++;
				else {
					intersection_n++;
					break;
				}
			}
		}
		for (int B = 0; B < (1<<coverage[i]); B++) {
			// DP(i, B) = Cost(i,B) + min{DP(i-1, b)|b extends B∩F(i-1)}
			int R0 = 0, R1 = 0, S0 = 0, S1 = 0;
			for (int r = 0; r < coverage[i]; r++) {
				if ((B & (1 << r)) == 0) { // read r belongs to R
					R0 += allele_column[i][r] == 0 ?1 :0;
					R1 += allele_column[i][r] == 1 ?1 :0;
				} else { // read r belongs to S
					S0 += allele_column[i][r] == 0 ?1 :0;
					S1 += allele_column[i][r] == 1 ?1 :0;
				}
			}
			int cost = std::min(R0, R1) + std::min(S0, S1);

			if (intersection_n == 0) {
				// No intersection between the two columns
				dp_table[B] = cost;
			} else {
				// Index projection
				int b = 0;
				for (int j = 0, k = 0; j < idx_column[i-1].size() and k < idx_column[i].size(); ) {
					if ((B & (1 << k)) != 0) { k++; continue; }
					if (idx_column[i-1][j] < idx_column[i][k]) j++;
					else if (idx_column[i-1][j] > idx_column[i][k]) k++;
					else {
						b = b | (1 << j);
						j++; k++;
					}
				}
				dp_table[i * SHIFT_COLUMN + B] = cost + min_extend[(i-1) * SHIFT_COLUMN + b];
			}
			if (i < m-1) {
				// Compute min{DP(i, b) | b extends B∩F(i+1)} for next column
				int b = 0;
				for (int j = 0, k = 0; j < idx_column[i].size() and k < idx_column[i+1].size();) {
					if ((B & (1 << j)) != 0) { j++; continue; }
					if (idx_column[i][j] < idx_column[i+1][k]) j++;
					else if (idx_column[i][j] > idx_column[i+1][k]) k++;
					else {
						b = b | (1 << j);
						j++; k++;
					}
				}
				min_extend[i * SHIFT_COLUMN + b] = std::min(
					min_extend[i * SHIFT_COLUMN + b], dp_table[i * SHIFT_COLUMN + B]
				);
			}
		}
	}
	int min_dp_value = INT32_MAX;
	for (int B = 0; B < (1<<coverage[m-1]); B++) {
		min_dp_value = std::min(min_dp_value, dp_table[(m-1) * SHIFT_COLUMN + B]);
	}
	fprintf(stderr, "MEC = %d\n", min_dp_value);

	delete [] back_trace;
	delete [] min_extend;
	delete [] dp_table;
	delete [] idx_column;
	delete [] allele_column;
	delete [] coverage;
}

void SNP_DBG::cut_tangle_are(const Bubble &bubble, char *hap1) {
	// Before CUT, find a intermedia haplotype in greedy manner
	greedy_tangle_area(bubble);
	int s = bubble.l_bound, w = 0;
	for (w = 0; w < NODE_SHIFT; w++) {
		if (out_degree(s, w) == 1) {
			break;
		}
	}
	for (int i = 0; i < K; i++) {
		if ((w >> (K-1-i)) != 0) hap1[i] = '1';
		else hap1[i] = '0';
	}
	for (int i = bubble.l_bound + K; i <= bubble.r_bound; i++) {
		if (get_out_weight(s, w, 0) > 0) {
			hap1[i - bubble.l_bound] = '0';
			s++; w = next_node(w, 0);
		} else {
			hap1[i - bubble.l_bound] = '1';
			s++; w = next_node(w, 1);
		}
	}

	std::vector<Read_Allele> sub_row;
	for (auto r : bubble.read_set) sub_row.push_back(read_row[r]);
	max_cut.haplotyping_bubble(snp_column, sub_row, bubble.l_bound, bubble.r_bound + 1, hap1);

	for (int s = bubble.l_bound; s + K <= bubble.r_bound; s++) {
		for (int e = 0; e < EDGE_SHIFT; e++) {
			weight[EID(s, e)] = 0;
		}
	}
	uint32_t e = 0;
	for (int k = 0; k < K; k++) e = (e << 1U | (hap1[k] == '1'));
	for (int s = bubble.l_bound + K; s <= bubble.r_bound; s++) {
		e = (e << 1U | (hap1[s - bubble.l_bound] == '1')) & (uint32_t)EDGE_MASK;
		weight[EID(s - K, e)] = edge_cnt[EID(s - K, e)];
		weight[EID(s - K, e ^ EDGE_MASK)] = edge_cnt[EID(s - K, e ^ EDGE_MASK)];
	}
}

void SNP_DBG::solve_bubbles() {
	const uint8_t PHASE0 = 1U, PHASE1 = 2U, AMB_PHASE = 3U;
	auto *phase_record = new uint8_t[KMER_N * NODE_SHIFT];
	memset(phase_record, 0U, KMER_N * NODE_SHIFT * sizeof(uint8_t));
	for (int32_t s = 0; s < KMER_N - 1; s++) {
		for (int32_t e = 0; e < EDGE_SHIFT / 2; e++) {
			if (weight[EID(s, e)] == 0) continue;
			uint8_t prev = e >> K & 1;
			for (int k = 1; k <= K; k++) {
				uint8_t curr = e >> (K - k) & 1;
				phase_record[s + k] |= (prev == curr) ? PHASE0 : PHASE1;
				prev = curr;
			}
		}
	}

	std::vector<Bubble> all_bubbles;
	Bubble temp_bubble; int prev_amb = snp_column.size();
	for (int s = 1; s < snp_column.size(); s++) {
		if (phase_record[s] != AMB_PHASE) continue; // Phase between s-1 and s is undetermined.
		const auto &prev_snp = snp_column[s - 1], curr_snp = snp_column[s];
		int l_bound = INT32_MAX, r_bound = INT32_MIN;
		std::vector<int> read_set;
		for (int i = 0, j = 0; i < prev_snp.size() and j < curr_snp.size();) {
			if (prev_snp.rid[i] < curr_snp.rid[j]) i++;
			else if (prev_snp.rid[i] > curr_snp.rid[j]) j++;
			else {
				int rid = prev_snp.rid[i];
				l_bound = std::min(l_bound, read_row[rid].front().snp_idx);
				r_bound = std::max(r_bound, read_row[rid].back().snp_idx);
				read_set.push_back(rid);
				i++; j++;
			}
		}
		assert(!read_set.empty());
		if (l_bound > prev_amb) {
			all_bubbles.push_back(temp_bubble);
			temp_bubble = Bubble(); // reset the temporary variable for a new bubble area
		}
		for (auto rid : read_set) temp_bubble.read_set.insert(rid);
		temp_bubble.l_bound = std::min(temp_bubble.l_bound, l_bound);
		temp_bubble.r_bound = std::max(temp_bubble.r_bound, r_bound);
		prev_amb = s;
	}
	if (!temp_bubble.read_set.empty()) all_bubbles.push_back(temp_bubble);

	const int MAX_COMPLEXITY = 256; // Maximum number of paths in bubble to resolve
	auto *dp = new int[KMER_N * NODE_SHIFT]; // How many paths reach the end from the start of a bubble

	int overflow_n = 0, incomplete_n = 0, simple_n = 0;
	for (auto &bubble: all_bubbles) {
		assert(bubble.l_bound >= 0 and bubble.l_bound < snp_column.size());
		assert(bubble.r_bound >= 0 and bubble.r_bound < snp_column.size());
		for (int s = bubble.l_bound; s + K - 1 <= bubble.r_bound; s++) {
			for (int w = 0; w < NODE_SHIFT; w++) {
				// Connect unclosed path in brute way which might connect false edges.
				// And when graph is patched, the reads on bubble might change. However, the affect is negligible.
				if (s > bubble.l_bound and node_cnt[VID(s, w)] > 0 and in_degree(s, w) == 0) {
					if (node_cnt[VID(s-1, prev_node(w, 0))] > 0) {
						set_in_weight(s, w, 0, 1);
						bubble.connect_n++;
					}
					if (node_cnt[VID(s-1, prev_node(w, 1))] > 0) {
						set_in_weight(s, w, 1, 1);
						bubble.connect_n++;
					}
				}
				if (s + K - 1 < bubble.r_bound and node_cnt[VID(s, w)] > 0 and out_degree(s, w) == 0) {
					if (node_cnt[VID(s+1, next_node(w, 0))] > 0) {
						set_out_weight(s, w, 0, 1);
						bubble.connect_n++;
					}
					if (node_cnt[VID(s+1, next_node(w, 1))] > 0) {
						set_out_weight(s, w, 1, 1);
						bubble.connect_n++;
					}
				}
			}
		}
		// Graph has been patched, so the backup needs update.
		for (int s = bubble.l_bound; s + K <= bubble.r_bound; s++) {
			for (int e = 0; e < EDGE_SHIFT; e++) {
				edge_cnt[EID(s, e)] = weight[EID(s, e)];
			}
		}
		for (int w = 0; w < NODE_SHIFT; w++)
			dp[VID(bubble.l_bound, w)] = (node_cnt[VID(bubble.l_bound, w)] > 0);
		for (int s = bubble.l_bound + 1; s + K - 1 <= bubble.r_bound; s++) {
			for (int w = 0; w < NODE_SHIFT; w++) {
				dp[VID(s, w)] = 0;
				if (get_in_weight(s, w, 0) > 0)
					dp[VID(s, w)] += dp[VID(s - 1, prev_node(w, 0))];
				if (get_in_weight(s, w, 1) > 0)
					dp[VID(s, w)] += dp[VID(s - 1, prev_node(w, 1))];
				dp[VID(s, w)] = std::min(dp[VID(s, w)], MAX_COMPLEXITY * 2); // Multiply 2 because of the graph symmetry
			}
		}
		for (int s = bubble.l_bound; s + K - 1 <= bubble.r_bound; s++) {
			for (int w = 0; w < NODE_SHIFT/2; w++)
				assert(dp[VID(s, w)] == dp[VID(s, w ^ NODE_MASK)]);
		}
		for (int w = 0; w < NODE_SHIFT / 2; w++) bubble.path_n += dp[VID(bubble.r_bound-K+1, w)];
		if (bubble.path_n >= MAX_COMPLEXITY) overflow_n++;
		else if (bubble.path_n == 0) incomplete_n++;
		else simple_n++;
	}
	delete[] dp;

	fprintf(stderr, "    Find %ld bubbles (%d simple, %d incomplete and %d complex)\n",
	        all_bubbles.size(), simple_n, incomplete_n, overflow_n);

	for (auto &bubble : all_bubbles) {
		if (0 < bubble.path_n and bubble.path_n < MAX_COMPLEXITY) {
			solve_small_bubble(bubble);
		}
	}

	char *hap1 = new char[snp_column.size()];
	int tan_size = 0;
	for (const auto &bubble : all_bubbles) {
		if (0 == bubble.path_n or bubble.path_n >= MAX_COMPLEXITY) {
			greedy_tangle_area(bubble);
			for (int i = bubble.l_bound; i <= bubble.r_bound; i++) {
				hint_post[i] = 3;
			}
			tan_size += bubble.r_bound - bubble.l_bound + 1;
		}
	}
	delete [] hap1;

	for (int s = 0; s < KMER_N; s++) {
		for (int w = 0; w < NODE_SHIFT; w++) {
			assert(out_degree(s, w) != 2);
		}
	}
	delete[] phase_record;
	fprintf(stderr, "    Proportion of tangle area in DBG is %.2f %% \n", 100.0 * tan_size / snp_column.size());

	output_graph("4_bubble");
}

std::vector<Phased_Block> SNP_DBG::haplotype_block() {
	std::vector<Phased_Block> ret;
	for (int32_t s = 0; s < KMER_N; s++) {
		for (int32_t w = 0; w < NODE_SHIFT/2; w++) {
			if (in_degree(s, w) == 0 and out_degree(s, w) == 1) {
				Phased_Block block(s);
				for (int k = 0; k < K; k++) block.push_back(s + k, (w >> (K-1-k)) & 1);
				int ns = s, nw = w;
				while (true) {
					int w0 = get_out_weight(ns, nw, 0);
					int w1 = get_out_weight(ns, nw, 1);
					if (w0 == 0 and w1 == 0) break;
					assert(!(w0 > 0 and w1 > 0));
					block.push_back(ns + K, w0>0 ?0 :1);
					ns++; nw = next_node(nw, w0>0 ?0 :1);
				}
				ret.push_back(block);
			}
		}
	}
	fprintf(stderr, "    %ld haplotype blocks on DBG\n", ret.size());
	return ret;
}

void SNP_DBG::resolve() {
	construct_graph();
	trim_graph();
	remove_small_tips();
	solve_bubbles();
}

void SNP_DBG::destroy() {
	free(node_cnt);
	free(edge_cnt);
	free(weight);
}

void SNP_DBG::check_graph_symmetry() {
	for (int i = 0; i < KMER_N - 1; i++) {
		for (int edge = 0; edge < EDGE_SHIFT/2; edge++) {
			assert(weight[i * EDGE_SHIFT + edge] == weight[i * EDGE_SHIFT + (edge ^ EDGE_MASK)]);
		}
	}
}

void SNP_DBG::output_graph(const char *suffix) {
	if (og_prefix == nullptr) return;
	std::string fn = std::string(og_prefix) + "." + std::string(suffix);
	FILE *fp = fopen(fn.c_str(), "w"); if (fp == nullptr) return;
	fprintf(fp, "#K=%d\n", K);
	fprintf(fp, "#Format=index(u) node(u) node(v) weight(u) weight(v) weight(u->v)\n");
	for (int32_t p = 0; p < KMER_N - 2; p++) {
		for (int32_t e = 0; e < EDGE_SHIFT; e++) {
			int w = weight[p * EDGE_SHIFT + e];
			if (w == 0) continue;

			fprintf(fp, "%d\t", p); // index(u)
			for (int k = 0; k < K; k++) fprintf(fp, "%d", ((e >> (K-k)) & 1)); fprintf(fp, "\t"); // node(u)
			for (int k = 1; k <= K; k++) fprintf(fp, "%d", ((e >> (K-k)) & 1)); fprintf(fp, "\t"); // node(v)
			fprintf(fp, "%d\t", node_cnt[VID(p, (e >> 1) & NODE_MASK)]); // weight(u)
			fprintf(fp, "%d\t", node_cnt[VID(p + 1, e & NODE_MASK)]); // weight(v)
			fprintf(fp, "%d\n", weight[EID(p, e)]); // weight(u->v)
		}
	}
	fclose(fp);
}