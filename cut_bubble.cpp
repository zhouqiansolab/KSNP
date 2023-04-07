//
// Created by ixiaohu on 2023/2/28.
//

#include <fstream>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <queue>
#include <set>
#include "cut_bubble.h"

CUT_Bubble::CUT_Bubble() {
	for (int i = 0; i < LOG_SUM10_SIZE; i++) {
		log_sum10_lookup[i] = log10(1.0 + pow(10.0, -1.0 * i / LOG_SUM10_PCS));
	}
}

void CUT_Bubble::haplotyping_bubble(const std::vector<SNP> &snp_column, const std::vector<Read_Allele> &read_row, int L, int R, char *hap1) {
	variants.clear();
	fragments.clear();

	for (int i = L; i < R; i++) { variants.emplace_back(Variant()); } // Only for nodes

	const int QUAL_OFFSET = 33; // Phred-scaled quality value offset
	for (const auto &r : read_row) {
		Fragment fm;
		for (const auto &a : r) {
			if (a.allele == -1) continue;
			assert(L <= a.snp_idx and a.snp_idx < R); // Fragments must in bubble
			fm.allele_call.push_back(a.allele == 0 ?'0' :'1');
			fm.variant_idx.push_back(a.snp_idx - L);
		}
		for (int i = 0; i < fm.size(); i++) {
			int q = 50; // Set a same value because re-alignment provide no DP score
			fm.phred_quality.push_back(q - 33);
			double log_miscall = (q - QUAL_OFFSET) / -10.0;
			fm.prob_miscall.push_back(log_miscall);
			double prob = pow(10, log_miscall);
			fm.probability.push_back(prob);
			fm.prob_correct.push_back(log10(1 - prob));
		}
		fragments.push_back(fm);
	}

	build_read_variant_graph();
	max_likelihood_cut(hap1);
	for (int idx = 0; idx < variants.size(); idx++) {
		if (not variants[idx].flipped) {
			hap1[idx] = (hap1[idx] == '1') ? '0' : '1';
		}
	}
}

double CUT_Bubble::add_logs(double a, double b) {
	double max = std::max(a, b), min = std::min(a, b);
	return ((max-min) >= 15.7) ? max : max + log_sum10_lookup[(int)((max - min) * LOG_SUM10_PCS)];
}

void CUT_Bubble::build_read_variant_graph() {
	for (int r = 0; r < fragments.size(); r++) {
		const auto &f = fragments[r];
		for (int i = 0; i < f.size(); i++) {
			auto &v = variants[f.variant_idx[i]];
			v.read_idx.push_back(r);
			v.allele_idx.push_back(i);
		}
		for (int i = 1; i < f.size(); i++) {
			int node_u = f.variant_idx[i-1], node_v = f.variant_idx[i];
			char allele_u = f.allele_call[i-1], allele_v = f.allele_call[i];
			double prob_u = f.probability[i-1], prob_v = f.probability[i];
			variants[node_u].edges.emplace_back(Edge(node_v, r, allele_u, allele_v, prob_u, prob_v));
			variants[node_v].edges.emplace_back(Edge(node_u, r, allele_v, allele_u, prob_v, prob_u));
		}
	}
	for (auto &v : variants) std::sort(v.edges.begin(), v.edges.end());

	// All fragment and variants belong to a same non-trivial block
}

double CUT_Bubble::calc_fragment_score(const Fragment &f, const char *hap) {
	double prob_RH0 = 0, prob_RH1 = 0;
	for (int i = 0; i < f.size(); i++) {
		if (f.phred_quality[i] < MIN_ALLELE_QUAL) continue;
		double prob_miscall = f.prob_miscall[i]; // Log(p)
		double prob_correct = f.prob_correct[i]; // Log(1-p)
		assert(prob_miscall <= 0 and prob_correct <= 0);
		if (f.allele_call[i] == hap[f.variant_idx[i]]) {
			prob_RH0 += prob_correct; prob_RH1 += prob_miscall;
		} else {
			prob_RH0 += prob_miscall; prob_RH1 += prob_correct;
		}
	}
	return add_logs(prob_RH0, prob_RH1);
}

/** Calculate a edge weight for a given haplotype pair hap_u and hap_v */
static inline double calc_edge_weight(const Edge &k, char hap_u, char hap_v) {
	double p1 = k.prob_u * k.prob_v + (1 - k.prob_u) * (1 - k.prob_v);
	double p2 = k.prob_u * (1 - k.prob_v) + k.prob_v * (1 - k.prob_u);
	if (hap_u == hap_v && k.allele_u == k.allele_v) return log10(p1 / p2);
	else if (hap_u != hap_v && k.allele_u != k.allele_v) return log10(p1 / p2);
	else if (hap_u == hap_v && k.allele_u != k.allele_v) return log10(p2 / p1);
	else if (hap_u != hap_v && k.allele_u == k.allele_v) return log10(p2 / p1);
	else return 0;
}

void CUT_Bubble::init_node_score(const char *hap, int S, int T) {
	for (auto node : {S, T}) {
		if (node < 0) continue; // T might be minus
		for (int i = 0; i < variants[node].size(); i++) {
			auto &f = fragments[variants[node].read_idx[i]];
			int j = variants[node].allele_idx[i]; // The index to allele call on fragment with respect to the node
			if (f.phred_quality[j] < MIN_ALLELE_QUAL) continue;

			double prob_miscall = f.prob_miscall[j];
			double prob_correct = f.prob_correct[j];
			if (node == S) {
				if (hap[S] == f.allele_call[j]) {
					f.prob_RH0 += prob_correct; f.prob_RH1 += prob_miscall;
					f.prob_RC0 += prob_correct; f.prob_RC1 += prob_miscall;
				} else {
					f.prob_RH0 += prob_miscall; f.prob_RH1 += prob_correct;
					f.prob_RC0 += prob_miscall; f.prob_RC1 += prob_correct;
				}
			} else if (node == T) { // T is flipped(cut) compared to S
				if (hap[T] == f.allele_call[j]) {
					f.prob_RH0 += prob_miscall; f.prob_RH1 += prob_correct;
					f.prob_RC0 += prob_miscall; f.prob_RC1 += prob_correct;
				} else {
					f.prob_RH0 += prob_correct; f.prob_RH1 += prob_miscall;
					f.prob_RC0 += prob_correct; f.prob_RC1 += prob_miscall;
				}
			}

			// Update score of every node outside 2 shores(S,T) covered by the fragment
			for (j = 0; j < f.size(); j++) {
				if (f.phred_quality[j] < MIN_ALLELE_QUAL) continue;
				if (f.variant_idx[j] == S or f.variant_idx[j] == T) continue;
				prob_miscall = f.prob_miscall[j];
				prob_correct = f.prob_correct[j];
				double ph0 = f.prob_RH0, ph1 = f.prob_RH1;
				double pc0 = f.prob_RC0, pc1 = f.prob_RC1;
				if (hap[f.variant_idx[j]] == f.allele_call[j]) {
					ph0 += prob_correct; ph1 += prob_miscall; // Add node to S side of cut
					pc0 += prob_miscall; pc1 += prob_correct; // Add node to T side of cut
				} else {
					ph0 += prob_miscall; ph1 += prob_correct; // Allele mismatch so flip probability
					pc0 += prob_correct; pc1 += prob_miscall;
				}

				// old_ll is likelihood of fragment if node is added to S side
				double old_LL = add_logs(ph0, ph1);
				// new_ll is likelihood of fragment if node is added to T side
				double new_LL = add_logs(pc0, pc1);
				variants[f.variant_idx[j]].score += old_LL - new_LL;
			}
		}
	}
}

Max_Heap::Max_Heap(std::vector<Variant> &v, std::vector<int> &e):
		v(v), e(e) {
	for (int i = 0; i < e.size(); i++) v[e[i]].heap_location = i;
	for (int i = (int)e.size() / 2 - 1; i >= 0; i--) heapify(i);
}

void Max_Heap::heapify(int f) {
	int cl = f * 2 + 1;
	int cr = f * 2 + 2;
	int max_idx = f;
	if (cl < e.size() and fabs(v[e[cl]].score) > fabs(v[e[max_idx]].score)) {
		max_idx = cl;
	}
	if (cr < e.size() and fabs(v[e[cr]].score) > fabs(v[e[max_idx]].score)) {
		max_idx = cr;
	}
	if (max_idx != f) {
		std::swap(e[max_idx], e[f]);
		std::swap(v[e[max_idx]].heap_location, v[e[f]].heap_location);
		heapify(max_idx);
	}
}

int Max_Heap::top() {
	return e[0];
}

void Max_Heap::pop() {
	e[0] = e.back();
	e.resize(e.size() - 1);
	v[e[0]].heap_location = 0;
	heapify(0);
}

void Max_Heap::update(int i, double val) {
	if (fabs(v[e[i]].score) < fabs(val)) { // score increased, propagate upward
		v[e[i]].score = val;
		while (i > 0) {
			int f = (i - 1) / 2;
			if (fabs(v[e[f]].score) < fabs(v[e[i]].score)) {
				std::swap(e[f], e[i]);
				std::swap(v[e[f]].heap_location, v[e[i]].heap_location);
			} else break;
			i = f;
		}
	} else if (fabs(v[e[i]].score) > fabs(val)) { // score decreased, propagate downward
		v[e[i]].score = val;
		heapify(i);
	}
}

void CUT_Bubble::update_node_score(const char *hap, int S, int T, int a, Max_Heap &max_heap) {
	auto &add = variants[a];
	for (int i = 0; i < add.size(); i++) {
		auto &f = fragments[add.read_idx[i]];
		int j = add.allele_idx[i];
		if (f.phred_quality[j] < MIN_ALLELE_QUAL) continue;

		// Store previous fragment scores before updating
		double prev_H0 = f.prob_RH0, prev_H1 = f.prob_RH1;
		double prev_C0 = f.prob_RC0, prev_C1 = f.prob_RC1;
		double prob_miscall = f.prob_miscall[j];
		double prob_correct = f.prob_correct[j];

		if (add.parent == S) {
			// Node is added to S side, original and new likelihoods are updated identically
			if (hap[a] == f.allele_call[j]) {
				f.prob_RH0 += prob_correct; f.prob_RH1 += prob_miscall;
				f.prob_RC0 += prob_correct; f.prob_RC1 += prob_miscall;
			} else {
				f.prob_RH0 += prob_miscall; f.prob_RH1 += prob_correct;
				f.prob_RC0 += prob_miscall; f.prob_RC1 += prob_correct;
			}
		} else if (add.parent == T) {
			// Node is added to T side, likelihoods are flipped
			if (hap[a] == f.allele_call[j]) {
				f.prob_RH0 += prob_miscall; f.prob_RH1 += prob_correct;
				f.prob_RC0 += prob_miscall; f.prob_RC1 += prob_correct;
			} else {
				f.prob_RH0 += prob_correct; f.prob_RH1 += prob_miscall;
				f.prob_RC0 += prob_correct; f.prob_RC1 += prob_miscall;
			}
		}
		for (j = 0; j < f.size(); j++) {
			if (f.phred_quality[j] < MIN_ALLELE_QUAL) continue;
			int node = f.variant_idx[j];
			if (variants[node].parent == S or
			    variants[node].parent == T or node == a) continue;

			double new_score = variants[node].score;
			double curr_H0 = prev_H0, curr_H1 = prev_H1;
			double curr_C0 = prev_C0, curr_C1 = prev_C1;
			prob_miscall = f.prob_miscall[j];
			prob_correct = f.prob_correct[j];
			if (hap[node] == f.allele_call[j]) {
				curr_H0 += prob_correct; curr_H1 += prob_miscall;
				curr_C0 += prob_miscall; curr_C1 += prob_correct;
			} else {
				curr_H0 += prob_miscall; curr_H1 += prob_correct;
				curr_C0 += prob_correct; curr_C1 += prob_miscall;
			}
			// Subtract old score
			double old_LL = add_logs(curr_H0, curr_H1);
			double new_LL = add_logs(curr_C0, curr_C1);
			new_score -= (old_LL - new_LL);

			curr_H0 = f.prob_RH0; curr_H1 = f.prob_RH1;
			curr_C0 = f.prob_RC0; curr_C1 = f.prob_RC1;
			if (hap[node] == f.allele_call[j]) {
				curr_H0 += prob_correct; curr_H1 += prob_miscall;
				curr_C0 += prob_miscall; curr_C1 += prob_correct;
			} else {
				curr_H0 += prob_miscall; curr_H1 += prob_correct;
				curr_C0 += prob_correct; curr_C1 += prob_miscall;
			}
			// Add new delta LL
			old_LL = add_logs(curr_H0, curr_H1);
			new_LL = add_logs(curr_C0, curr_C1);
			new_score += (old_LL - new_LL);

			// Update max heap also the variant/node score
			max_heap.update(variants[node].heap_location, new_score);
		}
	}
}

static inline bool cmp(const Cut_Edge &a, const Cut_Edge &b) {
	if (fabs(a.w - b.w) > 0.00001) return a.w > b.w;
	else return a.u < b.u;
}

double CUT_Bubble::max_likelihood_cut(const char *hap) {
	// Update cut edge weight for the current haplotype
	int total_edges = 0;
	for (int idx = 0; idx < variants.size(); idx++) {
		auto &v = variants[idx];
		int prev = -1; v.cut_edges.clear();
		for (const auto &e : v.edges) {
			if (e.v != prev) {
				double w = calc_edge_weight(e, hap[idx], hap[e.v]);
				w /= (fragments[e.frag_idx].size() - 1); // scale of 1/(k-1) in HapCut2 paper
				v.cut_edges.emplace_back(Cut_Edge(idx, e.v, w));
				total_edges++;
			} else {
				double w = calc_edge_weight(e, hap[idx], hap[e.v]);
				w /= (fragments[e.frag_idx].size() - 1);
				v.cut_edges.back().w += w;
			}
			prev = e.v;
		}
	}
	// Mentioned in HapCUT2 paper, first k edges to choose
	int K = std::min(total_edges / 2, 5);

	// Find top K biggest edges in graph
	std::vector<Cut_Edge> topK_edges;
	double min_weight = 10000; int min_id = 0;
	for (const auto &v : variants) {
		for (const auto &e : v.cut_edges) {
			if (topK_edges.size() < K) { // Maintain the top-K weight edges
				topK_edges.push_back(e);
				if (e.w < min_weight) {
					min_weight = e.w;
					min_id = (int)topK_edges.size() - 1;
				}
			} else if (e.w > min_weight) { // Replace the minimum-weight edge
				topK_edges[min_id] = e;
				min_weight = 10000;
				for (int k = 0; k < topK_edges.size(); k++) {
					if (topK_edges[k].w < min_weight) {
						min_weight = topK_edges[k].w;
						min_id = k;
					}
				}
			}
		}
	}
	std::sort(topK_edges.begin(), topK_edges.end(), cmp);

	// Edge contraction algorithm: merge vertices until only two nodes left or weight of graph is negative
	int max_iteration_n = (int)variants.size() / 10;
	max_iteration_n = std::max(max_iteration_n, 1);
	max_iteration_n = std::min(max_iteration_n, MAX_CUT_ROUND);
	int iter_since_improved = 0;
	double best_cut_score = -10000;
	for (auto &v : variants) v.flipped = false;
	for (int it = 0; it < max_iteration_n + K; it++) { // At most N/10 + K iterations
		int S, T; // Edge (s,t) is to cut
		double r = -1;
		if (it < K) { // In first K iteration, determinately take nodes of top K edges as shores
			S = topK_edges[it].u;
			T = topK_edges[it].v;
		} else {
			if (drand48() < 0.5) {
				r = drand48();
				int n = (int)(r * total_edges - 0.0001);
				for (int idx = 0; idx < variants.size(); idx++) {
					if (n < variants[idx].cut_edges.size()) { S = idx; break; };
					n -= variants[idx].cut_edges.size();
				}
				T = variants[S].cut_edges[n].v;
				if (variants[S].cut_edges[n].w >= 1) continue;
			} else {
				// Find node with high MEC score, initialize as S shore
				r = drand48();
				int n = (int)(r * variants.size());
				if (n >= variants.size()) n--;
				S = n;
				T = -1;
			}
		}

		for (int idx = 0; idx < variants.size(); idx++) { // Reset node scores
			auto &v = variants[idx];
			v.parent = idx;
			v.score = 0;
		}

		for (auto &f : fragments) { // Reset fragment scores
			f.prob_RH0 = f.prob_RH1 = 0;
			f.prob_RC0 = f.prob_RC1 = 0;
		}
		init_node_score(hap, S, T);

		std::vector<int> heap_elements;
		for (int idx = 0; idx < variants.size(); idx++) {
			// s and t are not in max heap
			if (idx != S and idx != T) {
				heap_elements.push_back(idx);
			}
		}
		Max_Heap max_heap(variants, heap_elements);

		int N = (int)variants.size() - 1 - (T != -1);
		while (N > 0) { // More than two clusters, this loop is O(N^2)
			int top_id = max_heap.top(); max_heap.pop();
			auto &v = variants[top_id];
			if (v.score > 0) {
				v.parent = S;
			} else if (v.score < 0) {
				if (T < 0) T = top_id;
				v.parent = T;
			} else if (T < 0) {
				T = top_id;
			} else { // score is 0
				if (drand48() < 0.5) v.parent = S;
				else v.parent = T;
			}
			N--;
			update_node_score(hap, S, T, top_id, max_heap);
		}

		if (T == -1) continue; // Cut is empty, ignore it

		// Calculate current cut score
		double new_ll = 0, old_ll = 0;
		for (const auto &f : fragments) {
			new_ll += add_logs(f.prob_RC0, f.prob_RC1);
			old_ll += calc_fragment_score(f, hap);
		}
		double curr_cut_score = new_ll - old_ll;
		if (curr_cut_score > best_cut_score) {
			best_cut_score = curr_cut_score;
			for (auto &v : variants) {
				v.flipped = (v.parent == S);
			}
		} else iter_since_improved++;
		if (iter_since_improved > CUT_CONVERGE) break;
	}
}