//
// Created by ixiaohu on 2022/1/28.
//

#include <cassert>
#include <cmath>
#include <algorithm>
#include <cstring>

#include "snp_dbg.h"
#include "htslib/sam.h"

extern std::string global_chrom;

void SNP_DBG::construct_graph(const char *bam_fn, std::vector<SNP> &snps) const {
	samFile *bam_fp = sam_open(bam_fn, "r");
	bam_hdr_t *bam_header = sam_hdr_read(bam_fp);
	bam1_t *aln = bam_init1();
	int read_id = 0;
	while (sam_read1(bam_fp, bam_header, aln) >= 0) {
		read_id++;
		int ref_start = aln->core.pos + 1;
		int bs = binary_search(snps, ref_start);
		if (bs == -1) continue;

		int read_len = 0, ref_len = 0;
		const uint32_t *cigar_array = bam_get_cigar(aln);
		for (int i = 0; i < aln->core.n_cigar; i++) {
			char op_chr = bam_cigar_opchr(cigar_array[i]);
			int op_len = bam_cigar_oplen(cigar_array[i]);
			if (op_chr != 'D' && op_chr != 'H') read_len += op_len;
			if (op_chr != 'I' && op_chr != 'S' && op_chr != 'H') ref_len += op_len;
		}

		int cid = 0, que_pointer = 0, ref_pointer = ref_start;
		uint32_t node_kmer = 0, NODE_MASK = (1<<K) - 1;
		uint32_t edge_kmer = 0, EDGE_MARK = (1<<(K+1)) - 1;
		int count = 0;
		for (int i = bs; i < snps.size(); i++) {
			auto &s = snps[i];
			if (s.pos >= ref_start + ref_len) break;
			// Find the cigar interval overlapping the SNP.
			// The cigar intervals are consecutive.
			// For example, 12D2X193I= has cigar intervals [0,12), [12, 14), [14, 33), [33, 33).
			while (cid < aln->core.n_cigar) {
				char op_chr = bam_cigar_opchr(cigar_array[cid]);
				int op_len = bam_cigar_oplen(cigar_array[cid]);
				if (op_chr == 'H') { cid++; continue; }
				if (op_chr == 'I' || op_chr == 'S') {
					que_pointer += op_len;
					cid++;
					continue;
				}
				assert(ref_pointer <= s.pos);
				if (ref_pointer + op_len > s.pos) break; // The cigar interval overlapping the SNP
				ref_pointer += op_len;
				if (op_chr != 'D') que_pointer += op_len;
				cid++;
			}

			char op_chr = bam_cigar_opchr(cigar_array[cid]);
			int que_pos = op_chr == 'D' ?que_pointer - 1 :que_pointer + s.pos - ref_pointer;
			const uint8_t *encoded_seq = bam_get_seq(aln);
			char que_base = op_chr == 'D' ?'-' :seq_nt16_str[bam_seqi(encoded_seq, que_pos)];

			if (que_base == '-' || (que_base != s.ref && que_base != s.alt)) {
				node_kmer = 0; edge_kmer = 0; count = 0;
				s.add(read_id, 2);
				continue;
			}
			count++;
			uint32_t bit = que_base == s.ref ?0 :1;
			s.add(read_id, bit);
			node_kmer = ((node_kmer << 1) | bit) & NODE_MASK;
			edge_kmer = ((edge_kmer << 1) | bit) & EDGE_MARK;
			if (count >= K) {
				add_node(i-K+1, node_kmer); // (i-K+1, i-K+2, ..., i)
				add_node(i-K+1, (~node_kmer) & NODE_MASK);
			}
			if (count >= K + 1) {
				// u = (i-K, i-K+1, ..., i-1)
				// v = (i-K+1, i-K+2, ..., i)
				add_edge(i-K+1, edge_kmer);
				add_edge(i-K+1, (~edge_kmer) & EDGE_MARK);
			}
		}
	}
	bam_destroy1(aln);
	bam_hdr_destroy(bam_header);
	sam_close(bam_fp);

	int node_intv = (1<<K), edge_intv = (1<<(K+1));
	for (int i = 1; i+K-1 < snp_n; i++) {
		// We are considering the consecutive K+1 SNPs: {i-1, i, ..., i+K-2, i+K-1}
		// u = (i-1, i, ..., i+K-2)
		// v =      (i, ..., i+K-2, i+K-1)
		for (uint32_t j = 0; j < edge_intv; j++) {
			int node_u = (j >> 1) & (node_intv - 1), node_v = j & (node_intv - 1);
			int weight_edge = edge_cnt[i * edge_intv + j];
			int weight_u = node_cnt[(i-1) * node_intv + node_u];
			int weight_v = node_cnt[i * node_intv + node_v];
			if (weight_edge > 0) assert(weight_u > 0 && weight_v > 0);
			// No read covers the K+1 SNPs.
			if (weight_v == 0 || weight_u == 0) {
				assert(weight_edge == 0);
				continue;
			}
			// But need to consider one read covers SNPs[i-1, i+K-2), another read covers SNPs[i, i+K-1).
			weight[i * edge_intv + j] = (int)std::ceil(0.2 * (weight_u + weight_v)) + weight_edge;
		}
	}
	int nodes_n = 0, edges_n = 0;
	for (int i = 0; i+K-1 < snp_n; i++) {
		for (int j = 0; j < node_intv; j++) {
			if (node_cnt[i * node_intv + j] > 0) {
				nodes_n++;
				edges_n += out_degree(i, j);
			}
		}
	}
	std::cerr << "Construct graph done; " << edges_n << " edges; " <<
	          "averaged out-degree: " << 1.0 * edges_n / nodes_n << std::endl;
}

SNP_Block SNP_DBG::traversal_block(int p, int u) {
	const int LIMIT = 1000000;
	int path_len = 1;
	SNP_Block ret(p);
	for (int i = K-1; i >= 0; i--) {
		if ((u & (1<<i)) != 0) ret.add(1);
		else ret.add(0);
	}
	while (out_degree(p, u) == 1 && path_len <= LIMIT) {
		int v0 = ((u << 1) | 0) & ((1<<(K+1)) - 1);
		int v1 = ((u << 1) | 1) & ((1<<(K+1)) - 1);
		int b;
		if (weight[(p + 1) * (1<<(K+1)) + v0] != 0) {
			assert(weight[(p + 1) * (1<<(K+1)) + v1] == 0);
			b = 0;
			u = v0 & ((1<<K) - 1);
		} else {
			assert(weight[(p + 1) * (1<<(K+1)) + v1] != 0);
			b = 1;
			u = v1 & ((1<<K) - 1);
		}
		p++;
		ret.add(b);
		path_len++;
	}
	if (path_len == 1) {
		ret.anchor = -1;
		free(ret.bit);
	}
	return ret;
}

std::vector<SNP_Block> SNP_DBG::snp_haplotype() {
	for (int i = 0; i < snp_n; i++) {
		for (int j = 0; j < (1<<K); j++) {
			int in = in_degree(i, j);
			int out = out_degree(i, j);
			assert(in != 2 && out != 2);
		}
	}

	std::vector<SNP_Block> blocks;
	for (int i = 0; i < snp_n; i++) {
		bool vis[1<<K]; memset(vis, false, sizeof(vis));
		for (int j = 0; j < (1<<K); j++) {
			if (in_degree(i, j) == 0 && out_degree(i, j) == 1 && !vis[j]) {
				vis[j] = vis[(~j & ((1<<K)-1))] = true;
				auto b = traversal_block(i, j);
				if (b.anchor == -1) continue;
				blocks.push_back(b);
			}
		}
	}
	std::cerr << "By traversal the linear graph, obtain " << blocks.size() << " block pairs" << std::endl;
	return blocks;
}

bool SNP_DBG::filter_edge(int heaviest, int edge_weight) {
	if (heaviest == edge_weight) return false;
	if (heaviest >= 4 && edge_weight == 1) return true;
	if (heaviest >= 5 and edge_weight <= 2) return true;
	if (heaviest >= 2 * edge_weight and edge_weight >= 3 && edge_weight < 5) return true;
	if (heaviest >= 1.5 * edge_weight and edge_weight >= 5) return true;
	if (heaviest >= 1.5 * edge_weight and heaviest >= 10) return true;
	return false;
}

void SNP_DBG::trim_graph() {
	int trim_n = 0;
	int node_intv = (1<<K), edge_intv = (1<<(K+1));
	for (int i = 1; i+K-1 < snp_n; i++) {
		// Edge (i-1, i, ..., i+K-2)
		//           (i, ..., i+K-2, i+K-1)
		int heaviest[2] = {0, 0};
		for (uint32_t j = 0; j < edge_intv; j++) {
			if (weight[i * edge_intv + j] == 0) continue;
			int low = (j & 1);
			heaviest[low] = std::max(weight[i * edge_intv + j], heaviest[low]);
		}

		// Filtering the lighter edges out
		bool keep_node_u[1<<K]; memset(keep_node_u, false, sizeof(keep_node_u));
		bool keep_node_v[1<<K]; memset(keep_node_v, false, sizeof(keep_node_v));
		bool keep_edge[1<<(K+1)]; memset(keep_edge, false, sizeof(keep_edge));
		for (uint32_t j = 0; j < edge_intv; j++) {
			if (weight[i * edge_intv + j] == 0) continue;
			int low = (j & 1);
			if (!filter_edge(heaviest[low], weight[i * edge_intv + j])) {
				int node_u = ((j >> 1) & (node_intv - 1));
				int node_v = (j & (node_intv-1));
				keep_node_u[node_u] = true;
				keep_node_v[node_v] = true;
				keep_edge[j] = true;
			} else trim_n++;
		}

		for (uint32_t j = 0; j < edge_intv; j++) {
			if (weight[i * edge_intv + j] == 0) continue;
			int node_u = ((j >> 1) & (node_intv - 1));
			int node_v = (j & (node_intv-1));

			// Once a node has no outgoing edges, remove this node.
			// Or the corresponding incoming edges: (i-2, i-1, ..., i+K-3) -> (i-1, i, ..., i+K-2)
			if (!keep_node_u[node_u] && i-2 >= 0) {
				weight[(i-1) * edge_intv + ((1<<K)|node_u)] = 0;
				weight[(i-1) * edge_intv + ((0<<K)|node_u)] = 0;
			}

			// Once a node has no incoming edges,
			// remove the outgoing edges: (i, i+1, ..., i+K-1) -> (i+1, i+2, ..., i+K)
			if (keep_node_u[node_u] && !keep_node_v[node_v] && i+K < snp_n) {
				weight[(i+1) * edge_intv + ((node_v<<1)|0)] = 0;
				weight[(i+1) * edge_intv + ((node_v<<1)|1)] = 0;
			}

			if (!keep_edge[j]) weight[i * edge_intv + j] = 0;
		}
	}
	int nodes_n = 0, edges_n = 0;
	for (int i = 0; i+K-1 < snp_n; i++) {
		for (int j = 0; j < node_intv; j++) {
			if (node_cnt[i * node_intv + j] > 0) {
				nodes_n++;
				edges_n += out_degree(i, j);
			}
		}
	}
	std::cerr << "Trim " << trim_n << " edges; " <<
	          "averaged out-degree: " << 1.0 * edges_n / nodes_n << std::endl;
	check_graph_symmetry();
}

void SNP_DBG::remove_tip() {
	int found_tips_n = 0;
	int node_intv = (1 << K), edge_intv = (1 << (K+1));
	for (int i = snp_n-2; i >= 0; i--) {
		for (int node_u = 0; node_u < node_intv/2; node_u++) {
			if (out_degree(i, node_u) != 2) continue;

			int u0 = (node_u << 1 | 0) & (node_intv - 1), p0 = 1;
			while (i + p0 < snp_n - 1 && out_degree(i + p0, u0) == 1 && p0 <= 1000) {
				if (weight[(i + p0 + 1) * edge_intv + (u0 << 1 | 0)] != 0) {
					u0 = (u0 << 1 | 0) & (node_intv - 1);
				} else u0 = (u0 << 1 | 1) & (node_intv - 1);
				p0++;
			}

			int u1 = (node_u << 1 | 1) & (node_intv - 1), p1 = 1;
			while (i + p1 < snp_n - 1 && out_degree(i + p1, u1) == 1 && p1 <= 1000) {
				if (weight[(i + p1 + 1) * edge_intv + (u1 << 1 | 0)] != 0) {
					u1 = (u1 << 1 | 0) & (node_intv - 1);
				} else u1 = (u1 << 1 | 1) & (node_intv - 1);
				p1++;
			}

			int who_is_bad = -1;
			int deg0 = out_degree(i + p0, u0); // assert(deg0 != 1); // Don't worry snp_n-1, which is specially judged
			int deg1 = out_degree(i + p1, u1); // assert(deg1 != 1);

			if (deg0 == 0 && deg1 == 0) {
				if (p0 >= 2 * p1) who_is_bad = 1;
				else if (p1 >= 2 * p0) who_is_bad = 0;
			} else if (deg0 > 0 && deg1 == 0) {
				if (p0 >= p1) who_is_bad = 1;
			} else if (deg1 > 0 && deg0 == 0) {
				if (p1 >= p0) who_is_bad = 0;
			} // else: there is a bubble here

			if (who_is_bad == 0) {
				weight[(i+1) * edge_intv + (node_u << 1 | 0)] = 0;
				weight[(i+1) * edge_intv + ((~(node_u << 1 | 0)) & (edge_intv-1))] = 0;
				u0 = (node_u << 1 | 0) & (node_intv - 1);
				for (int j = 1; j < p0; j++) {
					if (weight[(i+j+1) * edge_intv + (u0 << 1 | 0)] != 0) {
						weight[(i+j+1) * edge_intv + (u0 << 1 | 0)] = 0;
						weight[(i+j+1) * edge_intv + ((~(u0 << 1 | 0)) & (edge_intv-1))] = 0;
						u0 = (u0 << 1 | 0) & (node_intv - 1);
					} else {
						weight[(i+j+1) * edge_intv + (u0 << 1 | 1)] = 0;
						weight[(i+j+1) * edge_intv + ((~(u0 << 1 | 1)) & (edge_intv-1))] = 0;
						u0 = (u0 << 1 | 1) & (node_intv - 1);
					}
				}
			} else if (who_is_bad == 1) {
				weight[(i+1) * edge_intv + (node_u << 1 | 1)] = 0;
				weight[(i+1) * edge_intv + ((~(node_u << 1 | 1)) & (edge_intv-1))] = 0;
				u1 = (node_u << 1 | 1) & (node_intv - 1);
				for (int j = 1; j < p1; j++) {
					if (weight[(i+j+1) * edge_intv + (u1 << 1 | 0)] != 0) {
						weight[(i+j+1) * edge_intv + (u1 << 1 | 0)] = 0;
						weight[(i+j+1) * edge_intv + ((~(u1 << 1 | 0)) & (edge_intv-1))] = 0;
						u1 = (u1 << 1 | 0) & (node_intv - 1);
					} else {
						weight[(i+j+1) * edge_intv + (u1 << 1 | 1)] = 0;
						weight[(i+j+1) * edge_intv + ((~(u1 << 1 | 1)) & (edge_intv-1))] = 0;
						u1 = (u1 << 1 | 1) & (node_intv - 1);
					}
				}
			}

			if (who_is_bad != -1) found_tips_n++;
		}
	}
	std::cerr << "Remove " << found_tips_n << " tips forwardly; ";
	check_graph_symmetry();

	// Actually, you can reverse the graph and let the code above to process incoming edge tips.
	// If I have time, I would make the code cleaner.
	found_tips_n = 0;
	for (int i = 1; i < snp_n; i++) {
		for (int node_v = 0; node_v < node_intv/2; node_v++) {
			if (in_degree(i, node_v) != 2) continue;

			int v0 = ((0 << K | node_v) >> 1) & (node_intv - 1), p0 = 1;
			while (i - p0 > 0 && in_degree(i - p0, v0) == 1 && p0 <= 1000) {
				if (weight[(i - p0) * edge_intv + (0 << K | v0)] != 0) {
					v0 = ((0 << K | v0) >> 1) & (node_intv - 1);
				} else v0 = ((1 << K | v0) >> 1) & (node_intv - 1);
				p0++;
			}

			int v1 = ((1 << K | node_v) >> 1) & (node_intv - 1), p1 = 1;
			while (i - p1 > 0 && in_degree(i - p1, v1) == 1 && p1 <= 1000) {
				if (weight[(i - p1) * edge_intv + (0 << K | v1)] != 0) {
					v1 = ((0 << K | v1) >> 1) & (node_intv - 1);
				} else v1 = ((1 << K | v1) >> 1) & (node_intv - 1);
				p1++;
			}

			int who_is_bad = -1;
			int deg0 = in_degree(i - p0, v0); // assert(deg0 != 1);
			int deg1 = in_degree(i - p1, v1); // assert(deg1 != 1);

			if (deg0 == 0 && deg1 == 0) {
				if (p0 >= 2 * p1) who_is_bad = 1;
				else if (p1 >= 2 * p0) who_is_bad = 0;
			} else if (deg0 > 0 && deg1 == 0) {
				if (p0 >= p1) who_is_bad = 1;
			} else if (deg1 > 0 && deg0 == 0) {
				if (p1 >= p0) who_is_bad = 0;
			} // else: there is a bubble here

			if (who_is_bad == 0) {
				weight[i * edge_intv + (0 << K | node_v)] = 0;
				weight[i * edge_intv + (1 << K | ((~node_v) & (node_intv-1)))] = 0;
				v0 = ((0 << K | node_v) >> 1) & (node_intv - 1);
				for (int j = 1; j < p0; j++) {
					if (weight[(i-j) * edge_intv + (0 << K | v0)] != 0) {
						weight[(i-j) * edge_intv + (0 << K | v0)] = 0;
						weight[(i-j) * edge_intv + (1 << K | ((~v0) & (node_intv-1)))] = 0;
						v0 = ((0 << K | v0) >> 1) & (node_intv - 1);
					} else {
						weight[(i-j) * edge_intv + (1 << K | v0)] = 0;
						weight[(i-j) * edge_intv + (0 << K | ((~v0) & (node_intv-1)))] = 0;
						v0 = ((1 << K | v0) >> 1) & (node_intv - 1);
					}
				}
			} else if (who_is_bad == 1) {
				weight[i * edge_intv + (1 << K | node_v)] = 0;
				weight[i * edge_intv + (0 << K | ((~node_v) & (node_intv-1)))] = 0;
				v1 = ((1 << K | node_v) >> 1) & (node_intv - 1);
				for (int j = 1; j < p1; j++) {
					if (weight[(i-j) * edge_intv + (0 << K | v1)] != 0) {
						weight[(i-j) * edge_intv + (0 << K | v1)] = 0;
						weight[(i-j) * edge_intv + (1 << K | ((~v1) & (node_intv-1)))] = 0;
						v1 = ((0 << K | v1) >> 1) & (node_intv - 1);
					} else {
						weight[(i-j) * edge_intv + (1 << K | v1)] = 0;
						weight[(i-j) * edge_intv + (0 << K | ((~v1) & (node_intv-1)))] = 0;
						v1 = ((1 << K | v1) >> 1) & (node_intv - 1);
					}
				}
			}

			if (who_is_bad != -1) found_tips_n++;
		}
	}
	std::cerr << "Remove " << found_tips_n << " tips backwardly" << std::endl;
	check_graph_symmetry();
}

void SNP_DBG::remove_bubble() {
	int node_intv = (1 << K), edge_intv = (1 << (K+1));
	int equal_n = 0;
	for (int i = 0; i < snp_n-1; i++) {
		for (int node_u = 0; node_u < node_intv/2; node_u++) { // Only enumerating a half of nodes
			int e0 = weight[(i+1) * edge_intv + (node_u << 1 | 0)];
			int e1 = weight[(i+1) * edge_intv + (node_u << 1 | 1)];
			if (e0 == 0 || e1 == 0) continue;
			int node_r = (~node_u) & (node_intv-1);
			if (e0 > e1) {
				weight[(i+1) * edge_intv + (node_u << 1 | 1)] = 0;
				weight[(i+1) * edge_intv + (node_r << 1 | 0)] = 0;
			} else if (e0 < e1) {
				weight[(i+1) * edge_intv + (node_u << 1 | 0)] = 0;
				weight[(i+1) * edge_intv + (node_r << 1 | 1)] = 0;
			} else {
				equal_n++;
				if (equal_n % 2 == 1) {
					weight[(i+1) * edge_intv + (node_u << 1 | 1)] = 0;
					weight[(i+1) * edge_intv + (node_r << 1 | 0)] = 0;
				} else {
					weight[(i+1) * edge_intv + (node_u << 1 | 0)] = 0;
					weight[(i+1) * edge_intv + (node_r << 1 | 1)] = 0;
				}
			}
		}
	}
	std::cerr << "Break all bubbles" << std::endl;
	check_graph_symmetry();
}

void SNP_DBG::remove_tip_again() {
	// When bubbles are removed, there is no node with two outgoing edges
	// In my opinion, you detect and break bubbles first, and then remove tips
	int found_tips_n = 0;
	int node_intv = (1 << K), edge_intv = (1 << (K+1));
	for (int i = 1; i < snp_n; i++) {
		for (int node_v = 0; node_v < node_intv/2; node_v++) {
			if (in_degree(i, node_v) != 2) continue;

			// I think it is probably a slow stage, because the single path is too long
			int v0 = ((0 << K | node_v) >> 1) & (node_intv - 1), p0 = 1;
			while (i - p0 > 0 && in_degree(i - p0, v0) == 1 && p0 <= 1000) {
				if (weight[(i - p0) * edge_intv + (0 << K | v0)] != 0) {
					v0 = ((0 << K | v0) >> 1) & (node_intv - 1);
				} else v0 = ((1 << K | v0) >> 1) & (node_intv - 1);
				p0++;
			}

			int v1 = ((1 << K | node_v) >> 1) & (node_intv - 1), p1 = 1;
			while (i - p1 > 0 && in_degree(i - p1, v1) == 1 && p1 <= 1000) {
				if (weight[(i - p1) * edge_intv + (0 << K | v1)] != 0) {
					v1 = ((0 << K | v1) >> 1) & (node_intv - 1);
				} else v1 = ((1 << K | v1) >> 1) & (node_intv - 1);
				p1++;
			}

			int who_is_bad;
			if (p0 > p1) who_is_bad = 1;
			else if (p0 < p1) who_is_bad = 0;
			else {
				int w0 = weight[i * edge_intv + (0 << K | node_v)];
				int w1 = weight[i * edge_intv + (1 << K | node_v)];
				if (w0 > w1) who_is_bad = 1;
				else if (w0 < w1) who_is_bad = 0;
				else {
					if (i % 2 == 1) who_is_bad = 1;
					else who_is_bad = 0;
				}
			}

			if (who_is_bad == 0) {
				weight[i * edge_intv + (0 << K | node_v)] = 0;
				weight[i * edge_intv + (1 << K | ((~node_v) & (node_intv-1)))] = 0;
				v0 = ((0 << K | node_v) >> 1) & (node_intv - 1);
				for (int j = 1; j < p0; j++) {
					if (weight[(i-j) * edge_intv + (0 << K | v0)] != 0) {
						weight[(i-j) * edge_intv + (0 << K | v0)] = 0;
						weight[(i-j) * edge_intv + (1 << K | ((~v0) & (node_intv-1)))] = 0;
						v0 = ((0 << K | v0) >> 1) & (node_intv - 1);
					} else {
						weight[(i-j) * edge_intv + (1 << K | v0)] = 0;
						weight[(i-j) * edge_intv + (0 << K | ((~v0) & (node_intv-1)))] = 0;
						v0 = ((1 << K | v0) >> 1) & (node_intv - 1);
					}
				}
			} else {
				weight[i * edge_intv + (1 << K | node_v)] = 0;
				weight[i * edge_intv + (0 << K | (~node_v & (node_intv-1)))] = 0;
				v1 = ((1 << K | node_v) >> 1) & (node_intv - 1);
				for (int j = 1; j < p1; j++) {
					if (weight[(i-j) * edge_intv + (0 << K | v1)] != 0) {
						weight[(i-j) * edge_intv + (0 << K | v1)] = 0;
						weight[(i-j) * edge_intv + (1 << K | (~v1 & (node_intv-1)))] = 0;
						v1 = ((0 << K | v1) >> 1) & (node_intv - 1);
					} else {
						weight[(i-j) * edge_intv + (1 << K | v1)] = 0;
						weight[(i-j) * edge_intv + (0 << K | (~v1 & (node_intv-1)))] = 0;
						v1 = ((1 << K | v1) >> 1) & (node_intv - 1);
					}
				}
			}

			found_tips_n += 2;
		}
	}
	std::cerr << "Remove " << found_tips_n << " tips again" << std::endl;
	check_graph_symmetry();
}

void SNP_DBG::check_graph_symmetry() {
	for (int i = 1; i < snp_n; i++) {
		for (int edge = 0; edge < (1 << K); edge++) {
			int edge_intv = (1<<(K+1));
			int edge_r = (~edge) & (edge_intv - 1);
			assert(weight[i * edge_intv + edge] == weight[i * edge_intv + edge_r]);
		}
	}
}