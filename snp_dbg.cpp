//
// Created by ixiaohu on 2022/1/28.
//

#include <cassert>
#include <cmath>
#include <algorithm>
#include <cstring>

#include "snp_dbg.h"

extern std::string global_chrom;

extern std::vector<SNP> global_snps; // For debug only

void SNP_DBG::construct_graph(std::vector<SNP> &snps) const {
	const int SAM_BUF_SIZE = 24 * 1024 * 1024;
	char *buf = new char[SAM_BUF_SIZE];
	Parsed_CIGAR cigar;
	int read_id = 0;
	while (true) {
		if (!fgets(buf, SAM_BUF_SIZE, stdin)) break;
		if (buf[0] == '@') continue;
		read_id++;

		// Parse SAM record line
		SAM_Record sam(buf);

		int bs = binary_search(snps, sam.ref_start);
		if (bs == -1) continue;

		cigar.parse(sam.cigar);
		int read_len = 0, ref_len = 0;
		for (int i = 0; i < cigar.n; i++) {
			const auto &p = std::make_pair(cigar.num[i], cigar.op[i]);
			if (p.second != 'D') read_len += p.first;
			if (p.second != 'I' && p.second != 'S') ref_len += p.first;
		}

		int cid = 0, que_pointer = 0, ref_pointer = sam.ref_start;
		int que_tag = 0;
		uint32_t node_kmer = 0, NODE_MASK = (1<<K) - 1;
		uint32_t edge_kmer = 0, EDGE_MARK = (1<<(K+1)) - 1;
		int count = 0;
		for (int i = bs; i < snps.size(); i++) {
			auto &s = snps[i];
			if (s.pos > sam.ref_start + ref_len) break;
			// Find the cigar interval overlapping the SNP.
			// The cigar intervals are consecutive.
			// For example, 12D2X193I= has cigar intervals [0,12), [12, 14), [14, 33), [33, 33).
			while (cid < cigar.n) {
				const auto &p = std::make_pair(cigar.num[cid], cigar.op[cid]);
				if (p.second == 'I' || p.second == 'S') {
					que_pointer += p.first;
					cid++;
					continue;
				}
				assert(ref_pointer <= s.pos);
				if (ref_pointer + p.first > s.pos) break; // The cigar interval overlapping the SNP
				ref_pointer += p.first;
				if (p.second != 'D') que_pointer += p.first;
				if ((sam.flag & 0x10) == 0) {
					if (p.second == '=' || p.second == 'X') que_tag = que_pointer + 1;
					else if (p.second == 'D') que_tag = que_pointer;
				} else {
					if (p.second == '=' || p.second == 'X') que_tag = read_len - 1 - que_pointer;
					else if (p.second == 'D') que_tag = que_pointer;
				}
				cid++;
			}
			char que_base;
			// FIXME: I think it is wrong, but I have to respect the original program
			if (s.pos == sam.ref_start + ref_len) {
				if (que_tag < read_len) que_base = que_tag == -1 ?sam.bases[read_len-1] :sam.bases[que_tag];
				else que_base = '-';
			} else {
				assert(cid < cigar.n);
				const auto &p = std::make_pair(cigar.num[cid], cigar.op[cid]);
				int que_pos = p.second == 'D' ?que_pointer - 1 :que_pointer + s.pos - ref_pointer;
				que_base = p.second == 'D' ?'-' :sam.bases[que_pos];
			}

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
	delete [] buf;
	cigar.destroy();

	int node_intv = (1<<K), edge_intv = (1<<(K+1));
	for (int i = K-1; i < snps.size() - K+1; i++) {
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
	std::cerr << "Construct graph done" << std::endl;
}

#include <sstream>
struct DBG_Node {
	std::string chrom;
	std::vector<int> pos;
	char *bases;
	DBG_Node(int K, const std::string &a);
	void output();
};

DBG_Node::DBG_Node(int K, const std::string &a) {
	int i;
	for (i = 1; i < a.length(); i++) {
		if (a[i] == '_') {
			i++;
			break;
		}
		chrom += a[i];
	}
	for (int k = 0; k < K; k++) {
		int num = 0;
		for (; i < a.length(); i++) {
			if (a[i] == '_') {
				pos.push_back(num);
				i++;
				break;
			}
			num *= 10;
			num += a[i] - '0';
		}
	}
	bases = (char*) malloc((K + 1) *  sizeof(char));
	for (int k = 0; k < K; k++) {
		bases[k] = a[i + k];
	}
	bases[K] = '\0';
}

void DBG_Node::output() {
	std::cerr << chrom << "\t";
	for (auto p : pos) std::cerr << p << "\t";
	std::cerr << bases << std::endl;
}

int parse_label(const std::string &a) {
	for (int i = 0; i < a.length(); i++) {
		if (a[i] == '"') {
			int ret = 0;
			for (int j = i+1; j < a.length(); j++) {
				if (a[j] == '"') return ret;
				ret *= 10;
				ret += a[j] - '0';
			}
		}
	}
	return -1;
}

void SNP_DBG::load_graph_from_stdin(const std::vector<SNP> &snps) const {
	std::string line;
	int cnt = 0;
	while (std::getline(std::cin, line)) {
		if (line.find('{') != -1 || line.find('}') != -1) continue;
		std::stringstream ss(line);
		std::string node_u; ss >> node_u;
		std::string dummy; ss >> dummy;
		std::string node_v; ss >> node_v;
		std::string label; ss >> label;
		std::string color; ss >> color;
		DBG_Node u(K, node_u);
		DBG_Node v(K, node_v);
		global_chrom = u.chrom;
		for (int i = 0; i < K-1; i++) {
			assert(u.pos[i+1] == v.pos[i]);
			assert(u.bases[i+1] == v.bases[i]);
		}

		// Locate the position of the edge on SNPs list
		int snp_id = binary_search(snps, u.pos[0]);
		assert(snp_id != -1 && snps[snp_id].pos == u.pos[0]);
		std::string edge_bases = u.bases; edge_bases += v.bases[K-1];

		// Encode the edge bases bitwise
		uint32_t encoded_bases = 0;
		for(int i = 0; i <= K; i++) {
			encoded_bases <<= 1;
			assert(edge_bases[i] == snps[snp_id + i].alt || edge_bases[i] == snps[snp_id + i].ref);
			uint32_t bit = (edge_bases[i] == snps[snp_id + i].alt) ?1 :0;
			encoded_bases |= bit;
		}
		weight[(snp_id + 1) * (1<<(K+1)) + encoded_bases] = parse_label(label);
		assert(parse_label(label) > 0);

		// Check the correctness of encoding
		std::string temp;
		for (int i = 0; i <= K; i++) {
			if ((encoded_bases & (1<<i)) != 0) temp += snps[snp_id + K-i].alt;
			else temp += snps[snp_id + K-i].ref;
		}
		std::reverse(temp.begin(), temp.end());
		assert(edge_bases == temp);
		cnt++;
	}
	std::cerr << "Load " << cnt << " edges from stdin" << std::endl;
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
			} else {
				trim_n++;
			}
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
	std::cerr << "Trim " << trim_n << " edges" << std::endl;
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
	for (int i = 0; i < snp_n-1; i++) {
		for (int node_u = 0; node_u < node_intv/2; node_u++) { // Only enumerating a half of nodes
			int e0 = weight[(i+1) * edge_intv + (node_u << 1 | 0)];
			int e1 = weight[(i+1) * edge_intv + (node_u << 1 | 1)];
			if (e0 == 0 || e1 == 0) continue;
			// Break bubbles in a relatively brute way
			int node_r = (~node_u) & (node_intv-1);
			if (e0 >= e1) {
				weight[(i+1) * edge_intv + (node_u << 1 | 0)] = 0;
				weight[(i+1) * edge_intv + (node_r << 1 | 1)] = 0;
			} else {
				weight[(i+1) * edge_intv + (node_u << 1 | 1)] = 0;
				weight[(i+1) * edge_intv + (node_r << 1 | 0)] = 0;
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
			int deg0 = in_degree(i - p0, v0); // assert(deg0 == 0);
			int deg1 = in_degree(i - p1, v1); // assert(deg1 == 0); // I found the case whose path length > 1000
			if (p0 > p1) who_is_bad = 1;
			else if (p0 < p1) who_is_bad = 0;
			else {
				int w0 = weight[i * edge_intv + (0 << K | node_v)];
				int w1 = weight[i * edge_intv + (1 << K | node_v)];
				if (w0 > w1) who_is_bad = 1;
				else if (w0 < w1) who_is_bad = 0;
				else {
					char base0 = global_snps[i-1].ref, base1 = global_snps[i-1].alt;
					if (base0 > base1) who_is_bad = 1;
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

void SNP_DBG::output(const std::vector<SNP> &snps) {
	int edge_intv = (1<<(K+1));
	int snp_buf[K+1];
	const std::string color[] = {"red", "blue", "green"};
	fprintf(stdout, "digraph 1 {\n");
	for (int i = 1; i+K-1 < snps.size(); i++) {
		for (uint32_t j = 0; j < edge_intv; j++) {
			for (int k = 0; k < K+1; k++) snp_buf[k] = ((j>>k) & 1);
			int depth = weight[i * edge_intv + j];
			if (depth == 0) continue;

			fprintf(stdout, "\"%s_", global_chrom.c_str());
			for (int k = 0; k < K; k++) fprintf(stdout, "%d_", snps[i-1+k].pos);
			// SNP_buf K,   K-1, ..., 1,     0
			// SNP     i-1, i,   ..., i+K-2, i+K-1
			for (int k = 0; k < K; k++) fprintf(stdout, "%c", snp_buf[K-k] == 1 ?snps[i-1+k].alt :snps[i-1+k].ref);
			fprintf(stdout, "\" -> ");

			fprintf(stdout, "\"%s_", global_chrom.c_str());
			for (int k = 0; k < K; k++) fprintf(stdout, "%d_", snps[i+k].pos);
			for (int k = 0; k < K; k++) fprintf(stdout, "%c", snp_buf[K-1-k] == 1 ?snps[i+k].alt :snps[i+k].ref);
			fprintf(stdout, "\" ");

			int cid;
			if (depth < 5) cid = 0;
			else if (depth < 11) cid = 1;
			else cid = 2;
			fprintf(stdout, "[label=\"%d\" color=\"%s\"]\n", depth, color[cid].c_str());
		}
	}
	fprintf(stdout, "}\n");
}

SAM_Record::SAM_Record(char *s) {
	int i = 0;
	que_name = s; for (; s[i]; i++) if (s[i] == '\t') { s[i] = '\0'; i++; break; }
	flag = 0; for (; s[i]; i++) if (s[i] == '\t') { i++; break;} else { flag *= 10; flag += s[i] - '0'; }
	ref_name = s + i; for (; s[i]; i++) if (s[i] == '\t') { s[i] = '\0'; i++; break; }
	ref_start = 0; for (; s[i]; i++) if (s[i] == '\t') { i++; break;} else { ref_start *= 10; ref_start += s[i] - '0'; }
	for (; s[i]; i++) if (s[i] == '\t') { i++; break; }
	cigar = s + i; for (; s[i]; i++) if (s[i] == '\t') { s[i] = '\0'; i++; break; }
	for (; s[i]; i++) if (s[i] == '\t') { i++; break; }
	for (; s[i]; i++) if (s[i] == '\t') { i++; break; }
	for (; s[i]; i++) if (s[i] == '\t') { i++; break; }
	bases = s + i; for (; s[i]; i++) if (s[i] == '\t') { s[i] = '\0'; break; }
}

void Parsed_CIGAR::parse(const char *cigar) {
	n = 0;
	int x = 0;
	for (int i = 0; cigar[i]; i++) {
		char c = cigar[i];
		if (c >= '0' && c <= '9') {
			x *= 10;
			x += c - '0';
		} else {
			if (n >= m) {
				m <<= 1;
				num = (int*) realloc(num, m * sizeof(int));
				op = (char*) realloc(op, m * sizeof(char));
			}
			num[n] = x;
			op[n] = c;
			n++;
			x = 0;
		}
	}
}