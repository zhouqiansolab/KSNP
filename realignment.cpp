//
// Created by ixiaohu on 2022/10/14.
//

#include <cstring>
#include <cassert>
#include <algorithm>
#include <cctype>
#include <string>

#include "realignment.h"

Realignment::Realignment(int l, const char *r) :gr_len(l), global_ref(r) {
	q_len = r_len = 0;
	memset(ALPHA_TABLE, 4, 256 * sizeof(uint8_t));
	ALPHA_TABLE['A'] = ALPHA_TABLE['a'] = 0;
	ALPHA_TABLE['C'] = ALPHA_TABLE['c'] = 1;
	ALPHA_TABLE['G'] = ALPHA_TABLE['g'] = 2;
	ALPHA_TABLE['T'] = ALPHA_TABLE['t'] = 3;
}

std::pair<int, int> Realignment::bit_vector_dp(const bam1_t *aln, int q_snp, int r_snp, char alt_allele) {
	// Extracted query sequence
	const uint8_t *enc_seq = bam_get_seq(aln);
	int que_l = std::max(q_snp - OVERHANG_LEN, 0); // Query interval [que_l, que_r)
	int que_r = std::min(q_snp + OVERHANG_LEN + 1, aln->core.l_qseq);
	q_len = que_r - que_l; // Rows of DP matrix
	for (int i = que_l; i < que_r; i++) {
		que[i-que_l+1] = ALPHA_TABLE[seq_nt16_str[bam_seqi(enc_seq, i)]];
	}

	// Extracted reference sequence
	int ref_l = std::max(r_snp - OVERHANG_LEN, 0);
	int ref_r = std::min(r_snp + OVERHANG_LEN + 1, gr_len);
	r_len = ref_r - ref_l; // Columns of DP matrix
	for (int i = ref_l; i < ref_r; i++) {
		alt[i-ref_l+1] = ref[i-ref_l+1] = ALPHA_TABLE[global_ref[i]];
	}
	alt[r_snp-ref_l+1] = ALPHA_TABLE[alt_allele];

	// Compute Peq[σ]
	for (int a = 0; a < ALPHABET_SIZE; a++) {
		peq[a] = 0u;
		for (int i = q_len; i >= 1; i--) {
			uint32_t bit = (a == que[i]) ?1u :0u;
			peq[a] <<= 1u;
			peq[a]  |= bit;
		}
	}

	uint32_t pv = ~0u, mv = 0u; // The first column vertical
	int ref_score = q_len; // Initial score
	const uint32_t HIGH_SET = (1u << (q_len-1));
	for (int j = 1; j <= r_len; j++) { // DP loop
		uint32_t eq = peq[ref[j]];
		uint32_t xv = eq | mv;
		uint32_t xh = (((eq & pv) + pv) ^ pv) | eq;

		uint32_t ph = mv | (~ (xh | pv));
		uint32_t mh = pv & xh;

		if ((ph & HIGH_SET) != 0) ++ref_score;
		else if ((mh & HIGH_SET) != 0) --ref_score;

		ph = (ph << 1u) | 1u; // This is different from Myers' paper; our first row is 1,2,3...
		mh <<= 1u;
		pv = mh | (~(xv | ph));
		mv = ph & xv;
	}

	pv = ~0u; mv = 0u; // Calculate for alternative sequence
	int alt_score = q_len;
	for (int j = 1; j <= r_len; j++) {
		uint32_t eq = peq[alt[j]];
		uint32_t xv = eq | mv;
		uint32_t xh = (((eq & pv) + pv) ^ pv) | eq;

		uint32_t ph = mv | (~ (xh | pv));
		uint32_t mh = pv & xh;

		if ((ph & HIGH_SET) != 0) ++alt_score;
		else if ((mh & HIGH_SET) != 0) --alt_score;

		ph = (ph << 1u) | 1u;
		mh <<= 1u;
		pv = mh | (~(xv | ph));
		mv = ph & xv;
	}
//	auto truth = edit_distance();
//	assert(ref_score == truth.first and alt_score == truth.second);
	return std::make_pair(ref_score, alt_score);
}

std::pair<int, int> Realignment::edit_distance() {
	// Initialized DP matrix
	int H[MATRIX_SIZE][MATRIX_SIZE];
	for (int j = 1; j <= r_len; j++) H[0][j] = j; // Fill the first row
	for (int i = 1; i <= q_len; i++) H[i][0] = i; // Fill the first column
	H[0][0] = 0;

	// Dynamic programming loop
	for (int i = 1; i <= q_len; i++) {
		for (int j = 1; j <= r_len; j++) {
			H[i][j] = std::min(H[i-1][j], H[i][j-1]) + 1;
			H[i][j] = std::min(H[i][j], H[i-1][j-1] + ((que[i] == ref[j]) ?0 :1));
		}
	}
	int ref_ed = H[q_len][r_len];

	for (int i = 1; i <= q_len; i++) {
		for (int j = 1; j <= r_len; j++) {
			H[i][j] = std::min(H[i-1][j], H[i][j-1]) + 1;
			H[i][j] = std::min(H[i][j], H[i-1][j-1] + ((que[i] == alt[j]) ?0 :1));
		}
	}
	int alt_ed = H[q_len][r_len];

	return std::make_pair(ref_ed, alt_ed);
}

std::pair<int, char> Realignment::detect_allele(int ql, const char *q, int tl, const char *t, int r_snp) {
	int H[MATRIX_SIZE][MATRIX_SIZE];
	for (int j = 1; j <= tl; j++) H[0][j] = j;
	for (int i = 1; i <= ql; i++) H[i][0] = i;
	H[0][0] = 0;

	// Dynamic programming loop
	for (int i = 1; i <= ql; i++) {
		for (int j = 1; j <= tl; j++) {
			H[i][j] = std::min(H[i-1][j], H[i][j-1]) + 1;
			H[i][j] = std::min(H[i][j], H[i-1][j-1] + ((q[i] == t[j]) ?0 :1));
		}
	}

	std::string u, d; // A pair of string describing the alignment CIGAR
	int n = ql, m = tl;
	while (n > 0 and m > 0) {
		int M = q[n] != t[m] ?1 :0;
		if (H[n][m] == H[n - 1][m - 1] + M) {
			d += q[n];
			u += t[m];
			n--; m--;
		} else if (H[n][m] == H[n-1][m] + 1) { // Insertion into reference
			d += q[n];
			u += '-';
			n--;
		} else if (H[n][m] == H[n][m-1] + 1) { // Deletion from reference
			d += '-';
			u += t[m];
			m--;
		} else std::abort();
	}
	while (n > 0) { u += '-'; d += q[n]; n--; }
	while (m > 0) { u += t[m]; d += '-'; m--; }
	std::reverse(u.begin(), u.end());
	std::reverse(d.begin(), d.end());

	int ref_l = std::max(r_snp - OVERHANG_LEN, 0);
	int cnt = 0; char support_allele = '*';
	for (int i = 0; i < u.length(); i++) {
		if (u[i] == '-') continue;
		if (ref_l + cnt == r_snp) {
			support_allele = d[i];
			break;
		}
		cnt++;
	}
	assert(support_allele != '*');
	return std::make_pair(H[ql][tl], support_allele);
}

std::pair<int, int> Realignment::affine_gap() {
	/** Global alignment using affine-gap penalty.
	 *  H(i,j) = max{H(i-1,j-1) + M(i,j), E(i,j), F(i,j)}
	 *  E(i,j) = max{H(i,j-1) - deletion_open, E(i,j-1)} - deletion_extend
	 *  F(i,j) = max{H(i-1,j) - insertion_open, F(i-1,j)} - insertion_extend
	 * */
	const int MATCH_SCORE      = 2; // Parameters from Minimap2 for ONT/CLR alignment
	const int MISMATCH_PENALTY = 4;
	const int GAP_OPEN         = 4;
	const int GAP_EXTEND       = 2;

	int H[MATRIX_SIZE][MATRIX_SIZE], F[MATRIX_SIZE][MATRIX_SIZE], E[MATRIX_SIZE][MATRIX_SIZE];

	const int INF = 100000000;
	for (int j = 0; j <= r_len; j++) H[0][j] = E[0][j] = F[0][j] = -INF;
	H[0][0] = 0; // Must start at the first position of reference
	for (int i = 1; i <= q_len; i++) {
		H[i][0] = E[i][0] = F[i][0] = -GAP_OPEN - i * GAP_EXTEND;
		E[i][0] = -INF;
	}

	// Dynamic programming loop
	for (int i = 1; i <= q_len; i++) {
		for (int j = 1; j <= r_len; j++) {
			E[i][j] = std::max(H[i][j-1] - GAP_OPEN, E[i][j-1]) - GAP_EXTEND;
			F[i][j] = std::max(H[i-1][j] - GAP_OPEN, F[i-1][j]) - GAP_EXTEND;
			int M = que[i] != ref[j] ?-MISMATCH_PENALTY :MATCH_SCORE;
			H[i][j] = std::max(H[i-1][j-1] + M, E[i][j]);
			H[i][j] = std::max(H[i][j], F[i][j]);
		}
	}
	int ref_score = H[q_len][r_len];

	for (int i = 1; i <= q_len; i++) {
		for (int j = 1; j <= r_len; j++) {
			E[i][j] = std::max(H[i][j-1] - GAP_OPEN, E[i][j-1]) - GAP_EXTEND;
			F[i][j] = std::max(H[i-1][j] - GAP_OPEN, F[i-1][j]) - GAP_EXTEND;
			int M = que[i] != alt[j] ?-MISMATCH_PENALTY :MATCH_SCORE;
			H[i][j] = std::max(H[i-1][j-1] + M, E[i][j]);
			H[i][j] = std::max(H[i][j], F[i][j]);
		}
	}
	int alt_score = H[q_len][r_len];

	return std::make_pair(ref_score, alt_score);
}

void Realignment::show_alignment(const bam1_t *aln, const char *r, int r_snp) const {
	std::string u, d;
	auto *q = bam_get_seq(aln);
	int read_p = 0, ref_p = aln->core.pos;
	const uint32_t *cigar_array = bam_get_cigar(aln);
	for (int i = 0; i < aln->core.n_cigar; i++) {
		char op_chr = bam_cigar_opchr(cigar_array[i]); // HSIDX=
		int op_len = bam_cigar_oplen(cigar_array[i]);
		if (op_chr == 'H') { // Hard clipping is already removed from read
			continue;
		} else if (op_chr == 'S') { // Soft clipping, like insertions
			for (int j = 0; j < op_len; j++) u += '-';
			for (int j = 0; j < op_len; j++) d += seq_nt16_str[bam_seqi(q, read_p + j)];
			read_p += op_len;
		} else if (op_chr == 'D') { // Deletions from reference
			for (int j = 0; j < op_len; j++) u += r[ref_p + j];
			for (int j = 0; j < op_len; j++) d += '-';
			ref_p += op_len;
		} else if (op_chr == 'I') { // Insertions into reference
			for (int j = 0; j < op_len; j++) u += '-';
			for (int j = 0; j < op_len; j++) d += seq_nt16_str[bam_seqi(q, read_p + j)];
			read_p += op_len;
		} else { // Matches and mismatches
			for (int j = 0; j < op_len; j++) u += r[ref_p + j];
			for (int j = 0; j < op_len; j++) d += seq_nt16_str[bam_seqi(q, read_p + j)];
			read_p += op_len; ref_p += op_len;
		}
	}

	// Focus on allele area
	int p = 0, center = -1;
	for (int i = 0; i < u.length(); i++) {
		if (u[i] == '-') continue;
		if (p == r_snp - aln->core.pos) {
			center = i;
			break;
		}
		p++;
	}
	int cnt = 0, intv_l = 0, intv_r = u.length();
	for (int i = center-1; i >= 0; i--) {
		if (u[i] == '-') continue;
		cnt++;
		if (cnt == OVERHANG_LEN) {
			intv_l = i;
			break;
		}
	}
	cnt = 0;
	for (int i = center+1; i < u.length(); i++) {
		if (u[i] == '-') continue;
		cnt++;
		if (cnt == OVERHANG_LEN) {
			intv_r = i + 1;
			break;
		}
	}
	for (int i = intv_l; i < center; i++) fprintf(stdout, "%c", std::toupper(u[i]));
	fprintf(stdout, " %c ", std::toupper(u[center]));
	for (int i = center+1; i < intv_r; i++) fprintf(stdout, "%c", std::toupper(u[i]));
	fprintf(stdout, "\n");

	for (int i = intv_l; i < center; i++) fprintf(stdout, "%c", d[i]);
	fprintf(stdout, " %c ", d[center]);
	for (int i = center+1; i < intv_r; i++) fprintf(stdout, "%c", d[i]);
	fprintf(stdout, "\n");
}

std::pair<int, int> edit_distance(const uint8_t *que, int que_pos, int q_len,
								  const char *ref, int ref_pos, int r_len, char alt_allele) {
	const int SEQ_EXTEND_LEN = 15; /** Take SNP as center, extract 15bp from read/reference forwardly and backwardly. */

	/** query position and reference position is already aligned (except to deletion errors) */
	const int MATRIX_ROW = SEQ_EXTEND_LEN * 2 + 5; // Buffer size for DP matrix
	const int MATRIX_COL = SEQ_EXTEND_LEN * 2 + 5;

	/** Extracted query sequence */
	int que_l = std::max(que_pos - SEQ_EXTEND_LEN, 0); // Query interval [que_l, que_r)
	int que_r = std::min(que_pos + SEQ_EXTEND_LEN + 1, q_len);
	int row_n = que_r - que_l; // Rows of DP matrix
	char q[MATRIX_ROW];
	for (int i = que_l; i < que_r; i++) q[i - que_l + 1] = seq_nt16_str[bam_seqi(que, i)];

	/** Extracted reference sequence */
	int ref_l = std::max(ref_pos - SEQ_EXTEND_LEN, 0);
	int ref_r = std::min(ref_pos + SEQ_EXTEND_LEN + 1, r_len);
	int col_m = ref_r - ref_l; // Columns of DP matrix
	char t[MATRIX_COL];
	for (int i = ref_l; i < ref_r; i++) t[i - ref_l + 1] = (char)std::toupper(ref[i]); // Reference has lowercase bases

	int H[MATRIX_ROW][MATRIX_COL];
	for (int j = 1; j <= col_m; j++) H[0][j] = j; /** Fill the first row. */
	for (int i = 1; i <= row_n; i++) H[i][0] = i; /** Fill the first column. */
	H[0][0] = 0;

	/** Dynamic programming loop */
	for (int i = 1; i <= row_n; i++) {
		for (int j = 1; j <= col_m; j++) {
			H[i][j] = std::min(H[i-1][j], H[i][j-1]) + 1;
			H[i][j] = std::min(H[i][j], H[i-1][j-1] + ((q[i] == t[j]) ?0 :1));
		}
	}
	int ref_ed = H[row_n][col_m];

	/** Replace reference allele with alternative allele */
	t[ref_pos - ref_l + 1] = (char)std::toupper(alt_allele);
	for (int i = 1; i <= row_n; i++) {
		for (int j = 1; j <= col_m; j++) {
			H[i][j] = std::min(H[i-1][j], H[i][j-1]) + 1;
			H[i][j] = std::min(H[i][j], H[i-1][j-1] + ((q[i] == t[j]) ?0 :1));
		}
	}
	int alt_ed = H[row_n][col_m];

	return std::make_pair(ref_ed, alt_ed);
}