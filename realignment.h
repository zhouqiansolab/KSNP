//
// Created by ixiaohu on 2022/10/14.
//

#ifndef KSNP_REALIGNMENT_H
#define KSNP_REALIGNMENT_H

#include <cstdio>
#include "htslib/sam.h"

class Realignment {
private:
	int gr_len; /** Total length of global reference sequence */
	const char *global_ref; /** Global reference sequence */

	/** Take SNP as center, extract 15bp from read/reference forwardly and backwardly. */
	static const int OVERHANG_LEN = 15;
	static const int MATRIX_SIZE = OVERHANG_LEN * 2 + 5;
	int q_len, r_len; /** Length of extracted query and reference sequence */
	uint8_t que[MATRIX_SIZE]; /** Extracted query sequence (index from 1) */
	uint8_t ref[MATRIX_SIZE], alt[MATRIX_SIZE]; /** Extracted reference/alternative sequence (1-based index) */

	/** Bit-vector algorithm */
	static const int ALPHABET_SIZE = 5; /** A, C, G, T and other. */
	uint8_t ALPHA_TABLE[256]; /** Alphabet table: A->0, C->1, G->2, T->3, Other->4 */
	uint32_t peq[ALPHABET_SIZE];
public:
	Realignment(int l, const char *r);

	/**
	 * Take the input SNP as center, extract a small faction of reference and query sequence.
	 * @param aln         Aligned read in BAM format
	 * @param q_snp       SNP position on query sequence, which is already aligned to r_snp
	 * @param r_snp       SNP position on reference sequence (0-based, note that SNP in VCF is 1-based)
	 * @param alt_allele  Alternative allele used to replace reference allele
	 * @return Pair of edit distance with ref and alt allele.
	 */
	std::pair<int, int> bit_vector_dp(const bam1_t *aln, int q_snp, int r_snp, char alt_allele);

	/**
	 * Compute edit distance between extracted query and reference/alternative sequence.
	 * Run realignment for alternative sequence to remove reference bias.
	 * Mode: Global Alignment.
	 * Objective function: Levenshtein Distance.
	 * @return Pair of edit distance with ref and alt allele.
	 */
	std::pair<int, int> edit_distance();

	/**
	 * @param q and @param t is the sequence pair to align
	 * @param ql and @param tl is alignment sequence length
	 * @param r_snp is SNP position on global reference sequence (0-based)
	 * @return a pair contains edit distance and updated allele on query read
	 */
	std::pair<int, char> detect_allele(int ql, const char *q, int tl, const char *t, int r_snp);

	/**
	 * Do NOT use affine gap penalty, it is worsen than edit distance alone,
	 * possibly only suitable for long sequence alignment.
	 */
	std::pair<int, int> affine_gap();

	/** Output a pair of string which describes the alignment */
	void show_alignment(const bam1_t *aln, const char *r, int r_snp) const;
};

std::pair<int, int> edit_distance(const uint8_t *que, int que_pos, int q_len,
								  const char *ref, int ref_pos, int r_len, char alt_allele);

#endif //KSNP_REALIGNMENT_H

