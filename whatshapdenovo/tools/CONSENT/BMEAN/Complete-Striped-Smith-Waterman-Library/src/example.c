/*	example.c
 *	This is a simple example to show you how to use the SSW C library.
 *	To run this example:
 *	1) gcc -Wall -lz ssw.c example.c
 *	2) ./a.out
 *	Created by Mengyao Zhao on 07/31/12.
 */

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include "ssw.h"

//	Print the BLAST like output.
static void ssw_write (const s_align* a,
			const char* ref_seq,
			const char* read_seq,
			const int8_t* table) {

	fprintf(stdout, "optimal_alignment_score: %d\tsub-optimal_alignment_score: %d\t", a->score1, a->score2);
	if (a->ref_begin1 + 1) fprintf(stdout, "target_begin: %d\t", a->ref_begin1 + 1);
	fprintf(stdout, "target_end: %d\t", a->ref_end1 + 1);
	if (a->read_begin1 + 1) fprintf(stdout, "query_begin: %d\t", a->read_begin1 + 1);
	fprintf(stdout, "query_end: %d\n\n", a->read_end1 + 1);
	if (a->cigar) {
		int32_t c = 0, left = 0, e = 0, qb = a->ref_begin1, pb = a->read_begin1;
		uint32_t i;
		while (e < a->cigarLen || left > 0) {
			int32_t count = 0;
			int32_t q = qb;
			int32_t p = pb;
			fprintf(stdout, "Target: %8d    ", q + 1);
			for (c = e; c < a->cigarLen; ++c) {
				char letter = cigar_int_to_op(a->cigar[c]);
				uint32_t length = cigar_int_to_len(a->cigar[c]);
				uint32_t l = (count == 0 && left > 0) ? left: length;
				for (i = 0; i < l; ++i) {
					if (letter == 'I') fprintf(stdout, "-");
					else {
						fprintf(stdout, "%c", *(ref_seq + q));
						++ q;
					}
					++ count;
					if (count == 60) goto step2;
				}
			}
step2:
			fprintf(stdout, "    %d\n                    ", q);
			q = qb;
			count = 0;
			for (c = e; c < a->cigarLen; ++c) {
				char letter = cigar_int_to_op(a->cigar[c]);
				uint32_t length = cigar_int_to_len(a->cigar[c]);
				uint32_t l = (count == 0 && left > 0) ? left: length;
				for (i = 0; i < l; ++i){
					if (letter == 'M') {
						if (table[(int)*(ref_seq + q)] == table[(int)*(read_seq + p)])fprintf(stdout, "|");
						else fprintf(stdout, "*");
						++q;
						++p;
					} else {
						fprintf(stdout, "*");
						if (letter == 'I') ++p;
						else ++q;
					}
					++ count;
					if (count == 60) {
						qb = q;
						goto step3;
					}
				}
			}
step3:
			p = pb;
			fprintf(stdout, "\nQuery:  %8d    ", p + 1);
			count = 0;
			for (c = e; c < a->cigarLen; ++c) {
				char letter = cigar_int_to_op(a->cigar[c]);
				uint32_t length = cigar_int_to_len(a->cigar[c]);
				uint32_t l = (count == 0 && left > 0) ? left: length;
				for (i = 0; i < l; ++i) {
					if (letter == 'D') fprintf(stdout, "-");
					else {
						fprintf(stdout, "%c", *(read_seq + p));
						++p;
					}
					++ count;
					if (count == 60) {
						pb = p;
						left = l - i - 1;
						e = (left == 0) ? (c + 1) : c;
						goto end;
					}
				}
			}
			e = c;
			left = 0;
end:
			fprintf(stdout, "    %d\n\n", p);
		}
	}
}

//	Align a pair of genome sequences.
int main (int argc, char * const argv[]) {
	int32_t l, m, k, match = 2, mismatch = 2, gap_open = 3, gap_extension = 1;	// default parameters for genome sequence alignment
	// reference sequence
	static const char ref_seq[40] = {'C', 'A', 'G', 'C', 'C', 'T', 'T', 'T', 'C', 'T', 'G', 'A', 'C', 'C', 'C', 'G', 'G', 'A', 'A', 'A', 'T',
						'C', 'A', 'A', 'A', 'A', 'T', 'A', 'G', 'G', 'C', 'A', 'C', 'A', 'A', 'C', 'A', 'A', 'A', '\0'};
	static const char read_seq[16] = {'C', 'T', 'G', 'A', 'G', 'C', 'C', 'G', 'G', 'T', 'A', 'A', 'A', 'T', 'C', '\0'};	// read sequence
	s_profile* profile;
	int8_t* num = (int8_t*)malloc(16);	// the read sequence represented in numbers
	int8_t* ref_num = (int8_t*)malloc(64);	// the read sequence represented in numbers
	s_align* result;

	/* This table is used to transform nucleotide letters into numbers. */
	static const int8_t nt_table[128] = {
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 4, 4, 4,  3, 0, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 4, 4, 4,  3, 0, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
	};

	// initialize scoring matrix for genome sequences
	//  A  C  G  T	N (or other ambiguous code)
	//  2 -2 -2 -2 	0	A
	// -2  2 -2 -2 	0	C
	// -2 -2  2 -2 	0	G
	// -2 -2 -2  2 	0	T
	//	0  0  0  0  0	N (or other ambiguous code)
	int8_t* mat = (int8_t*)calloc(25, sizeof(int8_t));
	for (l = k = 0; l < 4; ++l) {
		for (m = 0; m < 4; ++m) mat[k++] = l == m ? match : - mismatch;	/* weight_match : -weight_mismatch */
		mat[k++] = 0; // ambiguous base: no penalty
	}
	for (m = 0; m < 5; ++m) mat[k++] = 0;

	for (m = 0; m < 15; ++m) num[m] = nt_table[(int)read_seq[m]];
	profile = ssw_init(num, 15, mat, 5, 2);
	for (m = 0; m < 39; ++m) ref_num[m] = nt_table[(int)ref_seq[m]];

	// Only the 8 bit of the flag is setted. ssw_align will always return the best alignment beginning position and cigar.
	result = ssw_align (profile, ref_num, 39, gap_open, gap_extension, 1, 0, 0, 15);
	ssw_write(result, ref_seq, read_seq, nt_table);

	align_destroy(result);
	init_destroy(profile);
	free(mat);
	free(ref_num);
	free(num);
	return(0);
}
