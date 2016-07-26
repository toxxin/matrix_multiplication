#include <stdio.h>
#include <stdlib.h>
#include "matrix_mult_strassen.h"


void print_matrix(int *m, int n) {
	unsigned int row, col;

	for (row = 0; row < n; ++row) {
		for (col = 0; col < n; ++col) {
			printf("%d ", *(m + n*row + col));
		}
		printf("\n");
	}
	printf("\n");
}

static void sum(int *aMatrix, int *bMatrix, int *cMatrix, size_t n) {
	int i;
	for (i = 0; i < n*n; ++i)
		*(cMatrix + i) = *(aMatrix + i) + *(bMatrix + i);
}

/* Two-dimentional array represented in memory as a linear array, so
 * all we need is just sum each element in this array */
static void sum_quarter(subm *s1, subm *s2, int *sMatrix, size_t n) {
	int s1_i, s1_j, s2_i, s2_j, i, j;

	for (s1_i = s1->rs, s2_i = s2->rs, i = 0; i < n; ++i, ++s1_i, ++s2_i) {
		for (s1_j = s1->cs, s2_j = s2->cs, j = 0; j < n; ++j, ++s1_j, ++s2_j) {
			*(sMatrix + n*i + j) = *(s1->M + s1_i*n*2 + s1_j) + *(s2->M + s2_i*n*2 + s2_j);
		}
	}
}

static void matrix_sub(int *aMatrix, int *bMatrix, int *cMatrix, size_t n) {
	int i;
	for (i = 0; i < n*n; ++i)
		*(cMatrix + i) = *(aMatrix + i) - *(bMatrix + i);
}

static void sub_quarter(subm *s1, subm *s2, int *sMatrix, size_t n) {
	int s1_i, s1_j, s2_i, s2_j, i, j;

	for (s1_i = s1->rs, s2_i = s2->rs, i = 0; i < n; ++i, ++s1_i, ++s2_i) {
		for (s1_j = s1->cs, s2_j = s2->cs, j = 0; j < n; ++j, ++s1_j, ++s2_j) {
			*(sMatrix + n*i + j) = *(s1->M + s1_i*n*2 + s1_j) - *(s2->M + s2_i*n*2 + s2_j);
		}
	}
}

static void copy_matrix(subm *s, int *sMatrix, size_t n) {
	unsigned int si, sj, i, j;

	for (si = s->rs, i =0; i < n/2; ++i, ++si) {
		for (sj = s->cs, j=0; j < n/2; ++j, ++sj) {
			*(sMatrix + (n/2)*i + j) = *(s->M + si*n + sj);
		}
	}
}

static void merge(int *C11, int *C12, int *C21, int *C22, int *result, size_t n) {
	int i, j;

	for (i = 0; i < n/2; ++i) {
		for (j = 0; j < n/2; ++j) {
			*(result + i*n + j) = *(C11 + i*n/2 + j);
			*(result + i*n + n/2 + j) = *(C12 + i*n/2 + j);
			*(result + i*n + n*n/2 + j) = *(C21 + i*n/2 + j);
			*(result + i*n + n*n/2 + n/2 + j) = *(C22 + i*n/2 + j);
		}
	}
}

/* Naive multiplication for comparing result with others algorithms
 * and optimisation approaches. We know that we work with linear array
 * and we have [n x n] size for aMatrix and bMatrix and [n/2 x n/2] for
 * pMatrix */
void matrix_mult_naive(int *aMatrix, int *bMatrix, int *pMatrix, size_t n) {
	unsigned int row, col, inner;

	for (row = 0; row < n; ++row) {
		for (col = 0; col < n; ++col) {
			*(pMatrix + n*row + col) = 0;
			/* Multiply the row of A by the column of B to get the row, column of product. */
			for (inner = 0; inner < n; ++inner) {
				*(pMatrix + n*row + col) += (*(aMatrix + row*n + inner)) * (*(bMatrix + inner*n + col));
			}
		}
	}
}

/* Strassen's Matrix Multiplication Algorithm which implement 2^n dimentional
 * matrix multimplication. Assume that we have matrix with power of 2 size */
//void matrix_mult_strassen(int *aMatrix, int *bMatrix int *pMatrix, size_t n) {
void matrix_mult_strassen(subm *sub1, subm *sub2, int *pMatrix, size_t n) {

	if (n == 16) {
		matrix_mult_naive(sub1->M, sub2->M, pMatrix, n);
	} else {
		int M1[n/2][n/2];
		int M2[n/2][n/2];
		int M3[n/2][n/2];
		int M4[n/2][n/2];
		int M5[n/2][n/2];
		int M6[n/2][n/2];
		int M7[n/2][n/2];
		
		int C11[n/2][n/2];
		int C12[n/2][n/2];
		int C21[n/2][n/2];
		int C22[n/2][n/2];

		int tmp1[n/2][n/2];
		int tmp2[n/2][n/2];

		/* M1 calculation */
		subm s1 = { .rs = 0, .cs = 0, .M = sub1->M};
		subm s2 = { .rs = n/2, .cs = n/2, .M = sub1->M };
		sum_quarter(&s1, &s2, &tmp1[0][0], n/2);

		s1.rs = 0; s1.cs = 0; s1.M = sub2->M;
		s2.rs = n/2; s2.cs = n/2; s2.M = sub2->M;
		sum_quarter(&s1, &s2, &tmp2[0][0], n/2);

		s1.rs = 0; s1.cs = 0; s1.M = &tmp1[0][0];
		s2.rs = 0; s2.cs = 0; s2.M = &tmp2[0][0];
		matrix_mult_strassen(&s1, &s2, &M1[0][0], n/2);

		/* M2 calculation */
		s1.rs = n/2; s1.cs = 0; s1.M = sub1->M;
		s2.rs = n/2; s2.cs = n/2; s2.M = sub1->M;
		sum_quarter(&s1, &s2, &tmp1[0][0], n/2);

		s1.rs = 0; s1.cs = 0; s1.M = sub2->M;
		copy_matrix(&s1, &tmp2[0][0], n);

		s1.rs = 0; s1.cs = 0; s1.M = &tmp1[0][0];
		s2.rs = 0; s2.cs = 0; s2.M = &tmp2[0][0];
		matrix_mult_strassen(&s1, &s2, &M2[0][0], n/2);

		/* M3 calculation:  M3 := A11(B12 - B22) */
		s1.rs = 0; s1.cs = 0; s1.M = sub1->M;
		copy_matrix(&s1, &tmp1[0][0], n);

		s1.rs = 0; s1.cs = n/2; s1.M = sub2->M;
		s2.rs = n/2; s2.cs = n/2; s2.M = sub2->M;
		sub_quarter(&s1, &s2, &tmp2[0][0], n/2);

		s1.rs = 0; s1.cs = 0; s1.M = &tmp1[0][0];
		s2.rs = 0; s2.cs = 0; s2.M = &tmp2[0][0];
		matrix_mult_strassen(&s1, &s2, &M3[0][0], n/2);

		/* M4 calculation: M4 := A22(B21 - B11) */
		s1.rs = n/2; s1.cs = n/2; s1.M = sub1->M;
		copy_matrix(&s1, &tmp1[0][0], n);

		s1.rs = n/2; s1.cs = 0; s1.M = sub2->M;
		s2.rs = 0; s2.cs = 0; s2.M = sub2->M;
		sub_quarter(&s1, &s2, &tmp2[0][0], n/2);

		s1.rs = 0; s1.cs = 0; s1.M = &tmp1[0][0];
		s2.rs = 0; s2.cs = 0; s2.M = &tmp2[0][0];
		matrix_mult_strassen(&s1, &s2, &M4[0][0], n/2);

		/* M5 calculation: M5 := (A11 + A12)B22 */
		s1.rs = 0; s1.cs = 0; s1.M = sub1->M;
		s2.rs = 0; s2.cs = n/2; s2.M = sub1->M;
		sum_quarter(&s1, &s2, &tmp1[0][0], n/2);

		s1.rs = n/2; s1.cs = n/2; s1.M = sub2->M;
		copy_matrix(&s1, &tmp2[0][0], n);

		s1.rs = 0; s1.cs = 0; s1.M = &tmp1[0][0];
		s2.rs = 0; s2.cs = 0; s2.M = &tmp2[0][0];
		matrix_mult_strassen(&s1, &s2, &M5[0][0], n/2);

		/* M6 calculation: M6 := (A21 - A11)(B11 + B12) */
		s1.rs = n/2; s1.cs = 0; s1.M = sub1->M;
		s2.rs = 0; s2.cs = 0; s2.M = sub1->M;
		sub_quarter(&s1, &s2, &tmp1[0][0], n/2);

		s1.rs = 0; s1.cs = 0; s1.M = sub2->M;
		s2.rs = 0; s2.cs = n/2; s2.M = sub2->M;
		sum_quarter(&s1, &s2, &tmp2[0][0], n/2);

		s1.rs = 0; s1.cs = 0; s1.M = &tmp1[0][0];
		s2.rs = 0; s2.cs = 0; s2.M = &tmp2[0][0];
		matrix_mult_strassen(&s1, &s2, &M6[0][0], n/2);

		/* M7 calculation: M7 := (A12 - A22)(B21 + B22) */
		s1.rs = 0; s1.cs = n/2; s1.M = sub1->M;
		s2.rs = n/2; s2.cs = n/2; s2.M = sub1->M;
		sub_quarter(&s1, &s2, &tmp1[0][0], n/2);

		s1.rs = n/2; s1.cs = 0; s1.M = sub2->M;
		s2.rs = n/2; s2.cs = n/2; s2.M = sub2->M;
		sum_quarter(&s1, &s2, &tmp2[0][0], n/2);

		s1.rs = 0; s1.cs = 0; s1.M = &tmp1[0][0];
		s2.rs = 0; s2.cs = 0; s2.M = &tmp2[0][0];
		matrix_mult_strassen(&s1, &s2, &M7[0][0], n/2);

		/* C11 = M1 + M4 - M5 + M7 */
		sum(&M1[0][0], &M4[0][0], &tmp1[0][0], n/2);
		matrix_sub(&M7[0][0], &M5[0][0], &tmp2[0][0], n/2);
		sum(&tmp1[0][0], &tmp2[0][0], &C11[0][0], n/2);

		/* C12 = M3 + M5 */
		sum(&M3[0][0], &M5[0][0], &C12[0][0], n/2);

		/* C21 = M2 + M4 */
		sum(&M2[0][0], &M4[0][0], &C21[0][0], n/2);

		/* C22 = M1 - M2 + M3 + M6 */
		matrix_sub(&M1[0][0], &M2[0][0], &tmp1[0][0], n/2);

		sum(&M3[0][0], &M6[0][0], &tmp2[0][0], n/2);
		sum(&tmp1[0][0], &tmp2[0][0], &C22[0][0], n/2);

		merge(&C11[0][0], &C12[0][0], &C21[0][0], &C22[0][0], pMatrix, n);
	}
}

