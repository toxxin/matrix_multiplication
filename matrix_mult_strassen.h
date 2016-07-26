#ifndef MATRIX_MULT_STRASSEN_H
#define MATRIX_MULT_STRASSEN_H

typedef struct _subm {
	int rs;
	int cs;
	int *M;
} subm;

void print_matrix(int *m, int n);

void matrix_mult_naive(int *aMatrix, int *bMatrix, int *pMatrix, size_t n);
void matrix_mult_strassen(subm *sub1, subm *sub2, int *pMatrix, size_t n);

#endif // MATRIX_MULT_STRASSEN_H
