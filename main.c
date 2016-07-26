#include <stdio.h>
#include <stdlib.h>
#include "matrix_mult_strassen.h"


int main() {
	unsigned int n, i;
	srand(time(NULL));

	n = 16;

	int *aMatrix, *bMatrix, *pMatrix;
	aMatrix = (int *)calloc(n*n, sizeof(int)*n);
	bMatrix = (int *)calloc(n*n, sizeof(int)*n);

	if ((aMatrix != NULL) && (bMatrix != NULL)) {
		for(i = 0; i < n*n; ++i) {
			*(aMatrix + i) = rand() % 11 + (-5);
			*(bMatrix + i) = rand() % 11 + (-5);
		}
	} else {
		fprintf(stderr, "Cannot allocate memory for matrix with size: %d\n", n);
		exit(EXIT_FAILURE);
	}

	pMatrix = (int *)calloc(n*n, sizeof(int)*n);
	if (pMatrix == NULL) {
		fprintf(stderr, "Cannot allocate memory for matrix with size: %d\n", n);
		exit(EXIT_FAILURE);
	}

#if 0
	printf("MatrixA\n");
	//print_matrix(&aMatrix[0][0], n);
	print_matrix(aMatrix, n);
	printf("MatrixB\n");
	//print_matrix(&bMatrix[0][0], n);
	print_matrix(bMatrix, n);
#endif

	/* Matrix multiplication by naive algorithm */
	matrix_mult_naive(aMatrix, bMatrix, pMatrix, n);
	printf("Result with naive algorithm:\n");
	print_matrix(pMatrix, n);

	/* Matrix multiplication by Strassen algorithm */
	subm s1 = { .rs = 0, .cs = 0, .M = aMatrix };
	subm s2 = { .rs = 0, .cs = 0, .M = bMatrix };
	matrix_mult_strassen(&s1, &s2, pMatrix, n);
	printf("Result with Strassen algorithm:\n");
	print_matrix(pMatrix, n);

	free(aMatrix);
	free(bMatrix);
	free(pMatrix);

	return 0;
}
