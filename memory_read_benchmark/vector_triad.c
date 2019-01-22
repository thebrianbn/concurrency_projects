#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>

void get_walltime(double* wcTime) {

    struct timeval tp;
    gettimeofday(&tp, NULL);
    *wcTime = (double)(tp.tv_sec + tp.tv_usec/1000000.0);

}

void dummy(double a[], double b[], double c[], double d[]) {
	/* Function used to prevent the compiler from doing an obvious
	optimization */

	printf("Dummy function called.");

}

double vector_triad(int N, int R) {

	double *A, *B, *C, *D;
	double S, E, MFLOPS;

	// Dynamic allocation for arrays
	A = (double *)malloc(sizeof(double)*N);
	B = (double *)malloc(sizeof(double)*N);
	C = (double *)malloc(sizeof(double)*N);
	D = (double *)malloc(sizeof(double)*N);

	// Intitialize arrays
	for (int i = 1; i < N; i++) {
		A[i] = 0.0; B[i] = 1.0;
		C[i] = 2.0; D[i] = 3.0;
 	}

 	// Get starting timestamp
 	get_walltime(&S);

 	for (int j = 1; j < R; j++) {
 		for (int i = 1; i < N; i++) {
 			A[i] = B[i] + C[i] * D[i];
 		}
 		if (A[2] < 0) {
 			dummy(A, B, C, D);
 		}
 	}

 	// Get ending timestamp
 	get_walltime(&E);

 	MFLOPS = R * N * 2.0 / ((E - S) * 1000000.0);
 	
 	return MFLOPS;

}

void main() {
	/* Run vector triad program with different values for R and N */

	int N[18] = {10, 100, 250, 500, 1000, 5000, 10000, 25000, 50000, 100000, 250000, 500000, 750000, 1000000, 2500000, 5000000, 7500000, 10000000};
	int R[18] = {10000, 5000, 2500, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 100, 100, 100, 100, 100, 50};

	for (int i = 0; i < 18; i++) {
		double MFLOP = vector_triad(N[i], R[i]);
		printf("%f\n", MFLOP);
	}

}

