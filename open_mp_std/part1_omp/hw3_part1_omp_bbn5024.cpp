// Author: Brian Bao Nguyen

#include <iostream>
#include <bits/stdc++.h>
#include <math.h>
#include <omp.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>


void get_walltime(double* wcTime) {
	/* Calculate the execution wall-clock time. */

     struct timeval tp;
     gettimeofday(&tp, NULL);
     *wcTime = (double)(tp.tv_sec + tp.tv_usec/1000000.0);

}

double calculate_std(double *num_array, int N, int P) {
	/* Calculate the standard deviation of an array of floats. */

	double sum = 0;
	double squared_sum = 0;
	double std, mean, variance, S, E;

	get_walltime(&S);

	// parallel reduction for initial sum
	#pragma omp parallel for reduction (+: sum, squared_sum)
	for (int i = 0; i < N; i++) {
    	sum += num_array[i];
    	squared_sum += num_array[i] * num_array[i];
	}

	mean = sum / N;
	variance = (squared_sum / N) - (mean * mean);

	// take the square root of the variance to get standard deviation
	std = sqrt(variance);

	get_walltime(&E);

	printf("N=%d, P=%d, TIME=%f\n", N, P, E - S);
	printf("\tStandard Deviation: %f\n", std);

	return std;
}


int main() {
	/* Calculate the standard deviation of an array of size N parallely. */

	// initialize variables
	double S, E;
	double *A;

	// testing variables
	int N[9] = {10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000, 1000000000};
	int P[9] = {10, 10, 10, 10, 10, 10, 10, 10, 10};

	

	// test calculate_std functionality and show runtimes
	for (int i = 0; i < 9; i++){

		// set the number of threads to P
		omp_set_num_threads(P[i]);

		// dynamic allocation for array
		A = (double *)malloc(sizeof(double)*N[i]);

		// set values of arrays to random doubles
		for (int j = 0; j < N[i]; j++) {
			A[j] = random();
		}

		double std = calculate_std(A, N[i], P[i]);

		// free memory for array
		free(A);
	}
	
	

	return 0;
}
