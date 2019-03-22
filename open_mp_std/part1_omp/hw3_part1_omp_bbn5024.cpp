// Author: Brian Bao Nguyen

#include <iostream>
#include <bits/stdc++.h>
#include <math.h>
#include <omp.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>


#define N 1000000000
#define P 10

void get_walltime(double* wcTime) {
	/* Calculate the execution wall-clock time. */

     struct timeval tp;
     gettimeofday(&tp, NULL);
     *wcTime = (double)(tp.tv_sec + tp.tv_usec/1000000.0);

}

double calculate_std(double *num_array) {
	/* Calculate the standard deviation of an array of floats. */

	double sum = 0;
	double squared_sum = 0;
	double std, mean, variance;

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

	return std;
}


int main() {
	/* Calculate the standard deviation of an array of size N parallely. */

	// initialize variables
	double S, E;
	double *A;

	// set the number of threads
	omp_set_num_threads(P);

	// dynamic allocation for array
	A = (double *)malloc(sizeof(double)*N);

	// set values of arrays to random doubles
	for(int i = 0; i < N; i++) {
		A[i] = random();
	}

	// calculate std
	get_walltime(&S);
	double std = calculate_std(A);
	get_walltime(&E);

	// free memory for array
	free(A);

	// show standard deviation and execution time
	printf("%f\n", E - S);
	printf("%f\n", std);

	return 0;
}
