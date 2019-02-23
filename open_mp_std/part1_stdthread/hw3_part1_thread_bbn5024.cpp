// Author: Brian Bao Nguyen

#include <iostream>
#include <bits/stdc++.h>
#include <math.h>
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

double calculate_std(double *num_array, int N) {
	/* Calculate the standard deviation of an array of floats. */

	// initialize variables
	double sum = 0;
	double new_sum = 0;
	double std, diff, mean;

	// parallel reduction for initial sum
	for (int i = 0; i < N; i++) {
    	sum += num_array[i];
	}

	// calculate mean of initial array values
	mean = sum / N;

	// parallel reduction for sum of values subtracted by mean
	for (int i = 0; i < N; i++) {
		diff = num_array[i] - mean;
    	new_sum += diff * diff;
	}

	// take the square root of the new mean to get standard deviation
	std = sqrt(new_sum / N);

	return std;
}


int main() {
	/* Calculate the standard deviation of an array of size N parallely. */

	// initialize variables
	double S, E;
	double *A;
	int N = 1000000000;
	int P = 10;

	// dynamic allocation for array
	A = (double *)malloc(sizeof(double)*N);

	// set values of arrays to random doubles
	for(int i = 0; i < N; i++) {
		A[i] = random();
	}

	// calculate std
	get_walltime(&S);
	double std = calculate_std(A, N);
	get_walltime(&E);

	// free memory for array
	free(A);

	// show standard deviation and execution time
	printf("%f\n", E - S);
	printf("%f\n", std);

	return 0;
}