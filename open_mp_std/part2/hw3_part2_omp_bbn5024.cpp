// Author: Brian Bao Nguyen

#include <iostream>
#include <bits/stdc++.h>
#include <math.h>
#include <omp.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>


#define sigma 2

void get_walltime(double* wcTime) {
	/* Calculate the execution wall-clock time. */

     struct timeval tp;
     gettimeofday(&tp, NULL);
     *wcTime = (double)(tp.tv_sec + tp.tv_usec/1000000.0);

}

void over_sigma(std::vector<std::tuple<int, double>> &targets_array, double *num_array,
	double &std, double &mean, int N) {
	/* Return a vector of tuples containing indices and values in an array
	that is over a specified sigma value */

	double threshold_right = mean + (sigma * std);
	double threshold_left = mean - (sigma * std);

	#pragma omp parallel for
	for (int i = 0; i < N; i++) {
		if (num_array[i] <= threshold_left || num_array[i] >= threshold_right) {
			targets_array.push_back(std::make_tuple(i, num_array[i]));
		}
	}
}

double calculate_std(double *num_array, double &std_val, double &mean_val, int N) {
	/* Calculate the standard deviation of an array of floats. */

	double sum = 0;
	double squared_sum = 0;
	double variance;

	// parallel reduction for initial sum
	#pragma omp parallel for reduction (+: sum, squared_sum)
	for (int i = 0; i < N; i++) {
    	sum += num_array[i];
    	squared_sum += num_array[i] * num_array[i];
	}

	mean_val = sum / N;
	variance = (squared_sum / N) - (mean_val * mean_val);

	// take the square root of the variance to get standard deviation
	std_val = sqrt(variance);
}

using namespace std;
int main() {
	/* Calculate the standard deviation of an array of size N parallely. */

	// initialize variables
	double S, E, std, mean;
	double *A;
	std::vector<std::tuple<int, double>> targets;

	// testing variables
	int N[9] = {10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000, 1000000000};
	int P[9] = {10, 10, 10, 10, 10, 10, 10, 10, 10};

	for (int i = 0; i < 9; i++){

		std = 0, mean = 0;

		// set the number of threads
		omp_set_num_threads(P[i]);

		// dynamic allocation for array
		A = (double *)malloc(sizeof(double)*N[i]);

		// set values of arrays to random doubles
		for(int j = 0; j < N[i]; j++) {
			A[j] = random();
		}

		// calculate std
		calculate_std(A, std, mean, N[i]);
		
		// Return a vector of tuples of indices and values greater than sigma
		get_walltime(&S);
		over_sigma(targets, A, std, mean, N[i]);
		get_walltime(&E);

		// free memory for array
		free(A);

		// show execution time of greater_sigma
		printf("%f\n", E - S);
	}

	return 0;
}
