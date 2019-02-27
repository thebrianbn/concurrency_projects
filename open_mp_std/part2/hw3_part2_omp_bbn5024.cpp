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
#define sigma 1

void get_walltime(double* wcTime) {
	/* Calculate the execution wall-clock time. */

     struct timeval tp;
     gettimeofday(&tp, NULL);
     *wcTime = (double)(tp.tv_sec + tp.tv_usec/1000000.0);

}

void within_sigma(std::vector<std::tuple<int, double>> &targets_array, double *num_array,
	double &std, double &mean) {
	/* Return a vector of tuples containing indices and values in an array
	that is over a specified sigma value */

	double threshold_left = mean + (sigma * std);
	double threshold_right = mean - (sigma * std);

	#pragma omp parallel for
	for (int i = 0; i < N; i++) {
		if (num_array[i] >= threshold_left && num_array[i] <= threshold_right) {
			targets_array.push_back(std::make_tuple(i, num_array[i]));
		}
	}
}

void calculate_std(double *num_array, double &std_val, double &mean_val) {
	/* Calculate the standard deviation of an array of floats. */

	double sum = 0;
	double new_sum = 0;
	double diff;

	// parallel reduction for initial sum
	#pragma omp parallel for reduction (+: sum)
	for (int i = 0; i < N; i += 4) {
    	sum += num_array[i] + num_array[i+1] + num_array[i+2] + num_array[i+3];
	}

	// calculate mean of initial array values
	mean_val = sum / N;

	// parallel reduction for sum of values subtracted by mean
	#pragma omp parallel for reduction (+: new_sum)
	for (int i = 0; i < N; i++) {
		diff = num_array[i] - mean_val;
    	new_sum += diff * diff;
	}

	// take the square root of the new mean to get standard deviation
	std_val = sqrt(new_sum / N);
}

using namespace std;
int main() {
	/* Calculate the standard deviation of an array of size N parallely. */

	// initialize variables
	double S, E, std, mean;
	double *A;
	std::vector<std::tuple<int, double>> targets;

	// set the number of threads
	omp_set_num_threads(P);

	// dynamic allocation for array
	A = (double *)malloc(sizeof(double)*N);

	// set values of arrays to random doubles
	for(int i = 0; i < N; i++) {
		A[i] = random();
	}

	// calculate std
	calculate_std(A, std, mean);
	
	// Return a vector of tuples of indices and values greater than sigma
	get_walltime(&S);
	within_sigma(targets, A, std, mean);
	get_walltime(&E);

	// free memory for array
	free(A);

	// show execution time of greater_sigma
	printf("%f\n", E - S);

	return 0;
}
