// Author: Brian Bao Nguyen

#include <iostream>
#include <bits/stdc++.h>
#include <math.h>
#include <omp.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>


/*

1) (10 points) Write a serial implementation to calculate the standard-deviation of an array 
of randomly generated floating point values (or doubles!). Then thread it using std::threads 
 well as OpenMP (yes, two implementations!). Compare the performance results (over both N and P).
 Estimate the overhead costs for both parallelization methods.

Benchmark times for P = 10 and N = 1,000,000,000:

OpenMP: 0.20s

Std::threads: 0.20s

*/

void get_walltime(double* wcTime) {
	/* Calculate the execution wall-clock time. */

     struct timeval tp;
     gettimeofday(&tp, NULL);
     *wcTime = (double)(tp.tv_sec + tp.tv_usec/1000000.0);

}

double calculate_std(double *num_array, int N) {
	/* Calculate the standard deviation of an array of floats */

	double sum = 0;
	double new_sum = 0;
	double std, diff, mean;

	#pragma omp parallel for reduction (+: sum)
	for (int i = 0; i < N; i++) {
    	sum += num_array[i];
	}

	mean = sum / N;

	#pragma omp parallel for reduction (+: new_sum)
	for (int i = 0; i < N; i++) {
		diff = num_array[i] - mean;
    	new_sum += diff * diff;
	}

	std = sqrt(new_sum / N);

	return std;
}


int main() {

	double S, E;
	double *A;
	int N = 1000000000;
	int P = 10;

	omp_set_num_threads(P);

	A = (double *)malloc(sizeof(double)*N);

	for(int i = 0; i < N; i++) {
		A[i] = random();
	}

	get_walltime(&S);
	double std = calculate_std(A, N);
	get_walltime(&E);

	printf("%f\n", E - S);
	printf("%f\n", std);

	return 0;
}