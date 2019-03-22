// Author: Brian Bao Nguyen

#include <iostream>
#include <bits/stdc++.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <thread>
#include <chrono>
#include <stdio.h>


#define N 1000000000
#define P 10

struct MYPARAM{
	int i_start;
	int i_stop;
	double d_sum;
	double d_squared_sum;
    bool b_complete;
};

void get_walltime(double* wcTime) {
	/* Calculate the execution wall-clock time. */

     struct timeval tp;
     gettimeofday(&tp, NULL);
     *wcTime = (double)(tp.tv_sec + tp.tv_usec/1000000.0);

}

double calculate_std(struct MYPARAM *p_params, std::vector<double> &num_array) {
	/* Calculate the standard deviation of an array of floats. */

	double sum = 0;
	double squared_sum = 0;
	double std;

	for (int i = p_params->i_start; i < p_params->i_stop; i++) {
    	sum += num_array[i];
    	squared_sum += num_array[i] * num_array[i];
	}

	p_params->d_sum = sum;
	p_params->d_squared_sum = squared_sum;
}

int main() {
	/* Calculate the standard deviation of an array of size N parallely. */

	double all_sum = 0, squared_sum = 0, std = 0;
	double S, E, mean, variance;

	std::vector<double> A(N);

	// set values of arrays to random doubles
	for(int i = 0; i < N; i++) {
		A[i] = random();
	}

	// initialize structs for threads
	struct MYPARAM *p_params = new struct MYPARAM[P];

	// set start and end array indices
	for (int i = 0; i < P; i++) {
		p_params[i].i_start = i * (N/P);
		p_params[i].i_stop = (i + 1) * (N/P);
		p_params[i].d_sum = 0.0;
		p_params[i].d_squared_sum = 0.0;
		p_params[i].b_complete = false;
	}

	// initialize standard threads
	std::thread threads[P];

	get_walltime(&S);

	// perform initial calls for sums
	for (int i = 0; i < P; i++) {
		threads[i] = std::thread(calculate_std, &p_params[i], std::ref(A));
	}

	// join all threads for sums
	for (int i = 0; i < P; i++) {
		threads[i].join();
	}

	// sum all thread results
	for (int i = 0; i < P; i++) {
		all_sum += p_params[i].d_sum;
		squared_sum += p_params[i].d_squared_sum;
	}

	// calculate mean and variance of initial array
	mean = all_sum / N;
	variance = (squared_sum / N) - (mean * mean);

	// take the square root of the variance to get standard deviation
	std = sqrt(variance);

	get_walltime(&E);

	// show standard deviation and execution time
	printf("%f\n", E - S);
	printf("%f\n", std);

	delete[] p_params;
	return 0;
}