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

	double all_sum, squared_sum, std;
	double S, E, mean, variance;

	// testing variables
	//int N[9] = {10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000, 1000000000};
	//int P[9] = {10, 10, 10, 10, 10, 10, 10, 10, 10};

	int N[10] = {10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000};
	int P[10] = {10, 10, 10, 10, 10, 10, 10, 10, 10, 10};

	for (int i = 0; i < 10; i++) {

		all_sum = 0, squared_sum = 0, std = 0;
		std::vector<double> A(N[i]);

		// set values of arrays to random doubles
		for(int j = 0; j < N[i]; j++) {
			A[j] = random();
		}

		// initialize structs for threads
		struct MYPARAM *p_params = new struct MYPARAM[P[i]];

		// set start and end array indices
		for (int j = 0; j < P[i]; j++) {
			p_params[j].i_start = j * (N[i]/P[i]);
			p_params[j].i_stop = (j + 1) * (N[i]/P[i]);
			p_params[j].d_sum = 0.0;
			p_params[j].d_squared_sum = 0.0;
			p_params[j].b_complete = false;
		}

		// initialize standard threads
		std::thread threads[P[i]];

		get_walltime(&S);

		// perform initial calls for sums
		for (int j = 0; j < P[i]; j++) {
			threads[j] = std::thread(calculate_std, &p_params[j], std::ref(A));
		}

		// join all threads for sums
		for (int j = 0; j < P[i]; j++) {
			threads[j].join();
		}

		// sum all thread results
		for (int j = 0; j < P[i]; j++) {
			all_sum += p_params[j].d_sum;
			squared_sum += p_params[j].d_squared_sum;
		}

		// calculate mean and variance of initial array
		mean = all_sum / N[i];
		variance = (squared_sum / N[i]) - (mean * mean);

		// take the square root of the variance to get standard deviation
		std = sqrt(variance);

		get_walltime(&E);

		// show standard deviation and execution time
		printf("N=%d, P=%d, TIME=%f\n", N[i], P[i], E - S);
		printf("\tStandard Deviation: %f\n", std);

		delete[] p_params;
	}
	
	return 0;
}