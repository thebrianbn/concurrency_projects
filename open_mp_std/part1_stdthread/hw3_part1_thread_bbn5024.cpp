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
	double d_result;
    bool b_complete;
};

void get_walltime(double* wcTime) {
	/* Calculate the execution wall-clock time. */

     struct timeval tp;
     gettimeofday(&tp, NULL);
     *wcTime = (double)(tp.tv_sec + tp.tv_usec/1000000.0);

}

double calculate_std(double sum) {
	/* Calculate the standard deviation of an array of floats. */

	// take the square root of the new mean to get standard deviation
	double std = sqrt(sum / N);

	return std;
}

void calculate_sum(struct MYPARAM *p_params, std::vector<double> &num_array) {

	double sum = 0;

	for (int i = p_params->i_start; i < p_params->i_stop; i++) {
    	sum += num_array[i];
	}

	p_params->d_result = sum;
	p_params->b_complete = true;
}

void calculate_ssd(struct MYPARAM *p_params, std::vector<double> &num_array, double mean) {

	double sum = 0;
	double diff;

	for (int i = p_params->i_start; i < p_params->i_stop; i++) {
		diff = num_array[i] - mean;
		sum += diff * diff;
	}

	p_params->d_result = sum;
	p_params->b_complete = true;
}


int main() {
	/* Calculate the standard deviation of an array of size N parallely. */

	// initialize variables
	double all_sum = 0, ssd_sum = 0, initial_mean = 0, std = 0;
	double S, E;

	std::vector<double> A(N);

	// set values of arrays to random doubles
	for(int i = 0; i < N; i++) {
		A[i] = random();
	}

	get_walltime(&S);

	// initialize structs for threads
	struct MYPARAM *p_params = new struct MYPARAM[P];

	// set start and end array indices
	for (int i = 0; i < P; i++) {
		p_params[i].i_start = i * (N/P);
		p_params[i].i_stop = (i + 1) * (N/P);
		p_params[i].d_result = 0.0;	
		p_params[i].b_complete = false;
	}

	// initialize standard threads
	std::thread threads[P];

	// perform initial calls for initial mean
	for (int i = 0; i < P; i++) {
		threads[i] = std::thread(calculate_sum, std::ref(&p_params[i]), A);
	}

	// check for thread completion
	bool b_done = false;
	while (!b_done)
	{
		b_done = true;
		for (int i = 0; i < P; i++)
		{
			if (!p_params[i].b_complete)	
				b_done = false;
		}
	
		std::this_thread::sleep_for(std::chrono::milliseconds(10));
	}

	// join all threads for initial mean
	for (int i = 0; i < P; i++) {
		threads[i].join();
	}

	// sum all thread results, reset for next firing of threads
	for (int i = 0; i < P; i++) {
		all_sum += p_params[i].d_result;
		p_params[i].d_result = 0.0;
		p_params[i].b_complete = false;
	}

	// calculate mean of initial array
	initial_mean = all_sum / N;

	// re-fire threads for ssd calculation
	for (int i = 0; i < P; i++) {
		threads[i] = std::thread(calculate_ssd, std::ref(&p_params[i]), A, initial_mean);
	}

	// check for thread completion
	b_done = false;
	while (!b_done)
	{
		b_done = true;
		for (int i = 0; i < P; i++)
		{
			if (!p_params[i].b_complete)	
				b_done = false;
		}
	
		std::this_thread::sleep_for(std::chrono::milliseconds(10));
	}

	// join all threads for ssd
	for (int i = 0; i < P; i++) {
		threads[i].join();
	}

	// sum all thread results for ssd
	for (int i = 0; i < P; i++) {
		ssd_sum += p_params[i].d_result;
	}

	// perform final calculate for std, doesn't need parallelization
	std = calculate_std(ssd_sum);

	get_walltime(&E);

	// show standard deviation and execution time
	printf("%f\n", E - S);
	printf("%f\n", std);

	delete[] p_params;
	return 0;
}