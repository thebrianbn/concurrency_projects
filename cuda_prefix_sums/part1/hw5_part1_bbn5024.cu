// Author: Brian Nguyen

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#define N 64 // number of array elements
#define B 4  // number of elements in a block

__global__ void scan(float *g_odata, float *g_idata, int n);
__global__ void prescan(float *g_odata, float *g_idata, int n, float *g_sums);
__global__ void uniform_add(float *o_array, float *sum_array);
void scanCPU(float *f_out, float *f_in, int i_n);

double myDiffTime(struct timeval &start, struct timeval &end) {
	/* Calculate the time difference. */

	double d_start, d_end;
	d_start = (double)(start.tv_sec + start.tv_usec/1000000.0);
	d_end = (double)(end.tv_sec + end.tv_usec/1000000.0);
	return (d_end - d_start);
} 

int main() {
	/* Compare results between serial and parallel versions of the
	prefix-sums algorithm. */

	int grid_size = floor(N / B);
	int thread_size = B / 2;
	int grid_size2 = grid_size / B;

	// arrays to be used for initial, cpu-result, and gpu-result arrays
	// respectively.
	float a[N], c[N], g[N], sums[grid_size];
	timeval start, end;

	// temporary pointer arrays for computation
	float *dev_a, *dev_g, *dev_sums;
	int size = N * sizeof(float);
	int size_sums = grid_size * sizeof(float);
	int size_sums2 = grid_size2 * sizeof(float);

	double d_gpuTime, d_cpuTime;

	// initialize matrix a with random floats between 0 and 1000
	for (int i = 1; i <= N; i++) {
		a[i-1] = i;
	}

	// CPU version (serial) of prefix-sum
	gettimeofday(&start, NULL);
	scanCPU(c, a, N);
	gettimeofday(&end, NULL);
	d_cpuTime = myDiffTime(start, end);

	// START OF FIRST PRE-SCAN RUN

	// initialize a and b matrices here for CUDA
	cudaMalloc((void **) &dev_a, size);
	cudaMalloc((void **) &dev_g, size);
	cudaMalloc((void **) &dev_sums, size_sums);

	// GPU version (CUDA) of prefix-sum
	gettimeofday(&start, NULL);
	cudaMemcpy(dev_a, a, size, cudaMemcpyHostToDevice);

	// work-efficient scan for SUMS array
	prescan<<<grid_size, thread_size, B*sizeof(float)>>>(dev_g, dev_a, N, dev_sums);
	cudaDeviceSynchronize();
	cudaMemcpy(g, dev_g, size, cudaMemcpyDeviceToHost);
	cudaMemcpy(sums, dev_sums, size_sums, cudaMemcpyDeviceToHost);
	gettimeofday(&end, NULL);
	
	cudaFree(dev_a); cudaFree(dev_g); cudaFree(dev_sums);

	for (int j = 0; j < grid_size; j++) {
		printf("sums[%i] = %0.3f\n", j, sums[j]);
	}

	// START OF SECOND PRE-SCAN RUN

	float inc[grid_size], sums_inc[grid_size2], inc_final[grid_size2];
	float *dev_inc, *dev_sums_inc, *dev_inc_final, *dev_sums_input;

	cudaMalloc((void **) &dev_sums_input, size_sums);
	cudaMalloc((void **) &dev_inc, size_sums);
	cudaMalloc((void **) &dev_sums_inc, size_sums2);

	cudaMemcpy(dev_sums_input, sums, size_sums, cudaMemcpyHostToDevice);

	prescan<<<grid_size2, thread_size, B*sizeof(float)>>>(dev_inc, dev_sums_input, grid_size2, dev_sums_inc);
	cudaDeviceSynchronize();
	cudaMemcpy(inc, dev_inc, size_sums, cudaMemcpyDeviceToHost);
	cudaMemcpy(sums_inc, dev_sums_inc, size_sums2, cudaMemcpyDeviceToHost);

	cudaFree(dev_inc); cudaFree(dev_sums_inc); cudaFree(dev_sums_input);

	
	for (int j = 0; j < grid_size2; j++) {
		printf("inc[%i] = %0.3f\n", j, sums_inc[j]);
	}
	for (int j = 0; j < grid_size; j++) {
		printf("inc[%i] = %0.3f\n", j, inc[j]);
	}

	scanCPU(inc_final, sums_inc, size_sums2);

	// START OF UPDATING SUMS

	for (int j = 0; j < grid_size2; j++) {
		printf("pscan2[%i] = %0.3f\n", j, inc_final[j]);
	}

	float g2[grid_size];
	float *dev_g2;

	cudaMalloc((void **) &dev_g2, size);
	cudaMalloc((void **) &dev_inc_final, size_sums);

	cudaMemcpy(dev_inc_final, inc_final, size_sums, cudaMemcpyHostToDevice);
	cudaMemcpy(dev_g2, inc, size_sums, cudaMemcpyHostToDevice);

	uniform_add<<<grid_size2, thread_size, B*sizeof(float)>>>(dev_g2, dev_inc_final);
	cudaDeviceSynchronize();

	cudaMemcpy(g2, dev_g2, size, cudaMemcpyDeviceToHost);

	for (int j = 0; j < grid_size; j++) {
		printf("g2[%i] = %0.3f\n", j, g2[j]);
	}

	// START OF FINAL UPDATE TO FIRST PREFIX SCAN

	float g3[N], first_add[grid_size];
	float *dev_g3, *dev_first_add;

	cudaMalloc((void **) &dev_g3, size);
	cudaMalloc((void **) &dev_first_add, size_sums);

	cudaMemcpy(dev_first_add, g2, size_sums, cudaMemcpyHostToDevice);
	cudaMemcpy(dev_g3, g, size, cudaMemcpyHostToDevice);

	uniform_add<<<grid_size, thread_size, B*sizeof(float)>>>(dev_g3, dev_first_add);
	cudaDeviceSynchronize();

	cudaMemcpy(g3, dev_g3, size, cudaMemcpyDeviceToHost);

	// display results of the prefix-sum
	for (int i = 0; i < N; i++) {
		printf("c[%i] = %0.3f, g3[%i] = %0.3f\n", i, c[i], i, g3[i]);
		//if (c[i] != g[i])
		//{
		//	printf("Results do not match! c[%i]=%f, g[%i]=%f\n", i, c[i], i, g[i]);
		//	break;
		//}
	}
		
	printf("GPU Time for scan size %i: %f\n", N, d_gpuTime);
	printf("CPU Time for scan size %i: %f\n", N, d_cpuTime);

	return 0;
}


__global__ void scan(float *g_odata, float *g_idata, int n) {
	/* CUDA Naive Scan Algorithm (double buffered). */

	extern __shared__ float temp[]; // allocated on invocation
	int thid = threadIdx.x;
	int pout = 0, pin = 1;

	// Load input into shared memory.
	// This is exclusive scan, so shift right by one
	// and set first element to 0
	temp[thid] = (thid > 0) ? g_idata[thid-1] : 0;
	__syncthreads();
	for (int offset = 1; offset < n; offset *= 2) {
		pout = 1 - pout; // swap double buffer indices
		pin = 1 - pout;
		if (thid >= offset)
			temp[pout*n+thid] += temp[pin*n+thid - offset];
		else
			temp[pout*n+thid] = temp[pin*n+thid];

		__syncthreads();
	}
	g_odata[thid] = temp[pout*n+thid]; // write output
}


__global__ void prescan(float *g_odata, float *g_idata, int n, float *g_sums) {
	/* CUDA Work-Efficient Scan Algorithm. */

	extern  __shared__  float temp[]; // allocated on invocation 
	int thid = threadIdx.x;  // thread id of a thread in a block
	int gthid = (blockIdx.x * blockDim.x) + thid; // global thread id of grid
	int offset = 1;

	/*
	// for each thread in a block, put data into shared memory
	if (gthid > n) {
		// handle non-power of two arrays by padding elements in last block
		temp[2*thid] = 0;
		temp[2*thid+1] = 0;
	}
	else {
		// grab data from input array
		temp[2*thid] = g_idata[2*gthid];
		temp[2*thid+1] = g_idata[2*gthid+1];
	}
	*/
	temp[2*thid] = g_idata[2*gthid];
	temp[2*thid+1] = g_idata[2*gthid+1];

    // build sum in place up the tree 
	for (int d = B>>1; d > 0; d >>= 1) { 
        __syncthreads(); 
		if (thid < d) { 
			int ai = offset*(2*thid+1)-1; 
			int bi = offset*(2*thid+2)-1; 
		    	temp[bi] += temp[ai];         
  		}
  		offset *= 2; 
    } 

	if (thid == 0) { 
		g_sums[blockIdx.x] = temp[B - 1];
		temp[B - 1] = 0; 
	}

	// clear the last element 
	// traverse down tree & build scan
	for (int d = 1; d < B; d *= 2) { 
    	offset >>= 1; 
    	__syncthreads(); 
		if (thid < d) { 
			int ai = offset*(2*thid+1)-1; 
			int bi = offset*(2*thid+2)-1; 
			float t = temp[ai]; 
    		temp[ai] = temp[bi]; 
    		temp[bi] += t; 
    	} 
	} 
	__syncthreads(); 
	
	// write results to device memory 
	g_odata[2*gthid] = temp[2*thid]; 
	g_odata[2*gthid+1] = temp[2*thid+1]; 
}

__global__ void uniform_add(float *o_array, float *sum_array) {

	int bid = blockIdx.x;
	int gthid = (bid * blockDim.x) + threadIdx.x; // global thread id of grid

	o_array[2*gthid] = o_array[2*gthid] + sum_array[bid];
	o_array[2*gthid+1] = o_array[2*gthid+1] + sum_array[bid];
}


void scanCPU(float *f_out, float *f_in, int i_n) {
	/* Apply all-prefix sums to an array on the CPu
	without parallelization. */

	f_out[0] = 0;

	/* for each array element, the value is the previous sum
	plus the current array value */
	for (int i = 1; i < i_n; i++)
		f_out[i] = f_out[i-1] + f_in[i-1];

}
