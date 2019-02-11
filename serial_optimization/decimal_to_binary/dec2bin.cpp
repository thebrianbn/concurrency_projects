// Author: Brian Bao Nguyen

#include <iostream>
#include <bits/stdc++.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

void get_walltime(double* wcTime) {

     struct timeval tp;

     gettimeofday(&tp, NULL);

     *wcTime = (double)(tp.tv_sec + tp.tv_usec/1000000.0);

}
/****************************************************************
 *	dec2bin
 *	
 *	Populates a structure with a binary array representing the value.
 *
 *	ex: input = 0xDCBA
 *	output A[] = 0 1 0 1   1 1 0 1   0 0 1 1   1 0 1 1
 *	         A         B         C         D
 *
 */
void dec2bin(long int decimalValue, int *binaryArray, int numBits, int lookupTable[][32])
{
	int bp = numBits;
	long int N = decimalValue;
	if (N < 5000000000) {
		binaryArray = lookupTable[N];
	}
	else {
		while(N > 0) {
			binaryArray[bp] = N % 2;
			bp -= 1;
			N = N / 2;
		};
	}
	
}

// print the contents of an array
void printArray(int *array, int numBits)
{
	for (int ii = 0; ii < numBits; ii++)
		printf("%i ", array[numBits - ii - 1]);
	printf("\n");
}

int main(int argc, char **argv)
{
	int binArray[32];
	long int decimal;

	int lookupTable[5000000000][32];

	for (int i = 0; i < 5000000000; i++) {
		dec2bin(i, lookupTable[i], 32, lookupTable);
	}

	double d_S1, d_E1;

	int i_R = 10000000;

	srandom(716);
	
	get_walltime(&d_S1);
	for (uint32_t ii = 0; ii < i_R; ii++)
	{
		decimal = random();
		dec2bin(decimal, binArray, 32, lookupTable);
	}
	get_walltime(&d_E1);

	/*
	#ifdef DEBUG
	decimal = random();
	dec2bin(decimal, binArray, 32);
	*/

	printf("Time dec2bin: %f\n", d_E1-d_S1);

	return 0;
}
