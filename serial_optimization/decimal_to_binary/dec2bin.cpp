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
void dec2bin(long int decimalValue, int *binaryArray, int numBits, int lookupTable[][32]) {
	/* Convert decimal-base numbers to binary-base numbers. */

	int bp = numBits;
	long int N = decimalValue;

	// If N is within the range of the lookup table, use it
	if (N <= 2130000000) {
		binaryArray = lookupTable[N];
	}
	// Otherwise, calculate conversion
	else {
		while(N > 0) {
			// Store remainder in binary and update N
			binaryArray[bp] = N % 2;
			N = N / 2;

			// Move to next position
			bp -= 1;
		};
	}
	
}


void printArray(int *array, int numBits) {
	/* Print the contents of an array. */
	for (int ii = 0; ii < numBits; ii++)
		printf("%i ", array[numBits - ii - 1]);
	printf("\n");
}


int main(int argc, char **argv)
{
	int binArray[32];
	long int decimal;
	int lookupTable[2130000000][32];

	// Initialize lookup table
	for (int i = 0; i < 2130000000; i++) {
		dec2bin(i, lookupTable[i], 32, lookupTable);
	}

	double d_S1, d_E1;

	int i_R = 10000000;

	// Random seed
	srandom(716);
	
	// Start execution walltime
	get_walltime(&d_S1);

	// For i_R iterations, convert from dec to bin for a random long int
	for (uint32_t ii = 0; ii < i_R; ii++)
	{
		decimal = random();
		dec2bin(decimal, binArray, 32, lookupTable);
	}
	// End execution walltime
	get_walltime(&d_E1);

	#ifdef DEBUG
	decimal = random();
	dec2bin(decimal, binArray, 32);
	#endif

	printf("Time dec2bin: %f\n", d_E1-d_S1);

	return 0;
}
