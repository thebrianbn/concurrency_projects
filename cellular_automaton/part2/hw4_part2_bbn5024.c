// Author: Brian Nguyen

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/time.h>
#include <time.h>

static double get_walltime() {
    struct timeval tp;
    gettimeofday(&tp, NULL);
    return ((double) (tp.tv_sec) + 1e-6 * tp.tv_usec);
}

int main(int argc, char **argv) {

    if (argc != 3) {
        printf("%s <m> <k>\n", argv[0]);
        printf("Program for parallel Game of Life\n");
        printf("with 1D grid partitioning\n");
        printf("<m>: grid dimension (an mxm grid is created)\n");
        printf("<k>: number of time steps\n");
        printf("(initial pattern specified inside code)\n");
		exit(1);
    }

    unsigned long m, k;

    m = atol(argv[1]);

    k = atol(argv[2]);

    int *grid_current;
    int *grid_next;
    
    grid_current = (int *) malloc(m * m * sizeof(int));
	if (grid_current == NULL)
	{
		printf("Error allocating memory for grid_current!\n");
		exit(1);
	}

    grid_next = (int *) malloc(m * m * sizeof(int));
    if (grid_next == NULL)
   	{
		printf("Error allocating memory for grid_next!\n");
		exit(1);
	}

    int i, j, t;

    /* static initalization, so that we can verify output */
    /* using very simple initialization right now */
    /* this isn't a good check for parallel debugging */
    for (i=0; i<m; i++) {
        for (j=0; j<m; j++) {
            grid_current[i*m+j] = 0;
            grid_next[i*m+j] = 0;
        }
    }

    /* initializing some cells in the middle */
    grid_current[m*m/2 + m/2 + 0] = 1;
    grid_current[m*m/2 + m/2 + 1] = 1;
    grid_current[m*m/2 + m/2 + 2] = 1;
    grid_current[m*m/2 + m/2 + 3] = 1;

    double d_startTime = 0.0, d_endTime = 0.0;
	d_startTime = get_walltime();

    /* serial code */
    /* considering only internal cells */
    for (t=0; t<k; t++) {
        for (i=1; i<m-1; i++) {
            for (j=1; j<m-1; j++) {
                /* avoiding conditionals inside inner loop */
                int prev_state = grid_current[i*m+j];
                // avoid repeated calculations
                int top_index = (i-1)*m;
                int bottom_index = (i+1)*m;
                int center_index = i*m;
                int left_shift = j-1;
                int right_shift = j+1;

                /* avoiding conditionals inside inner loop */
                int prev_state = grid_current[i*m+j];
                int num_alive  = 
                                grid_current[center_index+left_shift] + 
                                grid_current[center_index+right_shift] + 
                                grid_current[top_index+left_shift] + 
                                grid_current[top_index+j] + 
                                grid_current[top_index+right_shift] + 
                                grid_current[bottom_index+left_shift] + 
                                grid_current[bottom_index+j  ] + 
                                grid_current[bottom_index+right_shift];

                grid_next[i*m+j] = prev_state * ((num_alive == 2) + (num_alive == 3)) + (1 - prev_state) * (num_alive == 3);
            }
        }
        /* swap current and next */
        int *grid_tmp  = grid_next;
        grid_next = grid_current;
        grid_current = grid_tmp;
    }

    d_endTime = get_walltime();

    /* Verify */
    int verify_failed = 0;
    for (i=0; i<m; i++) {
        for (j=0; j<m; j++) {
            /* Add verification code here */
        }
    }

    printf("Time taken: %3.3lf s.\n", d_endTime - d_startTime);
    printf("Performance: %3.3lf billion cell updates/s\n", 
                (1.0*m*m)*k/((d_endTime - d_startTime)*1e9));

    /* free memory */
    free(grid_current); free(grid_next);

    return 0;
}