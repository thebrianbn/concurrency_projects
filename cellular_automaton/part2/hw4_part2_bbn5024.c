// Author: Brian Nguyen

#include <mpi.h>
#include <math.h> 
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

    int taskid, numtasks, p;
    unsigned long m, k;
    MPI_Status status;

    if (argc != 3) {
        printf("%s <m> <k>\n", argv[0]);
        printf("Program for parallel Game of Life\n");
        printf("with 1D grid partitioning\n");
        printf("<m>: grid dimension (an mxm grid is created)\n");
        printf("<k>: number of time steps\n");
        printf("(initial pattern specified inside code)\n");
        exit(1);
    }

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);

    /* receive grid dimensions, generation amount, and worker amount from
    command line */
    m = atol(argv[1]);
    k = atol(argv[2]);
    p = numtasks - 1;

    int rows_per_worker = floor(m / numtasks);

    int *grid_current;
    int *grid_next;

    int i, j, t;

    return 0;
}
