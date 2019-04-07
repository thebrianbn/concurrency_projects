#include <mpi.h>
#include <math.h> 
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/time.h>
#include <time.h>

MPI_Status status;

static double get_walltime() {
    struct timeval tp;
    gettimeofday(&tp, NULL);
    return ((double) (tp.tv_sec) + 1e-6 * tp.tv_usec);
}

int main(int argc, char **argv) {

    int top_index, bottom_index, left_index, right_index;
    int taskid, numtasks, start_row, end_row;
    unsigned long m, k, p;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);

    if (argc != 3) {
        printf("%s <m> <k>\n", argv[0]);
        printf("Program for parallel Game of Life\n");
        printf("with 1D grid partitioning\n");
        printf("<m>: grid dimension (an mxm grid is created)\n");
        printf("<k>: number of time steps\n");
        printf("(initial pattern specified inside code)\n");
		exit(1);
    }

    /* receive grid dimensions, generation amount, and worker amount from
    command line */
    m = atol(argv[1]);
    k = atol(argv[2]);
    p = numtasks - 1;

    int rows_per_worker = floor(m / p);

    int *grid_current;
    int *grid_next;

    int i, j, t;

    printf("%d", taskid);
    return 0;

    // if last worker, handle any extra rows
    if (taskid == p) {

        grid_current = (int *) malloc((rows_per_worker + (m % p)) * m * sizeof(int));
        if (grid_current == NULL) {
            printf("Error allocating memory for grid_current!\n");
            exit(1);
        }

        grid_next = (int *) malloc((rows_per_worker + (m % p)) * m * sizeof(int));
        if (grid_next == NULL) {
            printf("Error allocating memory for grid_next!\n");
            exit(1);
        }

        for (i=0; i<rows_per_worker + (m % p); i++) {
            for (j=0; j<m; j++) {
                grid_current[i*m+j] = 0;
                grid_next[i*m+j] = 0;
            }
        }

    }
    else {

        grid_current = (int *) malloc(rows_per_worker * m * sizeof(int));
        if (grid_current == NULL) {
            printf("Error allocating memory for grid_current!\n");
            exit(1);
        }

        grid_next = (int *) malloc(rows_per_worker * m * sizeof(int));
        if (grid_next == NULL) {
            printf("Error allocating memory for grid_next!\n");
            exit(1);
        }

        for (i=0; i<rows_per_worker; i++) {
            for (j=0; j<m; j++) {
                grid_current[i*m+j] = 0;
                grid_next[i*m+j] = 0;
            }
        }
    }

    // make some cells alive
    if (taskid == floor(p/2)) {
        grid_current[m*rows_per_worker/2 + m/2 + 0] = 1;
        grid_current[m*rows_per_worker/2 + m/2 + 1] = 1;
        grid_current[m*rows_per_worker/2 + m/2 + 2] = 1;
        grid_current[m*rows_per_worker/2 + m/2 + 3] = 1;
    }

    int *send_top_row;
    int *send_bottom_row;
    int *recv_top_row;
    int *recv_bottom_row;
    int num_alive;

    // arrays used for sendrecv
    send_top_row = (int *) malloc(m * sizeof(int));
    send_bottom_row = (int *) malloc(m * sizeof(int));
    recv_top_row = (int *) malloc(m * sizeof(int));
    recv_bottom_row = (int *) malloc(m * sizeof(int));

    double d_startTime = 0.0, d_endTime = 0.0;
    d_startTime = get_walltime();
    
    /* for each generation, update the game of life board with p
    number of workers */
    for (t=0; t<k; t++) {
        num_alive = 0;

        for (j=0; j<m; j++) {
            send_bottom_row[j] = grid_current[(rows_per_worker - 1) * m + j];
        }
        for (j=0; j<m; j++) {
            send_top_row[j] = grid_current[j];
        }
        /* for first worker, no need to receive a top row,
           for last worker, no need to receive a bottom row,
           for all others, send and receive top and bottom */
        if (taskid == 0) {

            // send bottom row, receive top row of next worker
            MPI_Sendrecv(&send_bottom_row, m, MPI_INT, taskid+1, 1,
                &recv_top_row, m, MPI_INT, taskid+1, 2, MPI_COMM_WORLD,
                &status);

            for (i=1; i<rows_per_worker; i++) {
                for (j=1; j<m-1; j++) {
                    int prev_state = grid_current[i*m+j];

                    /* for the last row, update number of alive cells
                    from received top row, otherwise update with
                    assigned rows */
                    if (i == rows_per_worker - 1) {
                        num_alive += recv_top_row[j-1] + recv_top_row[j] +
                        recv_top_row[j+1];
                    }
                    else {
                        num_alive += grid_current[(i+1)*m+j-1] + 
                                     grid_current[(i+1)*m+j  ] + 
                                     grid_current[(i+1)*m+j+1];
                    }

                    // update for rows from assigned grid
                    num_alive  += 
                                grid_current[(i  )*m+j-1] + 
                                grid_current[(i  )*m+j+1] + 
                                grid_current[(i-1)*m+j-1] + 
                                grid_current[(i-1)*m+j  ] + 
                                grid_current[(i-1)*m+j+1];
                    
                    // update the cell in grid next to be alive or dead     
                    grid_next[i*m+j] = prev_state * ((num_alive == 2) + (num_alive == 3)) + (1 - prev_state) * (num_alive == 3);
                }
            }
        }
        else if (taskid == p) {

            // send top row, receive bottom row of previous worker
            MPI_Sendrecv(&send_top_row, m, MPI_INT, taskid-1, 1,
                &recv_bottom_row, m, MPI_INT, taskid-1, 2, MPI_COMM_WORLD,
                &status);

            for (i=0; i<rows_per_worker + (m % p); i++) {
                for (j=1; j<m-1; j++) {
                    int prev_state = grid_current[i*m+j];

                    /* for the first row, update number of alive cells
                    from received top row, otherwise update with
                    assigned rows */
                    if (i == 0) {
                        num_alive += recv_bottom_row[j-1] + recv_bottom_row[j]
                        + recv_bottom_row[j+1];
                    }
                    else {
                        num_alive += grid_current[(i-1)*m+j-1] + 
                                     grid_current[(i-1)*m+j  ] + 
                                     grid_current[(i-1)*m+j+1];
                    }
                    
                    // update for rows from assigned grid
                    num_alive  += 
                                grid_current[(i  )*m+j-1] + 
                                grid_current[(i  )*m+j+1] + 
                                grid_current[(i+1)*m+j-1] + 
                                grid_current[(i+1)*m+j  ] + 
                                grid_current[(i+1)*m+j+1];

                    // update the cell in grid next to be alive or dead
                    grid_next[i*m+j] = prev_state * ((num_alive == 2) + (num_alive == 3)) + (1 - prev_state) * (num_alive == 3);
                }
            }
        }
        else {

            /* send top row, receive bottom row of previous worker
               send bottom row, receive top row of next worker */
            MPI_Sendrecv(&send_top_row, m, MPI_INT, taskid-1, 1,
                &recv_bottom_row, m, MPI_INT, taskid-1, 2, MPI_COMM_WORLD,
                &status);
            MPI_Sendrecv(&send_bottom_row, m, MPI_INT, taskid+1, 1,
                &recv_top_row, m, MPI_INT, taskid+1, 2, MPI_COMM_WORLD,
                &status);

            for (i=0; i<rows_per_worker + (m % p); i++) {
                for (j=1; j<m-1; j++) {
                    int prev_state = grid_current[i*m+j];
                    
                    /* for the first row, update number of alive cells
                    from received top row, for the bottom row, update 
                    from received bottom row, otherwise update with
                    assigned rows */
                    if (i == 0) {
                        num_alive += recv_bottom_row[j-1] + recv_bottom_row[j] +
                        recv_bottom_row[j+1];
                        num_alive += grid_current[(i+1)*m+j-1] + 
                                      grid_current[(i+1)*m+j  ] + 
                                      grid_current[(i+1)*m+j+1];
                    }
                    else if (i == rows_per_worker - 1) {
                        num_alive += recv_top_row[j-1] + recv_top_row[j]
                        + recv_top_row[j+1];
                                      grid_current[(i-1)*m+j  ] + 
                                      grid_current[(i-1)*m+j+1];
                    }
                    else {
                        num_alive += grid_current[(i-1)*m+j-1] + 
                                      grid_current[(i-1)*m+j  ] + 
                                      grid_current[(i-1)*m+j+1] +
                                      grid_current[(i+1)*m+j-1] + 
                                      grid_current[(i+1)*m+j  ] + 
                                      grid_current[(i+1)*m+j+1];
                    }

                    // update for rows from assigned grid
                    num_alive  += 
                                grid_current[(i  )*m+j-1] + 
                                grid_current[(i  )*m+j+1];

                    // update the cell in grid next to be alive or dead
                    grid_next[i*m+j] = prev_state * ((num_alive == 2) + (num_alive == 3)) + (1 - prev_state) * (num_alive == 3);
                }
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
    free(send_top_row); free(send_bottom_row);
    free(recv_top_row); free(recv_bottom_row);

    MPI_Finalize();

    return 0;
}