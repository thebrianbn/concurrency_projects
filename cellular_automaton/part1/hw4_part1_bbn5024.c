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

int* serial_life(int m, int k) {

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

    /* initializing some cells */
    for (i=1; i<3; i++) {
        for (j=1; j<m-1; j++) {
            grid_current[i*m+j] = 1;
        }
    }

    double d_startTime = 0.0, d_endTime = 0.0;
    d_startTime = get_walltime();

    /* serial code */
    /* considering only internal cells */
    for (t=0; t<k; t++) {
        for (i=1; i<m-1; i++) {
            for (j=1; j<m-1; j++) {

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

    return grid_current;

    /* free memory */
    free(grid_current); free(grid_next);

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

    // if last worker, handle any extra rows
    if (taskid == p) {

        grid_current = (int *) malloc((rows_per_worker + (m % numtasks)) * m * sizeof(int));
        if (grid_current == NULL) {
            printf("Error allocating memory for grid_current!\n");
            MPI_Finalize();
            exit(1);
        }

        grid_next = (int *) malloc((rows_per_worker + (m % numtasks)) * m * sizeof(int));
        if (grid_next == NULL) {
            printf("Error allocating memory for grid_next!\n");
            MPI_Finalize();
            exit(1);
        }

        for (i=0; i<rows_per_worker + (m % numtasks); i++) {
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
            MPI_Finalize();
            exit(1);
        }

        grid_next = (int *) malloc(rows_per_worker * m * sizeof(int));
        if (grid_next == NULL) {
            printf("Error allocating memory for grid_next!\n");
            MPI_Finalize();
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
    if (taskid == 0) {
        for (i=1; i<3; i++) {
            for (j=1; j<m-1; j++) {
                grid_current[i*m+j] = 1;
            }
        }
    }

    int *recv_top_row;
    int *recv_bottom_row;
    int dest_up = taskid-1;
    int dest_down = taskid+1;

    // arrays used for sendrecv
    recv_top_row = (int *) malloc(m * sizeof(int));
    recv_bottom_row = (int *) malloc(m * sizeof(int));

    double d_startTime = 0.0, d_endTime = 0.0;
    d_startTime = get_walltime();
    
    /* for each generation, update the game of life board with p
    number of workers */
    for (t=0; t<k; t++) {

        int num_alive;
        int prev_state;

        /* for first worker, no need to receive a top row,
           for last worker, no need to receive a bottom row,
           for all others, send and receive top and bottom */
        if (taskid == 0) {

            // send bottom row, receive top row of next worker
            MPI_Sendrecv(&grid_current[m * (rows_per_worker-1)], m, MPI_INT, (taskid+1), 0,
                recv_top_row, m, MPI_INT, taskid+1, 0, MPI_COMM_WORLD,
                &status);

            // update with grid for top rows
            for (i=1; i<rows_per_worker-1; i++) {
                for (j=1; j<m-1; j++) {
                    prev_state = grid_current[i*m+j];

                    // avoid repeated calculations
                    int top_index = (i-1)*m;
                    int bottom_index = (i+1)*m;
                    int center_index = i*m;
                    int left_shift = j-1;
                    int right_shift = j+1;

                    num_alive = grid_current[center_index+left_shift] + 
                                grid_current[center_index+right_shift] + 
                                grid_current[top_index+left_shift] + 
                                grid_current[top_index+j] + 
                                grid_current[top_index+right_shift] + 
                                grid_current[bottom_index+left_shift] + 
                                grid_current[bottom_index+j  ] + 
                                grid_current[bottom_index+right_shift];
                                               
                    // update the cell in grid next to be alive or dead     
                    grid_next[i*m+j] = prev_state * ((num_alive == 2) + (num_alive == 3)) + (1 - prev_state) * (num_alive == 3);
                }
            }

            /* for the last row, update number of alive cells
                    from received top row */
            for (j=1; j<m-1; j++) {
                prev_state = grid_current[j];

                num_alive = recv_top_row[j-1] + recv_top_row[j]
                        + recv_top_row[j+1];

                num_alive += grid_current[(rows_per_worker-2)*m+j-1] + 
                             grid_current[(rows_per_worker-2)*m+j  ] + 
                             grid_current[(rows_per_worker-2)*m+j+1] +
                             grid_current[(rows_per_worker-1)*m+j-1] + 
                             grid_current[(rows_per_worker-1)*m+j+1];

                // update the cell in grid next to be alive or dead
                grid_next[(rows_per_worker-1)*m+j] = prev_state * ((num_alive == 2) + (num_alive == 3)) + (1 - prev_state) * (num_alive == 3);
            }
        }
        else if (taskid != p) {

            /* send top row, receive bottom row of previous worker
               send bottom row, receive top row of next worker */
            MPI_Sendrecv(grid_current, m, MPI_INT, taskid-1, 0,
                recv_bottom_row, m, MPI_INT, taskid-1, 0, MPI_COMM_WORLD,
                &status);
            MPI_Sendrecv(&grid_current[m * (rows_per_worker-1)], m, MPI_INT, taskid+1, 0,
                recv_top_row, m, MPI_INT, taskid+1, 0, MPI_COMM_WORLD,
                &status);

            /* for the first row, update number of alive cells
                    from received bottom row */
            for (j=1; j<m-1; j++) {
                prev_state = grid_current[j];

                num_alive = recv_bottom_row[j-1] + recv_bottom_row[j] +
                        recv_bottom_row[j+1];

                num_alive += grid_current[m+j-1] + 
                             grid_current[m+j  ] + 
                             grid_current[m+j+1] +
                             grid_current[j-1] + 
                             grid_current[j+1];

                // update the cell in grid next to be alive or dead
                grid_next[j] = prev_state * ((num_alive == 2) + (num_alive == 3)) + (1 - prev_state) * (num_alive == 3);
            }
            
            // update with grid for middle rows
            for (i=1; i<rows_per_worker-1; i++) {
                for (j=1; j<m-1; j++) {
                    prev_state = grid_current[i*m+j];

                    // avoid repeated calculations
                    int top_index = (i-1)*m;
                    int bottom_index = (i+1)*m;
                    int center_index = i*m;
                    int left_shift = j-1;
                    int right_shift = j+1;
                    
                    num_alive = grid_current[center_index+left_shift] + 
                                grid_current[center_index+right_shift] + 
                                grid_current[top_index+left_shift] + 
                                grid_current[top_index+j] + 
                                grid_current[top_index+right_shift] + 
                                grid_current[bottom_index+left_shift] + 
                                grid_current[bottom_index+j  ] + 
                                grid_current[bottom_index+right_shift];

                    // update the cell in grid next to be alive or dead
                    grid_next[i*m+j] = prev_state * ((num_alive == 2) + (num_alive == 3)) + (1 - prev_state) * (num_alive == 3);
            
                }
            }

            // for the bottom row, update from received top row
            for (j=1; j<m-1; j++) {
                prev_state = grid_current[j];

                num_alive = recv_top_row[j-1] + recv_top_row[j]
                        + recv_top_row[j+1];

                num_alive += grid_current[(rows_per_worker-2)*m+j-1] + 
                             grid_current[(rows_per_worker-2)*m+j  ] + 
                             grid_current[(rows_per_worker-2)*m+j+1] +
                             grid_current[(rows_per_worker-1)*m+j-1] + 
                             grid_current[(rows_per_worker-1)*m+j+1];

                // update the cell in grid next to be alive or dead
                grid_next[(rows_per_worker-1)*m+j] = prev_state * ((num_alive == 2) + (num_alive == 3)) + (1 - prev_state) * (num_alive == 3);
            }
        }
        else {

            // send top row, receive bottom row of previous worker
            MPI_Sendrecv(grid_current, m, MPI_INT, taskid-1, 0,
                recv_bottom_row, m, MPI_INT, taskid-1, 0, MPI_COMM_WORLD,
                &status);

            /* for the first row, update number of alive cells
                    from received bottom row, */
            for (j=1; j<m-1; j++) {
                prev_state = grid_current[j];

                num_alive = recv_bottom_row[j-1] + recv_bottom_row[j]
                        + recv_bottom_row[j+1];

                // update for rows from assigned grid
                num_alive  += 
                            grid_current[j-1] + 
                            grid_current[j+1] + 
                            grid_current[m+j-1] + 
                            grid_current[m+j  ] + 
                            grid_current[m+j+1];

                // update the cell in grid next to be alive or dead
                grid_next[j] = prev_state * ((num_alive == 2) + (num_alive == 3)) + (1 - prev_state) * (num_alive == 3);
            }

            // update with grid for rest of rows
            for (i=1; i<rows_per_worker + (m % numtasks)-1; i++) {
                for (j=1; j<m-1; j++) {
                    prev_state = grid_current[i*m+j];

                    // avoid repeated calculations
                    int top_index = (i-1)*m;
                    int bottom_index = (i+1)*m;
                    int center_index = i*m;
                    int left_shift = j-1;
                    int right_shift = j+1;

                    num_alive  = 
                                grid_current[center_index+left_shift] + 
                                grid_current[center_index+right_shift] + 
                                grid_current[top_index+left_shift] + 
                                grid_current[top_index+j] + 
                                grid_current[top_index+right_shift] + 
                                grid_current[bottom_index+left_shift] + 
                                grid_current[bottom_index+j  ] + 
                                grid_current[bottom_index+right_shift];

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

    double time_diff = d_endTime - d_startTime;

    printf("Time taken for worker %d: %3.3lf s.\n", taskid, time_diff);
    printf("Performance: %3.3lf billion cell updates/s\n", 
                (1.0*m*m)*k/((d_endTime - d_startTime)*1e9));

    // collect all grid results in worker 0, and create the final grid
    // other workers send their grids and times
    if (taskid == 0) {

        int *grid_whole;
        grid_whole = (int *) malloc(m * m * sizeof(int));

        int *last_grid;
        last_grid = (int *) malloc((rows_per_worker + (m % numtasks)) * m * sizeof(int));

        double *worker_times;
        worker_times = (double *) malloc(numtasks * sizeof(double));
        worker_times[0] = d_endTime - d_startTime;

        double worker_time_current;

        // put in data from worker 0
        for (i=0; i<rows_per_worker; i++) {
            for (j=0; j<m; j++) {
                int current_index = i*m+j;
                grid_whole[current_index] = grid_current[current_index];
            }
        }

        // put in data from middle workers
        int workerid;
        for (workerid=1; workerid < numtasks-1; workerid++) {

            MPI_Recv(grid_current, rows_per_worker*m, MPI_INT, workerid, 0,
                MPI_COMM_WORLD, &status);
            MPI_Recv(&worker_time_current, 1, MPI_DOUBLE, workerid, 0,
                MPI_COMM_WORLD, &status);

            for (i=0; i<rows_per_worker; i++) {
                for (j=0; j<m; j++) {
                    int current_index = i*m+j;
                    grid_whole[(rows_per_worker*workerid*m)+current_index] = grid_current[current_index];
                }
            }
            worker_times[workerid] = worker_time_current;
        }

        // put in data from last worker
        MPI_Recv(last_grid, (rows_per_worker + (m % p)) *m, MPI_INT, p, 0,
                MPI_COMM_WORLD, &status);
        MPI_Recv(&worker_time_current, 1, MPI_DOUBLE, p, 0, MPI_COMM_WORLD, &status);

        for (i=0; i<rows_per_worker + (m % numtasks); i++) {
            for (j=0; j<m; j++) {
                int current_index = i*m+j;
                grid_whole[(rows_per_worker*p*m)+current_index] = last_grid[current_index];
            }
        }
        worker_times[p] = worker_time_current;

        // retrieve serial implementation's board
        int *serial_board = serial_life(m, k);

        // compare the serial grid against the parallel grid
        int verify_failed = 0;
        for (i=0; i<m; i++) {
            for (j=0; j<m; j++) {
                int current_index = i*m+j;
                if (grid_whole[current_index] != serial_board[current_index]) {
                    verify_failed += 1;
                }
            }
        }

        printf("Failed count: %d\n", verify_failed);

        // free memeory
        free(last_grid); free(grid_whole);
        free(worker_times); free(serial_board);
    }
    else if (taskid == p) {
        MPI_Send(grid_current, (rows_per_worker + (m % p)) * m, MPI_INT, 0, 0, MPI_COMM_WORLD);
        MPI_Send(&time_diff, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
    else {
        MPI_Send(grid_current, rows_per_worker * m, MPI_INT, 0, 0, MPI_COMM_WORLD);
        MPI_Send(&time_diff, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }

    // free memory
    free(grid_current); free(grid_next);
    free(recv_top_row); free(recv_bottom_row);

    MPI_Finalize();

    return 0;
}
