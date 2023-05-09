// Sai Bhargav Mandavilli
// EEL6763

/* scatter_gather.c */

#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

/* ttype: type to use for representing time */
typedef double ttype;
ttype tdiff(struct timespec a, struct timespec b)
/* Find the time difference. */
{
  ttype dt = (( b.tv_sec - a.tv_sec ) + ( b.tv_nsec - a.tv_nsec ) / 1E9);
  return dt;
}

struct timespec now()
/* Return the current time. */
{
  struct timespec t;
  clock_gettime(CLOCK_REALTIME, &t);
  return t;
}

/* global variables */
struct timespec begin, end;
double timespent;
int rank,size; // for storing this process' rank, and the number of processes
int *sendcounts; 	// array describing how many elements to send to each process
int *displs; 		// array describing the displacements where each segment begins
int N;

/* function declaration */
void initialize_data(int N,int A[N][N]);
void distribute_data(int N,int A[N][N], int recv_buff[N][N]); // stores recv_buff for that rank
void mask_operation(int N, int recv_buff[N][N], int Ap[N][N], int updated_buff[N][N]); // stores updated_buff for that rank
void collect_results(int N, int updated_buff[N][N], int A[N][N]);

int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);
	int N = atof(argv[1]);
	int A[N][N], recv_buff[N][N], updated_buff[N][N], Ap[N][N];

	begin = now();
	
	initialize_data(N, A);
	distribute_data (N, A, recv_buff); // use scatterv
	mask_operation(N, recv_buff, Ap, updated_buff);
	collect_results(N, updated_buff, A); // use gatherv
	
	end = now();
	timespent = tdiff(begin, end);
	
	if(rank == 0)
		printf("time spent: %.8f sec\n", timespent); 
	MPI_Finalize();

	free(sendcounts);
	free(displs);
	
	return 0;
}

/* function definitions */
void initialize_data(int N, int A[N][N])
{
    srand(1);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if(0 == rank){
		printf("Matrix A:\n");
        for(int i=0; i<N; i++) {
            for(int j=0; j<N; j++) {
               A[i][j] = rand() % 256; // to limit the random numbers to be between 0 - 255
               printf("%d\t", A[i][j]);
            }
			printf("\n");
        }
    }
}

void distribute_data(int N, int A[N][N], int recv_buff[N][N])
{
    int rem;            // elements remaining after division among processes
    int sum=0; 		// Sum of counts. Used to calculate displacements
    int i, j, k=0; 
	int *temp; // temporary receive buffer
    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	for(i=0; i<N; i++)
        for(j=0; j<N; j++)
            recv_buff[i][j] = 0; // clearing out receive buffer for each rank

    rem = (N - 2) % size;
    sendcounts = malloc(sizeof(int)*size);
    displs = malloc(sizeof(int)*size);
    temp = malloc(sizeof(int)*size);

    // calculate send counts and displacements
    for(i=0; i<size; i++){
        temp[i] = (N-2)/size;
        if(rem > 0){
          temp[i]++;
          rem--;
        }

        displs[i] = sum * N;
        sum += temp[i];
    }

    for(i=0; i<size; i++){
		sendcounts[i] = (temp[i] + 2) * N;}

    // print calculated send counts and displacements for each process  
    if(0 == rank){
		printf("sendcounts and displs for scatterv()\n");
        for(i=0; i<size; i++){
            printf("sendcounts[%d] = %d\tdispls[%d] = %d\n", i, sendcounts[i], i, displs[i]);
        }        
    }
    
	// divide the data among processes as described by sendcounts and displs
	MPI_Scatterv(A, sendcounts, displs, MPI_INT, recv_buff, N*N, MPI_INT, 0, MPI_COMM_WORLD);

    printf("part of the matrix that were sent to %d \n", rank);
	for(i=0; i<N; i++){
		for(j=0; j<N; j++){
			printf("%d\t", recv_buff[i][j]);
		}
		printf("\n");
	}
	
	free(temp);
}

void mask_operation (int N, int recv_buff[N][N], int Ap[N][N], int updated_buff[N][N])
{
    int i,j;
	
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
    
    // apply filter to create processed matrix Ap
	for (i=1; i<((sendcounts[rank]/N) - 1); i++){
        for(j=1; j<N-1; j++){
			// covering all neighbor pixels
            Ap[i][j] = (
				recv_buff[i-1][j-1] + recv_buff[i-1][j] + recv_buff[i-1][j+1]
				+ recv_buff[i][j-1] + (recv_buff[i][j] * 2) + recv_buff[i][j+1]
				+ recv_buff[i+1][j-1] + recv_buff[i+1][j] + recv_buff[i+1][j+1]
				)/10;
        }
    }
    
    for(i=1; i<((sendcounts[rank]/N) - 1); i++){
        for(j=1; j<N-1; j++){
            recv_buff[i][j] = Ap[i][j];
        }
    }

    if(rank == 0){   
        printf("updated values for rank :%d\n",rank);
		for(i=0; i<((sendcounts[rank]/N) - 1); i++){
			for(j=0; j<N; j++){
				updated_buff[i][j] = recv_buff[i][j];
                printf("%d\t", updated_buff[i][j]);
			}	
            printf("\n");
		}
	}

    else if(rank == (size-1)){
        printf("updated values for rank:%d\n",rank);
		for(i=1; i<(sendcounts[rank]/N); i++){
			for(j=0; j<N; j++){
				updated_buff[i-1][j] = recv_buff[i][j];
                printf("%d\t", updated_buff[i-1][j]);
			}
            printf("\n");
		}

	}
    
	else{
        printf("updated values for rank :%d\n",rank);
		for(i=1; i<((sendcounts[rank]/N) - 1); i++){
			for(j=0; j<N; j++){
				updated_buff[i-1][j] = recv_buff[i][j];
                printf("%d\t",updated_buff[i-1][j]);
			}	
            printf("\n");
		}
	}
}

void collect_results (int N, int updated_buff[N][N], int A[N][N])
{   
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
	
    for(int i=0; i<size; i++){
		if(i==0)
			sendcounts[i] = sendcounts[i] - N;

		else if(i == (size-1)){
			sendcounts[i] = sendcounts[i] - N;
			displs[i] = displs[i] + N;
		}

		else{
			sendcounts[i] = sendcounts[i] - N * 2;
			displs[i] = displs[i] + N;
		}
	}
	
	if(0 == rank){
		printf("sendcounts and displs for gatherv()\n");
        for(int i=0; i<size; i++){
            printf("sendcounts[%d] = %d\tdispls[%d] = %d\n", i, sendcounts[i], i, displs[i]);
        }        
    }

    // Gather the data from all ranks
	MPI_Gatherv(updated_buff, sendcounts[rank], MPI_INT, A, sendcounts, displs, MPI_INT, 0, MPI_COMM_WORLD);
    
	if (0 == rank){
		printf ("\n");
		printf ("Updated final (total) matrix:\n");
		printf("******************************************************\n");
        for(int i=0; i<N; i++){
            for(int j=0; j<N; j++){
               printf ("%d\t", A[i][j]);
            }
			printf ("\n");
        }
		printf("******************************************************\n");
    }    
}