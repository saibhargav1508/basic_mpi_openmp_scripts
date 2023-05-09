// Sai Bhargav Mandavilli
// EEL6763

/* sed_hybrid.c */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <time.h>
#include <mpi.h>
#include <omp.h>

#define WIDTH 5000
#define HEIGHT 5000

/*   ttype: type to use for representing time */
typedef double ttype;
ttype tdiff(struct timespec a, struct timespec b)
/* Find the time difference. */
{
  ttype dt = ((b.tv_sec - a.tv_sec) + (b.tv_nsec - a.tv_nsec) / 1E9);
  return dt;
}

struct timespec now()
/* Return the current time. */
{
  struct timespec t;
  clock_gettime(CLOCK_REALTIME, &t);
  return t;
}

 // clock_t begin, end;
struct timespec begin, end;
double time_spent;

int main(int argc, char *argv[])
{
    int *edgeImage; // stores chunks of pixel values
    MPI_Init(&argc, &argv);
    FILE *imgInput, *imgOutput; // file pointers to store input and output image text
    int taskid, numtasks, I, J, K, L, SUM, previous_displ = 0, previous_displB = 0;
    int *num_rows, *processed_matrix;
    unsigned int X, Y;
    long sumX, sumY;
	static int scatter_temp = 1, gather_temp = 1;
	int *inputText = (int *)malloc((HEIGHT * WIDTH) * sizeof(int)); // array to store pixel values from input text file

	int GX[9] = {-1, 0, 1, -2, 0, 2, -1, 0, 1}; // GX sobel mask
	int GY[9] = {1, 2, 1, 0, 0, 0, -1, -2, -1}; // GY sobel mask

    // read pixel values from input text file
	imgInput = fopen("input.txt", "r");
    for(L = 0; L < HEIGHT; L++)
        for (K = 0; K < WIDTH; K++)
            fscanf(imgInput, "%d", &inputText[L*WIDTH + K]);
    fclose(imgInput);

    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);

    int N = HEIGHT;
    int averow, extra;
    averow = (N - 2) / numtasks;
    extra = (N - 2) % numtasks;
    
    int *sendcountsA = malloc(sizeof(int) * numtasks);
	int *sendcountsB = malloc(sizeof(int) * numtasks);
	int *displsA = malloc(sizeof(int) * numtasks);
    int *displsB = malloc(sizeof(int) * numtasks);
    num_rows = malloc(sizeof(int) * numtasks);
    
	for (I = 0; I < numtasks; I++){
        scatter_temp = (extra > 0) ? 1 : 0;
        num_rows[I] = (averow + scatter_temp + 2);
        sendcountsA[I] = num_rows[I] * N;
        displsA[I] = previous_displ; // updating displs for next iteration
        previous_displ += (averow + scatter_temp) * N;
        if (extra > 0)
            extra--;
    }

    int *rbuf = (int *)malloc(sendcountsA[taskid] * sizeof(int));

	begin = now(); // input is ready to be distributed
    MPI_Scatterv(inputText, sendcountsA, displsA, MPI_INT, rbuf, sendcountsA[taskid], MPI_INT, 0, MPI_COMM_WORLD);
    edgeImage = (int *)malloc((num_rows[taskid] * WIDTH) * sizeof(int));
	
	#pragma omp parallel for shared(SUM, sumX, sumY, sendcountsA, rbuf, num_rows, taskid, edgeImage) private(I, J, X, Y)
	for (Y = 0; Y <= num_rows[taskid] - 1; Y++){
		for (X = 0; X <= WIDTH - 1; X++){
			sumX = 0;
			sumY = 0;

			/* image boundaries */
			if (Y == 0 || Y == sendcountsA[taskid] - 1)
				SUM = 0;
			else if (X == 0 || X == WIDTH - 1)
				SUM = 0;

			// /* Convolution starts here */
			else{
				/*-------X GRADIENT APPROXIMATION------*/
				for (I = 0; I < 3; I++)
				{
					for (J = 0; J < 3; J++)
					{
						int t = 3 * I + J;
						sumX = sumX + (int)((*(rbuf + X + I + (Y + J) * WIDTH)) * GX[t]);
					}
				}
				/*-------Y GRADIENT APPROXIMATION-------*/
				for (I = 0; I < 3; I++)
				{
					for (J = 0; J < 3; J++)
					{
						int t = 3 * I + J;
						sumY = sumY + (int)((*(rbuf + X + I +(Y + J) * WIDTH)) *GY[t]);
					}
				}
				/*---GRADIENT MAGNITUDE APPROXIMATION (Myler p.218)----*/
				SUM = abs(sumX) + abs(sumY);
			}
			if(SUM > 255) // for 8-bit precision
				SUM = 255;
			else if(SUM < 0)
				SUM = 0;
			edgeImage[(Y*WIDTH) + X] = (unsigned int)(SUM);
		}
	}

	averow = (N - 2) / numtasks; // initialize rows and extra rows distribution again
    extra = (N - 2) % numtasks;
    
	for (I = 0; I < numtasks; I++){
        gather_temp = (extra > 0) ? 1 : 0;
        sendcountsB[I] = (averow + gather_temp) * N;
        displsB[I] = previous_displB;
        previous_displB += sendcountsB[I]; // updating displs for next iteration
        if(extra > 0)
            extra--;
    }
    
	int *worker_output = (int *)malloc(N * (N - 2) * sizeof(int)); // we have N matrix of N-2 rows each
    processed_matrix = (int *)malloc(N * N * sizeof(int));

    MPI_Gatherv(edgeImage, sendcountsB[taskid], MPI_INT, worker_output, sendcountsB, displsB, MPI_INT, 0, MPI_COMM_WORLD);
	end = now(); // result is gathered
	time_spent = tdiff(begin, end);

    if(0 == taskid){
		printf("total time for applying filter: %.8f sec\n", time_spent);
        imgOutput = fopen("processed_matrix.txt", "w+"); // output text file
        for(L=0; L<HEIGHT; L++){
            for(K=0; K<WIDTH; K++){
				if(L>0 && L < (WIDTH - 1)){
                    processed_matrix[L * WIDTH + K] = worker_output[(L-1) * WIDTH + K];
                    fprintf(imgOutput, "%d\t", processed_matrix[L * N + K]);
                }
                else if(L == 0 || L == (WIDTH - 1)){
                    processed_matrix[L * WIDTH + K] = 0;
                    fprintf(imgOutput, "%d\t", processed_matrix[L * WIDTH + K]);
                }
            }
            fprintf(imgOutput, "\n");
        }
		fclose(imgOutput);
    }
    MPI_Finalize();
    return 0;
}