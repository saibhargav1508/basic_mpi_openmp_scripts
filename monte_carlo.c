// Sai Bhargav Mandavilli
// EEL6763

/* monte_carlo.c */

#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define MASTER 0               /* taskid of first task */
#define FROM_MASTER 1          /* setting a message type */
#define FROM_WORKER 2          /* setting a message type */

/*   ttype: type to use for representing time */
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
int taskid, numtasks;
struct timespec begin, end;
double time_spent;

/* Function declarations */
void init_rand_seed(void);
double estimate_g(double lower_bound, double upper_bound, long long int N);
void collect_results(double *result);

/* main function */
int main(int argc, char **argv)
{
	double result = 0.0;
	MPI_Init(&argc, &argv);

	float lower_bound = atof(argv[1]);
	float upper_bound = atof(argv[2]);
	long long int N = atof(argv[3]);
	
	init_rand_seed(); // using srand()
	result = estimate_g(lower_bound, upper_bound, N);
	collect_results(&result);
	
	MPI_Finalize();
	return 0;
}

/* function definitions */
void init_rand_seed(void)
{
	// each rank should start with a different seed
	// only then, the next rand() will generate a different sequence for each rank
	MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
	srand((unsigned) taskid);
}

double estimate_g(double lower_bound, double upper_bound, long long int N)
{
	int min_tasks, remainder_tasks; /* used to determine rows sent to each worker */
	int iterations; /* determine how many iterations each rank has to perform */
	double x, y; // variables for computing function
	double partialSum = 0;
	
	MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
	MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
	
	min_tasks = N / numtasks;
	remainder_tasks = N % numtasks;
	
	iterations = (taskid <= remainder_tasks) ? min_tasks+1 : min_tasks;
	
	if(taskid == MASTER)
		begin = now();
	
	for(int i=0; i<iterations; i++)
	{
		x = ((double)rand() / (double)RAND_MAX) * (upper_bound - lower_bound) + lower_bound;
		y = ((double)(upper_bound - lower_bound) / (double)N) *  ((8 * sqrt(2*3.14)) / (exp((2 * x) * (2 * x))));
		partialSum += y;
	}
	printf("partial sum of task %d: %lf\n", taskid, partialSum);
	return partialSum;
}

void collect_results(double *result)
{
	double temp;
	int source;              /* task id of message source */
	MPI_Status status;
	
	MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
	MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
		
	if(taskid == MASTER)
	{
		for(source=1; source<numtasks; source++)
		{
			temp = 0;
			MPI_Recv(&temp, 1, MPI_DOUBLE, source, FROM_WORKER, MPI_COMM_WORLD, &status);
			*result += temp;
			
			printf("Received results from task %d: %lf\n",source, temp);
		}
		
		end = now();
		time_spent = tdiff(begin, end);
		
		/* print results */
		printf ("\n");
		printf("******************************************************\n");
		printf("Estimate of integral: %lf\n", *result);
		printf("Time taken for total estimation: %.8f sec\n", time_spent);
		printf("******************************************************\n");
	}
	
	else{
	MPI_Send(result, 1, MPI_DOUBLE, MASTER, FROM_WORKER, MPI_COMM_WORLD);}
}