# mpi_openmp_scripts
This repository contains MPI and OpenMP based parallel programs.

## Matrix multiple using MPI
Multiplies two matrices by dividing the rows of second matrix between MPI ranks. Synchronization between ranks was achieved using MPI_Send() and MPI_Recv().

## Monte-Carlo Integration with MPI
One method of numerically estimating integrals is by using Monte-Carlo simulation. Consider the following integral g(a,b) and its estimate h(x):
![image](https://github.com/saibhargav1508/basic_mpi_openmp_scripts/assets/20701792/4d6dce5d-3510-4eb2-b2bc-36bbeb37f0b2)

The numerical solution to g(a,b) can be estimated using a uniform random variable x that is evenly distributed over [a, b] (i.e., between a and b). The estimate h(x) will converge to the correct solution as the number of samples N grows. Since each sample is independent, the calculation can be easily parallelized (embarrassingly parallel). Using this method, write a short MPI program that will use N samples to calculate:

![image](https://github.com/saibhargav1508/basic_mpi_openmp_scripts/assets/20701792/2cf7c639-520d-44ee-aa1e-75c9756cb733)

For this problem you must provide code for three functions, init_rand_seed (), estimate_g(…) and collect_results(…), which will be used by the following main(…) function:
```C
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
  return
}
```

Other specifications:
- The total number of samples to calculate N (to be split among all MPI nodes), as well as the bounds of the integral a and b, will be provided through command-line arguments.
- Every MPI rank will generate its own random numbers. Ensure that each rank uses a different starting seed for its random number generator - this should be the only thing that occurs in , init_rand_seed ().

Then in the estimate_g function, you should use the rand function to generate the random numbers. Hint: While debugging print out the random numbers to verify that they are indeed random.

- Each rank should return a single value, which should be combined (i.e., sum) in the master rank in order to compute the final integral.
- Use only the following functions:
  -   MPI_Init
  -   MPI_Comm_rank
  -   MPI_Comm_Size
  -   MPI_Send
  -   MPI_Recv
  -   MPI_Finalize

## Sliding windows techniques for image processing
Commonly used for or image processing, a sliding window algorithm uses a “window” (a rectangular region of fixed width and height) that “slides” across an image, as shown in Figure 1. For each of instance of these windows, the value of the pixel of interest for that window is recalculated based on a filter mask applied to that pixel. In Figure 1, the pixel of interest is q2, q3, q4, and q5, respectively, for each of the windows. The mask is generally an equation(s) developed specifically for the filtering operation (e.g., a Sobel filter). The basic concept is to recalculate the value of each pixel of the entire picture based on the values of the current pixel the adjacent pixels, as the window slides across each pixel.

![image](https://github.com/saibhargav1508/basic_mpi_openmp_scripts/assets/20701792/ebeecbda-a4fd-45ac-9e67-93b27d68c72a)

For this lab, we will use a simple neighborhood weighted-averaging filter. The operator takes the weighted average of adjacent pixel values, as illustrated in Figure 2.

  e' = (a+b+c+d+2e+f+g+h+i)/10

![image](https://github.com/saibhargav1508/basic_mpi_openmp_scripts/assets/20701792/c46543dd-d466-4976-b0c3-d1cdb6e1ee08)

Lab specifications:
- Implement the mask filtering operation on an image, represented by an N x N integer matrix, using MPI on R ranks (varying N and R on different runs).
- The mask operations must be distributed as evenly as possible for all ranks, including the master rank.
- Use MPI_Scatterv() to distribute the initial image matrix and use MPI_Gatherv() to collect the processed matrix. These are the only two MPI statements you can use for data transfer (e.g., no MPI_Send or MPI_Recv).
- To keep program simple, it is not required to process the first and last rows and columns.
- The matrix A should be initialized in one of the nodes with synthetic data using the rand() function. To simplify grading (i.e., everyone has the same input array), please initialize with the seed of 1 (the default value) and limit the random numbers to be between 0 - 255.
- You will provide code for four functions (initialize_data, distribute_data, mask_operation, and collect_results) with the following prototypes (you can change the prototypes to fit your coding style):
```C
void initialize_data (int *A, int N);
int* distribute_data (int *A, int N); // returns recv_buff for that rank
int* mask_operation (int *recv_buff, int N); // returns updated_buff for that rank
void collect_results (int *updated_buff, int N, int *Ap); // *Ap is processed matrix
```

Use the following pseudocode for your main function (plus timing functions, etc.). N is size of the NXN Matrix. Depending on how you define function prototypes, you can change how you call them in the main() function.

```C
int main(int argc, char **argv)
{ 
  MPI_Init(&argc, &argv);
  int N = atof(argv[1]);
  initialize_data(A, N);
  temp1 = distribute_data (A, N); // use scatterv
  temp2 = mask_operation(temp1, N);
  collect_results(temp2, N, Ap); // use gatherv
  MPI_Finalize();
  return 0;
}
```

## Sobel Edge Detector – Hybrid Programming Model
Find a serial implementation of the Sobel Edge Filter and convert it to MPI+OpenMP model. The Sobel Edge filter converts an image as shown below -

![image](https://github.com/saibhargav1508/basic_mpi_openmp_scripts/assets/20701792/fe7b84af-523a-4114-b034-429b2436d1a0)

Matlab code to convert image to text and vice versa is given in this repository.
Reference serial code for the Sobel Filter is also given.
