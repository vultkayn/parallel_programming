
//-----------------------------------------------------------------------
// Serial program for solving the heat conduction problem 
// on a square using the Jacobi method. 
// Written by August Ernstsson 2015-2019
//-----------------------------------------------------------------------

#define _POSIX_C_SOURCE 199309L
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <omp.h>

#ifndef
#define NCOLS 2
#endif

double timediff(struct timespec *begin, struct timespec *end)
{
	double sec = 0.0, nsec = 0.0;
	if ((end->tv_nsec - begin->tv_nsec) < 0)
	{
		sec  = (double)(end->tv_sec  - begin->tv_sec  - 1);
		nsec = (double)(end->tv_nsec - begin->tv_nsec + 1000000000);
	} else
	{
		sec  = (double)(end->tv_sec  - begin->tv_sec );
		nsec = (double)(end->tv_nsec - begin->tv_nsec);
	}
	return sec + nsec / 1E9;
}

void printm(int n, double *M)
{
	for (int i = 0; i < n; i ++)
	{
		for (int j = 0; j < n; j ++)
		{
			printf("%f\t", *(M + n * i + j));
		}
		printf("\n");
	}
	printf("\n");
}


void arrcpy(double *dst, double *src, int len)
{
	for (int it = 0; it < len; it++)
		dst[it] = src[it];
}

// assumes we have a number of threads multiple of NCOLS
double compute_region (
	int thrd,
	int rectw,
	int recth,
	double **T,
	int n,
	double *top,
	double *bot,
	double *left,
	double *right)
{
	double error = -INFINITY;
	int beg_col = (thrd % NCOLS) * rectw; 
	int beg_row = (thrd / NCOLS) * recth; 
	double *tmp2 = bot;
	// Copy to temp buffers
	arrcpy(tmp1, &T[beg_row][beg_col + 1], rectw);

	// Loop for each of this thread's rows
	for (int i = beg_row + 1; i <= beg_row + recth; ++i)
	{
		arrcpy(tmp2, &T[i][beg_col + 1], rectw);
		
		// Apply the Jacobi algorithm to each element in this row
		double prev = T[i][beg_col];
		for (int j = beg_col + 1; j <= beg_col + rectw; ++j)
		{
			double next = (prev + T[i][j+1] + T[i+1][j] + tmp1[j-1]) / 4.0;
			prev = T[i][j];
			if (j == beg_col + 1) {
				left[i] = next;
			}
			else if (j == beg_col + recth) {
				right[i] = next;
			} else
				T[i][j] = next;
			error = fmax(error, fabs(tmp2[j-1] - next));
		}
		
		arrcpy(tmp1, tmp2, rectw);
		if (i == beg_row + 1) {
			arrcpy(top, tmp2, rectw);

		}
	}

	return error;
}

void write_future_to_T (
	int thrd,
	int rectw, int recth,
	double **T,
	double *top,
	double *bot,
	double *left,
	double *right)
{
	int beg_col = (thrd % NCOLS) * rectw;
	int beg_row = (thrd / NCOLS) * recth;

	arrcpy(&T[beg_row+1][beg_col + 1], top, rectw);
	arrcpy(&T[beg_row+recth][beg_col + 1], bot, rectw);

	for (int i = beg_row + 1; i <= beg_row + recth; ++i) {
		T[i][beg_col + 1] = left[i];
		T[i][beg_col + rectw] = right[i];
	}
}

void laplsolv(int n, int maxiter, double tol)
{
	double T[n+2][n+2];
	int k;
	
	struct timespec starttime, endtime;
	
	// Set boundary conditions and initial values for the unknowns
	for (int i = 0; i <= n+1; ++i)
	{
		for (int j = 0; j <= n+1; ++j)
		{
			if      (i == n+1)           T[i][j] = 2;
			else if (j == 0 || j == n+1) T[i][j] = 1;
			else                         T[i][j] = 0;
		}
	}
	
	clock_gettime(CLOCK_MONOTONIC, &starttime);
	
	#pragma omp parallel
	{
	printf("Hello from Thread %d of %d\n",
			omp_get_thread_num(),
			omp_get_num_threads() );
	}

	// Solve the linear system of equations using the Jacobi method
	int rectw = n / NCOLS;
	int recth = NCOLS * n / (nthrds);
	for (k = 0; k < maxiter; ++k)
	{
		double error = -INFINITY;
		#pragma omp parallel shared(error, T, n)
		{
			double top[rectw], bot[rectw], left[recth], right[recth];
			int local_error = compute_region (omp_get_thread_num(), rectw, recth, T, n, top, bot, left, right);
			#pragma omp barrier
			// all threads copy their saves to T
			write_future_to_T (omp_get_thread_num(), rectw, recth, T, top, bot, left, right);
			#pragma omp critical
			{
				error = (error > local_error) ? error : local_error;
			}
		}
		
		if (error < tol)
			break;
	}
	
	clock_gettime(CLOCK_MONOTONIC, &endtime);
	
	printf("Time: %f\n", timediff(&starttime, &endtime));
	printf("Number of iterations: %d\n", k);
	printf("Temperature of element T(1,1): %.17f\n", T[1][1]);
}


int main(int argc, char* argv[])
{
	if (argc != 4)
	{
		printf("Usage: %s [size] [maxiter] [tolerance] \n", argv[0]);
		exit(1);
	}
	
	int size = atoi(argv[1]);
	int maxiter = atoi(argv[2]);
	double tol = atof(argv[3]);
	
	printf("Size %d, max iter %d and tolerance %f.\n", size, maxiter, tol);
	laplsolv(size, maxiter, tol);
	return 0;
}