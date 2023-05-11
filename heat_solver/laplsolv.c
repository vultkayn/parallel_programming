
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
#include <assert.h>
#include <omp.h>

#ifndef NCOLS
#define NCOLS 2
#endif

double timediff(struct timespec *begin, struct timespec *end)
{
	double sec = 0.0, nsec = 0.0;
	if ((end->tv_nsec - begin->tv_nsec) < 0)
	{
		sec = (double)(end->tv_sec - begin->tv_sec - 1);
		nsec = (double)(end->tv_nsec - begin->tv_nsec + 1000000000);
	}
	else
	{
		sec = (double)(end->tv_sec - begin->tv_sec);
		nsec = (double)(end->tv_nsec - begin->tv_nsec);
	}
	return sec + nsec / 1E9;
}

void printm(int n, double *M)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
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

int compute_end_col(int thrd, int rectw, int n)
{
	int beg_col = (thrd % NCOLS) * rectw;
	if ((thrd % NCOLS) == NCOLS -1)
		return n;
	return beg_col + rectw;
}

int compute_end_row(int thrd, int nthrds, int recth, int n)
{
	int beg_row = (thrd / NCOLS) * recth;
	if (nthrds - NCOLS <= thrd && thrd <= nthrds-1)
		return n;
	return beg_row + recth;
}

// assumes we have a number of threads multiple of NCOLS
double compute_region(
	int thrd,
	int nthrds,
	int rectw,
	int recth,
	int n,
	double T[n + 2][n + 2],
	double *top,
	double *bot,
	double *left,
	double *right)
{
	double error = -INFINITY;
	int beg_col = (thrd % NCOLS) * rectw;
	int beg_row = (thrd / NCOLS) * recth;
	int end_row = compute_end_row(thrd, nthrds, recth, n);
	int end_col = compute_end_col(thrd, rectw, n);
	int width = end_col - beg_col;
	double tmp2[width];
	double tmp1[width];
	// Copy to temp buffers
	arrcpy(tmp1, &T[beg_row][beg_col + 1], width);

	// printf("thrd: %d, (rectw, recth): (%d,%d), beg_row: %d, beg_col: %d\n", thrd, rectw, recth, beg_row, beg_col);

	// Loop for each of this thread's rows
	for (int i = beg_row + 1; i <= end_row; ++i)
	{
		assert(i <= n && "i > n");
		assert(i >= 1 && "i < 0");
		arrcpy(tmp2, &T[i][beg_col + 1], width);

		// Apply the Jacobi algorithm to each element in this row
		double prev = T[i][beg_col];
		for (int j = beg_col + 1; j <= end_col; ++j)
		{
			assert(j <= n && "j > n");
			assert(j >= 1 && "j < 0");
			double next = (prev + T[i][j + 1] + T[i + 1][j] + tmp1[j - (beg_col + 1)]) / 4.0;
			prev = T[i][j];
			short border = 0;
			if (j == beg_col + 1)
			{
				left[i - (beg_row + 1)] = next;
				border = 1;
			}
			else if (j == end_col)
			{
				right[i - (beg_row + 1)] = next;
				border = 1;
			}

			if (i == beg_row + 1)
			{
				top[j - (beg_col + 1)] = next;
				border = 1;
			}
			else if (i == end_row)
			{
				bot[j - (beg_col + 1)] = next;
				border = 1;
			}

			if (!border)
				T[i][j] = next;

			error = fmax(error, fabs(tmp2[j - (beg_col + 1)] - next));
		}

		arrcpy(tmp1, tmp2, width);
	}

	return error;
}

void write_future_to_T(
	int thrd,
	int nthrds,
	int rectw, int recth,
	int n,
	double T[n + 2][n + 2],
	double *top,
	double *bot,
	double *left,
	double *right)
{
	int beg_col = (thrd % NCOLS) * rectw;
	int beg_row = (thrd / NCOLS) * recth;
	int end_row = compute_end_row(thrd, nthrds, recth, n);
	int end_col = compute_end_col(thrd, rectw, n);
	int width = end_col - beg_col;

	arrcpy(&T[beg_row + 1][beg_col + 1], top, width);
	arrcpy(&T[end_row][beg_col + 1], bot, width);

	for (int i = beg_row + 1; i <= end_row; ++i)
	{
		T[i][beg_col + 1] = left[i - (beg_row + 1)];
		T[i][end_col] = right[i - (beg_row + 1)];
	}
}

void laplsolv(int n, int maxiter, double tol)
{
	double T[n + 2][n + 2];
	int k;

	struct timespec starttime, endtime;

	// Set boundary conditions and initial values for the unknowns
	for (int i = 0; i <= n + 1; ++i)
	{
		for (int j = 0; j <= n + 1; ++j)
		{
			if (i == n + 1)
				T[i][j] = 2;
			else if (j == 0 || j == n + 1)
				T[i][j] = 1;
			else
				T[i][j] = 0;
		}
	}

	clock_gettime(CLOCK_MONOTONIC, &starttime);

#pragma omp parallel
	{
		printf("Hello from Thread %d of %d\n",
			   omp_get_thread_num(),
			   omp_get_num_threads());
	}

	// Solve the linear system of equations using the Jacobi method
	int rectw = n / NCOLS;
	for (k = 0; k < maxiter; ++k)
	{
		double error = -INFINITY;
#pragma omp parallel shared(error, T, n)
		{
			int recth = NCOLS * n / (omp_get_num_threads());
			double top[rectw], bot[rectw], left[recth], right[recth];
			double local_error = compute_region(omp_get_thread_num(), omp_get_num_threads(), rectw, recth, n, T, top, bot, left, right);
#pragma omp barrier
			// all threads copy their saves to T
			write_future_to_T(omp_get_thread_num(), omp_get_num_threads(), rectw, recth, n, T, top, bot, left, right);
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
	#ifdef DEBUG
	printm(n + 2, T);
	#endif
}

int main(int argc, char *argv[])
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