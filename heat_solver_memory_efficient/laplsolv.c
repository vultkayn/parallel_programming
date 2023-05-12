
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

double compute(int beg, int end, int n, double T[n + 2][n + 2])
{
	if (omp_get_thread_num() == 0)
	{
		beg = 1;
	}

	int length = end - beg + 1;
	double tmp1[length], tmp2[length], temp3, temp4, next;
	double error = -INFINITY;

	// Copy to temp buffers
	arrcpy(tmp1, &T[0][beg], length);

	for (int i = 1; i <= n; ++i)
	{
		temp3 = T[i][beg - 1];
		temp4 = T[i][end + 1];

		arrcpy(tmp2, &T[i][beg], length);

#pragma omp barrier

		double prev = T[i][beg - 1];
		for (int j = beg; j <= end; ++j)
		{
			if (j == beg)
			{
				next = (temp3 + T[i][j + 1] + T[i + 1][j] + tmp1[j - beg]) / 4.0;
			}
			else if (j == end)
			{
				next = (prev + temp4 + T[i + 1][j] + tmp1[j - beg]) / 4.0;
			}
			else
			{
				next = (prev + T[i][j + 1] + T[i + 1][j] + tmp1[j - beg]) / 4.0;
			}

			error = fmax(error, fabs(tmp2[j - beg] - next));
			prev = T[i][j];
			T[i][j] = next;
		}
		arrcpy(tmp1, tmp2, length);
	}
	return error;
}

void set_domain(int *partition, const int size)
{
	int step = floor(size / omp_get_num_threads());
	partition[0] = 1;
	partition[omp_get_num_threads()] = size;

	for (int k = 1; k < (omp_get_num_threads()); k++)
	{
		partition[k] = step * k;
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
		printf("Hello from Thread %d of %d\n", omp_get_thread_num(), omp_get_num_threads());
		// assert(n<omp_get_num_threads() && "n>=omp_get_num_threads()");
	}

	// Solve the linear system of equations using the Jacobi method
	for (k = 0; k < maxiter; ++k)
	{
		double error = -INFINITY;

#pragma omp parallel shared(error, T, n)
		{
			int split_region[omp_get_num_threads() + 1];
			set_domain(split_region, n);
			double local_error = compute(split_region[omp_get_thread_num()] + 1, split_region[omp_get_thread_num() + 1], n, T);
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
	// printm(n+2, &T[0][0]);
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