#include "thresfilter.h"

#ifndef NTHREADS
#define NTHREADS 1
#endif

unsigned waiting_threads = 0;
pthread_mutex_t waiting_mtx;
pthread_cond_t cond_var;

uint threshold = 0;

int thresfilter(const int xsize, const int ysize, pixel* src)
{
	uint i, nump;
	nump = xsize * ysize;
	
	int nextra_threads = NTHREADS;

	if(nextra_threads>1){
		printf("Filtering with %d threads.\n", NTHREADS);
	}
	else{
		puts("Sequential execution.");
	}

	pthread_t thread[nextra_threads];
	int* partition = (int*) malloc(sizeof(int)*(nextra_threads+1));
	Arguments *args = (Arguments*) malloc(sizeof(Arguments)*nextra_threads);
	uint* sum = (uint*)malloc(sizeof(uint)*nextra_threads);

	set_domain(partition,nump,nextra_threads);
	
	for(uint p=0;p<nextra_threads;p++)
	{
		args[p].sum = &sum[p];
		args[p].nump = nump;
		args[p].vmin = partition[p];
		args[p].vmax = partition[p+1];
		args[p].src = src;

		if(pthread_create(&thread[p], NULL, &threshold_value, (void*) &args[p]) != 0)
		{
			perror("Failed to create the threads.\n");
			return 1;
		}
	}	

	wait_threads(thread,nextra_threads);

	pthread_mutex_destroy(&waiting_mtx);
	pthread_cond_destroy(&cond_var);

	free(partition);
	free(args);
}

void* output(void* args)
{
	uint psum;
	for (int i = ((Arguments*)args)->vmin; i < ((Arguments*)args)->vmax; i++)
	{
		psum = (uint)(((Arguments*)args)->src[i].r) + (uint)(((Arguments*)args)->src[i].g) + (uint)(((Arguments*)args)->src[i].b);
		if (threshold > psum)
		{
			((Arguments*)args)->src[i].r = ((Arguments*)args)->src[i].g = ((Arguments*)args)->src[i].b = 0;
		}
		else
		{
			((Arguments*)args)->src[i].r = ((Arguments*)args)->src[i].g = ((Arguments*)args)->src[i].b = 255;
		}
	}
}

void* threshold_value(void* args)
{
	*(((Arguments*)args)->sum) = 0;

	for (int i = ((Arguments*)args)->vmin; i < ((Arguments*)args)->vmax; i++)
	{
		*(((Arguments*)args)->sum) += (uint)(((Arguments*)args)->src[i].r) + (uint)(((Arguments*)args)->src[i].g) + (uint)(((Arguments*)args)->src[i].b);
	}

	barrier((Arguments *) args);

	output(args);
}

void wait_threads(pthread_t* thread, const int nextra_threads)
{
	for (int p=0;p<nextra_threads;p++)
	{
		if(pthread_join(thread[p],NULL) != 0){
			exit(1);
		}
	}
}

void barrier(Arguments* args)
{
	pthread_mutex_lock(&waiting_mtx);
	waiting_threads++;

	if(waiting_threads != NTHREADS) {
		threshold += *args->sum;
		pthread_cond_wait(&cond_var, &waiting_mtx);
	}
	else { // last to arrive will be the only one to run the else
		threshold += *args->sum;
		threshold /= args->nump;
	}

	pthread_cond_signal(&cond_var);
	pthread_mutex_unlock(&waiting_mtx);
}

void set_domain(int* partition, const int size, const int nextra_threads)
{
	int step = floor(size/nextra_threads);
	partition[0] = 0;
	partition[nextra_threads] = size;

	for(int k=1;k<nextra_threads;k++)
	{ 
		partition[k] = step*k;
	}
}