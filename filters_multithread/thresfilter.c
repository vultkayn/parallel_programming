#include "thresfilter.h"

int thresfilter(const int xsize, const int ysize, pixel* src)
{
	uint i, nump;
	nump = xsize * ysize;
	
	
	int part = 6;

	if(part>1){
		printf("Filtering with multiple threads.\n");
	}
	else{
		printf("Sequential execution.\n");
	}

	pthread_t thread[part];
	int* partition = (int*) malloc(sizeof(int)*(part+1));
	Arguments *args = (Arguments*) malloc(sizeof(Arguments)*part);
	uint* sum = (uint*)malloc(sizeof(uint)*part);

	set_domain(partition,nump,part);
	
	for(int p=0;p<part;p++)
	{
		args[p].sum = &sum[p];
		args[p].nump = nump;
		args[p].vmin = partition[p];
		args[p].vmax = partition[p+1];
		args[p].src = src;
		args[p].threshold = 0;

		if(pthread_create(&thread[p], NULL, &threshold_computing, (void*) &args[p]) != 0)
		{
			perror("Failed to create the threads.\n");
			return 1;
		}
	}	

	wait_threads(thread,part);

	uint threshold = 0;
		
	for(int p=0;p<part;p++)
	{
		threshold += sum[p];
	}	

	threshold /= nump;
	
	for(int p=0;p<part;p++)
	{
		args[p].sum = &sum[p];
		args[p].nump = nump;
		args[p].vmin = partition[p];
		args[p].vmax = partition[p+1];
		args[p].threshold = threshold;

		// threshold_processing(&args[p]);

		if(pthread_create(&thread[p], NULL, &threshold_processing, (void*) &args[p]) != 0)
		{
			perror("Failed to create the threads.\n");
			return 1;
		}
	}

	wait_threads(thread,part);

	free(partition);
	free(args);
}

void* threshold_processing(void* args)
{
	uint psum;
	for (int i = 0; i < ((Arguments*)args)->nump; i++)
	{
		psum = (uint)(((Arguments*)args)->src[i].r) + (uint)(((Arguments*)args)->src[i].g) + (uint)(((Arguments*)args)->src[i].b);
		if (((Arguments*)args)->threshold > psum)
		{
			((Arguments*)args)->src[i].r = ((Arguments*)args)->src[i].g = ((Arguments*)args)->src[i].b = 0;
		}
		else
		{
			((Arguments*)args)->src[i].r = ((Arguments*)args)->src[i].g = ((Arguments*)args)->src[i].b = 255;
		}
	}
}

void* threshold_computing(void* args)
{
	*(((Arguments*)args)->sum) = 0;
	for (int i = ((Arguments*)args)->vmin; i < ((Arguments*)args)->vmax; i++)
	{
		*(((Arguments*)args)->sum) += (uint)(((Arguments*)args)->src[i].r) + (uint)(((Arguments*)args)->src[i].g) + (uint)(((Arguments*)args)->src[i].b);
	}
}

void wait_threads(pthread_t* thread, const int part)
{
	for (int p=0;p<part;p++)
	{
		if(pthread_join(thread[p],NULL) != 0){
			exit(1);
		}
	}
}

void set_domain(int* partition, const int size, const int part)
{
	int step = floor(size/part);
	partition[0] = 0;
	partition[part] = size;

	for(int k=1;k<part;k++)
	{ 
		partition[k] = step*k;
	}
}