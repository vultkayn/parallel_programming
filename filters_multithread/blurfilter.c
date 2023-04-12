/*
  File: blurfilter.c
  Implementation of blurfilter function.
 */
 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>

#include "blurfilter.h"
#include "ppmio.h"

pixel* pix(pixel* image, const int xx, const int yy, const int xsize)
{
	int off = xsize*yy + xx;
	return (image + off);
}

int blurfilter(const int xsize, const int ysize, pixel* src, const int radius, const double *w)
{
	int part = 6;

	if(part>1){
		printf("Filtering with multiple threads.\n");
	}
	else{
		printf("Sequential execution.\n");
	}
	
	pthread_t thread[part];
	pixel *dst = (pixel*) malloc(sizeof(pixel)*xsize*ysize);
	int* partition = (int*) malloc(sizeof(int)*(part+1));
	Arguments *args = (Arguments*) malloc(sizeof(Arguments)*part);

	set_domain(partition,ysize,part);

	for (int p=0;p<part;p++)
	{
		args[p].xsize = xsize;
		args[p].ysize = ysize;
		args[p].src = src;
		args[p].radius = radius;
		args[p].w = w;
		args[p].dst = dst;
		args[p].vmin = partition[p];
		args[p].vmax = partition[p+1];
		
		if(pthread_create(&thread[p], NULL, &row_processing, (void*) &args[p]) != 0)
		{
			perror("Failed to create the threads.\n");
			return 1;
		}	
	}

	wait_threads(thread,part);

	set_domain(partition,xsize,part);

	for (int p=0;p<part;p++)
	{
		args[p].xsize = xsize;
		args[p].ysize = ysize;
		args[p].src = src;
		args[p].radius = radius;
		args[p].w = w;
		args[p].dst = dst;
		args[p].vmin = partition[p];
		args[p].vmax = partition[p+1];
		
		if(pthread_create(&thread[p], NULL, &column_processing, (void*) &args[p]) != 0)
		{
			perror("Failed to create the threads.\n");
			return 1;
		}		
	}

	wait_threads(thread,part);

	free(partition);
	free(args);
	free(dst);

	return 0;
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

void* row_processing(void* args)
{
	int x2, wi;
	double r, g, b, n, wc;
	for (int y=(((Arguments*)args)->vmin); y<((Arguments*)args)->vmax; y++)
	{
		for (int x=0; x < ((Arguments*)args)->xsize; x++)
		{
			r = ((Arguments*)args)->w[0] * pix(((Arguments*)args)->src, x, y, ((Arguments*)args)->xsize)->r;
			g = ((Arguments*)args)->w[0] * pix(((Arguments*)args)->src, x, y, ((Arguments*)args)->xsize)->g;
			b = ((Arguments*)args)->w[0] * pix(((Arguments*)args)->src, x, y, ((Arguments*)args)->xsize)->b;
			n = ((Arguments*)args)->w[0];
			for ( wi=1; wi <= ((Arguments*)args)->radius; wi++)
			{
				wc = ((Arguments*)args)->w[wi];
				x2 = x - wi;
				if (x2 >= 0)
				{
					r += wc * pix(((Arguments*)args)->src, x2, y, ((Arguments*)args)->xsize)->r;
					g += wc * pix(((Arguments*)args)->src, x2, y, ((Arguments*)args)->xsize)->g;
					b += wc * pix(((Arguments*)args)->src, x2, y, ((Arguments*)args)->xsize)->b;
					n += wc;
				}
				x2 = x + wi;
				if (x2 < ((Arguments*)args)->xsize)
				{
					r += wc * pix(((Arguments*)args)->src, x2, y, ((Arguments*)args)->xsize)->r;
					g += wc * pix(((Arguments*)args)->src, x2, y, ((Arguments*)args)->xsize)->g;
					b += wc * pix(((Arguments*)args)->src, x2, y, ((Arguments*)args)->xsize)->b;
					n += wc;
				}
			}
			pix(((Arguments*)args)->dst, x, y, ((Arguments*)args)->xsize)->r = r/n;
			pix(((Arguments*)args)->dst, x, y, ((Arguments*)args)->xsize)->g = g/n;
			pix(((Arguments*)args)->dst, x, y, ((Arguments*)args)->xsize)->b = b/n;
		}	
	}	
	return NULL;
}

void* column_processing(void* args)
{
	int y2, wi;
	double r, g, b, n, wc;
	for (int x=((Arguments*)args)->vmin; x < ((Arguments*)args)->vmax; x++)
	{
		for (int y=0; y<((Arguments*)args)->ysize; y++)
		{
			r = ((Arguments*)args)->w[0] * pix(((Arguments*)args)->dst, x, y, ((Arguments*)args)->xsize)->r;
			g = ((Arguments*)args)->w[0] * pix(((Arguments*)args)->dst, x, y, ((Arguments*)args)->xsize)->g;
			b = ((Arguments*)args)->w[0] * pix(((Arguments*)args)->dst, x, y, ((Arguments*)args)->xsize)->b;
			n = ((Arguments*)args)->w[0];
			for ( wi=1; wi <= ((Arguments*)args)->radius; wi++)
			{
				wc = ((Arguments*)args)->w[wi];
				y2 = y - wi;
				if (y2 >= 0)
				{
					r += wc * pix(((Arguments*)args)->dst, x, y2, ((Arguments*)args)->xsize)->r;
					g += wc * pix(((Arguments*)args)->dst, x, y2, ((Arguments*)args)->xsize)->g;
					b += wc * pix(((Arguments*)args)->dst, x, y2, ((Arguments*)args)->xsize)->b;
					n += wc;
				}
				y2 = y + wi;
				if (y2 < ((Arguments*)args)->ysize)
				{
					r += wc * pix(((Arguments*)args)->dst, x, y2, ((Arguments*)args)->xsize)->r;
					g += wc * pix(((Arguments*)args)->dst, x, y2, ((Arguments*)args)->xsize)->g;
					b += wc * pix(((Arguments*)args)->dst, x, y2, ((Arguments*)args)->xsize)->b;
					n += wc;
				}
			}
			pix(((Arguments*)args)->src,x,y, ((Arguments*)args)->xsize)->r = r/n;
			pix(((Arguments*)args)->src,x,y, ((Arguments*)args)->xsize)->g = g/n;
			pix(((Arguments*)args)->src,x,y, ((Arguments*)args)->xsize)->b = b/n;
		}
	}		
	return NULL;
}
