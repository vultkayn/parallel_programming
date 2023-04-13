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

#ifndef NTHREADS
#define NTHREADS 1
#endif

unsigned waiting_threads = 0;
pthread_mutex_t waiting_mtx;
pthread_cond_t cond_var;

pixel* pix(pixel* image, const int xx, const int yy, const int xsize)
{
	int off = xsize*yy + xx;
	return (image + off);
}

int blurfilter(const int xsize, const int ysize, pixel* src, const int radius, const double *w)
{
	int nextra_threads = NTHREADS;
	pthread_cond_init(&cond_var, NULL);
	pthread_mutex_init(&waiting_mtx, NULL);


	if(nextra_threads>1){
		printf("Filtering with %d threads.\n", NTHREADS);
	}
	else{
		puts("Sequential execution.");
	}
	
	pthread_t thread[nextra_threads];
	pixel *dst = (pixel*) malloc(sizeof(pixel)*xsize*ysize);
	// thread i compute lines in range (partition[i], partition[i+1]-1)
	int partition[nextra_threads+1];
	Arguments *args = (Arguments*) malloc(sizeof(Arguments)*(nextra_threads+1));

	int p;
	for (p=0; p<nextra_threads; p++)
	{
		args[p].xsize = xsize;
		args[p].ysize = ysize;
		args[p].src = src;
		args[p].radius = radius;
		args[p].w = w;
		args[p].dst = dst;
		args[p].tid = p;
		
		if (pthread_create(&thread[p], NULL, &row_processing, (void*) &args[p]) != 0)
		{
			perror("Failed to create the threads.\n");
			return 1;
		}
	}

	for (int tid = 0; tid < NTHREADS; tid++)
		pthread_join(thread[tid], NULL);

	pthread_mutex_destroy(&waiting_mtx);
	pthread_cond_destroy(&cond_var);

	free(args);
	free(dst);

	return 0;
}

void barrier()
{
	pthread_mutex_lock(&waiting_mtx);
	waiting_threads++;
	if(waiting_threads != NTHREADS)
		pthread_cond_wait(&cond_var, &waiting_mtx);
	
	pthread_cond_signal(&cond_var);
	pthread_mutex_unlock(&waiting_mtx);
}


void* row_processing(void* data)
{
	int x2, wi;
	double r, g, b, n, wc;
	Arguments *args = (Arguments *) data;
	int step = floor(args->ysize/NTHREADS);
	args->vmin = step * args->tid;
	args->vmax = step * (args->tid+1);

	for (int y=(args->vmin); y<args->vmax; y++)
	{
		for (int x=0; x < args->xsize; x++)
		{
			r = args->w[0] * pix(args->src, x, y, args->xsize)->r;
			g = args->w[0] * pix(args->src, x, y, args->xsize)->g;
			b = args->w[0] * pix(args->src, x, y, args->xsize)->b;
			n = args->w[0];
			for ( wi=1; wi <= args->radius; wi++)
			{
				wc = args->w[wi];
				x2 = x - wi;
				if (x2 >= 0)
				{
					r += wc * pix(args->src, x2, y, args->xsize)->r;
					g += wc * pix(args->src, x2, y, args->xsize)->g;
					b += wc * pix(args->src, x2, y, args->xsize)->b;
					n += wc;
				}
				x2 = x + wi;
				if (x2 < args->xsize)
				{
					r += wc * pix(args->src, x2, y, args->xsize)->r;
					g += wc * pix(args->src, x2, y, args->xsize)->g;
					b += wc * pix(args->src, x2, y, args->xsize)->b;
					n += wc;
				}
			}
			pix(args->dst, x, y, args->xsize)->r = r/n;
			pix(args->dst, x, y, args->xsize)->g = g/n;
			pix(args->dst, x, y, args->xsize)->b = b/n;
		}	
	}

	barrier();

	column_processing(data);

	return NULL;
}

void* column_processing(void* data)
{
	int y2, wi;
	double r, g, b, n, wc;
	Arguments *args = (Arguments *) data;

	int step = floor(args->xsize/NTHREADS);
	args->vmin = step * args->tid;
	args->vmax = step * (args->tid+1);

	for (int x=args->vmin; x < args->vmax; x++)
	{
		for (int y=0; y< args->ysize; y++)
		{
			pixel *pxl = pix(args->dst, x, y, args->xsize);
			r = args->w[0] * pxl->r;
			g = args->w[0] * pxl->g;
			b = args->w[0] * pxl->b;
			n = args->w[0];
			for ( wi=1; wi <= args->radius; wi++)
			{
				wc = args->w[wi];
				y2 = y - wi;
				if (y2 >= 0)
				{
					pixel *pxl2 = pix(args->dst, x, y2, args->xsize);
					r += wc * pxl2->r;
					g += wc * pxl2->g;
					b += wc * pxl2->b;
					n += wc;
				}
				y2 = y + wi;
				if (y2 < args->ysize)
				{
					pixel *pxl3 = pix(args->dst, x, y2, args->xsize);
					r += wc * pxl3->r;
					g += wc * pxl3->g;
					b += wc * pxl3->b;
					n += wc;
				}
			}
			pix(args->src,x,y, args->xsize)->r = r/n;
			pix(args->src,x,y, args->xsize)->g = g/n;
			pix(args->src,x,y, args->xsize)->b = b/n;
		}
	}		
	return NULL;
}
