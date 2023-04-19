/*
  File: blurfilter.h
  Declaration of pixel structure and blurfilter function.
 */

#ifndef _BLURFILTER_H_
#define _BLURFILTER_H_

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h> 

/* NOTE: This structure must not be padded! */
typedef struct _pixel {
	unsigned char r,g,b;
} pixel;

void blurfilter_row(const int xsize, const int ysize, pixel* src, pixel* dst, const int radius, const double *w);
pixel* pix(pixel* image, const int xx, const int yy, const int xsize);
void rotate(pixel* output, pixel* input, const int xsize, const int ysize, const int clockwise);
void blurfilter(const int xsize, const int ysize, const int world_rank, int world_size, pixel *src, double *w, int radius);
void split_workload(int* chunksize, int* displs, const int size, const int split, const int world_size);

#endif
