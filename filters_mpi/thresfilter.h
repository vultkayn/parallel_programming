/*
  File: thresfilter.h
  Declaration of pixel structure and thresfilter function.
 */

#ifndef _THRESFILTER_H_
#define _THRESFILTER_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#define uint unsigned int 

/* NOTE: This structure must not be padded! */
typedef struct _pixel {
	unsigned char r,g,b;
} pixel;

void thresfilter(const int xsize, const int ysize, pixel* src, const int world_rank, const int world_size);
void output(const uint size, const uint sum, pixel* src);
void split_workload(int* chunksize, int* displs, const int size, const int split, const int world_size);
uint threshold_value(const uint size, pixel* src);

#endif
