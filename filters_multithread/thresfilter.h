/*
  File: thresfilter.h
  Declaration of pixel structure and thresfilter function.
 */

#ifndef _THRESFILTER_H_
#define _THRESFILTER_H_

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <math.h>

#define uint unsigned int 

/* NOTE: This structure must not be padded! */
typedef struct _pixel {
	unsigned char r,g,b;
} pixel;

struct Arguments {
  uint* sum;
  uint nump;
  int vmin;
  int vmax;
  pixel* src;
  uint threshold;
}; typedef struct Arguments Arguments;

int thresfilter(const int xsize, const int ysize, pixel* src);
void* threshold_processing(void* args);
void set_domain(int* partition, const int size, const int part);
void wait_threads(pthread_t* thread, const int part);
void* threshold_computing(void* args);

#endif
