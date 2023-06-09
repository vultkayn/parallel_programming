/*
  File: blurfilter.h
  Declaration of pixel structure and blurfilter function.
 */

#ifndef _BLURFILTER_H_
#define _BLURFILTER_H_

/* NOTE: This structure must not be padded! */
typedef struct _pixel {
	unsigned char r,g,b;
} pixel;

typedef struct Arguments {
  pixel* src;
  const double* w;
  pixel* dst;
  int xsize;
  int ysize;
  int vmin;
  int vmax;
  int radius;
  unsigned tid;
} Arguments;

int blurfilter(const int xsize, const int ysize, pixel* src, const int radius, const double *w);
void set_domain(int* partition, const int size);
void wait_threads(pthread_t* thread);
void* row_processing(void *args);
void* column_processing(void *args);

#endif
