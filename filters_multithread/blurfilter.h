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

struct Arguments {
  int xsize;
  int ysize;
  int vmin;
  int vmax;
  pixel* src;
  int radius;
  const double* w;
  pixel* dst;
}; typedef struct Arguments Arguments;


int blurfilter(const int xsize, const int ysize, pixel* src, const int radius, const double *w);
void set_domain(int* partition, const int size, const int part);
void wait_threads(pthread_t* thread, const int part);
void* row_processing(void *args);
void* column_processing(void *args);


//void line_xv(const int xsize, const int ysize, const int ymin, const int ymax, pixel* src, const int radius, const double *w, pixel *dst);



// void line_yv(const int xsize, const int ysize, const int ymin, const int ymax, pixel* src, const int radius, const double *w, pixel *dst);


#endif
