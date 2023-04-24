#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <mpi.h>

#include "ppmio.h"
#include "blurfilter.h"
#include "gaussw.h"

#define MAX_RAD 1000

struct img_info {
	unsigned radius;
	unsigned xsize;
	unsigned ysize;
};

int main (int argc, char ** argv)
{
	MPI_Init(NULL, NULL);

	int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	struct img_info img_info;
	MPI_Datatype MPI_IMG_INFO;
	MPI_Type_contiguous(3, MPI_UNSIGNED, &MPI_IMG_INFO);
	MPI_Type_commit(&MPI_IMG_INFO);
	// struct timespec stime, etime;
	double stime, etime;
	double w[MAX_RAD];

	pixel *src = NULL;
	if(world_rank==0)
	{
		src = (pixel*) malloc(sizeof(pixel) * MAX_PIXELS);
		int colmax;
		
		/* Take care of the arguments */

		if (argc != 4)
		{
			fprintf(stderr, "Usage: %s radius infile outfile\n", argv[0]);
			exit(1);
		}
		
		img_info.radius = atoi(argv[1]);
		if ((img_info.radius > MAX_RAD) || (img_info.radius < 1))
		{
			fprintf(stderr, "Radius (%d) must be greater than zero and less then %d\n", img_info.radius, MAX_RAD);
			exit(1);
		}

		/* Read file */

		if (read_ppm (argv[2], &img_info.xsize, &img_info.ysize, &colmax, (char *) src) != 0)
			exit(1);
		
		if (colmax > 255)
		{
			fprintf(stderr, "Too large maximum color-component value\n");
			exit(1);
		}

		printf("Has read the image, generating coefficients\n");
		
		/* Filtering */

		get_gauss_weights(img_info.radius, w);
		
		printf("Calling filter\n");

		clock_gettime(CLOCK_REALTIME, &stime);
		stime = MPI_Wtime();
	}

	MPI_Bcast(&img_info, 1, MPI_IMG_INFO, 0, MPI_COMM_WORLD);
	MPI_Bcast(w, (img_info.radius+1), MPI_DOUBLE, 0, MPI_COMM_WORLD);

	blurfilter(img_info.xsize, img_info.ysize, world_rank, world_size, src, w, img_info.radius);
	
	if(world_rank==0)
	{		
		// clock_gettime(CLOCK_REALTIME, &etime);
		etime = MPI_Wtime();

		// printf("Filtering took: %g secs\n", (etime.tv_sec  - stime.tv_sec) + 1e-9*(etime.tv_nsec - stime.tv_nsec)) ;
		printf("Filtering took: %g secs\n", (etime  - stime)) ;
	
		/* Write result */
		if(write_ppm(argv[3], img_info.xsize, img_info.ysize, (char *)src) != 0)
			exit(1);

		printf("Writing output file\n");
	}	

	free(src);

	MPI_Finalize();
}
