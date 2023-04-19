#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <mpi.h>

#include "ppmio.h"
#include "blurfilter.h"
#include "gaussw.h"

#define MAX_RAD 1000

int main (int argc, char ** argv)
{
	MPI_Init(NULL, NULL);

	int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	int radius, xsize, ysize;
	struct timespec stime, etime;
	double w[MAX_RAD];
	pixel *src = (pixel*) malloc(sizeof(pixel) * MAX_PIXELS);

	if(world_rank==0)
	{
		int colmax;
		
		/* Take care of the arguments */

		if (argc != 4)
		{
			fprintf(stderr, "Usage: %s radius infile outfile\n", argv[0]);
			exit(1);
		}
		
		radius = atoi(argv[1]);
		if ((radius > MAX_RAD) || (radius < 1))
		{
			fprintf(stderr, "Radius (%d) must be greater than zero and less then %d\n", radius, MAX_RAD);
			exit(1);
		}

		/* Read file */

		if (read_ppm (argv[2], &xsize, &ysize, &colmax, (char *) src) != 0)
			exit(1);
		
		if (colmax > 255)
		{
			fprintf(stderr, "Too large maximum color-component value\n");
			exit(1);
		}

		printf("Has read the image, generating coefficients\n");
		
		/* Filtering */

		get_gauss_weights(radius, w);
		
		printf("Calling filter\n");

		clock_gettime(CLOCK_REALTIME, &stime);
	}

	MPI_Bcast(&xsize,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&ysize,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&radius,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(w,(radius+1),MPI_DOUBLE,0,MPI_COMM_WORLD);

	blurfilter(xsize, ysize, world_rank, world_size, src, w, radius);
	
	if(world_rank==0)
	{		
		clock_gettime(CLOCK_REALTIME, &etime);
	
		printf("Filtering took: %g secs\n", (etime.tv_sec  - stime.tv_sec) + 1e-9*(etime.tv_nsec - stime.tv_nsec)) ;
	
		/* Write result */
		if(write_ppm(argv[3], xsize, ysize, (char *)src) != 0)
			exit(1);

		printf("Writing output file\n");
	}	

	free(src);

	MPI_Finalize();
}
