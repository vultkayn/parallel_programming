#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <mpi.h>

#include "ppmio.h"
#include "thresfilter.h"

int main (int argc, char ** argv)
{
	MPI_Init(NULL, NULL);

	int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	
	int xsize, ysize;
	pixel *src = (pixel*) malloc(sizeof(pixel) * MAX_PIXELS);
	struct timespec stime, etime;

	if(world_rank==0)
	{
		int colmax;
	
		/* Take care of the arguments */
		if (argc != 3)
		{
			fprintf(stderr, "Usage: %s infile outfile\n", argv[0]);
			exit(1);
		}
		
		/* Read file */
		if(read_ppm (argv[1], &xsize, &ysize, &colmax, (char *) src) != 0)
			exit(1);
		
		if (colmax > 255)
		{
			fprintf(stderr, "Too large maximum color-component value\n");
			exit(1);
		}
	
		printf("Has read the image, calling filter\n");

		clock_gettime(CLOCK_REALTIME, &stime);
	}

	MPI_Bcast(&xsize,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&ysize,1,MPI_INT,0,MPI_COMM_WORLD);

	thresfilter(xsize, ysize, src, world_rank, world_size);

	if(world_rank==0)
	{		
		clock_gettime(CLOCK_REALTIME, &etime);
	
		printf("Filtering took: %g secs\n", (etime.tv_sec  - stime.tv_sec) + 1e-9*(etime.tv_nsec - stime.tv_nsec)) ;
	
		/* Write result */
		if (write_ppm(argv[2], xsize, ysize, (char *)src) != 0)
			exit(1);

		printf("Writing output file\n");
	}	

	free(src);

	MPI_Finalize();
}
