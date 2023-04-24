#include "thresfilter.h"

void thresfilter(const int xsize, const int ysize, pixel* src, const int world_rank, const int world_size)
{
	MPI_Datatype MPI_PIXEL;
  	MPI_Type_contiguous(3, MPI_UNSIGNED_CHAR, &MPI_PIXEL);
  	MPI_Type_commit(&MPI_PIXEL);
	
	/* Row */
	
	int* displs = (int *)malloc(sizeof(int)*(world_size));
	int* chunksize = (int*) malloc(sizeof(int)*(world_size));
	split_workload(chunksize, displs, xsize, ysize, world_size);
	int send_size = chunksize[world_rank];
	uint sum, temp_sum;

	pixel* rcv = (pixel*) malloc(sizeof(pixel) * send_size);

	if(world_rank==0)
	{
		if(world_size>1){
			printf("Filtering with %d processes.\n", world_size);
		}
		else{
			puts("Sequential execution.");
		}
	}

	MPI_Scatterv(src, chunksize, displs, MPI_PIXEL, rcv, send_size, MPI_PIXEL, 0, MPI_COMM_WORLD);
	
	temp_sum = threshold_value(send_size, rcv);
	
	MPI_Reduce(&temp_sum, &sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

	if(world_rank==0)
	{
		sum /= (xsize*ysize);
	}

	MPI_Bcast(&sum, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	output(send_size, sum, rcv);

	MPI_Gatherv(rcv, send_size, MPI_PIXEL, src, chunksize, displs, MPI_PIXEL, 0, MPI_COMM_WORLD);

	free(displs);
	free(chunksize);
	free(rcv);
}

uint threshold_value(const uint size, pixel* src)
{
	uint i, sum_temp;
	for (i = 0; i < size; i++)
	{
		sum_temp += (uint)src[i].r + (uint)src[i].g + (uint)src[i].b;
	}
	return sum_temp;
}

void output(const uint size, const uint sum, pixel* src)
{
	uint i, psum;
	for (i = 0; i < size; i++)
	{
		psum = (uint)src[i].r + (uint)src[i].g + (uint)src[i].b;
		if (sum > psum)
		{
			src[i].r = src[i].g = src[i].b = 0;
		}
		else
		{
			src[i].r = src[i].g = src[i].b = 255;
		}
	}
}

void split_workload(int* chunksize, int* displs, const int size, const int split, const int world_size)
{
	int step = floor(split/world_size);

	displs[0] = 0; 

	for(int k=0;k<(world_size-1);k++)
	{ 
		chunksize[k] = step*size;
	}

	chunksize[world_size-1] = (split - step*(world_size-1))*size;

	for(int j=1;j<world_size;j++)
	{
		displs[j] = displs[j-1] + chunksize[j-1];
	}
}
