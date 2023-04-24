/*
  File: blurfilter.c
  Implementation of blurfilter function.
 */
 
#include "blurfilter.h"

pixel* pix(pixel* image, const int xx, const int yy, const int xsize)
{
	int off = xsize*yy + xx;
	return (image + off);
}

void blurfilter(const int xsize, const int ysize, const int world_rank, int world_size, pixel *src, double *w, int radius)
{	
	MPI_Datatype MPI_PIXEL;
  	MPI_Type_contiguous(3, MPI_UNSIGNED_CHAR, &MPI_PIXEL);
  	MPI_Type_commit(&MPI_PIXEL);
	
	/* Row */
	
	int* displs = (int *)malloc(sizeof(int)*(world_size));
	int* chunksize = (int*) malloc(sizeof(int)*(world_size));
	split_workload(chunksize, displs, xsize, ysize, world_size);
	int send_size = chunksize[world_rank];
	pixel* rcv = (pixel*) malloc(sizeof(pixel) * send_size);
	pixel* rcvb = (pixel*) malloc(sizeof(pixel) * send_size);

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
	blurfilter_row(xsize, (chunksize[world_rank]/xsize), rcv, rcvb, radius, w);
	MPI_Gatherv(rcvb, send_size, MPI_PIXEL, src, chunksize, displs, MPI_PIXEL, 0, MPI_COMM_WORLD);

	/* Column */

	split_workload(chunksize, displs, ysize, xsize, world_size);
	send_size = chunksize[world_rank];

	pixel* dst = (pixel*) malloc(sizeof(pixel) * xsize * ysize);
	pixel* rcv2 = (pixel*) malloc(sizeof(pixel) * send_size);
	pixel* rcvb2 = (pixel*) malloc(sizeof(pixel) * send_size);
	
	if(world_rank==0)
	{
		rotate(dst,src,xsize,ysize,0);
	}	

	MPI_Scatterv(dst, chunksize, displs, MPI_PIXEL, rcv2, send_size, MPI_PIXEL, 0, MPI_COMM_WORLD);
	blurfilter_row(ysize, (chunksize[world_rank]/ysize), rcv2, rcvb2, radius, w);
	MPI_Gatherv(rcvb2, send_size, MPI_PIXEL, dst, chunksize, displs, MPI_PIXEL, 0, MPI_COMM_WORLD);

	if(world_rank==0) 
	{
		rotate(src,dst,xsize,ysize,1);
	}

	free(displs);
	free(chunksize);
	free(rcv);
	free(rcvb);
	free(dst);
	free(rcv2);
	free(rcvb2);
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

void rotate(pixel* output, pixel* input, const int xsize, const int ysize, const int clockwise)
{
	int x2, y2, x0, y0;
	for(int yy=0; yy<ysize; yy++)
	{
		for(int xx=0; xx<xsize; xx++)
		{
			if(clockwise==0) // + 90° rotation.
			{
				x2 = yy;
				y2 = -xx + (xsize-1);

				pix(output,x2,y2,ysize)->r=pix(input,xx,yy,xsize)->r;
				pix(output,x2,y2,ysize)->g=pix(input,xx,yy,xsize)->g;
				pix(output,x2,y2,ysize)->b=pix(input,xx,yy,xsize)->b;
			}	
			
			if(clockwise==1) // - 90° rotation.
			{
				x2 = yy;
				y2 = -xx + (xsize-1);

				pix(output,xx,yy,xsize)->r=pix(input,x2,y2,ysize)->r;
				pix(output,xx,yy,xsize)->g=pix(input,x2,y2,ysize)->g;
				pix(output,xx,yy,xsize)->b=pix(input,x2,y2,ysize)->b;
			}	
			
			
		}
	}	
}

void blurfilter_row(const int xsize, const int ysize, pixel* src, pixel* dst, const int radius, const double *w)
{
	int x, y, x2, y2, wi;
	double r, g, b, n, wc;

	for (y=0; y<ysize; y++)
	{
		for (x=0; x<xsize; x++)
		{
			r = w[0] * pix(src, x, y, xsize)->r;
			g = w[0] * pix(src, x, y, xsize)->g;
			b = w[0] * pix(src, x, y, xsize)->b;
			n = w[0];
			for ( wi=1; wi <= radius; wi++)
			{
				wc = w[wi];
				x2 = x - wi;
				if (x2 >= 0)
				{
					r += wc * pix(src, x2, y, xsize)->r;
					g += wc * pix(src, x2, y, xsize)->g;
					b += wc * pix(src, x2, y, xsize)->b;
					n += wc;
				}
				x2 = x + wi;
				if (x2 < xsize)
				{
					r += wc * pix(src, x2, y, xsize)->r;
					g += wc * pix(src, x2, y, xsize)->g;
					b += wc * pix(src, x2, y, xsize)->b;
					n += wc;
				}
			}
			pix(dst,x,y, xsize)->r = r/n;
			pix(dst,x,y, xsize)->g = g/n;
			pix(dst,x,y, xsize)->b = b/n;
		}
	}
}
