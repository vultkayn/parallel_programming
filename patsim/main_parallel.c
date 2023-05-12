#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <mpi.h>

#include "coordinate.h"
#include "definitions.h"
#include "physics.h"

typedef struct Row Row;
struct Row
{
    pcord_t particles[INIT_NO_PARTICLES];
	int nb_particles;
    Row *next;
	Row *last;
};

typedef struct Row_list Row_list;
struct Row_list
{
    Row *first;
}; 

typedef struct Send_rows Send_rows;
struct Send_rows
{
    pcord_t particles_current[INIT_NO_PARTICLES];
	pcord_t particles_next[INIT_NO_PARTICLES];
	pcord_t particles_last[INIT_NO_PARTICLES];
	int size[3]; 
}; 

Row_list *init_link_list(int nb_new_particles)
{
    Row_list *row_list = malloc(sizeof(*row_list));
    Row *row = malloc(sizeof(*row));

    if (row_list == NULL || row == NULL)
    {
        exit(EXIT_FAILURE);
    }

	row->nb_particles = nb_new_particles;
    
	for(int i=0;i<nb_new_particles;i++)
	{
		row->particles[i].x = 0.0;
		row->particles[i].y = 0.0;
		row->particles[i].vx = 0.0;
		row->particles[i].vy = 0.0;
		row->particles[i].collide = false; 
	}

    row->next = NULL;
	row->last = NULL;
    row_list->first = row;

    return row_list;
}

void insert(Row_list *liste, int nb_new_particles)
{
    Row *new = malloc(sizeof(*new));
    if (liste == NULL || new == NULL)
    {
        exit(EXIT_FAILURE);
    }

    new->nb_particles = nb_new_particles;

	for(int i=0;i<nb_new_particles;i++)
	{
		new->particles[i].x = 0.0;
		new->particles[i].y = 0.0;
		new->particles[i].vx = 0.0;
		new->particles[i].vy = 0.0;
	} 

	Row *current = liste->first;

	while ((current->next != NULL))
	{
		current = current->next;
	} 

	new->next = NULL;
	new->last = current;
	current->next = new;
}


Row* get_row(Row_list *liste, int nb)
{
    if (liste == NULL)
    {
        exit(EXIT_FAILURE);
    }

    Row *current = liste->first;

	int i = 0;

	while ((i!=nb))
    {
		current = current->next;
		i++;
    }

    return current;
}

void print_particles_row(Row *current)
{
	if(current!=NULL)
	{
		for(int i=0;i<current->nb_particles;i++)
		{
			printf("Particule %d: x=%f, y=%f, Vx=%f, Vy=%f\n",i+1,current->particles[i].x,current->particles[i].y,current->particles[i].vx,current->particles[i].vy);
		}
	} 
	else{
		printf("Empty!\n");
	} 
	
}

float rand1(){
	return (float)( rand()/(float) RAND_MAX );
}

void init_collisions(bool *collisions, unsigned int max){
	for(unsigned int i=0;i<max;++i)
		collisions[i]=0;
}

int main(int argc, char** argv){

	// MPI initialization.
	MPI_Init(NULL, NULL);

	int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	int my_id;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_id);

	MPI_Status status;

	MPI_Datatype MPI_PARTICLES;
	MPI_Type_contiguous(4, MPI_DOUBLE, &MPI_PARTICLES);
	MPI_Type_commit(&MPI_PARTICLES);

	pcord_t last[INIT_NO_PARTICLES]; 
	pcord_t current[INIT_NO_PARTICLES]; 
	pcord_t next[INIT_NO_PARTICLES];

	int size[3];
	int step = 1;
	unsigned int time_stamp = 0, time_max;

	Row_list *link_list = NULL;
	Send_rows *tab = NULL;

	if(my_id==0)
	{
		tab = (Send_rows*) malloc(sizeof(Send_rows)*world_size);

		printf("\n");
		printf("\nWorld size = %d\n",world_size);
		
		float pressure = 0;

		// Parse arguments
		if(argc != 2) {
			fprintf(stderr, "Usage: %s simulation_time\n", argv[0]);
			fprintf(stderr, "For example: %s 10\n", argv[0]);
			exit(1);
		}

		time_max = atoi(argv[1]);

		/* Initialize */
		// 1. set the walls
		cord_t wall;
		wall.y0 = wall.x0 = 0;
		wall.x1 = BOX_HORIZ_SIZE;
		wall.y1 = BOX_VERT_SIZE;

		// 2. allocate particle bufer and initialize the particles
		srand( time(NULL) + 1234 );


		float r, a;

		int row_nb = 3*world_size;

		float domain[row_nb+1]; 

		for(int i=0;i<row_nb+1;i++)
		{
			domain[i] = i*(BOX_VERT_SIZE/row_nb);
			// printf("%d=%f\n",i,domain[i]);
		}	

		link_list = init_link_list(0);

		for(int i=0;i<row_nb+1;i++)
		{
			insert(link_list, 0);
		} 

		pcord_t temp;

		for(int i=0; i<INIT_NO_PARTICLES; i++){
			// initialize random position
			temp.x = wall.x0 + rand1()*BOX_HORIZ_SIZE;
			temp.y = wall.y0 + rand1()*BOX_VERT_SIZE;

			// initialize random velocity
			r = rand1()*MAX_INITIAL_VELOCITY;
			a = rand1()*2*PI;
			temp.vx = r*cos(a);
			temp.vy = r*sin(a);

			int row_num = ceil(temp.y/(BOX_VERT_SIZE/row_nb));

			get_row(link_list,row_num)->particles[get_row(link_list,row_num)->nb_particles].x = temp.x;
			get_row(link_list,row_num)->particles[get_row(link_list,row_num)->nb_particles].y = temp.y;
			get_row(link_list,row_num)->particles[get_row(link_list,row_num)->nb_particles].vx = temp.vx;
			get_row(link_list,row_num)->particles[get_row(link_list,row_num)->nb_particles].vy = temp.vy;
			get_row(link_list,row_num)->nb_particles++;
		}

		for(int i=0;i<row_nb;i++)
		{
			printf("Row %d:\n",i+1);
			print_particles_row(get_row(link_list, i+1));
		}


		
	}	

	// for (time_stamp=0; time_stamp<time_max; time_stamp++) { // for each time stamp
	for (time_stamp=0; time_stamp<1; time_stamp++) 
	{
		for(step=1;step<2;step++)
		{
			if(my_id == 0)
			{
				for(int m=0;m<world_size;m++)
				{
					tab[m].size[0] = get_row(link_list, (step+m*3))->last->nb_particles;
					tab[m].size[1] = get_row(link_list, (step+m*3))->nb_particles;
					tab[m].size[2] = get_row(link_list, (step+m*3))->next->nb_particles;

					// printf("ID = %d, Last = %d, Current = %d, Next = %d\n",m,tab[m].size[0],tab[m].size[1],tab[m].size[2]);

					for(int i=0;i<tab[m].size[0];i++)
					{
						tab[m].particles_last[i] = get_row(link_list, (step+m*3))->last->particles[i]; 
					}
					for(int i=0;i<tab[m].size[1];i++)
					{
						tab[m].particles_current[i] = get_row(link_list, (step+m*3))->particles[i]; 
					} 
					for(int i=0;i<tab[m].size[2];i++)
					{
						tab[m].particles_next[i] = get_row(link_list, (step+m*3))->next->particles[i]; 
					}

					MPI_Send(tab[m].size,3,MPI_INT,m,1000+m,MPI_COMM_WORLD);
					MPI_Send(tab[m].particles_last,tab[m].size[0],MPI_PARTICLES,m,m+0,MPI_COMM_WORLD);
					MPI_Send(tab[m].particles_current,tab[m].size[1],MPI_PARTICLES,m,m+1,MPI_COMM_WORLD);
					MPI_Send(tab[m].particles_next,tab[m].size[2],MPI_PARTICLES,m,m+2,MPI_COMM_WORLD);
				}
				
			}

			MPI_Recv(size,3,MPI_INT,0,1000+my_id,MPI_COMM_WORLD,&status);
			MPI_Recv(last,size[0],MPI_PARTICLES,0,0+my_id,MPI_COMM_WORLD,&status);
			MPI_Recv(current,size[1],MPI_PARTICLES,0,1+my_id,MPI_COMM_WORLD,&status);
			MPI_Recv(next,size[2],MPI_PARTICLES,0,2+my_id,MPI_COMM_WORLD,&status);

			printf("ID = %d, Last = %d, Current = %d, Next = %d\n",my_id,size[0],size[1],size[2]);

			// unsigned int p, pp;

			// for(p=0; p<size[1]; p++)
			// { 
			// 	if(current[p].collide==true) continue;

			// 	/* check for collisions */
			// 	for(pp=0; pp<size[0]; pp++)
			// 	{
			// 		if(last[pp].collide) continue;
			// 		float t=collide(&current[p], &last[pp]);
			// 		if(t!=-1)
			// 		{ // collision
			// 			current[p].collide = true;
			// 			last[p].collide = true;
			// 			interact(&current[p], &last[pp], t);
			// 			break; // only check collision of two particles
			// 		}
			// 	}

			// 	for(pp=0; pp<size[1]; pp++)
			// 	{
			// 		if(pp!=p)
			// 		{
			// 			if(current[pp].collide == true) continue;
			// 			float t=collide(&current[p], &current[pp]);
			// 			if(t!=-1)
			// 			{ // collision
			// 				current[p].collide = true;
			// 				current[pp].collide = true;
			// 				interact(&current[p], &current[pp], t);
			// 				break; // only check collision of two particles
			// 			}
			// 		}
			// 	}

			// 	for(pp=0; pp<size[2]; pp++)
			// 	{
			// 		if(next[pp].collide) continue;
			// 		float t=collide(&current[p], &next[pp]);
			// 		if(t!=-1)
			// 		{ // collision
			// 			current[p].collide = true;
			// 			next[p].collide = true;
			// 			interact(&current[p], &next[pp], t);
			// 			break; // only check collision of two particles
			// 		}
			// 	}


			// }

			// MPI_Send(size,3,MPI_INT,my_id,1000+my_id,MPI_COMM_WORLD);
			// MPI_Send(last,size[0],MPI_PARTICLES,0,my_id+0,MPI_COMM_WORLD);
			// MPI_Send(current,size[1],MPI_PARTICLES,0,my_id+1,MPI_COMM_WORLD);
			// MPI_Send(next,size[2],MPI_PARTICLES,0,my_id+2,MPI_COMM_WORLD);

			// if(my_id==0)
			// {
			// 	MPI_Recv(tab[my_id].size,3,MPI_INT,my_id,1000+my_id,MPI_COMM_WORLD,&status);
			// 	MPI_Recv(tab[my_id].particles_last,size[0],MPI_PARTICLES,my_id,my_id+0,MPI_COMM_WORLD,&status);
			// 	MPI_Recv(tab[my_id].particles_current,size[1],MPI_PARTICLES,my_id,my_id+1,MPI_COMM_WORLD,&status);
			// 	MPI_Recv(tab[my_id].particles_next,size[2],MPI_PARTICLES,my_id,my_id+2,MPI_COMM_WORLD,&status);
			// }



		}	



	}














	// unsigned int p, pp;

	/* Main loop */
	// for (time_stamp=0; time_stamp<time_max; time_stamp++) { // for each time stamp

	// 	init_collisions(collisions, INIT_NO_PARTICLES);

	// 	for(p=0; p<INIT_NO_PARTICLES; p++) { // for all particles
	// 		if(collisions[p]) continue;

	// 		/* check for collisions */
	// 		for(pp=p+1; pp<INIT_NO_PARTICLES; pp++){
	// 			if(collisions[pp]) continue;
	// 			float t=collide(&particles[p], &particles[pp]);
	// 			if(t!=-1){ // collision
	// 				collisions[p]=collisions[pp]=1;
	// 				interact(&particles[p], &particles[pp], t);
	// 				break; // only check collision of two particles
	// 			}
	// 		}
	// 	}

	// 	// move particles that has not collided with another
	// 	for(p=0; p<INIT_NO_PARTICLES; p++)
	// 		if(!collisions[p]){
	// 			feuler(&particles[p], 1);

	// 			/* check for wall interaction and add the momentum */
	// 			pressure += wall_collide(&particles[p], wall);
	// 		}
	// }

	if(my_id==0)
	{
		//printf("Average pressure = %f\n", pressure / (WALL_LENGTH*time_max));
		printf("\n");
	} 

	

	// free(particles);
	// free(collisions);

	
	MPI_Finalize();
	
	return 0;

}

