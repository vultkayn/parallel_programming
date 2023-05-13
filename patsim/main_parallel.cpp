#include <cstdlib>
#include <ctime>
#include <iostream>
#include <vector>
#include <list>
#include <cmath>
#include "mpi.h"

#include "coordinate.h"
#include "definitions.h"
#include "physics.h"

using Row = std::vector<pcord_t>; // BUG use vector or can smash the stack with array

namespace mine
{
	float rand1()
	{
		return (float)(rand() / (float)RAND_MAX);
	}
}

void generate_particules(Row &row, int rank, int world_size, float &top, float &bot, cord_t const &wall)
{
	float height = BOX_VERT_SIZE / (float)world_size;
	top = height * (float)rank;
	bot = top + height;
	for (int i = 0; i < INIT_NO_PARTICLES; i++)
	{
		// initialize random position
		pcord_t part;
		part.x = wall.x0 + mine::rand1() * BOX_HORIZ_SIZE;
		float dy = mine::rand1() * height;
		part.y = top + dy;
		// initialize random velocity
		float r = mine::rand1() * MAX_INITIAL_VELOCITY;
		float a = mine::rand1() * 2 * PI;
		part.vx = r * cos(a);
		part.vy = r * sin(a);
		row.push_back(part);

#ifdef DEBUG
		for (int i = 0; i < row_nb; i++)
		{
			std::cerr << "Row " << i + 1 << ":" << std::endl;
			std::cerr << tab[i + 1] << std::endl;
		}
#endif
	}
}

// P = world_size, W = BOX_HORIZ_SIZE, N = INIT_NO_PARTICLES
// return bottom-most band (minus some that did collide)
void collide_within_oneself(Row &row,													// = N / P
														std::vector<pcord_t *> &top_band, // = 0.5 * [50W / (W²/P)] * N / P = 50N / 2W
														std::vector<pcord_t *> &bot_band, // = 0.5 * [50W / (W²/P)] * N / P = 50N / 2W
														Row &leaving_top,									// = 0.5 * [50W / (W²/P)] * N / P = 50N / 2W
														Row &leaving_bot,									// = 0.5 * [50W / (W²/P)] * N / P = 50N / 2W
														std::vector<short> &invalids,			// = nb_elems(leaving_top + leaving_bot) = 50N / W
														float box_top_y, float box_bottom_y)
/*
TOTAL extra bytes:
= sizeof(short)*(50N/W) + 2*sizeof(pointer)*(50N / 2W) + 2*sizeof(pcord_t) * (50N / 2W)
*/
{
	for (int p = 0; p < row.size(); p++)
	{
		auto &part{row[p]};
		if (part.collide)
			continue;

		/* check for collisions */
		for (unsigned pp = p + 1; pp < row.size(); pp++)
		{
			if (row[pp].collide)
				continue;
			float t = collide(&part, &row[pp]);
			if (t != -1)
			{ // collision
				part.collide = row[pp].collide = true;
				interact(&part, &row[pp], t);
				break; // only check collision of two particles
			}
		}
		if (not part.collide)
		{
			float moved_y = part.y + part.vy;
			if (part.y >= box_bottom_y - 50.f)
			{
				if (moved_y > box_bottom_y)
				{
					leaving_bot.push_back(part);
					invalids.push_back(p);
				}
				else
					bot_band.push_back(&row[p]);
			}
			else if (part.y <= (box_top_y + 50.f))
			{
				if (moved_y < box_top_y)
				{
					leaving_top.push_back(part);
					invalids.push_back(p);
				}
				else
					top_band.push_back(&row[p]);
			}
		}
	}
}

void resolve_interbands_collision(Row &row, Row &from_other, std::vector<pcord_t *> band, std::vector<short> &invalids)
{
	int ninvalids = invalids.size();
	for (unsigned pp = 0; pp < from_other.size(); pp++)
	{
		for (pcord_t *part : band)
		{
			if (part->collide)
				continue;
			float t = collide(part, &from_other[pp]);
			if (t != -1)
			{ // collision
				part->collide = from_other[pp].collide = true;
				interact(part, &from_other[pp], t);
				break; // only check collision of two particles
			}
		}
		// insert incoming particle into current processor row.
		if (ninvalids)
		{
			row[invalids.back()] = from_other[pp];
			invalids.pop_back();
			ninvalids--;
		}
		else // No invalid room left
			row.push_back(from_other[pp]);
	}
}

float move_and_collide_with_wall(Row &row, cord_t const &wall)
{
	float pressure{0};
	// move particles that has not collided with another
	for (unsigned p = 0; p < row.size(); p++)
		if (not row[p].collide)
		{
			feuler(&row[p], 1);

			/* check for wall interaction and add the momentum */
			pressure += wall_collide(&row[p], wall);
		}
	return pressure;
}

std::ostream &operator<<(std::ostream &os, Row &row)
{
	for (auto &particle : row)
	{
		os << "(" << particle.x << "," << particle.y << "," << particle.vx << "," << particle.vy << ")  ";
	}
	return os;
}

int main(int argc, char **argv)
{

	// MPI initialization.
	MPI_Init(NULL, NULL);

	MPI_Status status;

	int BlockLengths[] = {
			4,
			1,
			1 // for alignment
	};
	MPI_Datatype Types[] = {
			MPI_FLOAT,
			MPI_CHAR,
			MPI_UB // upper bound
	};

	MPI_Aint Displacements[] = {
			offsetof(pcord_t, x),
			offsetof(pcord_t, y),
			offsetof(pcord_t, vx),
			offsetof(pcord_t, vy),
			sizeof(pcord_t)};

	MPI_Datatype MPI_PARTICLES;
	MPI_Type_create_struct(3, BlockLengths, Displacements, Types, &MPI_PARTICLES);
	MPI_Type_commit(&MPI_PARTICLES);

	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	int my_id;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_id);

	unsigned time_max = 0;
	if (my_id == 0)
	{
		std::cout << "World size = " << world_size << std::endl;
		// Parse arguments
		if (argc != 2)
		{
			std::cerr << "Usage: " << argv[0] << " simulation_time" << std::endl;
			std::cerr << "For example: " << argv[0] << " 10" << std::endl;
			exit(1);
		}
		time_max = std::atoi(argv[1]); // FIXME broadcast time_max.
	}

	MPI_Bcast(&time_max, 1, MPI_FLOAT, 0, MPI_COMM_WORLD); // send time_max to all

	/* Initialize */
	// 1. set the walls
	cord_t wall;
	wall.y0 = wall.x0 = 0;
	wall.x1 = BOX_HORIZ_SIZE;
	wall.y1 = BOX_VERT_SIZE;
	std::srand(std::time(NULL) + 1234);

	// 2. allocate particle bufer and initialize the particles

	Row row{};
	std::vector<short> invalids{};
	float top_y = 0;
	float bot_y = 0;
	float local_pressure{0};

	for (unsigned time_stamp = 0; time_stamp < time_max; time_stamp++)
	{
		if (time_stamp == 0)
		{
			generate_particules(row, my_id, world_size, top_y, bot_y, wall);
		}

		// vector of pointer is more performant than a list since don't have to allocate at every pointer.
		/* vector of pointer is also better than vector of indices, since it reduces by 1 the number of indirection.
		indeed: indirections for vector<short> v, accessing 'int index = v[i]'
		*/
		std::vector<pcord_t *> top_band{};
		std::vector<pcord_t *> bot_band{};
		{ // Scope leaving_top and leaving_bot so their storage is released as soon as possible
			Row leaving_bot{};
			Row leaving_top{};
			collide_within_oneself(row, top_band, bot_band, leaving_top, leaving_bot, invalids, top_y, bot_y);

			if (my_id < world_size - 1) // send particles leaving to processor below
			{
				unsigned target_rank = my_id + 1;
				unsigned short size = leaving_bot.size();
				MPI_Send(&size, 1, MPI_UNSIGNED_SHORT, target_rank, target_rank + 0, MPI_COMM_WORLD);
				MPI_Send(leaving_bot.data(), leaving_bot.size(), MPI_PARTICLES, target_rank, target_rank + 1, MPI_COMM_WORLD);
			}
			if (my_id > 0) // send particles leaving to processor above
			{
				unsigned target_rank = my_id - 1;
				unsigned short size = leaving_top.size();
				MPI_Send(&size, 1, MPI_UNSIGNED_SHORT, target_rank, target_rank + 0, MPI_COMM_WORLD);
				MPI_Send(leaving_top.data(), leaving_top.size(), MPI_PARTICLES, target_rank, target_rank + 1, MPI_COMM_WORLD);
			}
		} // leaving_bot and leaving_top outlived their usefulness, delete them.

		if (my_id > 0) // receive from processor above
		{
			unsigned short size;
			MPI_Recv(&size, 1, MPI_UNSIGNED_SHORT, my_id - 1, my_id + 0, MPI_COMM_WORLD, &status);
			Row from_above(size);
			MPI_Recv(from_above.data(), size, MPI_PARTICLES, my_id - 1, my_id + 1, MPI_COMM_WORLD, &status);
			resolve_interbands_collision(row, from_above, top_band, invalids);
		}
		if (my_id < world_size - 1) // receive from processor below
		{
			unsigned short size;
			MPI_Recv(&size, 1, MPI_UNSIGNED_SHORT, my_id - 1, my_id + 0, MPI_COMM_WORLD, &status);
			Row from_below(size);
			MPI_Recv(from_below.data(), size, MPI_PARTICLES, my_id - 1, my_id + 1, MPI_COMM_WORLD, &status);
			resolve_interbands_collision(row, from_below, bot_band, invalids);
		}

		local_pressure += move_and_collide_with_wall(row, wall); // move all remaining
	}

	float pressure;
	MPI_Reduce(&local_pressure, &pressure, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
	if (my_id == 0)
		std::cout << "Average pressure = " << pressure / (WALL_LENGTH * time_max) << '\n'
							<< std::endl;

	MPI_Finalize();

	return 0;
}
