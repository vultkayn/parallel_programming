#include <cstdlib>
#include <ctime>
#include <iostream>
#include <vector>
#include <unordered_set>
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

void generate_particles(std::list<pcord_t> &row, int rank, int world_size, float &top, float &bot, cord_t const &wall)
{
	float height = BOX_VERT_SIZE / (float)world_size;
	bot = height * (float)rank;
	top = bot + height;
	int i;
	for (i = 0; i < INIT_NO_PARTICLES / world_size; i++)
	{
		// initialize random position
		pcord_t part;
		part.x = wall.x0 + mine::rand1() * BOX_HORIZ_SIZE;
		float dy = mine::rand1() * height;
		part.y = bot + dy;
		// initialize random velocity
		float r = mine::rand1() * MAX_INITIAL_VELOCITY;
		float a = mine::rand1() * 2 * PI;
		part.vx = r * cos(a);
		part.vy = r * sin(a);
		part.collide = 0;
		row.push_back(part);
	}
}

// P = world_size, W = BOX_HORIZ_SIZE, N = INIT_NO_PARTICLES
// return bottom-most band (minus some that did collide)
void collide_within_oneself(std::list<pcord_t> &row,														 // = N / P
														std::vector<pcord_t *> &upper_band,			 // = 0.5 * [50W / (W²/P)] * N / P = 50N / 2W
														std::vector<pcord_t *> &lower_band,			 // = 0.5 * [50W / (W²/P)] * N / P = 50N / 2W
														Row &leaving_top,										 // = 0.5 * [50W / (W²/P)] * N / P = 50N / 2W
														float top_y, float bot_y,
														int rank, int world_size)
/*
TOTAL extra bytes:
= sizeof(int)*(50N/W) + 2*sizeof(pointer)*(50N / 2W) + 2*sizeof(pcord_t) * (50N / 2W)
*/
{
	for (auto it = row.begin(); it != row.end(); )
	{
		auto &part{*it};
		if (part.collide)
			continue;

		/* check for collisions */
		auto iit = it;
		for (iit++; pp < row.size(); pp++)
		{
			if (row[pp].collide or invalids.find(pp) != invalids.end())
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
			if (moved_y <= bot_y + 50.f or part.y <= bot_y + 50.f)
			{
				if (rank > 0)
					lower_band.push_back(p);
			}

			// priority to "moved up"
			else if (moved_y > top_y)
			{
				if (rank < world_size - 1)
				{
					leaving_top.push_back(part);
					invalids.insert(p);
				}
			}
			else if (moved_y >= top_y - 50.f or part.y>= top_y - 50.f)
			{
				if (rank < world_size - 1)
					upper_band.push_back(p);
			}
		}
	}
}

// Phase 1
Row resolve_lower_band_collision(Row &row, Row &from_below, std::vector<int> const &lower_band, std::unordered_set<int> &invalids, float bot_y)
{
	// COLLISION CHECK
	for (pcord_t part : from_below)
	{
		for (int p : lower_band)
		{
			if (row[p].collide)
				continue;
			float t = collide(&part, &row[p]);
			if (t != -1)
			{ // collision
				part.collide = row[p].collide = true;
				interact(&part, &row[p], t);
				break; // only check collision of two particles
			}
		}

		// REGISTER THE PARTICLE
		// ADD RECEIVED PART TO CURRENT ROW
		if (invalids.empty())
		{
			row.push_back(part);
		}
		else
		{
			auto it{invalids.begin()};
			row[*it] = part;
			invalids.erase(it);
		}
	}

	Row leaving_bot{};
	for (int p : lower_band)
	{
		float moved_y = row[p].y + row[p].vy;
		if (not row[p].collide and moved_y < bot_y) // FIXME "invalids" check should be overkill
		{
			leaving_bot.push_back(row[p]);
			invalids.insert(p);
		}
	}
	return leaving_bot;
}

// Phase 2
void resolve_upper_band_collision(Row &row, Row &from_above, std::vector<int> const &upper_band, std::unordered_set<int> &invalids)
{
	// COLLISION CHECK
	for (pcord_t part : from_above)
	{
		for (int p : upper_band)
		{
			if (row[p].collide)
				continue;
			float t = collide(&part, &row[p]);
			if (t != -1)
			{ // collision
				part.collide = row[p].collide = true;
				interact(&part, &row[p], t);
				break; // only check collision of two particles
			}
		}
		// ADD RECEIVED PART TO CURRENT ROW
		if (invalids.empty())
		{
			row.push_back(part);
		}
		else
		{
			auto it{invalids.begin()};
			row[*it] = part;
			invalids.erase(it);
		}
	}
	// END COLLISION CHECK
}

float move_and_collide_with_wall(Row &row, cord_t const &wall, std::unordered_set<int> const &invalids)
{
	float pressure{0};
	// move particles that has not collided with another
	unsigned p;
	for (p = 0; p < row.size(); p++)
		if (not row[p].collide and invalids.find(p) == invalids.end())
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
		time_max = std::atoi(argv[1]);
	}

	MPI_Bcast(&time_max, 1, MPI_FLOAT, 0, MPI_COMM_WORLD); // send time_max to all

	/* Initialize */
	// 1. set the walls
	cord_t wall;
	wall.y0 = wall.x0 = 0;
	wall.x1 = BOX_HORIZ_SIZE;
	wall.y1 = BOX_VERT_SIZE;
	std::srand(1234);

	// 2. allocate particle bufer and initialize the particles

	std::list<pcord_t> row(100 + INIT_NO_PARTICLES / world_size);
	std::unordered_set<int> invalids{};
	float top_y = 0;
	float bot_y = 0;
	float local_pressure{0};
	unsigned time_stamp;
	
	double stime, etime;
	if(my_id==0)
	{
	    stime = MPI_Wtime();
	}
	
	for (time_stamp = 0; time_stamp < time_max; time_stamp++)
	{
		if (time_stamp == 0)
		{
			generate_particles(row, my_id, world_size, top_y, bot_y, wall);
// #ifdef DEBUG
			std::cout << "proc#" << my_id << " (top, bot)=(" << top_y << ',' << bot_y << ')'
								<< " with particules=" << row.size() << std::endl;
// #endif
		}
		else // reset collisions
			for (auto &part : row)
				part.collide = false;
		// vector of pointer is more performant than a list since don't have to allocate at every pointer.
		/* vector of pointer is also better than vector of indices, since it reduces by 1 the number of indirection.
		indeed: indirections for vector<int> v, accessing 'int index = v[i]'
		*/
		std::vector<int> upper_band{};
		Row leaving_bot{};
		{ // Scope leaving_top and lower_band so their storage is released as soon as possible
			Row leaving_top{};
			std::vector<int> lower_band{};
			collide_within_oneself(row, upper_band, lower_band, leaving_top, invalids, top_y, bot_y, my_id, world_size);

			if (my_id < world_size - 1) // send particles leaving to processor above
			{
				unsigned target_rank = my_id + 1;
				int size = leaving_top.size();
				MPI_Request req1, req2;
				MPI_Isend(&size, 1, MPI_INT, target_rank, target_rank + 0, MPI_COMM_WORLD, &req1);
				MPI_Request_free(&req1);
				MPI_Isend(leaving_top.data(), leaving_top.size(), MPI_PARTICLES, target_rank, target_rank + 1, MPI_COMM_WORLD, &req2);
				MPI_Request_free(&req2);
			}
			if (my_id > 0) // recv particles from processor below
			{
				int size;
				MPI_Recv(&size, 1, MPI_INT, my_id - 1, my_id + 0, MPI_COMM_WORLD, &status);
				Row from_below(size);
				MPI_Recv(from_below.data(), size, MPI_PARTICLES, my_id - 1, my_id + 1, MPI_COMM_WORLD, &status);
				leaving_bot = resolve_lower_band_collision(row, from_below, lower_band, invalids, bot_y);
			}
		} // leaving_top and lower_band outlived their usefulness, delete them.

		if (my_id > 0) // send to processor below
		{
			unsigned target_rank = my_id - 1;
			int size = leaving_bot.size();
			MPI_Request req1, req2;
			MPI_Isend(&size, 1, MPI_INT, target_rank, target_rank + 2, MPI_COMM_WORLD, &req1);
			MPI_Request_free(&req1);
			MPI_Isend(leaving_bot.data(), leaving_bot.size(), MPI_PARTICLES, target_rank, target_rank + 3, MPI_COMM_WORLD, &req2);
			MPI_Request_free(&req2);
		}
		if (my_id < world_size - 1) // receive from processor above
		{
			int size;
			MPI_Recv(&size, 1, MPI_INT, my_id + 1, my_id + 2, MPI_COMM_WORLD, &status);
			Row from_above(size);
			MPI_Recv(from_above.data(), size, MPI_PARTICLES, my_id + 1, my_id + 3, MPI_COMM_WORLD, &status);
			resolve_upper_band_collision(row, from_above, upper_band, invalids);
		}

		local_pressure += move_and_collide_with_wall(row, wall, invalids); // move all remaining
	}

	float pressure{0};
	MPI_Reduce(&local_pressure, &pressure, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
	if (my_id == 0)
	{
		std::cout << "Average pressure = " << pressure / (WALL_LENGTH * time_max) << std::endl;
		
		etime = MPI_Wtime();
		std::cout << "Computation Time: " << (etime  - stime) << " secs" << std::endl;
	}
	MPI_Finalize();

	return 0;
}
