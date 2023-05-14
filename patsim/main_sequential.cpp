#include <cstdlib>
#include <ctime>
#include <iostream>
#include <cmath>
#include <vector>

#include "coordinate.h"
#include "definitions.h"
#include "physics.h"

// Feel free to change this program to facilitate parallelization.

namespace mine
{
	float rand1()
	{
		return (float)(rand() / (float)RAND_MAX);
	}
}

pcord_t generate_particules(int rank, int world_size, float &top, float &bot, cord_t const &wall)
{
	float height = BOX_VERT_SIZE / (float)world_size;
	top = height * (float)rank;
	bot = top + height;
	int i;
	for (i = 0; i < INIT_NO_PARTICLES; i++)
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
		part.collide = 0;
		return part;
	}
}

int main(int argc, char **argv)
{
	float pressure = 0;

	// parse arguments
	if (argc != 3)
	{
		std::cerr << "Usage: " << argv[0] << " simulation_time nb_procs" << std::endl;
		std::cerr << "For example: " << argv[0] << " 10 64" << std::endl;
		exit(1);
	}

	unsigned time_max{std::atoi(argv[1])};

	/* Initialize */
	// 1. set the walls
	cord_t wall;
	wall.y0 = wall.x0 = 0;
	wall.x1 = BOX_HORIZ_SIZE;
	wall.y1 = BOX_VERT_SIZE;

	const long TOTAL_NO_PARTICLES = INIT_NO_PARTICLES * std::atoi(argv[2]);
	// 2. allocate particle bufer and initialize the particles
	// std::vector<pcord_t> particles(TOTAL_NO_PARTICLES);
	std::vector<pcord_t> particles;
	std::vector<bool> collisions(TOTAL_NO_PARTICLES, false);

	std::srand(std::time(NULL) + 1234);
	
	struct timespec stime, etime;
    	clock_gettime(CLOCK_REALTIME, &stime);
	
	float top_y = 0;
	float bot_y = 0;

	for (int my_id=0;my_id<std::atoi(argv[2]);my_id++)
	{
	    for (int i=0;i<INIT_NO_PARTICLES;i++)
	    {
	        particles.push_back(generate_particules(my_id, std::atoi(argv[2]), top_y, bot_y, wall));
	    }
	}

//	for (auto &particle : particles)
//	{
//		// initialize random position
//		particle.x = wall.x0 + rand1() * BOX_HORIZ_SIZE;
//		particle.y = wall.y0 + rand1() * BOX_VERT_SIZE;
//
//		// initialize random velocity
//		float r = rand1() * MAX_INITIAL_VELOCITY;
//		float a = rand1() * 2 * PI;
//		particle.vx = r * std::cos(a);
//		particle.vy = r * std::sin(a);
//	}

	/* Main loop */
	for (unsigned time_stamp = 0; time_stamp < time_max; time_stamp++)
	{ // for each time stamp
		for (size_t i=0; i < collisions.size(); i++)
			collisions[i] = false;

		for (unsigned p = 0; p < TOTAL_NO_PARTICLES; p++)
		{ // for all particles
			if (collisions[p])
				continue;

			/* check for collisions */
			for (unsigned pp = p + 1; pp < TOTAL_NO_PARTICLES; pp++)
			{
				if (collisions[pp])
					continue;
				float t = collide(&particles[p], &particles[pp]);
				if (t != -1)
				{ // collision
					collisions[p] = collisions[pp] = 1;
					interact(&particles[p], &particles[pp], t);
					break; // only check collision of two particles
				}
			}
		}

		// move particles that has not collided with another
		for (unsigned p = 0; p < TOTAL_NO_PARTICLES; p++)
			if (!collisions[p])
			{
				feuler(&particles[p], 1);

				/* check for wall interaction and add the momentum */
				pressure += wall_collide(&particles[p], wall);
			}
	}

	std::cout << "Average pressure = " << pressure / (WALL_LENGTH * time_max) << std::endl;
	
	clock_gettime(CLOCK_REALTIME, &etime);
    	std::cout << "Computation Time: " << (etime.tv_sec  - stime.tv_sec) + 1e-9*(etime.tv_nsec - stime.tv_nsec) << " secs" << std::endl;

	return 0;
}
