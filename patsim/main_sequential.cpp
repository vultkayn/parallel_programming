#include <cstdlib>
#include <ctime>
#include <iostream>
#include <cmath>
#include <array>

#include "coordinate.h"
#include "definitions.h"
#include "physics.h"

// Feel free to change this program to facilitate parallelization.

float rand1()
{
	return (float)(rand() / (float)RAND_MAX);
}

int main(int argc, char **argv)
{
	float pressure = 0;

	// parse arguments
	if (argc != 2)
	{
		std::cerr << "Usage: " << argv[0] << " simulation_time" << std::endl;
		std::cerr << "For example: " << argv[0] << " 10" << std::endl;
		exit(1);
	}

	unsigned time_max{std::atoi(argv[1])};

	/* Initialize */
	// 1. set the walls
	cord_t wall;
	wall.y0 = wall.x0 = 0;
	wall.x1 = BOX_HORIZ_SIZE;
	wall.y1 = BOX_VERT_SIZE;

	// 2. allocate particle bufer and initialize the particles
	std::array<pcord_t, INIT_NO_PARTICLES> particles{};
	std::array<bool, INIT_NO_PARTICLES> collisions{};
	collisions.fill(false);

	std::srand(std::time(NULL) + 1234);

	for (auto &particle : particles)
	{
		// initialize random position
		particle.x = wall.x0 + rand1() * BOX_HORIZ_SIZE;
		particle.y = wall.y0 + rand1() * BOX_VERT_SIZE;

		// initialize random velocity
		float r = rand1() * MAX_INITIAL_VELOCITY;
		float a = rand1() * 2 * PI;
		particle.vx = r * std::cos(a);
		particle.vy = r * std::sin(a);
	}

	/* Main loop */
	for (unsigned time_stamp = 0; time_stamp < time_max; time_stamp++)
	{ // for each time stamp
		collisions.fill(false);

		for (unsigned p = 0; p < INIT_NO_PARTICLES; p++)
		{ // for all particles
			if (collisions[p])
				continue;

			/* check for collisions */
			for (unsigned pp = p + 1; pp < INIT_NO_PARTICLES; pp++)
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
		for (unsigned p = 0; p < INIT_NO_PARTICLES; p++)
			if (!collisions[p])
			{
				feuler(&particles[p], 1);

				/* check for wall interaction and add the momentum */
				pressure += wall_collide(&particles[p], wall);
			}
	}

	std::cout << "Average pressure = " << pressure / (WALL_LENGTH * time_max) << std::endl;

	return 0;
}