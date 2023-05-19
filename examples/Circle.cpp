/*
 * Circle.cpp
 * RVO2 Library
 *
 * Copyright 2008 University of North Carolina at Chapel Hill
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 * Please send all bug reports to <geom@cs.unc.edu>.
 *
 * The authors may be contacted via:
 *
 * Jur van den Berg, Stephen J. Guy, Jamie Snape, Ming C. Lin, Dinesh Manocha
 * Dept. of Computer Science
 * 201 S. Columbia St.
 * Frederick P. Brooks, Jr. Computer Science Bldg.
 * Chapel Hill, N.C. 27599-3175
 * United States of America
 *
 * <http://gamma.cs.unc.edu/RVO2/>
 */

/*
 * Example file showing a demo with 250 agents initially positioned evenly
 * distributed on a circle attempting to move to the antipodal position on the
 * circle.
 */

#ifndef RVO_OUTPUT_TIME_AND_POSITIONS
#define RVO_OUTPUT_TIME_AND_POSITIONS 1
#endif

#include <cmath>
#include <cstddef>
#include <vector>

#if RVO_OUTPUT_TIME_AND_POSITIONS
#include <iostream>
#endif
#include <iomanip>

#if _OPENMP
#include <omp.h>
#endif

#include <RVO.h>
#include<Definitions.h>
#include<time.h>


#ifndef M_PI
const float M_PI = 3.14159265358979323846f;
#endif

/* Store the goals of the agents. */
std::vector<RVO::Vector2> goals;

void setupScenario(RVO::RVOSimulator *sim)
{
	/* Specify the global time step of the simulation. */
	sim->setTimeStep(0.25f);

	float scale=2.5;
    int agent_num = int(250 / scale);
    float radius = 1.5 / scale;
    float map_size = 200 / scale;
    float max_speed = 2 / scale;
    float prefer_v = 1 / scale;
    float time_horizon = 16;
    int maxNeighbors_ = agent_num;
    float neighborDist_ = 15 / scale;


	/* Specify the default parameters for agents that are subsequently added. */
	sim->setAgentDefaults(neighborDist_, maxNeighbors_, time_horizon, time_horizon,radius, max_speed,prefer_v);

	/*
	 * Add agents, specifying their start position, and store their goals on the
	 * opposite side of the environment.
	 */
	for (size_t i = 0; i < agent_num; ++i) {
		sim->addAgent(map_size *
		              RVO::Vector2(std::cos(i * 2.0f * M_PI / agent_num),
		                           std::sin(i * 2.0f * M_PI / agent_num)));
		goals.push_back(-sim->getAgentPosition(i));
	}
}

#if RVO_OUTPUT_TIME_AND_POSITIONS
void updateVisualization(RVO::RVOSimulator *sim)
{
	/* Output the current global time. */
	std::cout << sim->getGlobalTime();

	/* Output the current position of all the agents. */
	for (size_t i = 0; i < sim->getNumAgents(); ++i) {
		std::cout << " " << sim->getAgentPosition(i);
	}

	std::cout << std::endl;
}
#endif

void setPreferredVelocities(RVO::RVOSimulator *sim)
{
	/*
	 * Set the preferred velocity to be a vector of unit magnitude (speed) in the
	 * direction of the goal.
	 */
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (int i = 0; i < static_cast<int>(sim->getNumAgents()); ++i) {
		RVO::Vector2 goalVector = (goals[i] - sim->getAgentPosition(i))/sim->getTimeStep();

		if (RVO::absSq(goalVector) >RVO::sqr(sim->getAgentPrefSpeed(i))) {
			goalVector =sim->getAgentPrefSpeed(i)* RVO::normalize(goalVector);
		}

		sim->setAgentPrefVelocity(i, goalVector);
	}
}

bool reachedGoal(RVO::RVOSimulator *sim)
{
	/* Check if all agents have reached their goals. */
	for (size_t i = 0; i < sim->getNumAgents(); ++i) {
		if (RVO::absSq(sim->getAgentPosition(i) - goals[i]) > 0) {
			return false;
		}
	}

	return true;
}

int main()
{
	/* Create a new simulator instance. */
	RVO::RVOSimulator *sim = new RVO::RVOSimulator();

	/* Set up the scenario. */
	setupScenario(sim);
	int step=0;
	float min_d=10000;
	float total_d=0;
	float v_d=0;

	double start_cpu = clock();


	/* Perform (and manipulate) the simulation. */
	do {
#if RVO_OUTPUT_TIME_AND_POSITIONS
		//updateVisualization(sim);
#endif
		setPreferredVelocities(sim);
		sim->doStep();
		step+=1;

		for (size_t i = 0; i < sim->getNumAgents(); ++i) {
			for (size_t j = i+1; j < sim->getNumAgents(); ++j){
				if (RVO::abs(sim->getAgentPosition(i) - sim->getAgentPosition(j))<min_d){
					min_d=RVO::abs(sim->getAgentPosition(i) - sim->getAgentPosition(j));
				}
			}
			total_d+=abs(sim->getAgentVelocity(i))*sim->getTimeStep();
			v_d+=std::abs(abs(sim->getAgentVelocity(i))-abs(sim->getAgentPrefVelocity(i)));

		}

		std::cout<<std::setw(2)<<step;
		std::cout<<' ';
		std::cout<<std::setw(7)<<min_d;
		std::cout<<' ';
		std::cout<<std::setw(7)<<total_d;
		std::cout<<' ';
		std::cout<<std::setw(7)<<v_d/step;
		std::cout<<' ';
		// std::cout<<std::setw(7)<<(clock()-start_cpu)/CLOCKS_PER_SEC;
		// std::cout<<' ';


		std::cout<<'\n';

	}
	while (!reachedGoal(sim));

	delete sim;

	return 0;
}

//2029 0.976279 16633.5 468.055 2.13553 
//gamma=0.001
//infeasible=5  2029 1.03195 16639.9 178.135 
//infeasible=10 1948 1.03935 16577.5  177.17 after improve dist//1948 1.03935 16577.5 154.997 
//infeasible=20 2025 1.03969 16659.3 199.029 
//gamma=0.01
//infeasible=10 1966  1.0774 16642.9 212.607 after improve dist 1966  1.0774 16642.9  1877.27 186.49 
//add the proj_dist when infeasible
//gamma=0.01 infeasible=10 2033 0.956561 16690.1  1315.8 220.948 
//max_neighbor 10 1965 1.12105 16585.9 3.02143 