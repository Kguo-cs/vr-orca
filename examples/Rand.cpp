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
#include <cstdlib>
#include <iomanip>
#include <Definitions.h>

#if RVO_OUTPUT_TIME_AND_POSITIONS
#include <iostream>
#endif

#if _OPENMP
#include <omp.h>
#endif
#include<fstream> 
#include <RVO.h>

#ifndef M_PI
const float M_PI = 3.14159265358979323846f;
#endif

/* Store the goals of the agents. */
std::vector<RVO::Vector2> starts;
std::vector<RVO::Vector2> goals;

bool fail(RVO::Vector2 s, std::vector<RVO::Vector2> list,float radius){
    for (size_t i=0;i<list.size();++i){
        if (RVO::absSq(s - list[i]) <radius*radius*16)
        {
            return true;
        }
    }
    return false;
}


void setupScenario(RVO::RVOSimulator *sim)
{
	/* Specify the global time step of the simulation. */
	sim->setTimeStep(0.25f);

	float scale=2.5;
    int agent_num = int(250 / scale);
    float radius = 1.5 / scale;
    float map_size = 15;
    float max_speed = 2 / scale;
    float prefer_v = 1 / scale;
    float time_horizon = 10;
    int maxNeighbors_ = agent_num;
    float neighborDist_ = 15 / scale;

	/* Specify the default parameters for agents that are subsequently added. */
	sim->setAgentDefaults(neighborDist_, maxNeighbors_,  time_horizon,time_horizon,radius, max_speed,prefer_v);
    //float neighborDist, size_t maxNeighbors, float timeHorizon, float timeHorizonObst, float radius, float maxSpeed, const Vector2 &velocity
	/*
	 * Add agents, specifying their start position, and store their goals on the
	 * opposite side of the environment.
	 */

	starts.clear();
	goals.clear();


    for (size_t i=0;i<agent_num;++i){
        RVO::Vector2 s=RVO::Vector2(-map_size+2*map_size*(float)std::rand()/RAND_MAX,-map_size+2*map_size*(float)std::rand()/RAND_MAX);

        while(fail(s,starts,radius)){
            s=RVO::Vector2(-map_size+2*map_size*(float)std::rand()/RAND_MAX,-map_size+2*map_size*(float)std::rand()/RAND_MAX);
        }

        RVO::Vector2 g=RVO::Vector2(-map_size+2*map_size*(float)std::rand()/RAND_MAX,-map_size+2*map_size*(float)std::rand()/RAND_MAX);

        while(fail(g,goals,radius)){
            g=RVO::Vector2(-map_size+2*map_size*(float)std::rand()/RAND_MAX,-map_size+2*map_size*(float)std::rand()/RAND_MAX);
        }

        starts.push_back(s);
		goals.push_back(g);
        sim->addAgent(s);
    }
}

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
	std::ofstream outfile;
	outfile.open("/home/guo/Desktop/all2.txt");

	for(int i=20;i<30;i++)
	{
		std::srand(1);

		int step_sum=0;
		float min_sum=0;
		float d_sum=0;
		float e_sum=0;
		double start_cpu = clock();

		for(int n=0;n<1000;n++){

			RVO::RVOSimulator *sim = new RVO::RVOSimulator();

			sim->WEIGHT=(i+1)*10;

			/* Set up the scenario. */
			setupScenario(sim);

			// if(n!=466)
			// {
			// 	delete sim;
			// 	continue;
			// }

			int step=0;
			float min_d=10000;
			float total_d=0;

			// outfile<<"\n";

			/* Perform (and manipulate) the simulation. */
			do {

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
					// outfile<<sim->getAgentPosition(i).x();
					// outfile<<" ";
					// outfile<<sim->getAgentPosition(i).y();
					// outfile<<" ";	
				}
				// std::cout<<std::setw(2)<<step;
				// std::cout<<' ';
				// //std::cout.precision(10);
				// std::cout<<std::setw(15)<<min_d;
				// std::cout<<' ';
				// std::cout<<std::setw(15)<<total_d;
				// std::cout<<' ';
				// std::cout<<std::setw(7)<<(clock()-start_cpu)/CLOCKS_PER_SEC;
				// std::cout<<' ';


				//std::cout<<'\n';

				// outfile<<"\n";
			}
			while (!reachedGoal(sim));

			std::cout<<std::setw(2)<<n;
			std::cout<<' ';
			std::cout<<std::setw(2)<<step;
			std::cout<<' ';
			std::cout<<std::setw(7)<<min_d;
			std::cout<<' ';
			std::cout<<std::setw(7)<<total_d;
			std::cout<<' ';
			std::cout<<std::setw(7)<<(clock()-start_cpu)/CLOCKS_PER_SEC;
			std::cout<<' ';

			std::cout<<'\n';

			step_sum+=step;
			d_sum+=total_d;
			min_sum+=min_d;
			e_sum+=sim->total_e;


			delete sim;
		
		}

		std::cout<<std::setw(2)<<i;
		std::cout<<' ';
		std::cout<<std::setw(2)<<step_sum;
		std::cout<<' ';
		std::cout<<std::setw(7)<<d_sum;
		std::cout<<' ';
		std::cout<<std::setw(7)<<min_sum;
		std::cout<<' ';
		std::cout<<std::setw(7)<<e_sum;
		std::cout<<std::endl;

		outfile<<step_sum;
		outfile<<" ";
		outfile<<d_sum;
		outfile<<" ";
		outfile<<min_sum;
		outfile<<" ";
		outfile<<e_sum;
		outfile<<" ";

		outfile<<"\n";
	}

	return 0;
}


