/*
 * Agent.cpp
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

#include "Agent.h"

#include "KdTree.h"
#include "Obstacle.h"
#include <iomanip>
#include <boost/boost/functional/hash.hpp>

namespace RVO {
	Agent::Agent(RVOSimulator *sim) : maxNeighbors_(0), maxSpeed_(0.0f), neighborDist_(0.0f), radius_(0.0f), sim_(sim), timeHorizon_(0.0f), timeHorizonObst_(0.0f), id_(0) { }

	void Agent::computeNeighbors()
	{

		agentNeighbors_.clear();

		if (maxNeighbors_ > 0) {
			float rangeSq = sqr(neighborDist_);
			sim_->kdTree_->computeAgentNeighbors(this, rangeSq);
		}
	}

	void Agent::computeORCALines(){
		initialNs_.clear();
		initialLinePoints_.clear();

		const float invTimeHorizon = 1.0f / timeHorizon_;

		/* Create agent ORCA lines. */
		for (size_t i = 0; i < agentNeighbors_.size(); ++i) {
			const Agent *const other = agentNeighbors_[i];

			const Vector2 relativePosition = other->position_ - position_;
			const Vector2 relativeVelocity = velocity_ - other->velocity_;
			const float distSq = absSq(relativePosition);
			const float combinedRadius = radius_ + other->radius_;
			const float combinedRadiusSq = sqr(combinedRadius);

			Vector2 n;
			Vector2 u;

			if (distSq > combinedRadiusSq) {
				/* No collision. */
				const Vector2 w = relativeVelocity - invTimeHorizon * relativePosition;
				/* Vector from cutoff center to relative velocity. */
				const float wLengthSq = absSq(w);

				const float dotProduct1 = w * relativePosition;

				if (dotProduct1 < 0.0f && sqr(dotProduct1) > combinedRadiusSq * wLengthSq) {
					/* Project on cut-off circle. */
					const float wLength = std::sqrt(wLengthSq);
					const Vector2 unitW = w / wLength;

					n = unitW;
					u = (combinedRadius * invTimeHorizon - wLength) * unitW;
				}
				else {
					/* Project on legs. */
					const float leg = std::sqrt(distSq - combinedRadiusSq);

					if (det(relativePosition, w) > 0.0f) {
						/* Project on left leg. */
						n = Vector2(relativePosition.x() * leg - relativePosition.y() * combinedRadius, relativePosition.x() * combinedRadius + relativePosition.y() * leg) / distSq;
					}
					else {
						/* Project on right leg. */
						n = -Vector2(relativePosition.x() * leg + relativePosition.y() * combinedRadius, -relativePosition.x() * combinedRadius + relativePosition.y() * leg) / distSq;
					}

					const float dotProduct2 = relativeVelocity * n;

					u = dotProduct2 * n - relativeVelocity;

					n=-perp(n);
				}
			}
			else {
				/* Collision. Project on cut-off circle of time timeStep. */
				const float invTimeStep = 1.0f / sim_->timeStep_;

				/* Vector from cutoff center to relative velocity. */
				const Vector2 w = relativeVelocity - invTimeStep * relativePosition;

				const float wLength = abs(w);
				const Vector2 unitW = w / wLength;

				n = unitW;
				u = (combinedRadius * invTimeStep - wLength) * unitW;
			}

			initialNs_.push_back(n);
			initialLinePoints_.push_back(velocity_ +0.5*u);
		}

		dict_.clear();

	}

	void Agent::circleDist(Vector2& linePoint, Vector2& lineDirection,float& leftDist, float& rightDist){

		float dotProduct = linePoint*perp(lineDirection);
		float discriminant = dotProduct * dotProduct + sqr(maxSpeed_) - linePoint*linePoint;
		float sqrtDiscriminant;

		if(discriminant < 0.0) // * Max speed circle fully invalidates line lineNo.* /-1.2397904931091386 6.136163913165748
			if(lineDirection*linePoint > 0)
			{	
				//("Max speed circle fully invalidates")
				sqrtDiscriminant = -INFINITY;
			}
			else
			{	//if the full circle satisfy, it will continue before.

				sqrtDiscriminant = INFINITY;
			}
		else{
			sqrtDiscriminant = sqrt(discriminant);
		}

		leftDist=-dotProduct - sqrtDiscriminant;
		rightDist=-dotProduct + sqrtDiscriminant;
	}

	void Agent::optimizeORCALines(Agent*& collider){

		float bestTotalProjDist = INFINITY;

		resetORCALines(collider->position_);
		computeColliderNum(collider->id_);
		collider->computeColliderNum(id_);

		for (int step=0;step<most_step;step++){
			computeGradient();

			if(totalProjDist>=bestTotalProjDist || totalProjDist<RVO_EPSILON )
			{
				break;
			}

			if(totalProjDist<bestTotalProjDist)
			{
				bestPoint_= linePoints_[colliderNum_];
				collider->bestPoint_=collider->linePoints_[collider->colliderNum_];
				bestTotalProjDist=totalProjDist;
			}

			for(int i_=0;i_<commonVisibles_.size();i_++)
			{
				Agent*& agent=commonVisibles_[i_];
				agent->changed_=false;

				for(int j_=0;j_<commonVisibles_.size();j_++){//agent[i_]=neighbor neighbor[j_]=agent
					
					if(i_==j_){
						continue;
					}
					Agent*& neighbor=commonVisibles_[j_];
					int i=agent->myNum_[j_] ;//i=0
					int j=neighbor->myNum_[i_] ;//

					float d_n = -GAMMA * ( agent->gradients_[i] -neighbor->gradients_[j]);

					if (d_n != 0){
						agent->changed_=true;
						agent->linePoints_[i] += d_n* agent->ns_[i];
					}
				}
			} 
		}

		collider->dict_[id_]=collider->bestPoint_;
		bestLinePoints_.push_back(bestPoint_);
	}

	void Agent::resetORCALines(Vector2& colliderPosition)
	{	

		commonVisibles_.clear();

		commonVisibles_.push_back(this);

		for(int i=0;i<agentNeighbors_.size();i++)
		{
			Agent* & agent=agentNeighbors_[i];

			if (absSq(position_+colliderPosition-2*agent->position_)<sqr(neighborDist_))
			// if(absSq(position_-agent->position_)<sqr(neighborDist_) && absSq(colliderPosition-agent->position_)<sqr(neighborDist_))
			{
				commonVisibles_.push_back(agent);
			}
		}
		
		for(int i=0;i<commonVisibles_.size();i++)
		{
			Agent* & agent=commonVisibles_[i];

			agent->linePoints_.clear();
			agent->ns_.clear();
			agent->visibleNeighbors_.clear();
			agent->orcaLines_.clear();

			for(int j=0;j<agent->agentNeighbors_.size();j++){
				if(absSq(agent->agentNeighbors_[j]->position_-colliderPosition)<sqr(neighborDist_) &&  
					absSq(agent->agentNeighbors_[j]->position_-position_)<sqr(neighborDist_)){
					agent->linePoints_.push_back(agent->initialLinePoints_[j]);
					agent->ns_.push_back(agent->initialNs_[j]);
					agent->visibleNeighbors_.push_back(j);
				
					Line line;

					line.point = agent->initialLinePoints_[j];
					line.direction=perp(agent->initialNs_[j]);
					agent->orcaLines_.push_back(line);

				}
			}

			agent->myNum_.clear();

			for(int j=0;j<commonVisibles_.size();j++){//agent[i_]=neighbor neighbor[j_]=agent
				if(i==j){
					agent->myNum_.push_back(-1);
					continue;
				}

				for(int k=0;k<agent->visibleNeighbors_.size();k++){
					if(agent->agentNeighbors_[agent->visibleNeighbors_[k]]==commonVisibles_[j]){
						agent->myNum_.push_back(k);
						break;
					}
				}
			}

			agent->changed_=true;

			agent->neighborNum_=agent->visibleNeighbors_.size();
		}

	}

	void Agent::computeColliderNum(int id){

		for(int i=0;i<visibleNeighbors_.size();i++)
		{
			if(agentNeighbors_[visibleNeighbors_[i]]->id_==id){
				colliderNum_=i;
				bestPoint_= linePoints_[colliderNum_];
				return;	
			}
		}
	}

	void Agent::computeGradient(){
		totalProjDist=0;
		
		for(int i=0;i<commonVisibles_.size();i++){
			if(commonVisibles_[i]->changed_==true){
				totalProjDist+=commonVisibles_[i]->computeAgentGradient();
			}
		}
	}

	float Agent::computeAgentGradient(){
		point = velocity_;
		gradients_.clear();

		for(int i=0;i<neighborNum_;i++){
			gradients_.push_back(0);
		}

		projNum = -1;

		for(int i=0;i<neighborNum_;i++){
			Vector2& line_direction = ns_[i];
			Vector2& line_point=linePoints_[i];

			if ((point - line_point)* line_direction >= 0){
				continue;
			}

			int leftNum = i;
			int rightNum =i;

			circleDist(line_point,line_direction,left_dist,right_dist);

			for(int j=0;j<i;j++){
				float sin=det(line_direction,ns_[j]);
				float dist=ns_[j]*( line_point - linePoints_[j])/ sin;

				if(sin>0){
					if(dist<right_dist)
					{
						rightNum=j;
						right_dist=dist;
					}
				}
				else
				{
					if(dist>left_dist)
					{
						leftNum=j;
						left_dist=dist;
					}
				}
			}

			if (left_dist > right_dist)//# all infeasible gap or one infeasible gap ,outside the circle
			{
				for (size_t j= 0;j < neighborNum_; ++j) {
					orcaLines_[j].point = linePoints_[j];
				}
				
				linearProgram3(orcaLines_, 0, i, maxSpeed_, point);
				infeasibleGradient();

				return abs(velocity_ - point)+sim_->WEIGHT*distance*distance;
			}

			if (left_dist>0){
				point = line_point + perp(line_direction) * left_dist;
				intersectNum = leftNum;
			}
			else if(right_dist<0){
				point = line_point + perp(line_direction) * right_dist;
            	intersectNum = rightNum;
			}
			else{
				point = line_point;
				intersectNum = -1;		
			}
		
        	projNum = i;
		}

		computeDistGradient();

		return abs(velocity_ - point);

	}

	void Agent::infeasibleGradient(){

		if(infeasibleLine1==-2){//only constrained by one line and the radius
			gradients_[infeasibleLine0]=sim_->WEIGHT*distance*2;
			return ;
		}

		if(infeasibleLine2==-2){
			gradients_[infeasibleLine0]= sim_->WEIGHT*distance*2;
			gradients_[infeasibleLine1]= sim_->WEIGHT*distance*2;
		}	
		else if (infeasibleLine2==-1)
		{
			gradients_[infeasibleLine0]=infeasibleGradient3(infeasibleLine0,infeasibleLine1,infeasibleLine2);
			gradients_[infeasibleLine1]=infeasibleGradient3(infeasibleLine1,infeasibleLine0,infeasibleLine2);

		}
		else
		{
			gradients_[infeasibleLine0]=infeasibleGradient3(infeasibleLine0,infeasibleLine1,infeasibleLine2);
			gradients_[infeasibleLine1]=infeasibleGradient3(infeasibleLine1,infeasibleLine0,infeasibleLine2);
			gradients_[infeasibleLine2]=infeasibleGradient3(infeasibleLine2,infeasibleLine0,infeasibleLine1);

		}
	}

	float  Agent::infeasibleGradient3(int i,int j,int k){

		Vector2& direction=ns_[i];

		Vector2& infeasibleDirection1=ns_[j];
		float cos1=infeasibleDirection1*direction;
		float sin1=det(direction,infeasibleDirection1);

		if (k==-1)
		{
			Vector2 infeasibleDirection2=-normalize(point);

			float cos2=infeasibleDirection2*direction;

			float sin2=det(direction,infeasibleDirection2);
			
			float base=sin1*(sin2*sin1+cos2+cos2*cos1);

			if(base==0){
				return sim_->WEIGHT*distance*2;
			}

			dx=(1+cos1)*cos2/base;//base==0??

			dy=(1+cos1)*sin2/base;
		
		}
		else
		{
			Vector2 infeasibleDirection2=ns_[k];

			float cos2=infeasibleDirection2*direction;

			float sin2=det(direction,infeasibleDirection2);

			float base=sin1*(1-cos2)-sin2*(1-cos1);//base=0??
			
			if(base==0){
				throw std::runtime_error("zero2");
			}

			dx=(cos1-cos2)/base;
			dy=(sin1-sin2)/base;
		
		}

		if(velocity_==point){
			return sim_->WEIGHT*(1-dy)*distance*2+abs(Vector2(dy*direction.x()+dx*direction.y(),dy*direction.y()-dx*direction.x()));
		}

		return sim_->WEIGHT*(1-dy)*distance*2-Vector2(dy*direction.x()+dx*direction.y(),dy*direction.y()-dx*direction.x())*normalize(velocity_-point);
	}

	void Agent::computeDistGradient(){

		if (projNum != -1)
		{
			if (intersectNum == -1){
				gradients_[projNum]= 1;
			}
			else if (velocity_!=point){
				Vector2& direction = ns_[projNum];
				// if(velocity_==point){

				// 	gradients_[projNum] = 1/std::abs(det(
				// 	ns_[intersectNum],
				// 	direction));
				// 	gradients_[intersectNum] = gradients_[projNum];
				// 	return ;
				// }

				Vector2 point_to_opt = normalize(velocity_ - point);

				if (intersectNum == projNum)//point on the circle
				{
					intersect_direction = -normalize(point);
				}
				else
				{
					intersect_direction = ns_[intersectNum];

					gradients_[intersectNum] = det(point_to_opt, direction)/ det(direction,intersect_direction);
				}
				

				gradients_[projNum] = det(point_to_opt, intersect_direction) / -det(direction,intersect_direction);
			}
		}
	}

	void Agent::computeNewVelocity()
	{
		orcaLines_.clear();

		/* Create agent ORCA lines. */
		for (size_t i = 0; i < agentNeighbors_.size(); ++i) {

			Line line;

			line.point = bestLinePoints_[i];
			line.direction=perp(initialNs_[i]);
			orcaLines_.push_back(line);
		}

		size_t lineFail = linearProgram2(orcaLines_, maxSpeed_, prefVelocity_, false, newVelocity_);

		if (lineFail < orcaLines_.size()) {
			linearProgram3(orcaLines_, 0, lineFail, maxSpeed_, newVelocity_);
		}
		
	}

	void Agent::insertAgentNeighbor( Agent* agent, float &rangeSq)
	{
		if (this != agent) {
			const float distSq = absSq(position_ - agent->position_);

			if (distSq < rangeSq) {
				// if (agentNeighbors_.size() < maxNeighbors_) {
				agentNeighbors_.push_back(agent);
				// }

				//visibleAgents_.push_back(agent);

				// size_t i = agentNeighbors_.size() - 1;

				// while (i != 0 && distSq < agentNeighbors_[i - 1].first) {
				// 	agentNeighbors_[i] = agentNeighbors_[i - 1];
				// 	--i;
				// }

				// agentNeighbors_[i] = std::make_pair(distSq, agent);

				// if (agentNeighbors_.size() == maxNeighbors_) {
				// 	rangeSq = agentNeighbors_.back().first;
				// }
			}
		}
	}

	void Agent::update()
	{
		sim_->total_e  +=abs(newVelocity_-velocity_);
		velocity_ = newVelocity_;
		position_ += velocity_ * sim_->timeStep_;
	}

	bool Agent::linearProgram1(const std::vector<Line> &lines, size_t lineNo, float radius, const Vector2 &optVelocity, bool directionOpt, Vector2 &result)
	{
		const float dotProduct = lines[lineNo].point * lines[lineNo].direction;
		const float discriminant = sqr(dotProduct) + sqr(radius) - absSq(lines[lineNo].point);

		infeasibleLine2=-2;

		if (discriminant <= 0.0f) {
			/* Max speed circle fully invalidates line lineNo. */
			return false;
		}

		const float sqrtDiscriminant = std::sqrt(discriminant);
		float tLeft = -dotProduct - sqrtDiscriminant;
		float tRight = -dotProduct + sqrtDiscriminant;

		int right_line=-1;
		int left_line=-1;

		for (size_t i = 0; i < lineNo; ++i) {
			const float denominator = det(lines[lineNo].direction, lines[i].direction);
			const float numerator = det(lines[i].direction, lines[lineNo].point - lines[i].point);

			if (std::fabs(denominator) <= RVO_EPSILON) {
				/* Lines lineNo and i are (almost) parallel. */
				if (numerator < 0.0f) {
					return false;
				}
				else {
					continue;
				}
			}

			const float t = numerator / denominator;

			if (denominator >= 0.0f) {
				/* Line i bounds line lineNo on the right. */
				//tRight = std::min(tRight, t);
				if(t<tRight){
					tRight=t;
					right_line=i;
				}
			}
			else {
				/* Line i bounds line lineNo on the left. */
				//tLeft = std::max(tLeft, t);
				if(t>tLeft){
					tLeft=t;
					left_line=i;
				}
			}

			if (tLeft > tRight) {
				return false;
			}
		}

		if (directionOpt) {
			/* Optimize direction. */
			if (optVelocity * lines[lineNo].direction > 0.0f) {
				/* Take right extreme. */
				result = lines[lineNo].point + tRight * lines[lineNo].direction;
				infeasibleLine2=right_line;
			}
			else {
				/* Take left extreme. */
				result = lines[lineNo].point + tLeft * lines[lineNo].direction;
				infeasibleLine2=left_line;
			}
		}
		else {
			/* Optimize closest point. */
			const float t = lines[lineNo].direction * (optVelocity - lines[lineNo].point);

			if (t < tLeft) {
				result = lines[lineNo].point + tLeft * lines[lineNo].direction;
			}
			else if (t > tRight) {
				result = lines[lineNo].point + tRight * lines[lineNo].direction;
			}
			else {
				result = lines[lineNo].point + t * lines[lineNo].direction;
			}
		}

		return true;
	}

	size_t Agent::linearProgram2(const std::vector<Line> &lines, float radius, const Vector2 &optVelocity, bool directionOpt, Vector2 &result)
	{
		if (directionOpt) {
			/*
			 * Optimize direction. Note that the optimization velocity is of unit
			 * length in this case.
			 */
			result = optVelocity * radius;
		}
		else if (absSq(optVelocity) > sqr(radius)) {
			/* Optimize closest point and outside circle. */
			result = normalize(optVelocity) * radius;
		}
		else {
			/* Optimize closest point and inside circle. */
			result = optVelocity;
		}

		infeasibleLine1=-2;

		for (size_t i = 0; i < lines.size(); ++i) {
			if (det(lines[i].direction, lines[i].point - result) > 0.0f) {
				/* Result does not satisfy constraint i. Compute new optimal result. */
				const Vector2 tempResult = result;
				
				infeasibleLine1=i;

				if (!linearProgram1(lines, i, radius, optVelocity, directionOpt, result)) {
					result = tempResult;
					return i;
				}
			}
		}

		return lines.size();
	}

	void Agent::linearProgram3( std::vector<Line> &lines, size_t numObstLines, size_t beginLine, float radius, Vector2 &result)
	{
		distance=0;

		for (size_t i = beginLine; i < lines.size(); ++i) {
			if (det(lines[i].direction, lines[i].point - result) > distance) {
				/* Result does not satisfy constraint of line i. */
				std::vector<Line> projLines(lines.begin(), lines.begin() + static_cast<ptrdiff_t>(numObstLines));

				infeasibleLine0=i;

				for (size_t j = numObstLines; j < i; ++j) {
					Line line;

					float determinant = det(lines[i].direction, lines[j].direction);

					if (std::fabs(determinant) <= RVO_EPSILON) {
						/* Line i and line j are parallel. */

						if (lines[i].direction * lines[j].direction > 0.0f) {
							/* Line i and line j point in the same direction. */
							line.point=Vector2(0,0);
							line.direction=Vector2(0,0);
							projLines.push_back(line);

							continue;
						}
						else {
							/* Line i and line j point in opposite direction. */
							line.point = 0.5f * (lines[i].point + lines[j].point);
						}
					}
					else {
						line.point = lines[i].point + (det(lines[j].direction, lines[i].point - lines[j].point) / determinant) * lines[i].direction;
					}

					line.direction = normalize(lines[j].direction - lines[i].direction);
					projLines.push_back(line);
				}

				const Vector2 tempResult = result;

				if (linearProgram2(projLines, radius, Vector2(-lines[i].direction.y(), lines[i].direction.x()), true, result) < projLines.size()) {
					/* This should in principle not happen.  The result is by definition
					 * already in the feasible region of this linear program. If it fails,
					 * it is due to small floating point error, and the current result is
					 * kept.
					 */
					result = tempResult;
				}

				distance = det(lines[i].direction, lines[i].point - result);
			}
		}
	}

}
