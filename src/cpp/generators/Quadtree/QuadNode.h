/*
 * QuadNode.h
 *
 *  Created on: 21.05.2014
 *      Author: Moritz v. Looz (moritz.looz-corswarem@kit.edu)
 */

#ifndef QUADNODE_H_
#define QUADNODE_H_

#include <vector>
#include <algorithm>
#include <assert.h>
#include "../../auxiliary/Log.h"
#include "../HyperbolicSpace.h"

using std::vector;

namespace NetworKit {

template <class T>
class QuadNode {
public:
	QuadNode() {
		leftAngle = 0;
		rightAngle = 2*M_PI;
		minR = 0;
		maxR = 10;//TODO: magic Number, careful.
		coshMinR = cosh(minR);
		sinhMinR = sinh(minR);
		coshMaxR = cosh(maxR);
		sinhMaxR = sinh(maxR);
		capacity = 20;
		isLeaf = true;
		minRegion = 0;
	}

	~QuadNode() {
		// TODO Auto-generated constructor stub
	}

	QuadNode(double leftAngle, double minR, double rightAngle, double maxR, unsigned capacity, double minDiameter) {
		this->leftAngle = leftAngle;
		this->minR = minR;
		coshMinR = cosh(minR);
		sinhMinR = sinh(minR);
		this->maxR = maxR;
		coshMaxR = cosh(maxR);
		sinhMaxR = sinh(maxR);
		this->rightAngle = rightAngle;
		this->capacity = 20;
		this->minRegion = minDiameter;
		isLeaf = true;
	}

	void addContent(T input, double angle, double R) {
		assert(this->responsible(angle, R));
		if (isLeaf) {
			if (content.size() + 1 < capacity ||  HyperbolicSpace::getDistance(leftAngle, minR, rightAngle, maxR) < minRegion) {
				content.push_back(input);
				angles.push_back(angle);
				coshradii.push_back(cosh(R));
				sinhradii.push_back(sinh(R));
				radii.push_back(R);
			} else {
				//heavy lifting: split up!
				double middleAngle = (rightAngle - leftAngle) / 2 + leftAngle;
				/**
				 * we want to make sure the space is evenly divided to obtain a balanced tree
				 * Simply halving the radius will cause a larger space for the outer Quadnode, resulting in an unbalanced tree
				 * The hyperbolic space grows with approximately e^R, so we try this.
				 */
				double rSpace = exp(maxR - minR);
				double middleR = log(rSpace/2) + minR;

				QuadNode southwest(leftAngle, minR, middleAngle, middleR, capacity, minRegion);
				QuadNode southeast(middleAngle, minR, rightAngle, middleR, capacity, minRegion);
				QuadNode northwest(leftAngle, middleR, middleAngle, maxR, capacity, minRegion);
				QuadNode northeast(middleAngle, middleR, rightAngle, maxR, capacity, minRegion);
				children = {southwest, southeast, northwest, northeast};

				isLeaf = false;
				for (uint i = 0; i < content.size(); i++) {
					this->addContent(content[i], angles[i], radii[i]);
				}
				content.clear();
				this->addContent(input, angle, R);
			}
		}
		else {
			assert(children.size() > 0);
			bool foundResponsibleChild = false;
			for (uint i = 0; i < children.size(); i++) {
				if (children[i].responsible(angle, R)) {
					assert(!foundResponsibleChild);//only one!
					children[i].addContent(input, angle, R);
					foundResponsibleChild = true;
				} else {
					//cout << "Not responsible for (" << angle << ", " << R << "). Borders are " << children[i].leftAngle << "-" << children[i].rightAngle << ", and " << children[i].minR << "-" << children[i].maxR << endl;
				}
			}
			assert(foundResponsibleChild);
		}
	}

	double distanceLowerBoundPrecached(double angle, double R, double coshR, double sinhR) {
		//return 0;//TODO: find proper lower Bound
		double nearestR = R;
		double coshNR, sinhNR;
		if (nearestR < minR) nearestR = minR;
		if (nearestR > maxR) nearestR = maxR;
		if (angle >= this->leftAngle && angle < this->rightAngle) {
			TRACE(this->leftAngle, " < ",  angle, " < " , this->rightAngle);
			return std::abs(R-nearestR);
		}

		if (nearestR == maxR) {
			coshNR = coshMaxR;
			sinhNR = sinhMaxR;
			TRACE(R, ">=", maxR);
		} else if (nearestR == minR) {
			coshNR = coshMinR;
			sinhNR = sinhMinR;
			TRACE(R, "<=", minR);
		} else {
			coshNR = cosh(nearestR);
			sinhNR = sinh(nearestR);
		}
		double leftDistance = HyperbolicSpace::getDistance(angle, R, this->leftAngle,  nearestR);
		double rightDistance = HyperbolicSpace::getDistance(angle, R, this->rightAngle, nearestR);
		TRACE("leftDistance:", leftDistance);
		TRACE("rightDistance:", rightDistance);
		return 0.89*std::min(leftDistance, rightDistance);
	}

	double distanceLowerBound(double angle, double R) {
		double coshR = cosh(R);
		double sinhR = sinh(R);
		return distanceLowerBoundPrecached(angle, R, coshR, sinhR);
	}

	double distanceUpperBound(double angle, double R) {
		double coshR = cosh(R);
		double sinhR = sinh(R);

		double leftLower = HyperbolicSpace::getDistance(angle, R, this->leftAngle, minR);
		double rightLower = HyperbolicSpace::getDistance(angle, R, this->rightAngle, minR);
		double leftUpper = HyperbolicSpace::getDistance(angle, R, this->leftAngle, maxR);
		double rightUpper = HyperbolicSpace::getDistance(angle, R, this->rightAngle, maxR);
		return std::max(std::max(leftLower, rightLower), std::max(leftUpper, rightUpper));
	}

	bool responsible(double angle, double R) {

		return (angle >= leftAngle && angle < rightAngle && R >= minR && R < maxR);
	}

	std::vector<T> getElements() {
		if (isLeaf) {
			return content;
		} else {
			vector<T> result;
			for (uint i = 0; i < children.size(); i++) {
				std::vector<T> subresult = children[i].getElements();
				result.insert(result.end(), subresult.begin(), subresult.end());
			}
			return result;
		}
	}

	std::vector<T> getCloseElementsPrecached(double angle, double R, double coshR, double sinhR, double maxDistance) {
		std::vector<T> result;
		if (isLeaf) {
			if (this->distanceLowerBoundPrecached(angle, R, coshR, sinhR) < maxDistance) {
				if (this->distanceUpperBound(angle, R) < maxDistance) {
					return content;
				}
				else {
					for (uint i = 0; i < content.size(); i++) {
						if (HyperbolicSpace::getDistance(angle, R, angles[i], radii[i]) < maxDistance) {
								result.push_back(content[i]);
							}
					}
				}
			}
		} else {
			for (uint i = 0; i < children.size(); i++) {
				QuadNode * child = &children[i];
				if (child->distanceLowerBoundPrecached(angle, R, coshR, sinhR) < maxDistance) {
					vector<T> subresult = child->getCloseElementsPrecached(angle, R, coshR, sinhR, maxDistance);
					result.insert(result.end(), subresult.begin(), subresult.end());
				}
			}
		}
		return result;
	}

	std::vector<T> getCloseElements(double angle, double R, double maxDistance) {
		double coshinput = cosh(R);
		double sinhinput = sinh(R);
		return getCloseElementsPrecached(angle, R, coshinput, sinhinput, maxDistance);
	}

	QuadNode<T> * getAppropriateLeaf(double angle, double R) {
		assert(this->responsible(angle, R));
		if (isLeaf) return this;
		else {
			for (uint i = 0; i < children.size(); i++) {
				bool foundResponsibleChild = false;
				if (children[i].responsible(angle, R)) {
					assert(foundResponsibleChild == false);
					foundResponsibleChild = true;
					return children[i].getAppropriateLeaf(angle, R);
				}
			}
			DEBUG("No responsible child for (", angle, ", ", R, ") found. Segfault imminent.");
		}
	}

	double getLeftAngle() {
		return leftAngle;
	}

	double getRightAngle() {
		return rightAngle;
	}

	double getMinR() {
		return minR;
	}

	double getMaxR() {
		return maxR;
	}

private:
	double leftAngle;
	double rightAngle;
	double minR;
	double posMinCircleCenter;
	double minCircleRadius;
	double coshMinR;
	double sinhMinR;
	double maxR;
	double posMaxCircleCenter;
	double maxCircleRadius;
	double coshMaxR;
	double sinhMaxR;
	unsigned capacity;
	double minRegion;//the minimal region a QuadNode should cover. If it is smaller, don't bother splitting up.
	std::vector<QuadNode> children;
	std::vector<T> content;
	std::vector<double> angles;
	std::vector<double> coshradii;
	std::vector<double> sinhradii;
	std::vector<double> radii;
	bool isLeaf;
};
}

#endif /* QUADNODE_H_ */
