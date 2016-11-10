/*
 * KDNode.h
 *
 *  Created on: 10.11.2016
 *      Author: moritzl
 *
 *      for now: hyperbolic geometry
 */

#ifndef KDNODE_H_
#define KDNODE_H_

#include "SpatialCell.h"
#include "../../geometric/HyperbolicSpace.h"


namespace NetworKit {

template <class T, bool poincare = false>
class KDNode: public NetworKit::SpatialCell<T> {
public:
	KDNode() = default;
	virtual ~KDNode() = default;

	KDNode(const Point<double> &minCoords, const Point<double> &maxCoords, count capacity=1000) {
		this->isLeaf = true;
		this->subTreeSize = 0;
		this->minCoords = minCoords;
		this->maxCoords = maxCoords;
		this->capacity = capacity;

		double leftAngle = this->minCoords[0];
		double minR = this->minCoords[1];
		double rightAngle = this->maxCoords[0];
		double maxR = this->maxCoords[1];

		assert(leftAngle >= 0);
		assert(rightAngle <= 2*M_PI);
		assert(leftAngle < rightAngle);

		assert(minR >= 0);
		assert(maxR > minR);
	}

	void split() override {
		assert(this->isLeaf);
		const index numPoints = this->positions.size();
		if (numPoints == 0) {
			throw std::runtime_error("Cannot split empty cell.");
		}

		const count dimension = this->minCoords.getDimensions();
		index mostSpreadDimension = dimension;
		double maximumSpread = 0;

		//find dimension in which to split
		for (index d = 0; d < dimension; d++) {
			double maxElement = this->minCoords[d] - 1;
			double minElement = this->maxCoords[d] + 1;
			for (index i = 0; i < numPoints; i++) {
				double pos = this->positions[i][d];
				if (pos < minElement) minElement = pos;
				if (pos > maxElement) maxElement = pos;
			}
			double spread = maxElement - minElement;
			if (spread > maximumSpread) {
				maximumSpread = spread;
				mostSpreadDimension = d;
			}
		}

		//find median
		vector<double> sorted(numPoints);
		for (index i = 0; i < numPoints; i++) {
			sorted[i] = this->positions[i][mostSpreadDimension];
		}

		std::sort(sorted.begin(), sorted.end());
		double middle = sorted[numPoints/2];
		assert(middle <= this->maxCoords[mostSpreadDimension]);
		assert(middle >= this->minCoords[mostSpreadDimension]);

		Point<double> newLower(this->minCoords);
		Point<double> newUpper(this->maxCoords);
		newLower[mostSpreadDimension] = middle;
		newUpper[mostSpreadDimension] = middle;

		std::shared_ptr<KDNode<T, poincare> > firstChild(new KDNode<T, poincare>(this->minCoords, newUpper, this->capacity));
		std::shared_ptr<KDNode<T, poincare> > secondChild(new KDNode<T, poincare>(newLower, this->maxCoords, this->capacity));

		this->children = {firstChild, secondChild};
		this->isLeaf = false;
	}

	std::pair<double, double> distances(const Point<double> &query) const override {
		return hyperbolicDistances(query);
	}

	double distance(const Point<double> &query, index k) const override {
		double result;
		Point<double> pos = this->positions[k];
		if (poincare) {
			result = HyperbolicSpace::poincareMetric(HyperbolicSpace::polarToCartesian(pos[0], pos[1]), HyperbolicSpace::polarToCartesian(query[0], query[1]));
		} else {
			result = HyperbolicSpace::nativeDistance(pos[0], pos[1], query[0], query[1]);
		}
		return result;
	}

	std::pair<double, double> hyperbolicDistances(Point<double> query) const {
		double phi = query[0];
		double r = query[1];

		double leftAngle = this->minCoords[0];
		double minR = this->minCoords[1];
		double rightAngle = this->maxCoords[0];
		double maxR = this->maxCoords[1];

		assert(leftAngle >= 0);
		assert(rightAngle <= 2*M_PI);
		assert(leftAngle < rightAngle);

		assert(minR >= 0);
		assert(maxR > minR);
		assert(phi >= 0);
		assert(phi < 2*M_PI);
		assert(r > 0);

		double minRHyper, maxRHyper, r_h;
		if (poincare) {
			minRHyper=HyperbolicSpace::EuclideanRadiusToHyperbolic(minR);
			maxRHyper=HyperbolicSpace::EuclideanRadiusToHyperbolic(maxR);
			r_h = HyperbolicSpace::EuclideanRadiusToHyperbolic(r);
		} else {
			minRHyper=minR;
			maxRHyper=maxR;
			r_h = r;
		}

		double coshr = cosh(r_h);
		double sinhr = sinh(r_h);
		double coshMinR = cosh(minRHyper);
		double coshMaxR = cosh(maxRHyper);
		double sinhMinR = sinh(minRHyper);
		double sinhMaxR = sinh(maxRHyper);
		double cosDiffLeft = cos(phi - leftAngle);
		double cosDiffRight = cos(phi - rightAngle);

		/**
		 * If the query point is not within the quadnode, the distance minimum is on the border.
		 * Need to check whether extremum is between corners:
		 */

		double coshMinDistance, coshMaxDistance;

		//Left border
		double lowerLeftDistance = coshMinR*coshr-sinhMinR*sinhr*cosDiffLeft;
		double upperLeftDistance = coshMaxR*coshr-sinhMaxR*sinhr*cosDiffLeft;
		if (this->responsible(query)) coshMinDistance = 1; //strictly speaking, this is wrong
		else coshMinDistance = min(lowerLeftDistance, upperLeftDistance);

		coshMaxDistance = max(lowerLeftDistance, upperLeftDistance);
		//double a = cosh(r_h);
		double b, extremum;

		if (phi != leftAngle) {
			b = sinhr*cosDiffLeft;
			extremum = log((coshr+b)/(coshr-b))/2;
			if (extremum < maxRHyper && extremum >= minRHyper) {
				double extremeDistance = cosh(extremum)*coshr-sinh(extremum)*sinhr*cosDiffLeft;
				coshMinDistance = min(coshMinDistance, extremeDistance);
				coshMaxDistance = max(coshMaxDistance, extremeDistance);
			}
		} else {
			coshMinDistance = 1;
		}

		/**
		 * cosh is a function from [0,\infty) to [1, \infty)
		 * Variables thus need
		 */
		assert(coshMaxDistance >= 1);
		assert(coshMinDistance >= 1);

		//Right border
		double lowerRightDistance = coshMinR*coshr-sinhMinR*sinhr*cosDiffRight;
		double upperRightDistance = coshMaxR*coshr-sinhMaxR*sinhr*cosDiffRight;
		coshMinDistance = min(coshMinDistance, lowerRightDistance);
		coshMinDistance = min(coshMinDistance, upperRightDistance);
		coshMaxDistance = max(coshMaxDistance, lowerRightDistance);
		coshMaxDistance = max(coshMaxDistance, upperRightDistance);

		assert(coshMaxDistance >= 1);
		assert(coshMinDistance >= 1);

		if (phi != rightAngle) {
			b = sinhr*cosDiffRight;
			extremum = log((coshr+b)/(coshr-b))/2;
			if (extremum < maxRHyper && extremum >= minRHyper) {
				double extremeDistance = cosh(extremum)*coshr-sinh(extremum)*sinhr*cosDiffRight;
				coshMinDistance = min(coshMinDistance, extremeDistance);
				coshMaxDistance = max(coshMaxDistance, extremeDistance);
			}
		} else {
			coshMinDistance = 1;
		}

		assert(coshMaxDistance >= 1);
		assert(coshMinDistance >= 1);

		//upper and lower borders
		if (phi >= leftAngle && phi < rightAngle) {
			double lower = cosh(abs(r_h-minRHyper));
			double upper = cosh(abs(r_h-maxRHyper));
			coshMinDistance = min(coshMinDistance, lower);
			coshMinDistance = min(coshMinDistance, upper);
			coshMaxDistance = max(coshMaxDistance, upper);
			coshMaxDistance = max(coshMaxDistance, lower);
		}

		assert(coshMaxDistance >= 1);
		assert(coshMinDistance >= 1);

		//again with mirrored phi
		double mirrorphi;
		if (phi >= M_PI) mirrorphi = phi - M_PI;
		else mirrorphi = phi + M_PI;
		if (mirrorphi >= leftAngle && mirrorphi < rightAngle) {
			double lower = coshMinR*coshr+sinhMinR*sinhr;
			double upper = coshMaxR*coshr+sinhMaxR*sinhr;
			coshMinDistance = min(coshMinDistance, lower);
			coshMinDistance = min(coshMinDistance, upper);
			coshMaxDistance = max(coshMaxDistance, upper);
			coshMaxDistance = max(coshMaxDistance, lower);
		}

		assert(coshMaxDistance >= 1);
		assert(coshMinDistance >= 1);

		double minDistance, maxDistance;
		minDistance = acosh(coshMinDistance);
		maxDistance = acosh(coshMaxDistance);
		assert(maxDistance >= 0);
		assert(minDistance >= 0);
		return std::pair<double, double>(minDistance, maxDistance);
	}

private:

};

} /* namespace NetworKit */
#endif /* KDNODE_H_ */
