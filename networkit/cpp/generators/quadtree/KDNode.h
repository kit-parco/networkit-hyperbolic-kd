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
		return this->hyperbolicPolarDistances(query, poincare);
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
};

} /* namespace NetworKit */
#endif /* KDNODE_H_ */
