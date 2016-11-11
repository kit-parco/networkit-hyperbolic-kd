/*
 * QuadNodePolarEuclid.h
 *
 *  Created on: 21.05.2014
 *      Author: Moritz v. Looz (moritz.looz-corswarem@kit.edu)
 *
 *  Note: This is similar enough to QuadNode.h that one could merge these two classes.
 */

#ifndef QUADNODECARTESIANEUCLID_H_
#define QUADNODECARTESIANEUCLID_H_

#include <vector>
#include <algorithm>
#include <functional>
#include <assert.h>
#include "SpatialCell.h"
#include "../../auxiliary/Log.h"
#include "../../geometric/HyperbolicSpace.h"

using std::vector;
using std::min;
using std::max;
using std::cos;

namespace NetworKit {

template <class T>
class QuadNodeCartesianEuclid : public NetworKit::SpatialCell<T> {
	friend class QuadTreeGTest;
private:
	static const long unsigned sanityNodeLimit = 10E15; //just assuming, for debug purposes, that this algorithm never runs on machines with more than 4 Petabyte RAM
	bool splitTheoretical;

public:
	virtual ~QuadNodeCartesianEuclid() = default;

	/**
	 * Construct a QuadNode for cartesian coordinates.
	 *
	 *
	 * @param lower Minimal coordinates of region
	 * @param upper Maximal coordinates of region (excluded)
	 * @param capacity Number of points a leaf cell can store before splitting
	 * @param splitTheoretical Whether to split in a theoretically optimal way or in a way to decrease measured running times
	 *
	 */
	QuadNodeCartesianEuclid(Point<double> lower = Point<double>({0.0, 0.0}), Point<double> upper = Point<double>({1.0, 1.0}), unsigned capacity = 1000, bool splitTheoretical = false) {
		this->minCoords = lower;
		this->maxCoords = upper;
		assert(this->maxCoords.getDimensions() == this->minCoords.getDimensions());
		this->capacity = capacity;
		this->splitTheoretical = splitTheoretical;
		this->ID = 0;
		this->subTreeSize = 0;
		this->isLeaf = true;
	}

	void split() override {
		assert(this->isLeaf);
		assert(this->children.size() == 0);
		const count dimension = this->minCoords.getDimensions();
		vector<double> middle(dimension);
		if (splitTheoretical) {
			//Euclidean space is distributed equally
			for (index d = 0; d < dimension; d++) {
				middle[d] = (this->minCoords[d] + this->maxCoords[d]) / 2;
			}
		} else {
			//median of points
			const count numPoints = this->positions.size();
			assert(numPoints > 0);//otherwise, why split?
			vector<vector<double> > sorted(dimension);
			for (index d = 0; d < dimension; d++) {
				sorted[d].resize(numPoints);
				for (index i = 0; i < numPoints; i++) {
					sorted[d][i] = this->positions[i][d];
				}
				std::sort(sorted[d].begin(), sorted[d].end());
				middle[d] = sorted[d][numPoints/2];//this will crash if no points are there!
				assert(middle[d] <= this->maxCoords[d]);
				assert(middle[d] >= this->minCoords[d]);
			}
		}
		count childCount = pow(2,dimension);
		for (index i = 0; i < childCount; i++) {
			vector<double> lowerValues(dimension);
			vector<double> upperValues(dimension);
			index bitCopy = i;
			for (index d = 0; d < dimension; d++) {
				if (bitCopy & 1) {
					lowerValues[d] = middle[d];
					upperValues[d] = this->maxCoords[d];
				} else {
					lowerValues[d] = this->minCoords[d];
					upperValues[d] = middle[d];
				}
				bitCopy = bitCopy >> 1;
			}
			std::shared_ptr<QuadNodeCartesianEuclid<T> > child(new QuadNodeCartesianEuclid<T>(Point<double>(lowerValues), Point<double>(upperValues), this->capacity, splitTheoretical));
			assert(child->isLeaf);
			this->children.push_back(child);
		}
		this->isLeaf = false;
	}

	/**
	 * @param query Position of the query point
	 */
	std::pair<double, double> EuclideanDistances(Point<double> query) const {
		/**
		 * If the query point is not within the quadnode, the distance minimum is on the border.
		 * Need to check whether extremum is between corners.
		 */
		double maxDistance = 0;
		double minDistance = std::numeric_limits<double>::max();
		const count dimension = this->minCoords.getDimensions();

		if (this->responsible(query)) minDistance = 0;

		auto updateMinMax = [&minDistance, &maxDistance, query](Point<double> pos){
			double extremalValue = pos.distance(query);
			maxDistance = std::max(extremalValue, maxDistance);
			minDistance = std::min(minDistance, extremalValue);
		};

		vector<double> closestValues(dimension);
		vector<double> farthestValues(dimension);

		for (index d = 0; d < dimension; d++) {
			if (std::abs(query[d] - this->minCoords.at(d)) < std::abs(query[d] - this->maxCoords.at(d))) {
				closestValues[d] = this->minCoords.at(d);
				farthestValues[d] = this->maxCoords.at(d);
			} else {
				farthestValues[d] = this->minCoords.at(d);
				closestValues[d] = this->maxCoords.at(d);
			}
			if (query[d] >= this->minCoords.at(d) && query[d] <= this->maxCoords.at(d)) {
				closestValues[d] = query[d];
			}
		}
		updateMinMax(Point<double>(closestValues));
		updateMinMax(Point<double>(farthestValues));

		assert(minDistance < query.length() + this->maxCoords.length());
		assert(minDistance < maxDistance);
		return std::pair<double, double>(minDistance, maxDistance);
	}

	virtual std::pair<double, double> distances(const Point<double> &query) const override {
		return EuclideanDistances(query);
	}

	virtual double distance(const Point<double> &query, index k) const override {
		return query.distance(this->positions[k]);
	}
};
}

#endif /* QUADNODE_H_ */
