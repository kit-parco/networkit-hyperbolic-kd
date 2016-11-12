/*
 * QuadNodePolarEuclid.h
 *
 *  Created on: 21.05.2014
 *      Author: Moritz v. Looz (moritz.looz-corswarem@kit.edu)
 *
 *  Note: This is similar enough to QuadNode.h that one could merge these two classes.
 */

#ifndef QUADNODEPOLAREUCLID_H_
#define QUADNODEPOLAREUCLID_H_

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
class QuadNodePolarEuclid : public NetworKit::SpatialCell<T> {
	friend class QuadTreeGTest;
private:
	static const long unsigned sanityNodeLimit = 10E15; //just assuming, for debug purposes, that this algorithm never runs on machines with more than 4 Petabyte RAM
	bool splitTheoretical;
	double balance;

public:
	/**
	 * Construct a QuadNode for polar coordinates.
	 *
	 *
	 */
	QuadNodePolarEuclid(Point<double> minCoords = {0,0}, Point<double> maxCoords = {2*M_PI, 1}, unsigned capacity = 1000, bool splitTheoretical = false, double balance = 0.5)
	: SpatialCell<T>(minCoords, maxCoords, capacity) {
		if (balance <= 0 || balance >= 1) throw std::runtime_error("Quadtree balance parameter must be between 0 and 1.");
		if (minCoords.getDimensions() != 2) throw std::runtime_error("Currently only supported for two dimensions");
		this->balance = balance;
		this->splitTheoretical = splitTheoretical;
	}

	void split() {
		assert(this->isLeaf);
		assert(this->minCoords.getDimensions() == 2);
		const double leftAngle = this->minCoords[0];
		const double rightAngle = this->maxCoords[0];
		const double minR = this->minCoords[1];
		const double maxR = this->maxCoords[1];
		//heavy lifting: split up!
		double middleAngle, middleR;
		if (splitTheoretical) {
			//Euclidean space is distributed equally
			middleAngle = (rightAngle - leftAngle) / 2 + leftAngle;
			middleR = pow(maxR*maxR*(1-balance)+minR*minR*balance, 0.5);
		} else {
			//median of points
			const index n = this->positions.size();
			vector<double> angles(n);
			vector<double> radii(n);
			for (index i = 0; i < n; i++) {
				angles[i] = this->positions[i][0];
				radii[i] = this->positions[i][1];
			}

			std::sort(angles.begin(), angles.end());
			middleAngle = angles[n/2];
			std::sort(radii.begin(), radii.end());
			middleR = radii[n/2];
		}
		assert(middleR < maxR);
		assert(middleR > minR);

		std::shared_ptr<QuadNodePolarEuclid<T> > southwest(new QuadNodePolarEuclid<T>({leftAngle, minR}, {middleAngle, middleR}, this->capacity, splitTheoretical, balance));
		std::shared_ptr<QuadNodePolarEuclid<T> > southeast(new QuadNodePolarEuclid<T>({middleAngle, minR}, {rightAngle, middleR}, this->capacity, splitTheoretical, balance));
		std::shared_ptr<QuadNodePolarEuclid<T> > northwest(new QuadNodePolarEuclid<T>({leftAngle, middleR}, {middleAngle, maxR}, this->capacity, splitTheoretical, balance));
		std::shared_ptr<QuadNodePolarEuclid<T> > northeast(new QuadNodePolarEuclid<T>({middleAngle, middleR}, {rightAngle, maxR}, this->capacity, splitTheoretical, balance));
		this->children = {southwest, southeast, northwest, northeast};
		this->isLeaf = false;
	}

	virtual std::pair<double, double> distances(const Point<double> &query) const override {
		return this->EuclideanPolarDistances(query);
	}

	virtual double distance(const Point<double> &query, index k) const override {
		return this->euclidDistancePolar(query[0], query[1], this->positions[k][0], this->positions[k][1]);
	}
};
}

#endif /* QUADNODE_H_ */
