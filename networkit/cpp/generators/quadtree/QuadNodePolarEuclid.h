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
	 * @param leftAngle Minimal angular coordinate of region, in radians from 0 to 2\pi
	 * @param rightAngle Maximal angular coordinate of region, in radians from 0 to 2\pi
	 * @param minR Minimal radial coordinate of region, between 0 and 1
	 * @param maxR Maximal radial coordinate of region, between 0 and 1
	 * @param capacity Number of points a leaf cell can store before splitting
	 * @param minDiameter Minimal diameter of a quadtree node. If the node is already smaller, don't split even if over capacity. Default is 0
	 * @param splitTheoretical Whether to split in a theoretically optimal way or in a way to decrease measured running times
	 * @param alpha dispersion Parameter of the point distribution. Only has an effect if theoretical split is true
	 * @param diagnostics Count how many necessary and unnecessary comparisons happen in leaf cells? Will cause race condition and false sharing in parallel use
	 *
	 */
	QuadNodePolarEuclid(Point<double> minCoords = {0,0}, Point<double> maxCoords = {2*M_PI, 1}, unsigned capacity = 1000, bool splitTheoretical = false, double balance = 0.5) {
		if (balance <= 0 || balance >= 1) throw std::runtime_error("Quadtree balance parameter must be between 0 and 1.");
		if (minCoords.getDimensions() != 2) throw std::runtime_error("Currently only supported for two dimensions");
		this->balance = balance;
		this->minCoords = minCoords;
		this->maxCoords = maxCoords;
		assert(this->maxCoords.getDimensions() == this->minCoords.getDimensions());
		this->capacity = capacity;
		this->splitTheoretical = splitTheoretical;
		this->ID = 0;
		this->subTreeSize = 0;
		this->isLeaf = true;
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
		return EuclideanDistances(query);
	}

	virtual double distance(const Point<double> &query, index k) const override {
		return euclidDistancePolar(query[0], query[1], this->positions[k][0], this->positions[k][1]);
	}

	static double euclidDistancePolar(double phi_a, double r_a, double phi_b, double r_b){
		return pow(r_a*r_a+r_b*r_b-2*r_a*r_b*cos(phi_a-phi_b), 0.5);
	}

	/**
	 * @param phi Angular coordinate of query point
	 * @param r_h radial coordinate of query point
	 */
	std::pair<double, double> EuclideanDistances(Point<double> query) const {
		const double phi = query[0];
		const double r = query[1];
		const double leftAngle = this->minCoords[0];
		const double rightAngle = this->maxCoords[0];
		const double minR = this->minCoords[1];
		const double maxR = this->maxCoords[1];
		/**
		 * If the query point is not within the quadnode, the distance minimum is on the border.
		 * Need to check whether extremum is between corners.
		 */
		double maxDistance = 0;
		double minDistance = std::numeric_limits<double>::max();

		if (this->responsible(query)) minDistance = 0;

		auto updateMinMax = [&minDistance, &maxDistance, phi, r](double phi_b, double r_b){
			double extremalValue = euclidDistancePolar(phi, r, phi_b, r_b);
			//assert(extremalValue <= r + r_b);
			maxDistance = std::max(extremalValue, maxDistance);
			minDistance = std::min(minDistance, extremalValue);
		};

		/**
		 * angular boundaries
		 */
		//left
		double extremum = r*cos(leftAngle - phi);
		if (extremum < maxR && extremum > minR) {
			updateMinMax(leftAngle, extremum);
		}

		//right
		extremum = r*cos(rightAngle - phi);
		if (extremum < maxR && extremum > minR) {
			updateMinMax(rightAngle, extremum);
		}


		/**
		 * radial boundaries.
		 */
		if (phi > leftAngle && phi < rightAngle) {
			updateMinMax(phi, maxR);
			updateMinMax(phi, minR);
		}
		if (phi + M_PI > leftAngle && phi + M_PI < rightAngle) {
			updateMinMax(phi + M_PI, maxR);
			updateMinMax(phi + M_PI, minR);
		}
		if (phi - M_PI > leftAngle && phi -M_PI < rightAngle) {
			updateMinMax(phi - M_PI, maxR);
			updateMinMax(phi - M_PI, minR);
		}

		/**
		 * corners
		 */
		updateMinMax(leftAngle, maxR);
		updateMinMax(rightAngle, maxR);
		updateMinMax(leftAngle, minR);
		updateMinMax(rightAngle, minR);

		//double shortCutGainMax = maxR + r - maxDistance;
		//assert(minDistance <= minR + r);
		//assert(maxDistance <= maxR + r);
		assert(minDistance < maxDistance);
		return std::pair<double, double>(minDistance, maxDistance);
	}
};
}

#endif /* QUADNODE_H_ */
