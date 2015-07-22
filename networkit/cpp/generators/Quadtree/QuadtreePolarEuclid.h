/*
 * Quadtree.h
 *
 *  Created on: 21.05.2014
 *      Author: Moritz v. Looz (moritz.looz-corswarem@kit.edu)
 */

#ifndef QUADTREE_H_
#define QUADTREE_H_

#include <vector>
#include <memory>
#include <cmath>
#include <omp.h>
#include <functional>
#include "QuadNodePolarEuclid.h"

namespace NetworKit {

template <class T>
class Quadtree {
	friend class QuadTreeGTest;
public:
	Quadtree() {
		root = QuadNodePolarEuclid<T>();
		this->maxRadius = 1;
	}

	/**
	 * @param maxR Radius of the managed area. Must be smaller than 1.
	 * @param theoreticalSplit If true, split cells to get the same area in each child cell. Default is false
	 * @param alpha dispersion Parameter of the point distribution. Only has an effect if theoretical split is true
	 * @param capacity How many points can inhabit a leaf cell before it is split up?
	 *
	 */
	Quadtree(double maxR,bool theoreticalSplit=false, double alpha=1, count capacity=1000, double balance = 0.5) {
		root = QuadNodePolarEuclid<T>(0, 0, 2*M_PI, maxR, capacity, 0,theoreticalSplit,alpha,balance);
		this->maxRadius = maxR;
	}

	Quadtree(const vector<double> &angles, const vector<double> &radii, const vector<T> &content, double r, bool theoreticalSplit=false, double alpha=1, count capacity=1000, double balance = 0.5) {
		const count n = angles.size();
		assert(angles.size() == radii.size());
		assert(radii.size() == content.size());
		root = QuadNodePolarEuclid<T>(0, 0, 2*M_PI, r, capacity, 0,theoreticalSplit,alpha,balance);
		maxRadius = r;
		for (index i = 0; i < n; i++) {
			assert(content[i] < n);
			root.addContent(content[i], angles[i], radii[i]);
		}
	}

	/**
	 * @param newcomer content to be added at point x
	 * @param angle angular coordinate of x
	 * @param R radial coordinate of x
	 */
	void addContent(T newcomer, double angle, double r) {
		root.addContent(newcomer, angle, r);
	}

	/**
	 * @param newcomer content to be removed at point x
	 * @param angle angular coordinate of x
	 * @param R radial coordinate of x
	 */
	bool removeContent(T toRemove, double angle, double r) {
		return root.removeContent(toRemove, angle, r);
	}

	/**
	 * Get all elements, regardless of position
	 *
	 * @return vector<T> of elements
	 */
	vector<T> getElements() const {
		return root.getElements();
	}

	void extractCoordinates(vector<double> &anglesContainer, vector<double> &radiiContainer) const {
		root.getCoordinates(anglesContainer, radiiContainer);
	}

	void getElementsInEuclideanCircle(const Point2D<double> circleCenter, const double radius, vector<T> &circleDenizens) const {
		root.getElementsInEuclideanCircle(circleCenter, radius, false, circleDenizens);
	}

	count getElementsProbabilistically(Point2D<double> euQuery, std::function<double(double)> prob, vector<T> &circleDenizens) {
		return root.getElementsProbabilistically(euQuery, prob, false, circleDenizens);
	}

	count getElementsProbabilistically(Point2D<double> euQuery, std::function<double(double)> prob, bool suppressLeft, vector<T> &circleDenizens) {
		return root.getElementsProbabilistically(euQuery, prob, suppressLeft, circleDenizens);
	}

	void recount() {
		root.recount();
	}

	count size() const {
		return root.size();
	}

	count height() const {
		return root.height();
	}

	count countLeaves() const {
		return root.countLeaves();
	}

	index indexSubtree(index nextID) {
		return root.indexSubtree(nextID);
	}

	index getCellID(double phi, double r) const {
		return root.getCellID(phi, r);
	}

	double getMaxRadius() const {
		return maxRadius;
	}

	void sortPointsInLeaves() {
		#pragma omp parallel
		{
			#pragma omp single nowait
			{
				root.sortPointsInLeaves();
			}
		}
	}

	void reindex() {
		#pragma omp parallel
		{
			#pragma omp single nowait
			{
				root.reindex(0);
			}
		}
	}

	/**
	 * trims the vectors used to hold the content in the leaf cells. Reduces memory usage, makes changes slower
	 */
	void trim() {
		root.trim();
	}

private:
	QuadNodePolarEuclid<T> root;
	double maxRadius;
};
}

#endif /* QUADTREE_H_ */
