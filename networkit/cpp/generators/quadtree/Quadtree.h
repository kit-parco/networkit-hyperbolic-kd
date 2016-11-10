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
#include "QuadNode.h"
#include "../../geometric/HyperbolicSpace.h"
#include "../../auxiliary/Parallel.h"

namespace NetworKit {

template <class T, bool poincare=true>
class Quadtree {
	friend class QuadTreeGTest;
public:
	Quadtree() {
		root = QuadNode<T, poincare>();
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
		root = QuadNode<T,poincare>({0, 0}, {2*M_PI, maxR}, capacity, theoreticalSplit,alpha,balance);
		this->maxRadius = maxR;
	}

	//	QuadNode(const Point<double> &minCoords, const Point<double> &maxCoords, count capacity=1000, bool splitTheoretical = false, double alpha = 1, double balance = 0.5) {


	Quadtree(const vector<double> &angles, const vector<double> &radii, const vector<T> &content, double stretch, bool theoreticalSplit=false, double alpha=1, count capacity=1000, double balance = 0.5) {
		const count n = angles.size();
		assert(angles.size() == radii.size());
		assert(radii.size() == content.size());
		double R = stretch*HyperbolicSpace::hyperbolicAreaToRadius(n);
		double r = HyperbolicSpace::hyperbolicRadiusToEuclidean(R);
		root = QuadNode<T,poincare>({0, 0}, {2*M_PI, r}, capacity, theoreticalSplit,alpha,balance);
		maxRadius = r;
		for (index i = 0; i < n; i++) {
			assert(content[i] < n);
			root.addContent(content[i], {angles[i], radii[i]});
		}
	}

	/**
	 * @param newcomer content to be added at point x
	 * @param angle angular coordinate of x
	 * @param R radial coordinate of x
	 */
	void addContent(T newcomer, double angle, double r) {
		root.addContent(newcomer, {angle, r});
	}

	/**
	 * @param newcomer content to be removed at point x
	 * @param angle angular coordinate of x
	 * @param R radial coordinate of x
	 */
	bool removeContent(T toRemove, double angle, double r) {
		return root.removeContent(toRemove, {angle, r});
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

	/**
	 * Get elements whose hyperbolic distance to the query point is less than the hyperbolic distance
	 *
	 *
	 * @param circleCenter Cartesian coordinates of the query circle's center
	 * @param hyperbolicRadius Radius of the query circle
	 */
	vector<T> getElementsInHyperbolicCircle(Point<double> circleCenter, double hyperbolicRadius) const {
		vector<T> circleDenizens;
		root.getElementsInCircle(circleCenter, hyperbolicRadius, circleDenizens);
		return circleDenizens;
	}

	void getElementsInHyperbolicCircle(const Point<double> circleCenter, const double hyperbolicRadius, vector<T> &circleDenizens) const {
		getElementsInHyperbolicCircle(circleCenter, hyperbolicRadius, circleDenizens);
	}

	count getElementsProbabilistically(Point<double> query, std::function<double(double)> prob, vector<T> &circleDenizens) {
		return root.getElementsProbabilistically(query, prob, circleDenizens);
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
		return root.getCellID({phi, r});
	}

	double getMaxRadius() const {
		return maxRadius;
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
	QuadNode<T, poincare> root;
	double maxRadius;
};
}

#endif /* QUADTREE_H_ */
