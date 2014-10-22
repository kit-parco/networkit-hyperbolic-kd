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
#include "QuadNode.h"
#include "../../geometric/HyperbolicSpace.h"

namespace NetworKit {

template <class T>
class Quadtree {
	friend class QuadTreeTest;
public:
	Quadtree() {
		root = QuadNode<T>();
		this->maxRadius = 1;
	}

	Quadtree(double maxR,bool theoreticalSplit=false, double alpha=1) {
		root = QuadNode<T>(0, 0, 2*M_PI, maxR, 1000, 0,theoreticalSplit,alpha);
		this->maxRadius = maxR;
	}

	virtual ~Quadtree() {

	}

	void addContent(T newcomer, double angle, double R) {
		root.addContent(newcomer, angle, R);
	}

	bool removeContent(T toRemove, double angle, double R) {
		return root.removeContent(toRemove, angle, R);
	}

	vector<T> getElements() {
		return root.getElements();
	}

	vector<T> getCloseElements(Point2D<double> query, double maxDistance) {
		Point2D<double> origin(0,0);
		vector<T> circleDenizens;
		Point2D<double> center;

		double minPhi, maxPhi, radius;
		HyperbolicSpace::getEuclideanCircle(query, maxDistance, center, radius);
		double minR = center.length() - radius;
		double maxR = center.length() + radius;
		//assert(maxR < 1);//this looks fishy
		if (maxR > 1) maxR = 1;
		if (minR < 0) {
			maxR = std::max(abs(minR), maxR);
			minR = 0;
			minPhi = 0;
			maxPhi = 2*M_PI;
		} else {
			double spread = asin(radius / center.length());
			double phi_c, r_c;
			HyperbolicSpace::cartesianToPolar(center, phi_c, r_c);
			minPhi = phi_c - spread;
			maxPhi = phi_c + spread;
			/**
			 * If the circle overlaps the 2\pi line, we have to make two separate calls and collect
			 */
		}

		/**
		 * get Elements in Circle
		 */

		bool wraparound = false;
		root.getElementsInEuclideanCircle(minPhi, maxPhi, minR, maxR, center, radius, circleDenizens);
		if (minPhi < 0) {
			root.getElementsInEuclideanCircle(2*M_PI+minPhi, 2*M_PI, minR, maxR, center, radius, circleDenizens);
			wraparound = true;
		}
		if (maxPhi > 2*M_PI) {
			root.getElementsInEuclideanCircle(0, maxPhi - 2*M_PI, minR, maxR, center, radius, circleDenizens);
			wraparound = true;
		}

		//we have sort(deg(v)) here! This is not good, but does not make the asymptotical complexity of O(deg(v) log n) worse.
		if (wraparound) {
			std::sort(circleDenizens.begin(), circleDenizens.end());
			auto newend = unique(circleDenizens.begin(), circleDenizens.end());
			circleDenizens.resize(newend - circleDenizens.begin());
		}

		/**
		 * return the elements
		 */
		return circleDenizens;
	}

	count size() {
		return root.size();
	}

	count height() {
		return root.height();
	}


private:
	QuadNode<T> root;
	double maxRadius;
};
}

#endif /* QUADTREE_H_ */
