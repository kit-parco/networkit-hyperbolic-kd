/*
 * KDTree.h
 *
 *  Created on: 10.11.2016
 *      Author: moritzl
 */

#ifndef KDTREE_H_
#define KDTREE_H_

#include <vector>

#include "SpatialTree.h"
#include "KDNode.h"

using std::vector;

namespace NetworKit {

template <class T, bool poincare = false>
class KDTree: public NetworKit::SpatialTree<T> {
public:
	KDTree() = default;
	virtual ~KDTree() = default;
	KDTree(const Point<double> &minCoords, const Point<double> &maxCoords, count capacity=1000) {
		root = KDNode<T, poincare>(minCoords, maxCoords, capacity);
	}

	void addContent(T content, const Point<double> &coords) {//why not template Point with the number of dimensions?
		root.addContent(content, coords);
	}

	/**
	 * @param newcomer content to be added at point x
	 * @param angle angular coordinate of x
	 * @param R radial coordinate of x
	 */
	void addContent(T content, double angle, double r) {
		root.addContent(content, {angle, r});
	}

	void getElementsInEuclideanCircle(const Point<double> query, const double radius, vector<T> &circleDenizens) const {
		root.getElementsInEuclideanCircle(query, radius, circleDenizens);
	}

	count getElementsProbabilistically(Point<double> query, std::function<double(double)> prob, vector<T> &circleDenizens) {
		return root.getElementsProbabilistically(query, prob, circleDenizens);
	}

	count size() const {
		return root.size();
	}

	count height() const {
		return root.height();
	}

	void trim() {
		root.trim();
	}

private:
	KDNode<T, poincare> root;

};

} /* namespace NetworKit */
#endif /* KDTREE_H_ */
