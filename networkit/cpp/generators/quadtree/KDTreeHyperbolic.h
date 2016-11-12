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
#include "KDNodeHyperbolic.h"

using std::vector;

namespace NetworKit {

template <class T, bool poincare = false>
class KDTreeHyperbolic: public NetworKit::SpatialTree<T> {
public:
	KDTreeHyperbolic() = default;
	virtual ~KDTreeHyperbolic() = default;
	KDTreeHyperbolic(const Point<double> &minCoords, const Point<double> &maxCoords, count capacity=1000) {
		this->root = std::shared_ptr<KDNodeHyperbolic<T, poincare> >(new KDNodeHyperbolic<T, poincare>(minCoords, maxCoords, capacity));
	}
};

} /* namespace NetworKit */
#endif /* KDTREE_H_ */
