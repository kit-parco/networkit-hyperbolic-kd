/*
 * KDTreeEuclidean.h
 *
 *  Created on: 12.11.2016
 *      Author: moritzl
 */

#ifndef KDTREEEUCLIDEAN_H_
#define KDTREEEUCLIDEAN_H_

#include <vector>

#include "SpatialTree.h"
#include "KDNodeEuclidean.h"

using std::vector;

namespace NetworKit {

template <class T, bool cartesian=true>
class KDTreeEuclidean: public NetworKit::SpatialTree<T> {
public:
	KDTreeEuclidean() = default;
	virtual ~KDTreeEuclidean() = default;
	KDTreeEuclidean(const Point<double> &minCoords, const Point<double> &maxCoords, count capacity=1000) {
		this->root = std::shared_ptr<KDNodeEuclidean<T, cartesian> >(new KDNodeEuclidean<T, cartesian>(minCoords, maxCoords, capacity));
	}
};

} /* namespace NetworKit */
#endif /* KDTREEEUCLIDEAN_H_ */
