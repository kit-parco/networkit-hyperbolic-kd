/*
 * NeighborhoodDistanceIndex.h
 *
 *  Created on: 24.06.2013
 *      Author: cls
 */

#ifndef NEIGHBORHOODDISTANCEINDEX_H_
#define NEIGHBORHOODDISTANCEINDEX_H_

#include "LinkPredictor.h"
#include <math.h>
#include <algorithm>

namespace NetworKit {

/**
 * @ingroup linkprediction
 *
 * Assigns a distance value to pairs of nodes according to the
 * overlap of their neighborhoods.
 */
class NeighborhoodDistanceIndex : public LinkPredictor {

public:

	NeighborhoodDistanceIndex(const Graph& G);

	virtual void preprocess();

	virtual double runImpl(node u, node v);
};

} /* namespace NetworKit */
#endif /* NEIGHBORHOODDISTANCEINDEX_H_ */
