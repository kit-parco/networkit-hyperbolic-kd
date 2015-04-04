/*
 * NeighborhoodDistanceIndex.h
 *
 *  Created on: 24.06.2013
 *      Authors: cls, Kolja Esders
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
private:
  virtual double runImpl(node u, node v) override;

public:
  using LinkPredictor::LinkPredictor;
  
};

} /* namespace NetworKit */
#endif /* NEIGHBORHOODDISTANCEINDEX_H_ */
