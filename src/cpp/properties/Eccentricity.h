/*
 * Eccentricity.h
 *
 *  Created on: 19.02.2014
 *      Author: cls
 */

#ifndef ECCENTRICITY_H_
#define ECCENTRICITY_H_

#include "../graph/Graph.h"

namespace NetworKit {

class Eccentricity {

public:

	/**
	 * TODO: documentation
	 */
	static std::pair<edgeweight, edgeweight> getValue(const IGraph& G, node u);
};

} /* namespace NetworKit */

#endif /* ECCENTRICITY_H_ */
