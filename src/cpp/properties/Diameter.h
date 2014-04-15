/*
 * Diameter.h
 *
 *  Created on: 19.02.2014
 *      Author: Daniel Hoske, Christian Staudt
 */

#ifndef DIAMETER_H_
#define DIAMETER_H_

#include "../graph/Graph.h"

namespace NetworKit {

class Diameter {

public:

	/**
	 * Estimates a range for the diameter of @a G. Based on the algorithm suggested in
	 * C. Magnien, M. Latapy, M. Habib: Fast Computation of Empirically Tight Bounds for
	 * the Diameter of Massive Graphs. Journal of Experimental Algorithmics, Volume 13, Feb 2009.
	 *
	 * @return Pair of lower and upper bound for diameter.
	 */
	static std::pair<edgeweight, edgeweight> estimatedDiameterRange(const Graph& G, double error);

	/** @return exact diameter of the graph @a G */
	static edgeweight exactDiameter(const Graph& G);


	/** @return a 2-approximation of the vertex diameter (unweighted diameter) of @a G,
				the maximum over all components in case the graph is disconnected
	 */
	static count estimatedVertexDiameter(const Graph& G);
};

} /* namespace NetworKit */

#endif /* DIAMETER_H_ */
