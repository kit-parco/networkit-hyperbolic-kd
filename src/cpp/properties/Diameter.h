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
	static std::pair<edgeweight, edgeweight> estimatedDiameterRange(const IGraph& G, double error);

	/** @return exact diameter of the graph @a G */
	static edgeweight exactDiameter(const IGraph& G);


	/** @return a 2-approximation of the vertex diameter (unweighted diameter) of @a G.

		@param[in]	samples		One sample is enough if the graph is connected. If there 
								are multiple connected components, then the number of samples
								must be chosen so that the probability of sampling the component
								with the largest diameter ist high. 
	 */
	static edgeweight estimatedVertexDiameter(const IGraph& G, count samples);


	/** @return a 2-approximation of the vertex diameter (unweighted diameter) of @a G.
			Considers each connected component and returns the maximum diameter.
	 */
	static edgeweight estimatedVertexDiameterPedantic(const IGraph& G);
};

} /* namespace NetworKit */

#endif /* DIAMETER_H_ */
