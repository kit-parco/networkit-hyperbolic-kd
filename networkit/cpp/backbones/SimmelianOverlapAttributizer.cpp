/*
 * SimmelianOverlapAttributizer.cpp
 *
 *  Created on: 22.07.2014
 *      Author: Gerd Lindner
 */

#include "SimmelianOverlapAttributizer.h"
#include <limits>

namespace NetworKit {

SimmelianOverlapAttributizer::SimmelianOverlapAttributizer(count maxRank) :
		maxRank(maxRank) {}

std::vector<double> SimmelianOverlapAttributizer::getAttribute(const Graph& graph, const std::vector<count>& attribute) {
	std::vector<RankedNeighbors> neighbors = getRankedNeighborhood(graph, attribute);
	std::vector<double> overlapAttribute(graph.upperEdgeIdBound(), 0.0);

	graph.parallelForEdges([&](node u, node v, edgeid eid) {
		Redundancy redundancy = getOverlap(u, v, neighbors, maxRank);

		overlapAttribute[eid] = (double) redundancy.overlap;
	});

	return overlapAttribute;
}

} /* namespace NetworKit */
