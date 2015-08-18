/*
 * LocalDegreeScore.cpp
 *
 *  Created on: 28.08.2014
 *      Author: Gerd Lindner
 */

#include "LocalDegreeScore.h"
#include "LocalSimilarityAttributizer.h"

namespace NetworKit {

LocalDegreeScore::LocalDegreeScore(const Graph& graph) : graph(graph) {
}

std::vector<double> LocalDegreeScore::scores() {
	return scoreData;
}

double LocalDegreeScore::score(node u, node v) {
	throw std::runtime_error("Not implemented: Use scores() instead.");
}

double LocalDegreeScore::score(edgeid eid) {
	throw std::runtime_error("Not implemented: Use scores() instead.");
}

void LocalDegreeScore::run() {
	if (!graph.hasEdgeIds()) {
		throw std::runtime_error("edges have not been indexed - call indexEdges first");
	}

	scoreData.resize(graph.upperEdgeIdBound(), 0.0);

	graph.balancedParallelForNodes([&](node i) {
		count d = graph.degree(i);

		/**
		 *  The top d^e edges (sorted by degree)
		 * are to be kept in the backbone */

		std::vector<AttributizedEdge<count>> neighbors;
		graph.forNeighborsOf(i, [&](node _i, node j, edgeid eid) {
			if (graph.degree(j) > d)
				neighbors.push_back(AttributizedEdge<count>(i, j, eid, graph.degree(j)));
		});
		std::sort(neighbors.begin(), neighbors.end());

		count rank = 1;

		/**
		 * By convention, we want to the edges with highest "similarity" or "cohesion" to have values close to 1,
		 * so we invert the range.
		 */

		#pragma omp critical
		for (auto neighborEdge : neighbors) {
			edgeid eid = neighborEdge.eid;

			double e = 1.0; // If the node has only one neighbor, the edge should be kept anyway.
			if (d > 1)
				e = 1.0 - (log(rank) / log(d));

			scoreData[eid] = std::max(e, scoreData[eid]);
			rank++;
		}

	});

	hasRun = true;
}

} /* namespace NetworKit */
