/*
 * LocalSimilarityAttributizer.cpp
 *
 *  Created on: 26.07.2014
 *      Author: Gerd Lindner
 */

#include "LocalSimilarityAttributizer.h"
#include <math.h> //log
#include <set>
#include "../auxiliary/Log.h"

namespace NetworKit {

LocalSimilarityAttributizer::LocalSimilarityAttributizer() {}

EdgeAttribute LocalSimilarityAttributizer::getAttribute(const Graph& graph, const EdgeAttribute& attribute) {
	/*
	 * For each edge, we calculate the minimum required sparsification exponent e
	 * such that the edge is contained in the backbone.
	 */

	EdgeAttribute sparsificationExp(1.0);

	graph.forNodes([&](node i) {
		count d = graph.degree(i);

		/*
		 * The top d^e edges (sorted by similarity in descending order)
		 * are to be kept in the backbone.
		 */

		std::vector<AttributizedEdge> neighbors;
		graph.forNeighborsOf(i, [&](node j) {
			double sim = getSimilarity(graph, i, j);
			//ERROR("sim: ", sim, "\n");
			neighbors.push_back(AttributizedEdge(i, j, sim));
		});
		std::sort(neighbors.begin(), neighbors.end(), greater());

		count rank = 1;
		for(std::vector<AttributizedEdge>::iterator it = neighbors.begin(); it != neighbors.end(); ++it) {
			uEdge edgeKey = uEdge(it->ego, it->alter);

			double e = 0.0;
			if (d > 1)
				e = log(rank) / log(d);

			//ERROR("e: ", e, ",",rank,",",d, "\n");
			sparsificationExp.set(edgeKey, std::min(e, sparsificationExp[edgeKey]));
			rank++;
		}

	});

	return sparsificationExp;
}

/**
 * Returns the similarity between two nodes.
 */
double LocalSimilarityAttributizer::getSimilarity(const Graph& graph, node u, node v) {
	//Use the jaccard measure as similarity measure.
	/* TODO: The following implementation might be quite inefficient....
	 * can the implementation in community/JaccardMeasure be used?
	 */
	std::set<node> uNeighbors;
	graph.forNeighborsOf(u, [&](node n) {
		uNeighbors.insert(n);
	});

	count inUnion = graph.degree(u);
	count inIntersection = 0;

	graph.forNeighborsOf(v, [&](node n) {
		if (uNeighbors.erase(n))
			inIntersection++;
		else
			inUnion++;
	});

	return (double) inIntersection / (double) inUnion;
}

} /* namespace NetworKit */
