/*
 * MultiscaleBackbone.cpp
 *
 *  Created on: 20.06.2014
 *      Author: Gerd Lindner
 */

#include "MultiscaleBackbone.h"

namespace NetworKit {

MultiscaleBackbone::MultiscaleBackbone(double a) : alpha(a) {}

Graph MultiscaleBackbone::calculate(const Graph& graph) {
	Graph backboneGraph = cloneNodes(graph, true);

	//The following vector is used for the _local_ normalization of edgeweights.
	//We use a global vector for performance reasons.
	std::vector<edgeweight> normalizedWeights(graph.upperNodeIdBound());

	graph.forNodes([&](node u) {
		count k = graph.degree(u);

		//Normalize edgeweights of N(u)
		edgeweight sum = 0.0;
		graph.forNeighborsOf(u, [&](node v) {
			sum += graph.weight(u, v);
		});
		graph.forNeighborsOf(u, [&](node v) {
			normalizedWeights[v] = graph.weight(u, v) / sum;
		});

		//Filter edges by probability
		graph.forNeighborsOf(u, [&](node v) {
			//In case d(u) == 1 and d(v) > 1: ignore u
			if (k > 1 || graph.degree(v) == 1) {
				edgeweight p = normalizedWeights[v];
				double probability = getProbability(k, p);

				if (probability < alpha) {
					backboneGraph.setWeight(u, v, graph.weight(u, v));
				}
			}
		});
	});

	return backboneGraph;
}

/**
 * Returns the probability that a node of the given
 * degree has an edge of the given weight.
 *
 * The null hypothesis is the following: the normalized weights of the
 * edges connected to a node of degree k are uniformly distributed.
 */
double MultiscaleBackbone::getProbability(count degree, edgeweight normalizedWeight) {
	return 1 - (1 - pow(1 - normalizedWeight, degree - 1));
}

} /* namespace NetworKit */
