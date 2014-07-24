/*
 * MultiscaleAttributizer.cpp
 *
 *  Created on: 20.06.2014
 *      Author: Gerd Lindner
 */

#include "MultiscaleAttributizer.h"

namespace NetworKit {

MultiscaleAttributizer::MultiscaleAttributizer() {}

EdgeAttribute MultiscaleAttributizer::getAttribute(const Graph& graph, const EdgeAttribute& attribute) {
	//The following vector is used for the _local_ normalization of edgeweights.
	//We use a global vector for performance reasons.
	std::vector<edgeweight> normalizedWeights(graph.upperNodeIdBound());

	EdgeAttribute multiscaleAttribute;

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

				uEdge key = uEdge(u, v);
				if (multiscaleAttribute.contains(key))
					multiscaleAttribute.set(key, std::min(multiscaleAttribute[key], probability));
				else
					multiscaleAttribute.set(key, probability);
			}
		});
	});

	return multiscaleAttribute;
}

/**
 * Returns the probability that a node of the given
 * degree has an edge of the given weight.
 *
 * The null hypothesis is the following: the normalized weights of the
 * edges connected to a node of degree k are uniformly distributed.
 */
double MultiscaleAttributizer::getProbability(count degree, edgeweight normalizedWeight) {
	return 1 - (1 - pow(1 - normalizedWeight, degree - 1));
}

} /* namespace NetworKit */
