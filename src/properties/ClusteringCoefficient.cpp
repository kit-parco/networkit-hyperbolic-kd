/*
 * ClusteringCoefficient.cpp
 *
 *  Created on: 08.04.2013
 *      Author: cls
 */

#include "ClusteringCoefficient.h"

namespace NetworKit {

ClusteringCoefficient::ClusteringCoefficient() {
	// TODO Auto-generated constructor stub

}

ClusteringCoefficient::~ClusteringCoefficient() {
	// TODO Auto-generated destructor stub
}

double ClusteringCoefficient::calculate(Graph& G) {

	count n = G.numberOfNodes();
	std::vector<count> numerator(n); //
	std::vector<count> denominator(n); // $\deg(u) \cdot ( \deg(u) - 1 )$
	std::vector<double> coefficient(n); // $c(u) := \frac{2 \cdot |E(N(u))| }{\deg(u) \cdot ( \deg(u) - 1)}$


	G.parallelForNodes([&](node u){
		count edgeCount = 0;
		G.forEdgesOf(u, [&](node u, node v) {
			G.forEdgesOf(v, [&](node v, node w){
				if (G.hasEdge(u, w)) {
					edgeCount += 1;
				}
			});
		});

		numerator[u] = edgeCount; // factor 2 is omitted because each edge has been counted twice
	});

	G.parallelForNodes([&](node u){
		denominator[u] = G.degree(u) * (G.degree(u) - 1);
	});

	G.parallelForNodes([&](node u){
		coefficient[u] = numerator[u] / denominator[u];
	});

	double cc = G.parallelSumForNodes([&](node u){
		return coefficient[u];
	});

	cc = (1.0 / n) * cc;

	return cc;
}

} /* namespace NetworKit */
