/*
 * GraphStructuralRandMeasure.cpp
 *
 *  Created on: 16.04.2013
 *      Author: cls
 */

#include "GraphStructuralRandMeasure.h"

namespace NetworKit {

GraphStructuralRandMeasure::GraphStructuralRandMeasure() {
	// TODO Auto-generated constructor stub

}

GraphStructuralRandMeasure::~GraphStructuralRandMeasure() {
	// TODO Auto-generated destructor stub
}

double GraphStructuralRandMeasure::getDissimilarity(Graph& G, Partition& first,
		Partition& second) {
	count m = G.numberOfEdges();
	assert (m > 0);

	count e11 = 0; 	// number of connected node pairs for which clusterings aggree
	count e00 = 0;	// number of connected node pairs for which clusterings disagree

	G.forEdges([&](node u, node v){
		if ((first[u] == first[v]) && (second[u] == second[v])) {
			e11 += 1;
		} else if ((first[u] != first[v]) && (second[u] != second[v])) {
			e00 += 1;
		}
	});

	double rand = 1 - (e11 + e00) * (1.0 / m);

	// assert range [0, 1]
	assert (rand <= 1.0);
	assert (rand >= 0.0);
	return rand;
}

/*double GraphStructuralRandMeasure::getDissimilarity(Graph& G, Clustering& zeta1, Graph& G2, Clustering& zeta2) {
	// TODO:
}*/

} /* namespace NetworKit */
