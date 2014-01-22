/*
 * JaccardMeasure2.cpp
 *
 *  Created on: 19.01.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include "JaccardMeasure2.h"

namespace NetworKit {

JaccardMeasure2::JaccardMeasure2() {
	// TODO Auto-generated constructor stub

}

JaccardMeasure2::~JaccardMeasure2() {
	// TODO Auto-generated destructor stub
}

double JaccardMeasure2::getDissimilarity(Graph& G, Partition& first,
		Partition& second) {

	int64_t n = G.numberOfNodes();
	assert (n > 0);

	int64_t s11 = 0; 	// number of node pairs for which clusterings aggree
	int64_t s00 = 0;	// number of node pairs for which clusterings disagree

	G.forNodePairs([&](node u, node v){
		if ((first[u] == first[v]) && (second[u] == second[v])) {
			s11 += 1;
		} else if ((first[u] != first[v]) && (second[u] != second[v])) {
			s00 += 1;
		}
	});

	double jaccard;
	double divisor = n * (n - 1) - 2 * s00;
	if (divisor > 0) {
		jaccard = 1 - ((2 * s11) / divisor);
	} else {
		jaccard = 0;
	}

	// assert range [0, 1]
	assert (jaccard <= 1.0);
	assert (jaccard >= 0.0);
	return jaccard;

}

} /* namespace NetworKit */
