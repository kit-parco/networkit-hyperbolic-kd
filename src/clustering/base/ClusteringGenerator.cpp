/*
 * ClusteringGenerator.cpp
 *
 *  Created on: 10.12.2012
 *      Author: cls
 */

#include "ClusteringGenerator.h"

namespace EnsembleClustering {

ClusteringGenerator::ClusteringGenerator() {
	// TODO Auto-generated constructor stub

}

ClusteringGenerator::~ClusteringGenerator() {
	// TODO Auto-generated destructor stub
}

Clustering ClusteringGenerator::makeSingletonClustering(Graph& G) {
	int64_t n = G.numberOfNodes();
	Clustering zeta(n);
	zeta.allToSingletons();
	return zeta;
}

Clustering ClusteringGenerator::makeOneClustering(Graph& G) {
	int64_t n = G.numberOfNodes();
	Clustering zeta(n);
	cluster one = zeta.addCluster();
	G.forNodes([&](node v){
		zeta.addToCluster(one, v);
	});
	return zeta;
}

Clustering ClusteringGenerator::makeRandomClustering(Graph& G, int k) {
	// random number generator
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(1, k);
	//
	int64_t n = G.numberOfNodes();
	Clustering zeta(n);

	for (int64_t i = 0; i < k; ++i) {
		zeta.addCluster();
	}

	G.parallelForNodes([&](node v){
		cluster c = dis(gen);
		zeta.addToCluster(c, v);
	});

	assert (zeta.isProper(G));
	return zeta;
}

} /* namespace EnsembleClustering */
