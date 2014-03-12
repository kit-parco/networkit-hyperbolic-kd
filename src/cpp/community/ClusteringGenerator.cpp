/*
 * ClusteringGenerator.cpp
 *
 *  Created on: 10.12.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include "ClusteringGenerator.h"
#include "../auxiliary/Random.h"
#include "GraphClusteringTools.h"

namespace NetworKit {

ClusteringGenerator::ClusteringGenerator() {

}

ClusteringGenerator::~ClusteringGenerator() {

}

Partition ClusteringGenerator::makeSingletonClustering(Graph& G) {
	count n = G.upperNodeIdBound();
	Partition zeta(n);
	zeta.allToSingletons();
	return zeta;
}

Partition ClusteringGenerator::makeOneClustering(Graph& G) {
	count n = G.upperNodeIdBound();
	Partition zeta(n);
	zeta.allToOnePartition();
	return zeta;
}

Partition ClusteringGenerator::makeRandomClustering(Graph& G, count k) {
	count n = G.upperNodeIdBound();
	Partition zeta(n);

	zeta.setUpperBound(k);

	G.forNodes([&](node v) { //parallel
		index c = Aux::Random::integer(k-1);
		zeta.addToSubset(c, v);
	});

	assert (zeta.numberOfSubsets() == k);
	assert (GraphClusteringTools::isProperClustering(G, zeta));
	return zeta;
}

Partition ClusteringGenerator::makeContinuousBalancedClustering(Graph& G, count k) {
	count n = G.upperNodeIdBound(); 
	Partition clustering(n);
	clustering.setUpperBound(k);

	std::vector<count> blockSize(k, 0);

	// compute block sizes
	for (index block = 0; block < k; ++block) {
		blockSize[block] = n / k + (n % k > block);
	}

	// compute prefix sums of block sizes
	for (index block = 1; block < k; ++block) {
		blockSize[block] += blockSize[block-1];
	}

	// fill clustering according to blocks
	node v = 0;
	for (index block = 0; block < k; ++block) {
		while (v < blockSize[block]) {
			clustering.addToSubset(block,v);
			++v;
		}
	}

	return clustering;
}

} /* namespace NetworKit */
