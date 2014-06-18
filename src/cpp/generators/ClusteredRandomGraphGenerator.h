/*
 * ClusteredRandomGraphGenerator.h
 *
 *  Created on: 28.02.2014
 *      Author: cls
 */

#ifndef CLUSTEREDRANDOMGRAPHGENERATOR_H_
#define CLUSTEREDRANDOMGRAPHGENERATOR_H_

#include "StaticGraphGenerator.h"

namespace NetworKit {

/**
 * The ClusteredRandomGraphGenerator class is used to create a clustered random graph.
 * The number of nodes and the number of edges are adjustable as well as the probabilities
 * for intra-cluster and inter-cluster edges.
 */
class ClusteredRandomGraphGenerator: public NetworKit::StaticGraphGenerator {
public:
	/**
	 * Creates a clustered random graph:
	 *
	 * @param[in]	n	number of nodes
	 * @param[in]	k	number of clusters
	 * @param[in]	pin		intra-cluster edge probability
	 * @param[in]	pout	inter-cluster edge probability
	 */
	ClusteredRandomGraphGenerator(count n, count k, double pin, double pout);

	/**
	 * Generates a clustered random graph with the properties given in the constructor.
	 * @return The generated graph.
	 */
	Graph generate() override;

private:

	count n;
	count k;
	double pin;
	double pout;
};

} /* namespace NetworKit */

#endif /* CLUSTEREDRANDOMGRAPHGENERATOR_H_ */
