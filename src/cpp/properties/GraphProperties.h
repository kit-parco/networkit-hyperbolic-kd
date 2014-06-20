/*
 * GraphProperties.h
 *
 *  Created on: 03.06.2013
 *      Author: cls
 */

#ifndef GRAPHPROPERTIES_H_
#define GRAPHPROPERTIES_H_


#include "../graph/Graph.h"
#include "../io/METISGraphReader.h"

namespace NetworKit {

/**
 * DEPRECATED: Implement algorithms in their own classes.
 * 
 * Collection of methods for basic network properties.
 */
class GraphProperties {
protected:
	/**
	 * @return Degree assortativity of the graph @a G.
	 * Based on Eq. (4) in Newman: Assortative mixing in networks.
	 * URL: http://arxiv.org/pdf/cond-mat/0205405.pdf.
	 * A similar description of this algorithm can be found in
	 * Newman: Networks. An Introduction. Chapter 8.7.
	 */
	static double degreeAssortativitySlower(const Graph& G, bool useWeights = false);


public:
	GraphProperties();
	virtual ~GraphProperties();

	static std::vector<count> degreeDistribution(const Graph& G);

	static std::vector<unsigned int> degreeSequence(const Graph& G); // TODO: revert to count when cython issue fixed


	/**
	 * The local clustering coefficient for a node is the number of edges among its
	 * neighbors divided by the number of possible edges.
	 * 	For an undirected graph where N(v) does not include v itself:
	 * 		$c_v := \frac{2 \cdot |E(N(v))| }{\deg(v) ( \deg(v) - 1 )}$
	 *
	 *
	 * @param[in]	G	the graph
	 * @param[out]		node -> local clustering coefficient
	 */
	static std::vector<double> localClusteringCoefficients(const Graph& G);


	/**
	 * The average local clustering coefficient for the graph.
	 * 		$\frac{1}{n} \cdot \sum_{v \in V} c_v$
	 *
	 * @param[in]	G	the graph
	 */
	static double averageLocalClusteringCoefficient(const Graph& G);

	static std::vector<double> localClusteringCoefficientPerDegree(const Graph& G);

	static std::pair<count, count> minMaxDegree(const Graph& G);

	static double averageDegree(const Graph& G);

	/**
	 * @return Degree assortativity of the graph @a G.
	 * Based on description in
	 * Newman: Networks. An Introduction. Chapter 8.7.
	 */
	static double degreeAssortativity(const Graph& G, bool useWeights = false);
};

} /* namespace NetworKit */
#endif /* GRAPHPROPERTIES_H_ */
