/*
 * JaccardDistance.h
 *
 *  Created on: 17.11.2014
 *      Author: Michael Hamann, Gerd Lindner
 */

#ifndef JACCARDDISTANCE_H_
#define JACCARDDISTANCE_H_

#include "NodeDistance.h"
#include "../graph/Graph.h"
#include "../auxiliary/Timer.h"


namespace NetworKit {

/**
 * @ingroup distmeasures
 * Jaccard distance assigns a distance value to pairs of nodes
 * according to the similarity of their neighborhoods. Note that we define the JaccardDistance as 1-JaccardSimilarity.
 */
class JaccardDistance: public NetworKit::NodeDistance {

public:

	/**
	 * @param G The graph.
	 * @param triangles Edge attribute containing the number of triangles each edge is contained in.
	 */
	JaccardDistance(const Graph& G, const std::vector<int>& triangles);

	/**
	 * REQ: Needs to be called before getEdgeAttribute delivers meaningful results.
	 */
	virtual void preprocess();

	/**
	 * Returns the Jaccard distance between node @a u and node @a v.
	 * @return Jaccard distance between the two nodes.
	 */
	virtual double distance(node u, node v);

	/**
	 * Returns the Jaccard distances between all connected nodes.
	 * @return Vector containing the Jaccard distances between all connected pairs of nodes.
	 */
	std::vector<double> getEdgeAttribute();

protected:
	const std::vector<int>& triangles;
	std::vector<double> jDistance; //result vector

	inline double getJaccardDistance(count degU, count degV, count t);

};

} /* namespace NetworKit */
#endif /* JACCARDDISTANCE_H_ */
