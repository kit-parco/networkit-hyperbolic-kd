/*
 * Matching.h
 *
 *  Created on: 03.12.2012
 *      Author: cls
 */

#ifndef MATCHING_H_
#define MATCHING_H_

# include "../aux/log.h"
#include "../graph/Graph.h"
#include "../graph/NodeMap.h"

namespace EnsembleClustering {

class Matching : public NodeMap<node> {


public:

	/**
	 * Construct new matching.
	 *
	 * @param[in]	n 	maximum number of nodes
	 */
	Matching(int64_t n);

	/**
	 * Destructor.
	 */
	virtual ~Matching();


	/**
	 * Set two nodes as eachothers matching
	 * partners.
	 *
	 *
	 */
	void match(const node& u, const node& v);


	/**
	 * Reset the two nodes to unmatched.
	 */
	void unmatch(const node& u, const node& v);


	/**
	 * Check if node is matched
	 *
	 * @param[in]	u 	a node
	 * @param[out]		true if u is matched
	 */
	bool isMatched(const node& u) const;


	/**
	 * Check if the two nodes are matched.
	 *
	 */
	bool areMatched(const node& u, const node& v) const;

	/**
	 * Check whether this is a proper matching
	 * in the graph, i.e. no two edges are adjacent.
	 *
	 *
	 * @paramt[in]	G	a graph
	 * @param[out]		true if this is a proper matching
	 */
	bool isProper(Graph& G) const;


	/**  copy semantics **/


	/**
	 * Assignment operator.
	 */
	Matching& operator=(const Matching& from);


	/**
	 * Properly copy this object.
	 */
	 void clone(const Matching& from);


	/**
	 * Properly destruct this object.
	 */
	void dispose();



};

} /* namespace EnsembleClustering */
#endif /* MATCHING_H_ */
