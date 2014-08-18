/*
 * SSSP.h
 *
 *  Created on: 15.04.2014
 *      Author: cls
 */

#ifndef SSSP_H_
#define SSSP_H_

#include <set>
#include <stack>

#include "Graph.h"

namespace NetworKit {

/**
 * @ingroup graph
 * Abstract base class for single-source shortest path algorithms.
 */
class SSSP {

public:

	/**
	 * Creates the SSSP class for @a G and source @a s.
	 *
	 * @param G The graph.
	 * @param s The source node.
	 */
	SSSP(const Graph& G, node s, bool storePaths=true, bool storeStack=false);

	/** Computes the shortest paths from the source to all other nodes. */
	virtual void run(node t = none) = 0;

	/**
	 * Returns a vector of weighted distances from the source node, i.e. the
 	 * length of the shortest path from the source node to any other node.
 	 *
 	 * @return The weighted distances from the source node to any other node in the graph.
	 */
	virtual std::vector<edgeweight> getDistances() const;

	/**
	 * Returns the distance from the source node to @a t.
	 * @param  t Target node.
	 * @return The distance from source to target node @a t.
	 */
	edgeweight distance(node t) const;

	/**
	 * Returns the number of shortest paths between the source node and @a t.
	 * @param  t Target node.
	 * @return The number of shortest paths between source and @a t.
	 */
	count numberOfPaths(node t) const;

	/**
	 * Returns the predecessor nodes of @a t on all shortest paths from source to @a t.
	 * @param t Target node.
	 * @return The predecessors of @a t on all shortest paths from source to @a t.
	 */
	std::vector<node> getPredecessors(node t) const;

	/**
	 * Returns a shortest path from source to @a t and an empty path if source and @a t are not connected.
	 *
	 * @param t Target node.
	 * @param forward If @c true (default) the path is directed from source to @a t, otherwise the path is reversed.
	 * @return A shortest path from source to @a t or an empty path.
	 */
	virtual std::vector<node> getPath(node t, bool forward=true) const;

	/**
	 * Returns all shortest paths from source to @a t and an empty set if source and @a t are not connected.
	 *
	 * @param t Target node.
	 * @param forward If @c true (default) the path is directed from source to @a t, otherwise the path is reversed.
	 * @return All shortest paths from source node to target node @a t.
	 */
	virtual std::set<std::vector<node> > getPaths(node t, bool forward=true) const;


	/**
	* Returns a stack of nodes ordered in decreasing distance from the source
	*
	* @return stack of nodes
	*/
	virtual std::stack<node> getStack() const;

protected:

	const Graph& G;
	const node source;
	std::vector<edgeweight> distances;
	std::vector<std::vector<node> > previous; // predecessors on shortest path
	std::vector<count> npaths;

	std::stack<node> stack;

	bool storePaths;		//!< if true, paths are reconstructable and the number of paths is stored
	bool storeStack;		//!< if true, store a stack of nodes ordered in decreasing distance from the source
};

inline edgeweight SSSP::distance(node t) const {
	return distances[t];
}

inline count SSSP::numberOfPaths(node t) const {
	if (! storePaths) {
		throw std::runtime_error("number of paths have not been stored");
	}
	return npaths[t];
}

inline std::vector<node> SSSP::getPredecessors(node t) const {
	if (! storePaths) {
		throw std::runtime_error("predecessors have not been stored");
	}
	return previous[t];
}

} /* namespace NetworKit */

#endif /* SSSP_H_ */
