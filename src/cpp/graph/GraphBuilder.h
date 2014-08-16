/*
 * GraphBuilder.h
 *
 *  Created on: 15.07.2014
 *      Author: Marvin Ritter (marvin.ritter@gmail.com)
 */

#ifndef GRAPH_BUILDER_H
#define GRAPH_BUILDER_H

#include <vector>

#include "../Globals.h"
#include "Graph.h"

namespace NetworKit {

class GraphBuilder {
protected:
	count n; //!< current number of nodes

	bool weighted; //!< true if the graph will be weighted, false otherwise
	bool directed; //!< true if the graph will be directed, false otherwise

	std::vector< std::vector<node> > halfEdges;
	std::vector< std::vector<edgeweight> > halfEdgeWeights;

	index indexHalfEdgeArray(node u, node v) const;

public:
	GraphBuilder(count n = 0, bool weighted = false, bool directed = false);

	/**
	 * Returns <code>true</code> if this graph supports edge weights other than 1.0.
	 * @return <code>true</code> if this graph supports edge weights other than 1.0.
	 */
	bool isWeighted() const { return weighted; }

	/**
	 * Return @c true if this graph supports directed edges.
	 * @return @c true if this graph supports directed edges.
	 */
	bool isDirected() const { return directed; }

	/**
	 * Return <code>true</code> if graph contains no nodes.
	 * @return <code>true</code> if graph contains no nodes.
	 */
	bool isEmpty() const { return n == 0; }

	/**
	 * Return the number of nodes in the graph.
	 * @return The number of nodes.
	 */
	count numberOfNodes() const { return n; }

 	/**
	 * Get an upper bound for the node ids in the graph.
	 * @return An upper bound for the node ids.
	 */
	index upperNodeIdBound() const { return n; }

	/**
	 * Add a new node to the graph and return it.
	 * @return The new node.
	 */
	node addNode();
	
	/**
	 * Insert an edge between the nodes @a u and @a v. If the graph is weighted you can optionally
	 * set a weight for this edge. The default weight is 1.0.
	 * @param u Endpoint of edge.
	 * @param v Endpoint of edge.
	 * @param weight Optional edge weight.
	 */
	void addEdge(node u, node v, edgeweight ew = defaultEdgeWeight);

	/**
	 * Set the weight of an edge. If the edge does not exist,
	 * it will be inserted.
	 *
	 * @param[in]	u	endpoint of edge
	 * @param[in]	v	endpoint of edge
	 * @param[in]	weight	edge weight
	 */
	void setWeight(node u, node v, edgeweight ew);

	/**
	 * Increase the weight of an edge. If the edge does not exist,
	 * it will be inserted.
	 *
	 * @param[in]	u	endpoint of edge
	 * @param[in]	v	endpoint of edge
	 * @param[in]	weight	edge weight
	 */
	void increaseWeight(node u, node v, edgeweight ew);

	/**
	 * Generates a Graph instance. The graph builder will be reseted at the end.
	 */
	Graph toGraph(bool parallel = true) { return parallel ? toGraphParallel() : toGraphSequential(); }

	/**
	 * Iterate over all nodes of the graph and call @a handle (lambda closure).
	 *
	 * @param handle Takes parameter <code>(node)</code>.
	 */
	template<typename L> void forNodes(L handle) const;

	/**
	 * Iterate randomly over all nodes of the graph and call @a handle (lambda closure).
	 *
	 * @param handle Takes parameter <code>(node)</code>.
	 */
	template<typename L> void parallelForNodes(L handle) const;

	/**
	 * Iterate over all undirected pairs of nodes and call @a handle (lambda closure).
	 *
	 * @param handle Takes parameters <code>(node, node)</code>.
	 */
	template<typename L> void forNodePairs(L handle) const;


	/**
	 * Iterate over all undirected pairs of nodes in parallel and call @a handle (lambda closure).
	 *
	 * @param handle Takes parameters <code>(node, node)</code>.
	 */
	template<typename L> void parallelForNodePairs(L handle) const;

private:
	Graph toGraphParallel();
	Graph toGraphSequential();
};

template<typename L>
void GraphBuilder::forNodes(L handle) const {
	for (node v = 0; v < n; v++) {
		handle(v);
	}
}

template<typename L>
void GraphBuilder::parallelForNodes(L handle) const {
	#pragma omp parallel for schedule(dynamic)
	for (node v = 0; v < n; v++) {
		handle(v);
	}
}

template<typename L>
void GraphBuilder::forNodePairs(L handle) const {
	for (node u = 0; u < n; u++) {
		for (node v = u + 1; v < n; v++) {
			handle(u, v);
		}
	}
}

template<typename L>
void GraphBuilder::parallelForNodePairs(L handle) const {
	#pragma omp parallel for schedule(dynamic)
	for (node u = 0; u < n; u++) {
		for (node v = u + 1; v < n; v++) {
			handle(u, v);
		}
	}
}

} /* namespace NetworKit */

#endif /* GRAPH_BUILDER_H */
