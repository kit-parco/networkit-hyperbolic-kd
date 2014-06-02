/*
 * Graph.h
 *
 *  Created on: 04.02.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu), Henning Meyerhenke (henning.meyerhenke@kit.edu)
 */

#include <vector>
#include <algorithm>
#include <stdexcept>

#include "AbstractGraph.h"
#include "Coordinates.h"
#include "../viz/Point.h"



#include "BasicGraph.h"
namespace NetworKit {
	using Graph = WeightedGraph_T;
}
#define GRAPH_H_

#ifndef GRAPH_H_
#define GRAPH_H_

#include <vector>
#include <algorithm>
#include <stdexcept>

#include "AbstractGraph.h"
#include "Coordinates.h"
#include "../viz/Point.h"

namespace NetworKit {

/**
 * An undirected graph (with optional weights) and parallel iterator methods.
 */
// TODO: add final to class
class Graph : public AbstractGraph {

protected:

	// per node data
	std::vector<count> deg; //!< degree of each node (size of neighborhood)

	// per edge data
	std::vector<std::vector<node> > adja; //!< neighbors/adjacencies
	std::vector<std::vector<edgeweight> > eweights; //!< edge weights

	// user-defined edge attributes
	std::vector<std::vector<std::vector<double> > > edgeMaps_double; // contains edge maps (u, v) -> double
	std::vector<double> edgeAttrDefaults_double; // stores default value for edgeMaps_double[i] at index i

	/**
	 * Return the index of v in the adjacency array of u.
	 */
	index find(node u, node v) const;

public:

	/** GRAPH INTERFACE **/

	/** 
	 * Create a graph of n nodes.
	 */
	Graph(count n = 0, bool weighted = false);

	/* move assignments and more currently removed because a problems with gcc 4.7.1 on phipute1 */
	// Graph(const Graph& other) = default;
	// Graph(Graph&& other) = default;
	// Graph& operator=(Graph&& other) = default;
	// Graph& operator=(const Graph& other) = default;

	/**
	 * Calculate an approximation of the memory used by this graph. Only memory increasing with the
	 * number of edges or nodes of this graph is taken into account. 
	 */
	count getMemoryUsage() const ;

	/**
	 * Try to save some memory by shrinking internal data structures of the graph. Only run this
	 * once you finished editing the graph. Otherwise it will cause unnecessary reallocation of
	 * memory. 
	 */
	void shrinkToFit();

	/** Only to be used from Cython */
	void stealFrom(Graph& input);


	/** NODE MODIFIERS **/

	/**
	 * Add a new node to the graph and return it.
	 */
	node addNode() override;

	/**
	 * Add a new node to the graph with coordinates @a x and @y and return it.
	 */
	node addNode(float x, float y) override;


	/** NODE PROPERTIES **/
 
	/**
	 * @return true if the node is isolated (= degree is 0)
	 */
	bool isIsolated(node v) const;

	/**
	 * Return the number of neighbors for node v.
	 */
	count degree(node v) const;

	/**
	 * @return Weighted degree of @a v.
	 */
	edgeweight weightedDegree(node v) const;

	/**
	 * @return Random neighbor of @a v. None if degree is zero.
	 */
	node randomNeighbor(node v) const;


	/** EDGE MODIFIERS **/

	/**
	 * Insert an undirected edge between two nodes.
	 */
	void addEdge(node u, node v, edgeweight weight = defaultEdgeWeight);

	/**
	 * Remove undirected edge between two nodes.
	 */
	void removeEdge(node u, node v);

	/**
	 * Check if undirected edge {u,v} exists in G
	 *
	 */
	bool hasEdge(node u, node v) const;

	/**
	 * Merges edge {u,v} to become a supernode. Edges to u and v are
	 * rewired, multiple edges merged and their weights added.
	 * The vertex weights of @a u and @a v are added.
	 * A self-loop is only created if @a discardSelfLoop is set to false.
	 *
	 * @return New node that has been created if u != v. Otherwise none.
	 */
	node mergeEdge(node u, node v, bool discardSelfLoop = true);


	/** GLOBAL PROPERTIES **/

	/** 
	 * Return true if this graph supports directed edges.
	 */
	bool isDirected() const;


	/** EDGE ATTRIBUTES **/

	/**
	 * Return edge weight.
	 *
	 * Return 0 if edge does not exist.
	 */
	edgeweight weight(node u, node v) const;

	/**
	 * Set the weight of an edge. If the edge does not exist,
	 * it will be inserted.
	 *
	 * @param[in]	u	endpoint of edge
	 * @param[in]	v	endpoint of edge
	 * @param[in]	weight	edge weight
	 */
	void setWeight(node u, node v, edgeweight w);

	/**
	 * Increase the weight of an edge. If the edge does not exist,
	 * it will be inserted.
	 *
	 * @param[in]	u	endpoint of edge
	 * @param[in]	v	endpoint of edge
	 * @param[in]	weight	edge weight
	 */
	void increaseWeight(node u, node v, edgeweight w);

	/**
	 * Add new edge map for an attribute of type double.
	 */
	int addEdgeAttribute_double(double defaultValue);

	/**
	 * @return attribute of type double for an edge.
	 *
	 * @param[in]	u	node
	 * @param[in]	v	node
	 * @param[in]	attrId	attribute id
	 */
	double attribute_double(node u, node v, int attrId) const;

	/**
	 * Set edge attribute of type double If the edge does not exist,
	 * it will be inserted.
	 *
	 * @param[in]	u	endpoint of edge
	 * @param[in]	v	endpoint of edge
	 * @param[in]	attr	double edge attribute
	 */
	void setAttribute_double(node u, node v, int attrId, double attr);


	/** EDGE ITERATORS **/

	/**
	 * Iterate over all edges of the graph and call f (lambda closure).
	 */
	void forEdges(FEdge f) const;

	/**
	 * Iterate in parallel over all edges of the graph and call f (lambda closure).
	 */
	virtual void parallelForEdges(FEdge f) const;

	/**
	 * Iterate over all edges of the graph and call f (lambda closure).
	 *
	 * Handler takes arguments (u, v, w) where u and v are the nodes of the edge and w is its weight.
	 */
	virtual void forWeightedEdges(FEdgeWeighted f) const;

	/**
	 * Iterate over all edges of the graph and call f (lambda closure).
	 *
	 * Handler takes arguments (u, v, w) where u and v are the nodes of the edge and w is its weight.
	 */
	virtual void parallelForWeightedEdges(FEdgeWeighted f) const;

	/**
	 * Iterate over all edges of the const graph and call f (lambda closure).
	 *
	 *	@param[in]	attrId		attribute id
	 *	@param[in]	f 		takes arguments (u, v, a) where a is an edge attribute of edge {u, v}
	 *
	 */
	virtual void forEdgesWithAttribute_double(int attrId, FEdgeWeighted f) const;


	/** NEIGHBORHOOD ITERATORS **/

	/**
	 * Iterate over all neighbors of a node and call f (lamdba closure).
	 */
	void forNeighborsOf(node u, FNode f) const;

	/**
	 * Iterate over all edge weights of a node and call f (lamdba closure).
	 */
	virtual void forWeightedNeighborsOf(node u, FNodeWeighted f) const;

	/**
	 * Iterate over all incident edges of a node and call f (lamdba closure).
	 */
	virtual void forEdgesOf(node u, FEdge f) const;

	/**
	 * Iterate over all incident edges of a node and call f (lamdba closure).
	 *
	 * Handle takes parameters (u, v, w) where w is the edge weight.
	 *
	 */
	virtual void forWeightedEdgesOf(node u, FEdgeWeighted f) const;


	/** REDUCTION ITERATORS **/

	/**
	 * Iterate in parallel over all edges and sum (reduce +) the values returned by the f
	 */
	virtual double parallelSumForWeightedEdges(FEdgeWeightedSum f) const;

};

} /* namespace NetworKit */


/** EDGE ITERATORS **/

inline void NetworKit::Graph::forEdges(NetworKit::FEdge f) const {
	for (node u = 0; u < z; ++u) {
		for (node v : this->adja[u]) {
			if (u >= v) { // {u, v} instead of (u, v); if v == none, u > v is not fulfilled
				f(u, v);
			}
		}
	}
}

inline void NetworKit::Graph::parallelForEdges(NetworKit::FEdge f) const {
	#pragma omp parallel for
	for (node u = 0; u < z; ++u) {
		for (node v : this->adja[u]) {
			if (u >= v) { // {u, v} instead of (u, v); if v == none, u > v is not fulfilled
				f(u, v);
			}
		}
	}
}


inline void NetworKit::Graph::forWeightedEdges(NetworKit::FEdgeWeighted f) const {
	for (node u = 0; u < z; ++u) {
		for (index vi = 0; vi < adja[u].size(); ++vi) {
			node v = this->adja[u][vi];
			if (u >= v) { // {u, v} instead of (u, v); if v == none, u > v is not fulfilled
				if (weighted) {
					edgeweight w = this->eweights[u][vi];
					f(u, v, w);
				} else {
					f(u, v, defaultEdgeWeight);
				}
			}
		}
	}
}

inline void NetworKit::Graph::parallelForWeightedEdges(NetworKit::FEdgeWeighted f) const {
	#pragma omp parallel for
	for (node u = 0; u < z; ++u) {
		for (index vi = 0; vi < adja[u].size(); ++vi) {
			node v = this->adja[u][vi];
			if (u >= v) { // {u, v} instead of (u, v); if v == none, u > v is not fulfilled
				if (weighted) {
					edgeweight w = this->eweights[u][vi];
					f(u, v, w);
				} else {
					f(u, v, defaultEdgeWeight);
				}
			}
		}
	}
}

inline void NetworKit::Graph::forEdgesWithAttribute_double(int attrId, NetworKit::FEdgeWeighted f) const {
	std::vector<std::vector<double> > edgeMap = this->edgeMaps_double[attrId];
	for (node u = 0; u < z; ++u) {
		for (index vi = 0; vi < (index) adja[u].size(); ++vi) {
			node v = this->adja[u][vi];
			double attr = edgeMap[u][vi];
			if (u >= v) { // {u, v} instead of (u, v); if v == none, u > v is not fulfilled
				f(u, v, attr);
			}
		}
	}
}


/** NEIGHBORHOOD ITERATORS **/

inline void NetworKit::Graph::forNeighborsOf(node u, NetworKit::FNode f) const {
	for (node v : this->adja[u]) {
		if (v != none) {
			f(v);
		}
	}
}

inline void NetworKit::Graph::forWeightedNeighborsOf(node u, NetworKit::FNodeWeighted f) const {
	if (weighted) {
		for (index i = 0; i < (index) adja[u].size(); i++) {
			node v = adja[u][i];
			if (v != none) {
				edgeweight ew = eweights[u][i];
				f(v, ew);
				assert(ew == weight(u, v));
			}
		}
	} else {
		for (index i = 0; i < (index) adja[u].size(); i++) {
			node v = adja[u][i];
			if (v != none) {
				f(v, defaultEdgeWeight);
			}
		}
	}
}

inline void NetworKit::Graph::forEdgesOf(node u, NetworKit::FEdge f) const {
	for (node v : this->adja[u]) {
		if (v != none) {
			f(u, v);
		}
	}
}

inline void NetworKit::Graph::forWeightedEdgesOf(node u, NetworKit::FEdgeWeighted f) const {
	for (index i = 0; i < adja[u].size(); ++i) {
		node v = adja[u][i];
		if (v != none) {
			if (weighted) {
				edgeweight w = this->eweights[u][i];
				f(u, v, w);
			} else {
				f(u, v, defaultEdgeWeight);
			}
		}
	}
}

/** REDUCTION ITERATORS **/

inline double NetworKit::Graph::parallelSumForWeightedEdges(NetworKit::FEdgeWeightedSum f) const {
	double sum = 0.0;
	#pragma omp parallel for reduction(+:sum)
	for (node u = 0; u < z; ++u) {
		for (index i = 0; i < this->adja[u].size(); i++) {
			node v = this->adja[u][i];
			edgeweight ew = this->eweights[u][i];
			if (u >= v) { // {u, v} instead of (u, v); if v == none, u > v is not fulfilled
				sum += f(u, v, ew);
			}
		}
	}
	return sum;
}

#endif /* GRAPH_H_ */
