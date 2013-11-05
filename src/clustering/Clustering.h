/*
 * Clustering.h
 *
 *  Created on: 31.10.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef CLUSTERING_H_
#define CLUSTERING_H_



#include "../graph/NodeMap.h"
#include <cassert>



namespace NetworKit {

typedef index cluster;	//!< cluster is represented as a 0-based index

class Clustering: public NodeMap<cluster> {

protected:

	cluster nextCluster;	//!< next free cluster id for new cluster
	std::string name;
	cluster upperIdBound;	//!< upper bound for cluster ids

	inline cluster getNextCluster() {
		// TODO: performance - is this a bottleneck?
		cluster c = this->nextCluster;
		this->nextCluster++;
		this->upperIdBound++;
		return c;
	}

	/**
	 * Check if clustering can hold a valid entry for the node because
	 * it is in the range mapped.
	 */
	bool isInRange(node v);

public:

	/*
	 * Construct an empty clustering.
	 */
	Clustering();

	// FIXME: the constructor is outdated and will lead to an error if nodes from the graph have been deleted
	/**
	 * Construct new clustering.
	 *
	 * @param[in]	z	upper node id bound
	 */
	Clustering(index z);

	virtual ~Clustering();


	/**
	 * Set a human-readable identifier for the instance.
	 */
	void setName(std::string name);


	/**
	 * Get the human-readable identifier
	 */
	std::string getName() const;

	/**
	 *  Index operator.
	 *
	 *  @param[in]	u	a node
	 */
	inline cluster& operator [](const node& u) {
		return this->data[u];
	}
	/**
	 * Index operator for const instances of this class.
	 *
	 * @param[in]	u 	a node
	 */
	inline const cluster& operator [](const node& u) const {
		return this->data[u];
	}

	/**
	 * Return the cluster (id) in which a node
	 * is contained.
	 */
	inline cluster clusterOf(node u) const {
		assert (u < this->numberOfNodes());
		return this->data[u];
	}

	/**
	 * Call this before assigning nodes to cluster ids.
	 */
	cluster addCluster();

	/**
	 * Add a (previously unassigned) node to a cluster
	 */
	void addToCluster(cluster c, node u);

	/**
	 * Move a (previously assigned) node to a cluster.
	 */
	void moveToCluster(cluster c, node u);

	/**
	 * Creates a singleton cluster containing the node.
	 */
	void toSingleton(node u);


	/**
	 * Assigns every node to a singleton cluster.
	 * Cluster id is equal to node id.
	 */
	void allToSingletons();

	/**
	 * Assigns the nodes from both clusters to a new cluster.
	 */
	void mergeClusters(cluster c, cluster d);

	/**
	 * Check whether this clustering is a proper clustering of
	 * the graph, i.e. a disjoint partition of the whole node set.
	 *
	 */
	bool isProper(Graph& G);


	/**
	 * Check if clustering is a 1-clustering,
	 * i.e. every node is assigned to the same cluster.
	 */
	bool isOneClustering(Graph& G);


	/**
	 * Check if clustering is a singleton clustering,
	 * i.e. every node is assigned to a different cluster.
	 */
	bool isSingletonClustering(Graph& G);

	/**
	 * Get the current number of clusters in this clustering.
	 */
	count numberOfClusters() const;




	/**
	 * Return an upper bound for the cluster ids that have been assigned.
	 *
	 * (This is the maximum id + 1.)
	 */
	cluster upperBound() const;

	/**
	 * Return a lower bound for the cluster ids that have been assigned.
	 */
	cluster lowerBound() const;

	/**
	 * Set an upper bound for the cluster ids.
	 */
	void setUpperBound(cluster id);


	/**
	 * Change cluster IDs to be consecutive, starting at 0.
	 */
	void compact();


	/**
	 * Check if clustering assigns a valid cluster to the node.
	 */
	bool contains(node v);


	/**
	 * Check if two nodes belong to the same cluster
	 */
	bool inSameCluster(node u, node v);


	/**
	 * @return Quotient of size of largest cluster and ceil(average cluster size).
	 */
	float getImbalance();


	/**
	 * Check if this clustering equals another clustering (with respect to a graph).
	 * Criterion for equality:
	 *
	 * 		$$\zeta_1(G) = \zeta_2(G) \iff  \for \{u, v\} \in E: \zeta_1(u) = \zeta_1(v) \implies \zeta_2(u) = \zeta_2(v) \and \zeta_1(u) \neq \zeta_1(v) \implies \zeta_2(u) \neq \zeta_2(v) $$
	 */
	bool equals(Clustering& other, Graph& G);


	/**
	 * Iterate over all entries (node, cluster) and execute callback function (lambda closure).
	 */
	template<typename Callback> void forEntries(Callback func);


	/**
	 * Iterate over all entries (node, cluster) in parallel and execute callback function (lambda closure).
	 */
	template<typename Callback> void parallelForEntries(Callback handle);

	/**
	 * Iterate over all entries (node, cluster) and execute callback function (lambda closure).
	 */
	template<typename Callback> void forEntries(Callback func) const;


	/**
	 * Iterate over all entries (node, cluster) in parallel and execute callback function (lambda closure).
	 */
	template<typename Callback> void parallelForEntries(Callback handle) const;


	/**
	 * Get a list of cluster sizes. Indices do not necessarily correspond to cluster ids.
	 */
	std::vector<count> clusterSizes() const;


	/**
	 * Get a map from cluster id to size of the cluster.
	 */
	std::map<cluster, count> clusterSizeMap() const;

	/**
	 * Append a node.
	 */
	void append(node u);


	/**
	 * Get the members of a specific cluster.
	 */
	std::vector<node> getMembers(const cluster C) const;


	/**
	 * Compute communication graph (sometimes called quotient graph) induced by
	 * graph @a graph and this clustering/partition.
	 */
	Graph communicationGraph(const Graph& graph);


	/**
	 * @return Sum of weights of edges that connect @a v with nodes in cluster @a cid.
	 */
	edgeweight weightedDegreeWithCluster(const Graph& graph, node v, cluster cid) const;
};

} /* namespace NetworKit */

template<typename Callback>
inline void NetworKit::Clustering::forEntries(Callback func) {
	for (node v = 0; v < this->n; v += 1) {
		cluster c = (*this)[v];
		func(v, c);
	}

}

template<typename Callback>
inline void NetworKit::Clustering::parallelForEntries(
		Callback handle) {
	#pragma omp parallel for
	for (node v = 0; v < this->n; v += 1) {
		cluster c = (*this)[v];
		handle(v, c);
	}
}


template<typename Callback>
inline void NetworKit::Clustering::forEntries(Callback func) const {
	for (node v = 0; v < this->n; v += 1) {
		cluster c = (*this)[v];
		func(v, c);
	}

}

template<typename Callback>
inline void NetworKit::Clustering::parallelForEntries(
		Callback handle) const {
	#pragma omp parallel for
	for (node v = 0; v < this->n; v += 1) {
		cluster c = (*this)[v];
		handle(v, c);
	}
}


#endif /* CLUSTERING_H_ */
