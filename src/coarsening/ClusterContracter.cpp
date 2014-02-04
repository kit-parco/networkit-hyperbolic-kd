/*
 * ClusterContracter.cpp
 *
 *  Created on: 30.10.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include "ClusterContracter.h"

namespace NetworKit {

ClusterContracter::ClusterContracter() {
	// TODO Auto-generated constructor stub

}

ClusterContracter::~ClusterContracter() {
	// TODO Auto-generated destructor stub
}

std::pair<Graph, NodeMap<node> > ClusterContracter::run(Graph& G, Partition& zeta) {

	Graph Gcon(0); // empty graph
	Gcon.markAsWeighted(); // Gcon will be a weighted graph

	IndexMap<index, node> clusterToSuperNode(zeta.upperBound(), none); // there is one supernode for each cluster

	// populate map cluster -> supernode
	G.forNodes([&](node v){
		index c = zeta[v];//zeta.clusterOf(v);Cluster
		if (! clusterToSuperNode.hasBeenSet(c)) {
			clusterToSuperNode[c] = Gcon.addNode(); // TODO: probably does not scale well, think about allocating ranges of nodes
		}
	});


	int64_t n = G.numberOfNodes();
	NodeMap<node> nodeToSuperNode(n);

	// set entries node -> supernode
	G.parallelForNodes([&](node v){
		nodeToSuperNode[v] = clusterToSuperNode[zeta.subsetOf(v)];
	});


	// iterate over edges of G and create edges in Gcon or update edge and node weights in Gcon
	G.forWeightedEdges([&](node u, node v, edgeweight ew) {
		node su = nodeToSuperNode[u];
		node sv = nodeToSuperNode[v];
		TRACE("edge (", su, ", ", sv, ")");
		if (zeta.subsetOf(u) == zeta.subsetOf(v)) {
			// add edge weight to supernode (self-loop) weight
			Gcon.setWeight(su, su, Gcon.weight(su, su) + ew);
		} else {
			// add edge weight to weight between two supernodes (or insert edge)
			Gcon.setWeight(su, sv, Gcon.weight(su, sv) + ew);
		}
	}); // TODO: parallel?

	return std::make_pair(Gcon, nodeToSuperNode);

}

}
