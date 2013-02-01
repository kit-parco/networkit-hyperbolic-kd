/*
 * Generator.cpp
 *
 *  Created on: 05.12.2012
 *      Author: cls
 */

#include "GraphGenerator.h"

namespace EnsembleClustering {

GraphGenerator::GraphGenerator() {
	// TODO Auto-generated constructor stub

}

GraphGenerator::~GraphGenerator() {
	// TODO Auto-generated destructor stub
}


// TODO: parallel? is insertEdge thread safe?


Graph GraphGenerator::makeErdosRenyiGraph(int64_t n, double p) {
	Aux::RandomProbability randP;
	Graph G(n);
	for (node u = 1; u <= n; ++u) {
		for (node v = u + 1; v <= n; ++v) {
			if (randP.generate() <= p) {
				G.insertEdge(u, v);
			}
		}
	}
	return G;
}

Graph GraphGenerator::makeRandomGraph(int64_t n, double p) {
	return this->makeErdosRenyiGraph(n, p);	// alias
}

Graph GraphGenerator::makeCircularGraph(int64_t n) {
	// TODO: modernize
	Graph G(n);
	for (int i = 0; i < n; ++i) {
		G.insertEdge(i + 1, ((i+1) % n) + 1);
	}
	return G;
}

Graph GraphGenerator::makeCompleteGraph(int64_t n) {
	// TODO: modernize
	Aux::RandomProbability randP;
	Graph G(n);
	for (node u = 1; u <= n; ++u) {
		for (node v = u + 1; v <= n; ++v) {
			G.insertEdge(u, v);
		}
	}
	return G;
}



Graph GraphGenerator::makeClusteredRandomGraph(int64_t n, int64_t k, double pin, double pout) {
	assert(pin >= pout);

	Graph G(n);
	Aux::RandomProbability randP;
	Aux::RandomInteger randInt(1, k);
	// assign nodes evenly to clusters
	Clustering zeta(n);
	G.forallNodes([&](node v){
		cluster c = randInt.generate();
		zeta.addToCluster(c, v);
	});

	assert (zeta.numberOfClusters() == k);

	for (node u = 1; u <= n; ++u) {
		for (node v = u + 1; v <= n; ++v) {
			if (zeta.clusterOf(u) == zeta.clusterOf(v)) {
				if (randP.generate() <= pin) {
					G.insertEdge(u, v);
				}
			} else {
				if (randP.generate() <= pout) {
					G.insertEdge(u, v);
				}
			}
		}
	}

	return G;
}

std::pair<Graph, Clustering> GraphGenerator::makeClusteredRandomGraphWithReferenceClustering(
		int64_t n, int64_t k, double pin, double pout) {
	assert(pin >= pout);

	Graph G(n);
	Aux::RandomProbability randP;
	Aux::RandomInteger randInt(1, k);
	// assign nodes evenly to clusters
	Clustering zeta(n);
	G.forallNodes([&](node v){
		cluster c = randInt.generate();
		zeta.addToCluster(c, v);
	});

	assert (zeta.numberOfClusters() == k);

	for (node u = 1; u <= n; ++u) {
		for (node v = u + 1; v <= n; ++v) {
			if (zeta.clusterOf(u) == zeta.clusterOf(v)) {
				if (randP.generate() <= pin) {
					G.insertEdge(u, v);
				}
			} else {
				if (randP.generate() <= pout) {
					G.insertEdge(u, v);
				}
			}
		}
	}

	return std::make_pair(G, zeta);
}

Graph GraphGenerator::makeClusteredRandomGraph(Clustering& zeta, double pin,
		double pout) {
	assert (pin >= pout);

	int64_t n = zeta.numberOfNodes();
	Graph G(n);

	Aux::RandomProbability randP;
	G.forallNodePairs([&](node u, node v){
		if (zeta.inSameCluster(u, v)) {
			if (randP.generate() <= pin) {
				G.insertEdge(u, v);
			}
		} else {
			if (randP.generate() <= pout) {
				G.insertEdge(u, v);
			}
		}
	});

	return G;
}

Graph GraphGenerator::makeBarabasiAlbertGraph(int64_t n, int64_t k) {

	Graph G(n);

	// all nodes need to have at least degree 1 - create a path
	for (node v = 1; v < n; ++v) {
		G.insertEdge(v, (v + 1));
	}

	int64_t m = G.numberOfEdges(); // number of edges
	int64_t r = 0;

	G.forallNodes([&](node u) {
		TRACE("connecting node " << u);
		for (int64_t i = 0; i < k; ++i) { // for all k new edges
			TRACE("2m = " << 2 * m);
			Aux::RandomInteger randInt(0, 2*m);	// TODO: n * k instantiations of RandomInteger are inefficient because random device reads from /dev/random
			r = randInt.generate();
			TRACE("r = " << r);
			for (node v = 1; v <= n; ++v) {
				if (r <= G.degree(v)) {
					// select v
					G.insertEdge(u, v);
					TRACE("inserting edge (" << u << "," << v << ")");
					m += 1;
					r = 0;
					break;
				} else {
					TRACE("skipping node " << v);
				}
				r -= G.degree(v);
			}
		}
	});


	return G;
}

} /* namespace EnsembleClustering */
