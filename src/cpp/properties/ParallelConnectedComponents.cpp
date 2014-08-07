/*
 * ConnectedComponents.cpp
 *
 *  Created on: Dec 16, 2013
 *      Author: cls
 */

#include <set>

#include "ParallelConnectedComponents.h"
#include "../structures/Partition.h"
#include "../coarsening/PartitionCoarsening.h"
#include "../auxiliary/Log.h"

namespace NetworKit {

ParallelConnectedComponents::ParallelConnectedComponents(const Graph& G, bool coarsening) : G(G), coarsening(coarsening) {

}


void ParallelConnectedComponents::run() {
	if (G.isDirected()) {
		throw std::runtime_error("algorithm does not accept directed graphs");
	}

	// calculate connected components by label propagation
	count z = G.upperNodeIdBound();

	DEBUG("initializing labels");
	component = Partition(z);
	component.allToSingletons();

	DEBUG("initializing active nodes");
	const char INACTIVE = 0;
	const char ACTIVE = 1;
	std::vector<char> activeNodes(z); // record if node must be processed
	std::vector<char> nextActiveNodes(z, ACTIVE); // for next iteration

	DEBUG("main loop");
//	count numActive = 0; // for debugging purposes only
	count numIterations = 0;
	bool change = true;
	// only 8 iterations when coarsening is on, otherwise till no more changes happened
	while (change && (!coarsening || numIterations < 8)) {
//		TRACE("label propagation iteration");
		activeNodes.swap(nextActiveNodes);
		nextActiveNodes.assign(z, INACTIVE);
		change = false;
//		numActive = 0;
		G.balancedParallelForNodes([&](node u) {
			if (activeNodes[u] == ACTIVE) {
//				++numActive;
				// get smallest
				index smallest = component[u];
				G.forNeighborsOf(u, [&](node v) {
					smallest = std::min(smallest, component[v]);
				});

				if (component[u] != smallest) {
					component.moveToSubset(smallest, u);
					change = true;
					G.forNeighborsOf(u, [&](node v) {
						// only nodes that do not have the smallest component label
						// will see a new component label because of the change of u,
						// only they need to be activated
						if (component[v] != smallest) {
							nextActiveNodes[v] = ACTIVE;
						}
					});
				}
			}
		});
//		TRACE("num active: ", numActive);
		++numIterations;
	}
	if (coarsening && numIterations == 8) { // TODO: externalize constant
		// coarsen and make recursive call
		ParallelPartitionCoarsening con;
		std::pair<Graph, std::vector<node> > coarse = con.run(G, component);
		ParallelConnectedComponents cc(coarse.first);
		cc.run();

		// apply to current graph
		G.parallelForNodes([&](node u) {
			component[u] = cc.componentOfNode(coarse.second[u]);
		});
	}
}

void ParallelConnectedComponents::runSequential() {
	if (G.isDirected()) {
		throw std::runtime_error("algorithm does not accept directed graphs");
	}
	// calculate connected components by label propagation
	count z = G.upperNodeIdBound();
	DEBUG("initializing labels");
	component = Partition(z);
	component.allToSingletons();
	DEBUG("initializing active nodes");
	std::vector<bool> activeNodes(z, true); // record if node must be processed

	DEBUG("main loop");
	// count numActive = 0; // for debugging purposes only
	count numIterations = 0;
	bool change = true;
	// only 8 iterations when coarsening is on, otherwise till no more changes happened
	while (change && (!coarsening || numIterations < 8)) {
		// TRACE("label propagation iteration");
		change = false;
		// numActive = 0;
		G.forNodes([&](node u) {
			if (activeNodes[u]) {
				// ++numActive;
				// get smallest
				index smallest = component[u];
				G.forNeighborsOf(u, [&](node v) {
					smallest = std::min(smallest, component[v]);
				});

				if (component[u] != smallest) {
					component.moveToSubset(smallest, u);
					change = true;
					G.forNeighborsOf(u, [&](node v) {
						// only nodes that do not have the smallest component label
						// will see a new component label because of the change of u,
						// only they need to be activated
						if (component[v] != smallest) {
							activeNodes[v] = true;
						}
					});
				} else {
					activeNodes[u] = false; // current node becomes inactive
				}
			}
		});
		// TRACE("num active: ", numActive);
		++numIterations;
	}
	if (coarsening && numIterations == 8) { // TODO: externalize constant
		// coarsen and make recursive call
		PartitionCoarsening con;
		std::pair<Graph, std::vector<node> > coarse = con.run(G, component);
		ParallelConnectedComponents cc(coarse.first);
		cc.run();
		// apply to current graph
		G.forNodes([&](node u) {
			component[u] = cc.componentOfNode(coarse.second[u]);
		});
	}
}


Partition ParallelConnectedComponents::getPartition() {
	return this->component;
}

count ParallelConnectedComponents::numberOfComponents() {
	return this->component.numberOfSubsets();
}

count ParallelConnectedComponents::componentOfNode(node u) {
	assert (component[u] != none);
	return component[u];
}

}
