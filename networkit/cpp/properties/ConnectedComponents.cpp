/*
 * ConnectedComponents.cpp
 *
 *  Created on: Dec 16, 2013
 *      Author: cls
 */

#include <set>

#include "ConnectedComponents.h"
#include "../structures/Partition.h"
#include "../auxiliary/Log.h"

namespace NetworKit {

ConnectedComponents::ConnectedComponents(const Graph& G) : G(G), hasRun(false) {
	if (G.isDirected()) {
		throw std::runtime_error("Error, connected components of directed graphs cannot be computed, use StronglyConnectedComponents for them.");
	}
}

void ConnectedComponents::runImpl() {
	DEBUG("initializing labels");
	component = Partition(G.upperNodeIdBound(), none);
	numComponents = 0;

	// perform breadth-first searches
	G.forNodes([&](node u) {
		assureRunning();
		if (component[u] == none) {
			component.setUpperBound(numComponents+1);
			index c = numComponents;
			G.BFSfrom(u, [&](node v) {
				component[v] = c;
			});
			assert (component[u] != none);
			++numComponents;
		}
	});

	hasRun = true;
}


Partition ConnectedComponents::getPartition() {
	if (!hasRun) throw std::runtime_error("run method has not been called");
	return this->component;
}


std::map<index, count> ConnectedComponents::getComponentSizes() {
	if (!hasRun) throw std::runtime_error("run method has not been called");
	return this->component.subsetSizeMap();
}

}
