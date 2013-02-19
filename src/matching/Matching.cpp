/*
 * Matching.cpp
 *
 *  Created on: 03.12.2012
 *      Author: cls
 */

#include "Matching.h"

namespace EnsembleClustering {

Matching::Matching(int64_t n) :
		NodeMap<node>(n, 0) {
	// initialize each node's matching partner to itself
	for (index i = 0; i < n; ++i) {
		this->data[i] = i;
	}
}

Matching::~Matching() {
}

bool Matching::isMatched(const node& u) const {
	return (this->data[u] != u);
}

bool Matching::isProper(Graph& G) const {
	/**
	 * The content of this data structure represents a matching iff
	 * 	(for all v in V: M[v] = v or M[M[v]] = v) and
	 * 	(for all (u,v) in M): (u,v) in E
	 *
	 */
	bool sym = true;
	// check if entries are symmetric
	for (node v = 0; v < G.numberOfNodes(); ++v) {
		sym = (((*this)[v] == v) || ((*this)[(*this)[v]] == v));
		if (!sym) {
			DEBUG("node " << v << " is not symmetrically matched");
			return false;
		}
	}

	bool inGraph = true;
	// check if every pair exists as an edge
	for (node v = 0; v < G.numberOfNodes(); ++v) {
		node w = (*this)[v];
		if (v != w) {
			inGraph = G.hasEdge(v, w);
			if (!inGraph) {
				DEBUG("matched pair (" << v << "," << w << ") is not an edge");
				return false;
			}
		}
	}

	return (sym && inGraph);
}

void Matching::match(const node& u, const node& v) {
	(*this)[u] = v;
	(*this)[v] = u;
}

void Matching::unmatch(const node& u, const node& v) {
	(*this)[u] = u;
	(*this)[v] = v;
}

bool Matching::areMatched(const node& u, const node& v) const {
	return ((*this)[u] == v);
}

count Matching::matchingSize() const {
	count size = 0;
	for (index i = 0; i < n; ++i) { // TODO: parallel
		if (isMatched(i)) {
			++size;
		}
	}
	return size / 2;
}

index Matching::mate(node v) const {
	if (isMatched(v)) {
		return data[v];
	}
	else return none;
}

}
 /* namespace EnsembleClustering */
