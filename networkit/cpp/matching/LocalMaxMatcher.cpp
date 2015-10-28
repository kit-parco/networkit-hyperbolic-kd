/*
 * ParallelMatcher.cpp
 *
 *  Created on: 05.12.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include "LocalMaxMatcher.h"

namespace NetworKit {

LocalMaxMatcher::LocalMaxMatcher(const Graph& graph): Matcher(graph)
{
	//TODO: check for self loops here, issue a warning
}

// TODO: update to new edge attribute system
// TODO: make local max matching parallel


Matching LocalMaxMatcher::run() {
	int64_t z = G.upperNodeIdBound();
	count E = G.numberOfEdges();
	Matching M(z);

	// put edges into array of triples
	struct MyEdge {
		node s; // source
		node t; // target
		edgeweight w; // weight
	};

	std::vector<MyEdge> edges(E);
	index e = 0;
	G.forEdges([&](node u, node v, edgeweight w) {
		edges[e].s = u;
		edges[e].t = v;
		edges[e].w = w + Aux::Random::real(1e-6);
		++e;
	});

	// candidates records mating candidates
	std::vector<MyEdge> candidates(z);
	G.parallelForNodes([&](node u) {
		candidates[u].w = (edgeweight) 0;
		candidates[u].t = u; // itself as mating partner => unmatched
	});

	while (E > 0) {
		// for each edge find out if it is locally maximum
		for (auto edge: edges) {
			if (edge.w > candidates[edge.s].w && edge.w > candidates[edge.t].w) {
				candidates[edge.s].t = edge.t;
				candidates[edge.s].w = edge.w;
				candidates[edge.t].t = edge.s;
				candidates[edge.t].w = edge.w;
			}
		}

		// check if candidates agree to match; if so, then match them
		for (auto edge: edges) {
			node u = edge.s;
			node v = edge.t;
			if (candidates[u].t == v && candidates[v].t == u) {
				// both nodes agree
				M.match(u, v);
			}
		}

		// create remaining "graph" by selecting remaining edges (as triples)
		// adjust candidates
		std::vector<MyEdge> newEdges;
		for (auto edge: edges) {
			if (! M.isMatched(edge.s) && ! M.isMatched(edge.t) && edge.s != edge.t) {
				newEdges.push_back(edge);
				candidates[edge.s].w = (edgeweight) 0;
				candidates[edge.t].w = (edgeweight) 0;
			}
		}
		edges = newEdges;
		E = edges.size();
	}

	return M;
}

} /* namespace NetworKit */
