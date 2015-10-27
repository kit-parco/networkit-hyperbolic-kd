/*
* ApproxBetweenness2.cpp
*
*  Created on: 13.06.2014
*      Author: Christian Staudt, Elisabetta Bergamini
*/


#include "ApproxBetweenness2.h"
#include "../graph/BFS.h"
#include "../graph/Dijkstra.h"
#include "../graph/SSSP.h"
#include "../auxiliary/SignalHandling.h"
#include "../auxiliary/Parallelism.h"


#include <memory>
#include <omp.h>

namespace NetworKit {

ApproxBetweenness2::ApproxBetweenness2(const Graph& G, count nSamples, bool normalized, bool parallel) : Centrality(G, normalized), nSamples(nSamples), parallel(parallel) {
}

void ApproxBetweenness2::run() {
	hasRun = false;

	Aux::SignalHandler handler;

	//std::vector<node> sampledNodes = G.nodes();
	std::vector<node> sampledNodes;

	// sample nodes
	for (count i = 0; i <= nSamples; ++i) {
	 	sampledNodes.push_back(G.randomNode());
	}


	// thread-local scores for efficient parallelism
	count maxThreads = omp_get_max_threads();
	if (!parallel) maxThreads = 1;
	std::vector<std::vector<double> > scorePerThread(maxThreads, std::vector<double>(G.upperNodeIdBound()));


	auto computeDependencies = [&](node s){
		// run single-source shortest path algorithm
		std::unique_ptr<SSSP> sssp;
		if (G.isWeighted()) {
			sssp.reset(new Dijkstra(G, s, true, true));
		} else {
			sssp.reset(new BFS(G, s, true, true));
		}
		if (!handler.isRunning()) return;
		sssp->run();
		if (!handler.isRunning()) return;


		// create stack of nodes in non-decreasing order of distance
		std::vector<node> stack = sssp->getStack(true);

		// compute dependencies and add the contributions to the centrality score
		std::vector<double> dependency(G.upperNodeIdBound(), 0.0);
		for (auto it = stack.rbegin(); it != stack.rend(); ++it) {
			node t = *it;
			if (t == s){
				continue;
			}
			for (node p : sssp->getPredecessors(t)) {
				// TODO: make weighting factor configurable

				// workaround for integer overflow in large graphs
				bigfloat tmp = sssp->numberOfPaths(p) / sssp->numberOfPaths(t);
				double weight;
				tmp.ToDouble(weight);

				dependency[p] += (double(sssp->distance(p)) / sssp->distance(t)) * weight * (1 + dependency[t]);
			}
			scorePerThread[omp_get_thread_num()][t] += dependency[t];
		}
	};


	if (parallel) {
		#pragma omp parallel for
		for (index i = 0; i < sampledNodes.size(); ++i) {
			computeDependencies(sampledNodes[i]);
		}

		scoreData = std::vector<double>(G.upperNodeIdBound(), 0.0);

		// add up all thread-local values
		for (const auto &local : scorePerThread) {
			G.parallelForNodes([&](node v){
				scoreData[v] += local[v];
			});
		}
	} else {
		for (auto u : sampledNodes) {
			computeDependencies(u);
		}

		scoreData.swap(scorePerThread[0]);
	}


	const count n = G.numberOfNodes();
	const count pairs = (n-2) * (n-1);

	// extrapolate
	G.parallelForNodes([&](node u) {
		scoreData[u] = scoreData[u] * (2 * n / double(nSamples));

		if (normalized) {
			// divide by the number of possible pairs
			scoreData[u] = scoreData[u] / pairs;
		}
	});

	handler.assureRunning();
	hasRun = true;
}


} /* namespace NetworKit */
