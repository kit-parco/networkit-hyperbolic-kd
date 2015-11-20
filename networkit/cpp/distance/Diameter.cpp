/*
 * Diameter.cpp
 *
 *  Created on: 19.02.2014
 *      Author: Daniel Hoske, Christian Staudt
 */

#include <numeric>

#include "Diameter.h"
#include "Eccentricity.h"
#include "../graph/BFS.h"
#include "../graph/Dijkstra.h"
#include "../components/ConnectedComponents.h"
#include "../structures/Partition.h"
#include "../graph/BFS.h"
#include "../auxiliary/Parallel.h"

namespace NetworKit {


edgeweight Diameter::exactDiameter(const Graph& G) {
	using namespace std;

	Aux::SignalHandler handler;

	edgeweight diameter = 0.0;

	if (! G.isWeighted()) {
		std::tie(diameter, std::ignore) = estimatedDiameterRange(G, 0);
	} else {
		 G.forNodes([&](node v) {
			handler.assureRunning();
		 	Dijkstra dijkstra(G, v);
		 	dijkstra.run();
		 	auto distances = dijkstra.getDistances();
		 	for (auto distance : distances) {
		 		if (diameter < distance) {
		 			diameter = distance;
		 		}
		 	}
//			DEBUG("ecc(", v, "): ", *std::max_element(distances.begin(), distances.end()), " of ", distances);
		 });
	}

	if (diameter == std::numeric_limits<edgeweight>::max()) {
		throw std::runtime_error("Graph not connected - diameter is infinite");
	}
	return diameter;
}





std::pair<edgeweight, edgeweight> Diameter::estimatedDiameterRange(const NetworKit::Graph &G, double error) {
	// TODO: make abortable with ctrl+c using SignalHandling code
	if (G.isDirected()) {
		throw std::runtime_error("Error, the diameter of directed graphs cannot be computed yet.");
	}

	Aux::SignalHandler handler;
	/*
	 * This is an implementation of a slightly modified version of the exactSumSweep algorithm as presented in
	 * Fast diameter and radius BFS-based computation in (weakly connected) real-world graphs: With an application to the six degrees of separation games
	 * by Michele Borassi, Pierluigi Crescenzi, Michel Habib, Walter A. Kosters, Andrea Marino, Frank W. Takes
	 * http://www.sciencedirect.com/science/article/pii/S0304397515001644
	 */

	std::vector<count> sum(G.upperNodeIdBound(), 0);
	std::vector<count> eccLowerBound(G.upperNodeIdBound(), 0), eccUpperBound(G.upperNodeIdBound(), std::numeric_limits<count>::max());

	#pragma omp parallel for
	for (node u = 0; u < G.upperNodeIdBound(); ++u) {
		if (!G.hasNode(u)) {
			eccUpperBound[u] = 0;
		}
	}

	std::vector<count> distances(G.upperNodeIdBound(), 0);
	std::vector<bool> finished(G.upperNodeIdBound(), false);
	count k = 4;


	ConnectedComponents comp(G);
	comp.run();
	count numberOfComponents = comp.numberOfComponents();
	std::vector<node> startNodes(numberOfComponents, 0), maxDist(numberOfComponents, 0);
	std::vector<node> firstDeg2Node(numberOfComponents, none);
	std::vector<node> distFirst(numberOfComponents, 0);
	std::vector<node> ecc(numberOfComponents, 0);


	auto updateBounds = [&]() {
		G.parallelForNodes([&](node u) {
			if (finished[u]) return;

			auto c = comp.componentOfNode(u);

			if (distances[u] <= distFirst[c]) {
				eccUpperBound[u] = std::max(distances[u], ecc[c] - distances[u]);
				eccLowerBound[u] = eccUpperBound[u];
				finished[u] = true;
			} else {
				eccUpperBound[u] = std::min(distances[u] + ecc[c] - 2 * distFirst[c], eccUpperBound[u]);
				eccLowerBound[u] = std::max(eccLowerBound[u], distances[u]);
				finished[u] = (eccUpperBound[u] == eccLowerBound[u]);
			}
		});

		ecc.clear();
		ecc.resize(numberOfComponents, 0);
		distFirst.clear();
		distFirst.resize(numberOfComponents, 0);
	};

	auto diameterBounds = [&]() {
		auto maxExact = *Aux::Parallel::max_element(eccLowerBound.begin(), eccLowerBound.end());
		auto maxPotential = *Aux::Parallel::max_element(eccUpperBound.begin(), eccUpperBound.end());
		return std::make_pair(maxExact, maxPotential);
	};

	auto visitNode = [&](node v, count dist) {
		sum[v] += dist;
		distances[v] = dist;

		index c = comp.componentOfNode(v);
		ecc[c] = std::max(dist, ecc[c]);
		if (firstDeg2Node[c] == none && G.degree(v) > 1) {
			firstDeg2Node[c] = v;
			distFirst[c] = dist;
		}
	};

	for (index i = 0; i < k; ++i) {
		handler.assureRunning();
		if (i == 0) {
			std::vector<count> minDeg(numberOfComponents, G.numberOfNodes());

			// for each component, find the node with the minimum degreee and add it as start node
			G.forNodes([&](node v) {
				count d = G.degree(v);
				count c = comp.componentOfNode(v);
				if (d < minDeg[c]) {
					startNodes[c] = v;
					minDeg[c] = d;
				}
			});
		}

		G.BFSfrom(startNodes, [&](node v, count dist) {
			visitNode(v, dist);
			index c = comp.componentOfNode(v);
			if (sum[v] >= maxDist[c]) {
				maxDist[c] = sum[v];
				startNodes[c] = v;
			}
		});

		updateBounds();
	}


	handler.assureRunning();

	std::vector<count> minDist(numberOfComponents, G.numberOfNodes());
	G.forNodes([&](node u) {
		auto c = comp.componentOfNode(u);
		if (sum[u] < minDist[c]) {
			minDist[c] = sum[u];
			startNodes[c] = u;
		}
	});

	handler.assureRunning();

	G.BFSfrom(startNodes, visitNode);
	updateBounds();

	count lb, ub;
	std::tie(lb, ub) = diameterBounds();

	for (index i = 0; i < G.numberOfNodes() && ub > (lb + error*lb); ++i) {
		handler.assureRunning();
		startNodes.clear();
		startNodes.resize(numberOfComponents, none);

		if ((i % 2) == 0) {
			G.forNodes([&](node u) {
				auto c = comp.componentOfNode(u);
				if (startNodes[c] == none || std::tie(eccUpperBound[u], sum[u]) > std::tie(eccUpperBound[startNodes[c]], sum[startNodes[c]])) {
					startNodes[c] = u;
				}
			});
		} else {
			G.forNodes([&](node u) {
				auto c = comp.componentOfNode(u);
				if (startNodes[c] == none || std::tie(eccLowerBound[u], sum[u]) < std::tie(eccLowerBound[startNodes[c]], sum[startNodes[c]])) {
					startNodes[c] = u;
				}
			});
		}

		handler.assureRunning();

		G.BFSfrom(startNodes, visitNode);

		handler.assureRunning();

		updateBounds();

		std::tie(lb, ub) = diameterBounds();
	}

	return {lb, ub};
}


edgeweight Diameter::estimatedVertexDiameter(const Graph& G, count samples) {

	edgeweight infDist = std::numeric_limits<edgeweight>::max();

	// TODO: consider weights

	auto estimateFrom = [&](node v) -> count {
		BFS bfs(G, v);
		bfs.run();
		auto distances = bfs.getDistances();

		// get two largest path lengths
		edgeweight maxD = 0;
		edgeweight maxD2 = 0; // second largest distance
		for (auto d : distances) {
			if ((d != infDist) && (d >= maxD)) {
				maxD2 = maxD;
				maxD = d;
			}
		}

		edgeweight dMax = maxD + maxD2;
		count vd = (count) dMax + 1; 	// count the nodes, not the edges
		return vd;
	};

	edgeweight vdMax = 0;
	#pragma omp parallel for
	for (count i = 0; i < samples; ++i) {
		node u = G.randomNode();
		edgeweight vd = estimateFrom(u);
		DEBUG("sampled vertex diameter from node ", u, ": ", vd);
		#pragma omp critical
		{
			if (vd > vdMax) {
				vdMax = vd;
			}
		}
	}

	return vdMax;

}

edgeweight Diameter::estimatedVertexDiameterPedantic(const Graph& G) {
	count vd = 0;
	if (!G.isWeighted()) {
		std::vector<bool> visited(G.upperNodeIdBound(), false);
		// perform breadth-first searches
		G.forNodes([&](node u) {
			if (visited[u] == false) {
				count maxDist = 0, maxDist2 = 0;
				G.BFSfrom(u, [&](node v, count dist) {
					visited[v] = true;
					if (dist > maxDist) {
						maxDist2 = maxDist;
						maxDist = dist;
					}
					else if (dist > maxDist2) {
						maxDist2 = dist;
					}
				});
				if (maxDist + maxDist2 > vd) {
					vd = maxDist + maxDist2;
				}
				assert (visited[u] == true);
			}
		});
		vd ++; //we need the number of nodes, not the number of edges
	}
	else {
		ConnectedComponents cc(G);
		DEBUG("finding connected components");
		cc.run();
		INFO("Number of components ", cc.numberOfComponents());
		DEBUG("Estimating size of the largest component");
		std::map<count, count> sizes = cc.getComponentSizes();
		count largest_comp_size = 0;
		for(auto it = sizes.cbegin(); it != sizes.cend(); ++it) {
			DEBUG(it->second);
			if (it->second > largest_comp_size) {
				largest_comp_size = it->second;
			}
		}
		INFO("Largest component size: ", largest_comp_size);
		vd = largest_comp_size;
	}
	return vd;
}

} /* namespace NetworKit */
