/*
 * ForestFireScore.cpp
 *
 *  Created on: 26.08.2014
 *      Author: Gerd Lindner
 */

#include "ForestFireScore.h"
#include <limits>
#include <set>
#include <queue>
#include "../auxiliary/Log.h"

namespace NetworKit {

ForestFireScore::ForestFireScore(const Graph& G, double pf, double targetBurntRatio): EdgeScore<double>(G), pf(pf), targetBurntRatio(targetBurntRatio) {}

void ForestFireScore::run() {
	if (!G.hasEdgeIds()) {
		throw std::runtime_error("edges have not been indexed - call indexEdges first");
	}

	std::vector<count> burnt (G.upperEdgeIdBound(), 0);
	count edgesBurnt = 0;

	while (edgesBurnt < targetBurntRatio * G.numberOfEdges()) {
		//Start a new fire
		std::queue<node> activeNodes;
		std::vector<bool> visited (G.upperNodeIdBound(), false);
		activeNodes.push(G.randomNode());

		auto forwardNeighbors = [&](node u) {
			std::vector<std::pair<node, edgeid>> validEdges;
			G.forNeighborsOf(u, [&](node, node x, edgeid eid){
				if (! visited[x]) {
					validEdges.emplace_back(x, eid);
				}
			});
			return validEdges;
		};

		while (! activeNodes.empty()) {
			node v = activeNodes.front();
			activeNodes.pop();

			std::vector<std::pair<node, edgeid>> validNeighbors = forwardNeighbors(v);
			while (true) {
				double q = Aux::Random::real(1.0);
				if (q > pf || validNeighbors.empty()) {
					break;
				}
				count index = Aux::Random::integer(validNeighbors.size() - 1);

				{ // mark node as visited, burn edge
					node x;
					edgeid eid;
					std::tie(x, eid) = validNeighbors[index];
					activeNodes.push(x);
					burnt[eid]++;
					edgesBurnt++;
					visited[x] = true;
				}

				validNeighbors[index] = validNeighbors.back();
				validNeighbors.pop_back();
			}
		}
	}

	std::vector<double> burntNormalized (G.upperEdgeIdBound(), 0.0);
	double maxv = (double) *std::max_element(std::begin(burnt), std::end(burnt));

	if (maxv > 0) {
		count idx = 0;
		for (auto& b : burnt) {
			burntNormalized[idx++] = b / maxv;
		}
	}

	scoreData = std::move(burntNormalized);
	hasRun = true;
}

double ForestFireScore::score(node u, node v) {
	throw std::runtime_error("Not implemented: Use scores() instead.");
}

double ForestFireScore::score(edgeid eid) {
	throw std::runtime_error("Not implemented: Use scores() instead.");
}

} /* namespace NetworKit */
