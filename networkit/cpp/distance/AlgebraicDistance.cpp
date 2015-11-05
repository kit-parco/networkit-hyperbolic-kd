/*
 * AlgebraicDistance.cpp
 *
 *  Created on: 03.11.2015
 *      Author: Henning Meyerhenke, Christian Staudt
 */


#include "AlgebraicDistance.h"

#include "../auxiliary/Timer.h"


namespace NetworKit {

AlgebraicDistance::AlgebraicDistance(const Graph& G, count numberSystems, count numberIterations, double omega, index norm) : NodeDistance(G), numSystems(numberSystems), numIters(numberIterations), omega(omega), norm(norm) {
	if ((omega < 0.0) || (omega > 1.0)) throw std::invalid_argument("omega must be in [0,1]");
}

void AlgebraicDistance::randomInit() {
	// allocate space for loads
	loads.resize(numSystems*G.upperNodeIdBound());

	#pragma omp parallel for
	for (index i = 0; i < loads.size(); ++i) {
		loads[i] = Aux::Random::real();
	}
}

void AlgebraicDistance::preprocess() {
	Aux::Timer running1;
	running1.start();
	// random init
	randomInit();

	// main loop

	{
		std::vector<double> oldLoads(loads.size());

		for (index iter = 0; iter < numIters; ++iter) {
			// store previous iteration
			loads.swap(oldLoads);

			G.balancedParallelForNodes([&](node u) {
				std::vector<double> val(numSystems, 0.0);

				double weightedDeg = 0;
				// step 1
				G.forNeighborsOf(u, [&](node v, edgeweight weight) {
					for (index i = 0; i < numSystems; ++i) {
						val[i] += weight * oldLoads[v*numSystems + i];
					}

					weightedDeg += weight;
				});

				for (index i = 0; i < numSystems; ++i) {
					val[i] /= weightedDeg;

					// step 2
					loads[u*numSystems + i] = (1 - omega) * oldLoads[u*numSystems + i] + omega * val[i];
				}
			});
		}
	}

	// calculate edge scores
	if (!G.hasEdgeIds()) {
		throw std::runtime_error("edges have not been indexed - call indexEdges first");
	}

	edgeScores.resize(G.upperEdgeIdBound(), none);

	G.parallelForEdges([&](node u, node v, edgeid eid) {
		edgeScores[eid] = distance(u, v);
	});

	running1.stop();
	INFO("elapsed millisecs for AD preprocessing: ", running1.elapsedMilliseconds(), "\n");
}

double AlgebraicDistance::distance(node u, node v) {
	if (loads.size() == 0) {
		throw std::runtime_error("Call preprocess() first.");
	}
	double result = 0.0;

	if (norm == MAX_NORM) {
		for (index sys = 0; sys < numSystems; ++sys) {
			double absDiff = fabs(loads[u*numSystems + sys] - loads[u*numSystems + sys]);
			if (absDiff > result) {
				result = absDiff;
			}
		}
	} else {
		for (index sys = 0; sys < numSystems; ++sys) {
			double absDiff = fabs(loads[u*numSystems + sys] - loads[v*numSystems + sys]);
			result += pow(absDiff, norm);
		}
		result = pow(result, 1.0 / (double) norm);
	}

	return std::isnan(result) ? 0 : result;
}


std::vector<double> AlgebraicDistance::getEdgeAttribute() {
	return edgeScores;
}



} /* namespace NetworKit */