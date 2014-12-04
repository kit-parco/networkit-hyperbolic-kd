/*
 * AlgebraicDistance.cpp
 *
 *  Created on: 19.06.2013
 *      Author: cls
 */

#include "../auxiliary/Random.h"
#include "../auxiliary/Log.h"

#include "AlgebraicDistance.h"

namespace NetworKit {

AlgebraicDistance::AlgebraicDistance(const Graph& G, count numberSystems, count numberIterations, double omega, index norm) : NodeDistance(G), numSystems(numberSystems), numIters(numberIterations), omega(omega), norm(norm) {

}

void AlgebraicDistance::preprocess() {
	Aux::Timer running1;
	running1.start();
	// random init
	randomInit();

	// main loop
	for (index iter = 0; iter < numIters; ++iter) {
		// store previous iteration
		std::vector<std::vector<double> > oldLoads = loads;

		for (index sys = 0; sys < numSystems; ++sys) {
			G.forNodes([&](node u) {
				double val = 0.0;

				// step 1
				G.forNeighborsOf(u, [&](node v, edgeweight weight) {
					val += weight * oldLoads[sys][v];
				});
				val /= G.weightedDegree(u);

				// step 2
				loads[sys][u] = (1 - omega) * oldLoads[sys][u] + omega * val;
			});
		}
	}
	running1.stop();
	INFO("elapsed millisecs for AD preprocessing: ", running1.elapsedMilliseconds(), "\n");
}

double AlgebraicDistance::distance(node u, node v) {
	double result = 0.0;

	if (norm == MAX_NORM) { // maximum norm
		for (index sys = 0; sys < numSystems; ++sys) {
			double absDiff = fabs(loads[sys][u] - loads[sys][v]);
			if (absDiff > result) {
				result = absDiff;
			}
		}
	}
	else {
		for (index sys = 0; sys < numSystems; ++sys) {
			double absDiff = fabs(loads[sys][u] - loads[sys][v]);
			result += pow(absDiff, norm);
		}
		result = pow(result, 1.0 / (double) norm);
	}

	return result;
}

std::vector<double> AlgebraicDistance::getEdgeAttribute() {
	INFO("attributizing edges");
	std::vector<double> ad(G.upperEdgeIdBound(), 0.0);

	G.parallelForEdges([&](node u, node v, edgeid uv){
		ad[uv] = distance(u, v);
	});

	return ad;
}

void AlgebraicDistance::randomInit() {
	count n = G.numberOfNodes();

	// allocate space for loads
	loads.resize(numSystems);
	for (index i = 0; i < numSystems; ++i) {
		loads[i].resize(n);
	}

	for (index i = 0; i < numSystems; ++i) {
		G.forNodes([&](node v) {
			loads[i][v] = Aux::Random::real();
		});
	}
}

} /* namespace NetworKit */
