/*
 * EigenvectorCentrality.cpp
 *
 *  Created on: 19.03.2014
 *      Author: Henning
 */

#include "EigenvectorCentrality.h"
#include "../auxiliary/NumericTools.h"

namespace NetworKit {

EigenvectorCentrality::EigenvectorCentrality(const Graph& G, double tol):
		Centrality(G, true), tol(tol)
{

}

void EigenvectorCentrality::runImpl() {
	count z = G.upperNodeIdBound();
	std::vector<double> values(z, 1.0);
	scoreData = values;

	// do not execute algorithm on directed graphs since this is error prone
	// and can yield misleading results (wrong metric, not implementation fault!)
	if (G.isDirected()) {
		return;
	}

	double length = 0.0;
	double oldLength = 0.0;

	double NEAR_ZERO = 1e-16;

	auto converged([&](double val, double other) {
		// compute residual
		return (Aux::NumericTools::equal(val, other, tol));
	});

	do {
		oldLength = length;

		// iterate matrix-vector product
		G.parallelForNodes([&](node u) {
			values[u] = 0.0;
			G.forNeighborsOf(u, [&](node v) {
				values[u] += G.weight(u, v) * scoreData[v];
			});
		});

//		// set everything very small to zero
//		G.parallelForNodes([&](node u) {
//			if (values[u] < NEAR_ZERO) {
//				values[u] = 0.0;
//			}
//		});

		// normalize values
		length = 0.0;
		length = G.parallelSumForNodes([&](node u) {
			return (values[u] * values[u]);
		});
		length = sqrt(length);

//		TRACE("length: ", length);
//		TRACE(values);

		assert(! Aux::NumericTools::equal(length, NEAR_ZERO));
		G.parallelForNodes([&](node u) {
			values[u] /= length;
		});

//		TRACE(values);

		scoreData = values;
		assureRunning();
	} while (! converged(length, oldLength));

	// check sign and correct if necessary
	if (scoreData[0] < 0) {
		G.parallelForNodes([&](node u) {
			scoreData[u] = fabs(scoreData[u]);
		});
	}

	hasRun = true;
}

} /* namespace NetworKit */
