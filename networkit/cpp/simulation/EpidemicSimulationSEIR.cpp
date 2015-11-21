/*
 * EpidemicSimulationSEIR.h
 *
 *  Created on: 20.11.2015
 *      Author: Christian Staudt
 */

#include "EpidemicSimulationSEIR.h"
#include "../auxiliary/Random.h"

namespace NetworKit {

EpidemicSimulationSEIR::EpidemicSimulationSEIR(const Graph& G, count tMax, double transP, count eTime, count iTime) : Algorithm(), G(G), tMax(tMax), transP(transP), eTime(eTime), iTime(iTime)  {
}

void EpidemicSimulationSEIR::run() {

	typedef EpidemicSimulationSEIR::State State;

	index t = 0;

	//initialize state and timestamp arrays
	state.resize(G.upperNodeIdBound(), State::U);
	timestamp.resize(G.upperNodeIdBound(), none);

	auto setState = [&](node v, State X){
		state[v] = X;
		timestamp[v] = t;
	};

	// initialize nodes to Susceptible
	G.parallelForNodes([&](node v) {
		setState(v, State::S);
	});

	// contact may expose susceptible node to infection
	auto contact = [&](node v) {
		if ((state[v] == State::S) && (Aux::Random::probability() <= transP)) {
			setState(v, State::E);
		}
	};


	auto sweep = [&](node u) {
		if (state[u] == State::S) {
			// do nothing
		} else if (state[u] == State::E) {
			// exposed nodes become infectious after time
			if ((t - timestamp[u]) >= eTime) {
				setState(u, State::I);
			}
		} else if (state[u] == State::I) {
			// contact neighbors of infectious node
			G.forNeighborsOf(u, [&](node v){
				contact(v);
			});
			// infectious nodes become removed after time
			if ((t - timestamp[u]) >= iTime) {
				setState(u, State::R);
			}
		} else if (state[u] == State::R) {
			// do nothing
		} else if (state[u] == State::U) {
			throw std::runtime_error("node in undefined state encountered - should not happen");
		} else {
			throw std::runtime_error("else branch taken - should not happen");
		}
	};

	while (t < tMax) {
		G.forNodes(sweep);
		t += 1;
	}

	hasRun = true;
}


} /* namespace NetworKit */
