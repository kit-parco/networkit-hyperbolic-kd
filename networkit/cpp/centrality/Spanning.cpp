/*
 * Spanning.cpp
 *
 *  Created on: 29.07.2015
 *      Author: henningm
 */


#include "Spanning.h"
#include "../auxiliary/Log.h"
#include "../auxiliary/Timer.h"
#include "../spanning/RandomSpanningTree.h"
#include "../spanning/PseudoRandomSpanningTree.h"

#include <fstream>
#include <sstream>

namespace NetworKit {

Spanning::Spanning(const Graph& G, double tol): Centrality(G), tol(tol), lamg(1e-5) {
	// prepare LAMG
	CSRMatrix matrix = CSRMatrix::graphLaplacian(G);
	Aux::Timer t;
	t.start();
	lamg.setupConnected(matrix);
	t.stop();

	setupTime = t.elapsedMilliseconds();

	DEBUG("done setting up Spanning");
}

void Spanning::run() {
	count n = G.numberOfNodes();
	scoreData.clear();
	scoreData.resize(G.numberOfEdges(), 0.0);

	// set up solution vector and status
	Vector solution(n);

	Vector rhs(n, 0.0);
	Vector zeroVector(n, 0.0);

	// solve for each edge
	G.forEdges([&](node u, node v, edgeid e) {
		// set up right-hand side
		rhs[u] = +1.0;
		rhs[v] = -1.0;
		TRACE("before solve for ", u, " and ", v);

		solution = zeroVector;

		lamg.solve(rhs, solution);
		double diff = solution[u] - solution[v];
		scoreData[e] = fabs(diff); // TODO: check unweighted, fix weighted case, fix edge IDs!
		rhs[u] = 0.0;
		rhs[v] = 0.0;
	});

	hasRun = true;
}

uint64_t Spanning::getSetupTime() const {
	return setupTime;
}

void Spanning::runApproximation() {
	const count n = G.numberOfNodes();
	const count m = G.numberOfEdges();
	double epsilon2 = tol * tol;
	const count k = ceil(log(n)) / epsilon2;
	double randTab[3] = {1/sqrt(k), -1/sqrt(k)};
	Vector solution(n);
	scoreData.clear();
	scoreData.resize(m, 0.0);

	for (index i = 0; i < k; ++i) {
		Vector rhs(n, 0.0);

		// rhs(v) = \sum_e=1 ^m q(e) * B(e, v)
		//        = +/- q(e)
		G.forEdges([&](node u, node v) {
			double r = randTab[Aux::Random::integer(1)];

			if (u < v) {
				rhs[u] += r;
				rhs[v] -= r;
			}
			else {
				rhs[u] -= r;
				rhs[v] += r;
			}
		});

		lamg.solve(rhs, solution);

		G.forEdges([&](node u, node v, edgeid e) {
			double diff = solution[u] - solution[v];
			scoreData[e] += diff * diff; // TODO: fix weighted case!
		});
	}

	hasRun = true;
}

void Spanning::runParallelApproximation() {
	const count n = G.numberOfNodes();
	const count m = G.numberOfEdges();
	double epsilon2 = tol * tol;
	const count k = ceil(log(n)) / epsilon2;
	double randTab[3] = {1/sqrt(k), -1/sqrt(k)};
	std::vector<Vector> solutions(k, Vector(n));
	std::vector<Vector> rhs(k, Vector(n));
	scoreData.clear();
	scoreData.resize(m, 0.0);

#pragma omp parallel for
	for (index i = 0; i < k; ++i) {
		// rhs(v) = \sum_e=1 ^m q(e) * B(e, v)
		//        = +/- q(e)
		G.forEdges([&](node u, node v) {
			double r = randTab[Aux::Random::integer(1)];

			if (u < v) {
				rhs[i][u] += r;
				rhs[i][v] -= r;
			}
			else {
				rhs[i][u] -= r;
				rhs[i][v] += r;
			}
		});
	}

	lamg.parallelSolve(rhs, solutions);

#pragma omp parallel for
	for (index i = 0; i < k; ++i) {
		G.forEdges([&](node u, node v, edgeid e) {
			double diff = solutions[i][u] - solutions[i][v];
			scoreData[e] += diff * diff; // TODO: fix weighted case!
		});
	}

	hasRun = true;
}

void Spanning::runTreeApproximation() {
	const count reps = 500;
	scoreData.clear();
	scoreData.resize(G.numberOfEdges(), 0.0);

	RandomSpanningTree rst(G);

	for (index i = 0; i < reps; ++i) {
		rst.run();
		Graph tree = rst.getTree();
		G.forEdges([&](node u, node v, edgeid e) {
			scoreData[e] += tree.hasEdge(u, v);
		});
	}

	G.forEdges([&](node u, node v, edgeid e) {
		scoreData[e] /= reps;
	});

	hasRun = true;
}

uint64_t Spanning::runApproximationAndWriteVectors(const std::string &graphPath) {
	Aux::Timer t;
	const count n = G.numberOfNodes();
	const count m = G.numberOfEdges();
	const double epsilon2 = tol * tol;
	const count k = ceil(log(n)) / epsilon2;
	double randTab[3] = {1/sqrt(k), -1/sqrt(k)};
	Vector solution(n);
	scoreData.clear();
	scoreData.resize(m, 0.0);

	t.start();
	for (index i = 0; i < k; ++i) {
		Vector rhs(n, 0.0);

		// rhs(v) = \sum_e=1 ^m q(e) * B(e, v)
		//        = +/- q(e)
		G.forEdges([&](node u, node v) {
			double r = randTab[Aux::Random::integer(1)];
			if (u < v) {
				rhs[u] += r;
				rhs[v] -= r;
			}
			else {
				rhs[u] -= r;
				rhs[v] += r;
			}
		});

		lamg.solve(rhs, solution);

		G.forEdges([&](node u, node v, edgeid e) {
			double diff = solution[u] - solution[v];
			scoreData[e] += diff * diff; // TODO: fix weighted case!
		});
	}
	t.stop();
	hasRun = true;

	return t.elapsedMilliseconds();
}

void Spanning::runPseudoTreeApproximation() {
	const count reps = 500;
	scoreData.clear();
	scoreData.resize(G.numberOfEdges(), 0.0);

	PseudoRandomSpanningTree rst(G);

	for (index i = 0; i < reps; ++i) {
		rst.run();
		Graph tree = rst.getTree();
		G.forEdges([&](node u, node v, edgeid e) {
			scoreData[e] += tree.hasEdge(u, v);
		});
	}

	G.forEdges([&](node u, node v, edgeid e) {
		scoreData[e] /= reps;
	});

	hasRun = true;
}

double Spanning::runForEdge(node u, node v) {
	count n = G.numberOfNodes();

	// set up solution vector and status
	Vector solution(n, 0.0);
	Vector rhs(n, 0.0);

	// set up right-hand side
	rhs[u] = +1.0;
	rhs[v] = -1.0;
	TRACE("before solve for ", u, " and ", v);

	lamg.solve(rhs, solution);
	return fabs(solution[u] - solution[v]); // TODO: fix weighted case, fix edge IDs!
}


}
