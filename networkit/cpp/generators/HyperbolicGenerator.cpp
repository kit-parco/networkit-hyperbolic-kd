/*
 * HyperbolicGenerator.cpp
 *
 *  Created on: 20.05.2014
 *      Author: Moritz v. Looz (moritz.looz-corswarem@kit.edu)
 * This generator is described in the publication
 * "Generating Random Hyperbolic Graphs in Subquadratic Time" by
 * Moritz von Looz, Henning Meyerhenke and Roman Prutkin, presented at ISAAC 2015.
 *
 * The model of hyperbolic random graphs is presented in
 * "Hyperbolic geometry of complex networks. Physical Review E, 82:036106, Sep 2010." by
 *   Dmitri Krioukov, Fragkiskos Papadopoulos, Maksim Kitsak, Amin Vahdat, and Marian Boguna
 *
 */

#include <cstdlib>
#include <random>
#include <math.h>
#include <assert.h>
#include <omp.h>
#include <algorithm>

#include "../graph/GraphBuilder.h"
#include "HyperbolicGenerator.h"
#include "../auxiliary/Random.h"
#include "../auxiliary/Parallel.h"

namespace NetworKit {

/**
 * Construct a generator for n nodes and m edges
 */
HyperbolicGenerator::HyperbolicGenerator(count n, double avgDegree, double plexp, double T) {
	nodeCount = n;
	if (plexp < 2) throw std::runtime_error("Exponent of power-law degree distribution must be > 2");
	if (T < 0 || T == 1) throw std::runtime_error("Temperature must be non-negative and not 1.");//Really necessary? Graphs with T=1 can be generated, only their degree is not controllable
	if (avgDegree > n) throw std::runtime_error("Average Degree must be at most n");
	if (T < 1) {
		alpha = 0.5*(plexp-1);
	} else {
		alpha = 0.5*(plexp-1)/T;
	}

	double R = HyperbolicSpace::hyperbolicAreaToRadius(n);
	double targetR = HyperbolicSpace::getTargetRadius(n, n*avgDegree/2, alpha, T);// 2*log(8*n / (M_PI*(m/n)*2));
	stretch = targetR / R;
	factor = 1;
	temperature=T;
	initialize();
}

void HyperbolicGenerator::initialize() {
	if (temperature == 0) {
		capacity = 1000;
	} else {
		capacity = 10;
	}
	theoreticalSplit = false;
	threadtimers.resize(omp_get_max_threads());
	balance = 0.5;
}

Graph HyperbolicGenerator::generate() {
	return generate(nodeCount, factor, alpha, stretch, temperature);
}

Graph HyperbolicGenerator::generate(count n, double distanceFactor, double alpha, double stretchradius, double T) {
	double R = stretchradius*HyperbolicSpace::hyperbolicAreaToRadius(n);
	assert(R > 0);
	vector<double> angles(n);
	vector<double> radii(n);
	double r = HyperbolicSpace::hyperbolicRadiusToEuclidean(R);
	assert(r > 0);
	assert(r < 1);

	//sample points randomly
	HyperbolicSpace::fillPoints(angles, radii, stretchradius, alpha);
	vector<index> permutation(n);

	index p = 0;
	std::generate(permutation.begin(), permutation.end(), [&p](){return p++;});

	//can probably be parallelized easily, but doesn't bring much benefit
	Aux::Parallel::sort(permutation.begin(), permutation.end(), [&angles,&radii](index i, index j){return angles[i] < angles[j] || (angles[i] == angles[j] && radii[i] < radii[j]);});

	vector<double> anglecopy(n);
	vector<double> radiicopy(n);

	#pragma omp parallel for
	for (index j = 0; j < n; j++) {
		anglecopy[j] = angles[permutation[j]];
		radiicopy[j] = radii[permutation[j]];
	}

	INFO("Generated Points");
	return generate(anglecopy, radiicopy, r, R*distanceFactor, T);
}

Graph HyperbolicGenerator::generate(const vector<double> &angles, const vector<double> &radii, double r, double thresholdDistance, double T) {
	Aux::Timer timer;
	timer.start();
	index n = angles.size();
	assert(radii.size() == n);
	assert(r > 0);
	assert(r < 1);
	assert(alpha > 0);
	Quadtree<index> quad(r, theoreticalSplit, alpha, capacity, balance);

	//initialize a graph builder for n nodes and an undirected, unweighted graph with direct swap
	for (index i = 0; i < n; i++) {
		assert(radii[i] < r);
		quad.addContent(i, angles[i], radii[i]);
	}

	quad.trim();
	quad.recount();
	assert(quad.size() == n);
	timer.stop();
	INFO("Filled Quadtree, took ", timer.elapsedMilliseconds(), " milliseconds.");

	return generate(angles, radii, quad, thresholdDistance, T);
}

Graph HyperbolicGenerator::generateExternal(const vector<double> &angles, const vector<double> &radii, double k, double gamma, double T) {
	count n = angles.size();
	assert(angles.size() == radii.size());
	vector<double> radiiPoincare(n);
	double targetR = HyperbolicSpace::getTargetRadius(n, n*k/2, (gamma-1)/2, T, 0.001);
	for (index i = 0; i < n; i++) {
		assert(angles[i] > 0);
		assert(angles[i] <= 2*M_PI);
		assert(radii[i] >= 0);
		if (radii[i] > targetR) {
			DEBUG("Coordinate radii[",i, "] = ", radii[i],  " > ", targetR, " = targetR");
			targetR = std::nextafter(radii[i], std::numeric_limits<double>::max());
		}

		assert(radii[i] <= targetR);
		radiiPoincare[i] = HyperbolicSpace::hyperbolicRadiusToEuclidean(radii[i]);
	}

	double r = HyperbolicSpace::hyperbolicRadiusToEuclidean(targetR);

	for (double radius : radiiPoincare) {
		if (r <= radius) r = std::nextafter(radius, std::numeric_limits<double>::max());
	}

	return generate(angles, radiiPoincare, r, targetR, T);
}

Graph HyperbolicGenerator::generateCold(const vector<double> &angles, const vector<double> &radii, const Quadtree<index> &quad, double thresholdDistance) {
	index n = angles.size();
	assert(radii.size() == n);
	assert(quad.size() == n);

	for (index i = 0; i < n; i++) {
		assert(radii[i] < quad.getMaxRadius());
	}

	for (index elem : quad.getElements()) {
		assert(elem < n);
	}

	Aux::Timer timer;
	timer.start();
	vector<double> empty;
	GraphBuilder result(n, false, false);
	bool suppressLeft = directSwap ? false : std::is_sorted(angles.begin(), angles.end());//relying on lazy evaluation here

	#pragma omp parallel
	{
		index id = omp_get_thread_num();
		threadtimers[id].start();
		#pragma omp for schedule(guided) nowait
		for (index i = 0; i < n; i++) {
			//get neighbours for node i
			count expectedDegree = (4/M_PI)*n*exp(-HyperbolicSpace::EuclideanRadiusToHyperbolic(radii[i])/2);//TODO: adapt for alpha!=1
			vector<index> near;
			near.reserve(expectedDegree*1.1);
			quad.getElementsInHyperbolicCircle(HyperbolicSpace::polarToCartesian(angles[i], radii[i]), thresholdDistance, suppressLeft, near);
			//count realDegree = near.size();
			//std::swap(expectedDegree, realDegree);//dummy statement for debugging
			if (directSwap) {
				auto newend = std::remove(near.begin(), near.end(), i); //no self loops!
				if (newend != near.end()) {
					assert(newend+1 == near.end());
					assert(*(newend)==i);
					near.pop_back();//std::remove doesn't remove element but swaps it to the end
				}
				result.swapNeighborhood(i, near, empty, false);
			} else {
				for (index j : near) {
					if (j > i) result.addHalfEdge(i,j);
				}
			}

		}
		threadtimers[id].stop();
	}

	timer.stop();
	INFO("Generating Edges took ", timer.elapsedMilliseconds(), " milliseconds.");
	return result.toGraph(true);
}

Graph HyperbolicGenerator::generate(const vector<double> &angles, const vector<double> &radii, Quadtree<index> &quad, double thresholdDistance, double T) {
	if (T < 0) throw std::runtime_error("Temperature cannot be negative.");
	if (T == 0) return generateCold(angles, radii, quad, thresholdDistance);
	assert(T > 0);
	count n = angles.size();
	assert(radii.size() == n);
	assert(quad.size() == n);
	assert(thresholdDistance > 0);

	bool anglesSorted = std::is_sorted(angles.begin(), angles.end());

	//now define lambda
	double beta = 1/T;
	assert(beta == beta);
	auto edgeProb = [beta, thresholdDistance](double distance) -> double {return 1 / (exp(beta*(distance-thresholdDistance)/2)+1);};

	//just for now, write out coordinates
	for (index i = 0; i < n; i++) {
		auto cartesian = HyperbolicSpace::polarToCartesian(angles[i], HyperbolicSpace::EuclideanRadiusToHyperbolic(radii[i]));
	}

	//get Graph
	GraphBuilder result(n, false, false);//no direct swap with probabilistic graphs, sorry
	count totalCandidates = 0;
	#pragma omp parallel for
	for (index i = 0; i < n; i++) {
		vector<index> near;
		totalCandidates += quad.getElementsProbabilistically(HyperbolicSpace::polarToCartesian(angles[i], radii[i]), edgeProb, anglesSorted, near);
		for (index j : near) {
			if (j >= n) ERROR("Node ", j, " prospective neighbour of ", i, " does not actually exist. Oops.");
			if (j > i) {
				result.addHalfEdge(i, j);
			}
		}

	}
	DEBUG("Candidates tested: ", totalCandidates);
	return result.toGraph(true, true);

}
}
