/*
 * HyperbolicGenerator.cpp
 *
 *  Created on: 20.05.2014
 *      Author: Moritz v. Looz (moritz.looz-corswarem@kit.edu)
 * This generator is a simplified version of the model presented in
 * "Hyperbolic geometry of complex networks. Physical Review E, 82:036106, Sep 2010." by
 *   Dmitri Krioukov, Fragkiskos Papadopoulos, Maksim Kitsak, Amin Vahdat, and Marian Boguna
 *
 */

#include <cstdlib>
#include <random>
#include <math.h>
#include <assert.h>
#include <omp.h>

#include "../graph/GraphBuilder.h"
#include "HyperbolicGenerator.h"
#include "Quadtree/Quadtree.h"
#include "../auxiliary/Random.h"
#include "../auxiliary/ProgressMeter.h"


namespace NetworKit {


HyperbolicGenerator::HyperbolicGenerator() {
	stretch = 1;
	alpha = 1;
	factor = 1;
	nodeCount = 10000;
	initialize();
}
/**
 * Construct a generator for n nodes and the parameters, t, alpha and s.
 */
HyperbolicGenerator::HyperbolicGenerator(count n, double distanceFactor, double alpha, double stretchradius) {
	nodeCount = n;
	stretch = stretchradius;
	factor = distanceFactor;
	this->alpha = alpha;
	initialize();
}
/**
 * Construct a generator for n nodes and m edges
 */
HyperbolicGenerator::HyperbolicGenerator(count n, count m) {
	nodeCount = n;
	double R = HyperbolicSpace::hyperbolicAreaToRadius(n);
	double targetR = 2*log(8*n / (M_PI*(m/n)*2));
	stretch = targetR / R;
	factor = 1;
	alpha = 1;
	initialize();
}

void HyperbolicGenerator::initialize() {
	capacity = 1000;
	theoreticalSplit = false;
	threadtimers.resize(omp_get_max_threads());
}

Graph HyperbolicGenerator::generate() {
	return generate(nodeCount, factor, alpha, stretch);
}

Graph HyperbolicGenerator::generate(count n, double distanceFactor, double alpha, double stretchradius) {
	double R = stretchradius*HyperbolicSpace::hyperbolicAreaToRadius(n);
	vector<double> angles(n);
	vector<double> radii(n);
	double r = HyperbolicSpace::hyperbolicRadiusToEuclidean(R);
	//sample points randomly

	HyperbolicSpace::fillPoints(angles, radii, stretchradius, alpha);
	vector<index> permutation(n);

	index p = 0;
	std::generate(permutation.begin(), permutation.end(), [&p](){return p++;});

	//can probably be parallelized easily, but doesn't bring much benefit
	std::sort(permutation.begin(), permutation.end(), [&angles,&radii](index i, index j){return angles[i] < angles[j] || (angles[i] == angles[j] && radii[i] < radii[j]);});

	vector<double> anglecopy(n);
	vector<double> radiicopy(n);

	#pragma omp parallel for
	for (index j = 0; j < n; j++) {
		anglecopy[j] = angles[permutation[j]];
		radiicopy[j] = radii[permutation[j]];
	}

	INFO("Generated Points");
	return generate(anglecopy, radiicopy, r, R*distanceFactor);
}

double HyperbolicGenerator::expectedNumberOfEdges(count n, double stretch) {
	double R = stretch*HyperbolicSpace::hyperbolicAreaToRadius(n);
	return (8 / M_PI) * n * exp(-R/2)*(n/2);
}

Graph HyperbolicGenerator::generate(const vector<double> &angles, const vector<double> &radii, double R, double thresholdDistance) {
	Aux::Timer timer;
	timer.start();
	index n = angles.size();
	assert(radii.size() == n);
	Quadtree<index> quad(R, theoreticalSplit, alpha, capacity);

	//initialize a graph builder for n nodes and an undirected, unweighted graph with direct swap
	GraphBuilder result(n, false, false, true);
	for (index i = 0; i < n; i++) {
		assert(radii[i] < R);
		quad.addContent(i, angles[i], radii[i]);
	}
	quad.trim();
	timer.stop();
	INFO("Filled Quadtree, took ", timer.elapsedMilliseconds(), " milliseconds.");
	timer.start();

	Aux::ProgressMeter progress(n, 1000);
	#pragma omp parallel
	{
		index id = omp_get_thread_num();
		threadtimers[id].start();
		#pragma omp for schedule(guided) nowait
		for (index i = 0; i < n; i++) {
			//get neighbours for node i
			vector<index> near = quad.getElementsInHyperbolicCircle(HyperbolicSpace::polarToCartesian(angles[i], radii[i]), thresholdDistance);
			for (index j : near) {
				if (i != j) {
					//we add half-edges from both directions at the same time. Due to the symmetry of distances, the correct edges will be formed in parallel without a need for deduplication
					result.addEdge(i,j);
				}
			}
			if (i % 1000 == 0) {
				#pragma omp critical (progress)
				{
					progress.signal(i);
				}
			}
		}
		threadtimers[id].stop();
	}

	timer.stop();
	INFO("Generating Edges took ", timer.elapsedMilliseconds(), " milliseconds.");
	return result.toGraph(true);
}
}
