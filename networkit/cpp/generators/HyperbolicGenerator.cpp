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

#include "../graph/GraphBuilder.h"
#include "HyperbolicGenerator.h"
#include "Quadtree/Quadtree.h"
#include "../auxiliary/Random.h"
#include "../auxiliary/ProgressMeter.h"
#include "../auxiliary/Timer.h"

namespace NetworKit {


HyperbolicGenerator::HyperbolicGenerator() {
	stretch = 1;
	alpha = 1;
	factor = 1;
	nodeCount = 10000;
}
/**
 * Construct a generator for n nodes and the parameters, t, alpha and s.
 */
HyperbolicGenerator::HyperbolicGenerator(count n, double distanceFactor, double alpha, double stretchradius) {
	nodeCount = n;
	stretch = stretchradius;
	factor = distanceFactor;
	this->alpha = alpha;
}
/**
 * Construct a generator for n nodes and m edges
 */
HyperbolicGenerator::HyperbolicGenerator(count n, count m) {
	nodeCount = n;
	double R = acosh((double)n/(2*M_PI)+1);
	double targetR = 2*log(8*n / (M_PI*(m/n)*2));
	stretch = targetR / R;
	factor = 1;
	alpha = 1;
}

HyperbolicGenerator::~HyperbolicGenerator() {
}


Graph HyperbolicGenerator::generate() {
	return generate(nodeCount, factor, alpha, stretch);
}

Graph HyperbolicGenerator::generate(count n, double distanceFactor, double alpha, double stretchradius) {
	double R = stretchradius*acosh((double)n/(2*M_PI)+1);
	vector<double> angles(n);
	vector<double> radii(n);
	double r = HyperbolicSpace::hyperbolicRadiusToEuclidean(R);
	//sample points randomly
	HyperbolicSpace::fillPoints(&angles, &radii, stretchradius, alpha);
	INFO("Generated Points");
	return generate(&angles, &radii, r, R*distanceFactor);
}

double HyperbolicGenerator::expectedNumberOfEdges(count n, double stretch) {
	double R = stretch*acosh((double)n/(2*M_PI)+1);
	return (8 / M_PI) * n * exp(-R/2)*(n/2);
}

std::map<index, Point<float> > HyperbolicGenerator::getCoordinates(vector<double> &angles, vector<double> &radii) {
	count n = angles.size();
	assert(radii.size() == n);
	std::map<index, Point<float> > result;
	for (index i = 0; i < angles.size(); i++) {
		Point2D<double> coord = HyperbolicSpace::polarToCartesian(angles[i], radii[i]);
		Point<float> temp(coord[0], coord[1]);
		result.emplace(i, temp);
	}
	return result;
}

Graph HyperbolicGenerator::generate(vector<double> * angles, vector<double> * radii, double R, double thresholdDistance) {
	Aux::Timer timer;
	timer.start();
	index n = angles->size();
	assert(radii->size() == n);
	Quadtree<index> quad(R);

	//initialize a graph builder for n nodes and an undirected, unweighted graph with direct swap
	GraphBuilder result(n, false, false, true);
	for (index i = 0; i < n; i++) {
		assert(radii->at(i) < R);
		quad.addContent(i, angles->at(i), radii->at(i));
	}
	timer.stop();
	INFO("Filled Quadtree, took ", timer.elapsedMilliseconds(), " milliseconds.");
	timer.start();

	Aux::ProgressMeter progress(n, 1000);
	#pragma omp parallel for schedule(dynamic, 1000)
	for (index i = 0; i < n; i++) {
		//get neighbours for node i
		vector<index> near = quad.getCloseElements(HyperbolicSpace::polarToCartesian(angles->at(i), radii->at(i)), thresholdDistance);
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

	timer.stop();
	INFO("Generated Graph, took ", timer.elapsedMilliseconds(), " milliseconds.");
	return result.toGraph(true);
}
}
