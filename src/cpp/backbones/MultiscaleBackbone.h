/*
 * MultiscaleBackbone.h
 *
 *  Created on: 20.06.2014
 *      Author: Gerd Lindner
 */

#ifndef MULTISCALEBACKBONE_H_
#define MULTISCALEBACKBONE_H_

#include "BackboneCalculator.h"
#include "gtest/gtest_prod.h"

namespace NetworKit {

/** 
 * Calculates the multiscale backbone for a given input graph.
 *
 * See "Extracting the multiscale backbone of complex weighted networks" by Serrano et al.
 */
class MultiscaleBackbone : public BackboneCalculator {

public:

	/**
	 * Creates a new instance of the Multiscale Backbone calculator.
	 * @param alpha 		filter parameter
	 */
	MultiscaleBackbone(double alpha);

	Graph calculate(const Graph& graph, const edgeAttribute& attribute);

private:
	//Calculation parameters
	double alpha;

	//Private helper functions
	double getProbability(count degree, edgeweight normalizedWeight);

	FRIEND_TEST(MultiscaleBackboneGTest, testSimpleMultiscaleBackbone);
};

}
/* namespace NetworKit */
#endif /* SIMMELIANBACKBONE_H_ */
