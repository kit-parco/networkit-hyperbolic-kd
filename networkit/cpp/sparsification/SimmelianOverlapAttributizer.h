/*
 * SimmelianOverlapAttributizer.h
 *
 *  Created on: 22.07.2014
 *      Author: Gerd Lindner
 */

#ifndef SIMMELIANOVERLAPATTRIBUTIZER_H_
#define SIMMELIANOVERLAPATTRIBUTIZER_H_

#include "SimmelianAttributizer.h"
#include <set>

namespace NetworKit {

/**
 * Calculates the Simmelian backbone (paramaetric variant) for a given input graph.
 */
class SimmelianOverlapAttributizer : public SimmelianAttributizer {

public:

	/**
	 * Creates a new instance of the parametric variant of the Simmelian Backbone calculator.
	 * @param maxRank 		the maximum rank that is considered for overlap calculation
	 */
	SimmelianOverlapAttributizer(const Graph& graph, const std::vector<int>& triangles, count maxRank);
	
	virtual std::vector<double> getAttribute() override;

private:
	count maxRank;
};

}
/* namespace NetworKit */
#endif /* SIMMELIANOVERLAPATTRIBUTIZER_H_ */
