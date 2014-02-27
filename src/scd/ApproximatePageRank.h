/*
 * ApproximatePageRank.h
 *
 *  Created on: 26.02.2014
 *      Author: Henning
 */

#ifndef APPROXIMATEPAGERANK_H_
#define APPROXIMATEPAGERANK_H_

#include <vector>
#include "../graph/Graph.h"

namespace NetworKit {

class ApproximatePageRank {
protected:
	const Graph& G;
	double alpha;
	double eps;

	void push(node u, node seed, std::vector<double>& pr, std::vector<double>& residual);

public:
	ApproximatePageRank(Graph& g, double alpha, double epsilon);
	virtual ~ApproximatePageRank();

	std::vector<double> run(node seed);
};

} /* namespace NetworKit */
#endif /* APPROXIMATEPAGERANK_H_ */
