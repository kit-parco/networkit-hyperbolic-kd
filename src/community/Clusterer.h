/*
 * Clusterer.h
 *
 *  Created on: 30.10.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef CLUSTERER_H_
#define CLUSTERER_H_

#include "../graph/Graph.h"
#include "../clustering/Clustering.h"

namespace EnsembleClustering {


class Clusterer {
public:

	Clusterer();

	virtual ~Clusterer();

	virtual Clustering run(Graph& G) = 0;

	/**
	 * @return string representation of algorithm and parameters.
	 */
	virtual std::string toString() const;
};

} /* namespace EnsembleClustering */
#endif /* CLUSTERER_H_ */
