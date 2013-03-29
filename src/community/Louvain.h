/*
 * Louvain.h
 *
 *  Created on: 25.02.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu), Henning Meyerhenke (henning.meyerhenke@kit.edu)
 */

#ifndef LOUVAIN_H_
#define LOUVAIN_H_

#include "Clusterer.h"
#include "../coarsening/ClusterContracter.h"
#include "../coarsening/ClusteringProjector.h"
#include "../base/IndexMap.h"
#include "../independentset/Luby.h"
#include "../Globals.h"

#include "omp.h"

namespace NetworKit {

class Louvain: public NetworKit::Clusterer {

protected:
	bool anyChange;	//!< indicates whether any change was made to the clustering in the last pass over the nodes
	std::string parallelism; //!< switch for the kind of parallelization strategy to use


public:

	std::string VERSION;	// algorithm version number - increment in constructor for significant changes to the implementation

	/**
	 * @param[in]	par		parallelization strategy
	 */
	Louvain(std::string par="none");

	virtual ~Louvain();

	virtual Clustering pass(Graph& G);

	virtual Clustering run(Graph& G);

	/**
	 * @return string representation of algorithm and parameters.
	 */
	virtual std::string toString() const;
};

} /* namespace NetworKit */
#endif /* LOUVAIN_H_ */
