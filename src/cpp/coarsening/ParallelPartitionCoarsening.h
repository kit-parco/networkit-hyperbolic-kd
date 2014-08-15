/*
 * ParallelPartitionCoarsening.h
 *
 *  Created on: 03.07.2014
 *      Author: cls
 */

#ifndef PARALLELPARTITIONCOARSENING_H_
#define PARALLELPARTITIONCOARSENING_H_

#include "../Globals.h"
#include "GraphCoarsening.h"
#include "../structures/Partition.h"

namespace NetworKit {

/**
 * @ingroup coarsening
 */
class ParallelPartitionCoarsening: public NetworKit::GraphCoarsening {

public:

	virtual std::pair<Graph, std::vector<node> > run(const Graph& G, const Partition& zeta);


};

} /* namespace NetworKit */

#endif /* PARALLELPARTITIONCOARSENING_H_ */
