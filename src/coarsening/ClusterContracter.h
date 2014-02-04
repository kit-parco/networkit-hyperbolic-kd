/*
 * ClusterContracter.h
 *
 *  Created on: 30.10.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef CLUSTERCONTRACTER_H_
#define CLUSTERCONTRACTER_H_


#include "Contracter.h"
#include "../structures/Partition.h"
#include "../graph/NodeMap.h"

namespace NetworKit {

class ClusterContracter: public Contracter {

public:

	ClusterContracter();

	virtual ~ClusterContracter();

	virtual std::pair<Graph, NodeMap<node> > run(Graph& G, Partition& zeta);




};


} // namespace

#endif /* CLUSTERCONTRACTER_H_ */
