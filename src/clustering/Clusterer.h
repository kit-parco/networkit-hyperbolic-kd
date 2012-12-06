/*
 * Clusterer.h
 *
 *  Created on: 30.10.2012
 *      Author: cls
 */

#ifndef CLUSTERER_H_
#define CLUSTERER_H_

#include "../graph/Graph.h"

namespace EnsembleClustering {

// TODO: import
class Clustering;

class Clusterer {
public:

	Clusterer();

	virtual ~Clusterer();

	virtual Clustering& run(Graph& G) = 0;
};

} /* namespace EnsembleClustering */
#endif /* CLUSTERER_H_ */
