/*
 * HashingOverlapper.h
 *
 *  Created on: 31.01.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef HASHINGOVERLAPPER_H_
#define HASHINGOVERLAPPER_H_

#include <functional>

#include "Overlapper.h"

namespace NetworKit {

/**
 * @ingroup overlap
 * Determines the overlap of multiple partitions by hashing partition identifiers.
 */
class HashingOverlapper: public NetworKit::Overlapper {

public:

	/** Default destructor */
	virtual ~HashingOverlapper() = default;

	virtual Partition run(Graph& G, std::vector<Partition>& clusterings);
};

} /* namespace NetworKit */
#endif /* HASHINGOVERLAPPER_H_ */
