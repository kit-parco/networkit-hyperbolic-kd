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

class HashingOverlapper: public NetworKit::Overlapper {

public:

	HashingOverlapper();

	virtual ~HashingOverlapper();

	virtual Clustering run(Graph& G, std::vector<Clustering>& clusterings);
};

} /* namespace NetworKit */
#endif /* HASHINGOVERLAPPER_H_ */
