/*
 * CNM.h
 *
 *  Created on: Jun 10, 2013
 *      Author: michi
 */

#ifndef CNM_H_
#define CNM_H_

#include "CommunityDetectionAlgorithm.h"

namespace NetworKit {

/**
 * Clustering algorithm due to Clauset, Newman and Moore.
 * Probably not the fastest possible implementation, but it already uses a priority queue
 * and local updates.
 */
class CNM : public NetworKit::CommunityDetectionAlgorithm {
public:

	Partition run(Graph &graph) override;

	std::string toString() const override {
		return "CNM";
	}
};

}

#endif /* CNM_H_ */
