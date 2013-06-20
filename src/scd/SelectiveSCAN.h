/*
 * SelectiveSCAN.h
 *
 *  Created on: 14.06.2013
 *      Author: cls
 */

#ifndef SELECTIVESCAN_H_
#define SELECTIVESCAN_H_

#include "SelectiveCommunityDetector.h"

namespace NetworKit {

class SelectiveSCAN: public NetworKit::SelectiveCommunityDetector {

public:
	double epsilon;
	double mu;

public:

	SelectiveSCAN();

	virtual ~SelectiveSCAN();

	virtual std::unordered_map<node, std::unordered_set<node> > run(Graph& G, std::unordered_set<node> seeds);

protected:

	virtual double nodeNodeSimilarity (node u, node v, Graph& G);

	virtual std::pair<bool,std::vector<node>> isCore(node u, Graph& G);
};

} /* namespace NetworKit */
#endif /* SELECTIVESCAN_H_ */
