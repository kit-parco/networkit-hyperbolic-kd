/*
 * HyperbolicGraphReader.h
 *
 *  Created on: 29.11.2016
 *      Author: moritzl
 */

#ifndef HYPERBOLICGRAPHREADER_H_
#define HYPERBOLICGRAPHREADER_H_

#include <vector>
#include <string>

#include "../graph/Graph.h"

namespace NetworKit {

class HyperbolicGraphReader {
public:
	static void readGraph(std::string filenamePrefix, Graph &G, std::vector<double> &angleOutput, std::vector<double> &radiiOutput);
};

} /* namespace NetworKit */
#endif /* HYPERBOLICGRAPHREADER_H_ */
