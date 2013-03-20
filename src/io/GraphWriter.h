/*
 * GraphWriter.h
 *
 *  Created on: 30.01.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef GRAPHWRITER_H_
#define GRAPHWRITER_H_

#include "../graph/Graph.h"

namespace NetworKit {

class GraphWriter {
public:

	GraphWriter();

	virtual ~GraphWriter();

	virtual void write(Graph& G, std::string path) = 0;
};

} /* namespace NetworKit */
#endif /* GRAPHWRITER_H_ */
