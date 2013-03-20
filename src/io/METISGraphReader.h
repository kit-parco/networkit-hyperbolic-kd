/*
 * METISGraphReader.h
 *
 *  Created on: 17.01.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef METISGRAPHREADER_H_
#define METISGRAPHREADER_H_

#include "GraphReader.h"

#include "METISParser.h"
#include "../aux/StringTools.h"

namespace NetworKit {

class METISGraphReader: public NetworKit::GraphReader {

public:

	METISGraphReader();

	virtual ~METISGraphReader();

	virtual Graph read(std::string path);
};

} /* namespace NetworKit */
#endif /* METISGRAPHREADER_H_ */
