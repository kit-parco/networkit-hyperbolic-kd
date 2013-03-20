/*
 * ClusteringReader.h
 *
 *  Created on: 15.02.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef CLUSTERINGREADER_H_
#define CLUSTERINGREADER_H_

#include <fstream>

#include "../clustering/Clustering.h"
#include "../graph/Graph.h"

namespace NetworKit {

class ClusteringReader {

public:

	ClusteringReader();

	virtual ~ClusteringReader();

	/**
	 * Read a clustering from a file. File format:
	 * 		line n contains cluster id of node (n - 1)
	 *
	 * @param[in]	path	path to file
	 */
	virtual Clustering read(std::string path);
};

} /* namespace NetworKit */
#endif /* CLUSTERINGREADER_H_ */
