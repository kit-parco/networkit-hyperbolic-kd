/*
 * ClusteringWriter.cpp
 *
 *  Created on: 22.01.2013
 *      Author: cls
 */

#include "ClusteringWriter.h"

namespace EnsembleClustering {

ClusteringWriter::ClusteringWriter() {
	// TODO Auto-generated constructor stub

}

ClusteringWriter::~ClusteringWriter() {
	// TODO Auto-generated destructor stub
}

void ClusteringWriter::write(Clustering& zeta, std::string path) {
	std::ofstream file;
	file.open(path.c_str());

	zeta.forEntries([&](node v, cluster c){
		file << v << " " << c << std::endl;
	});

	file.close();
}

} /* namespace EnsembleClustering */
