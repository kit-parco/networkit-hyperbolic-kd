/*
 * GraphIO.cpp
 *
 *  Created on: 09.01.2013
 *      Author: cls
 */

#include "GraphIO.h"

namespace EnsembleClustering {

GraphIO::GraphIO() {
	// TODO Auto-generated constructor stub

}

GraphIO::~GraphIO() {
	// TODO Auto-generated destructor stub
}

void GraphIO::writeEdgeList(Graph& G, std::string path) {

	std::ofstream file;
	file.open(path.c_str());

	G.forEdges([&](node v, node w) {
		file << v << " " << w << std::endl;
	});

	file.close();
	INFO("wrote graph to file: " << path);

}

void GraphIO::writeAdjacencyList(Graph& G, std::string path) {
	std::ofstream file;
	file.open(path.c_str());

	G.forNodes([&](node v) {
		file << v;
		G.forNeighborsOf(v, [&](node x) {
			file << " " << x;
		});
		file << std::endl;
	});
	INFO("wrote graph to file: " << path);
}

} /* namespace EnsembleClustering */
