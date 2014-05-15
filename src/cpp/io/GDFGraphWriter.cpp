/*
 * GDFGraphWriter.cpp
 *
 *  Created on: 27.10.2013
 *      Author: Stefan Bertsch
 */

#include "GDFGraphWriter.h"

namespace NetworKit {

void GDFGraphWriter::write(Graph& G, const std::string& path) {
	this->write(G, G.isWeighted(), path);
}

void GDFGraphWriter::writeGeneric(Graph& G, bool weighted, const std::string& path, count dim) {
	std::ofstream file(path);
	assert(file.good());

//	int64_t n = G.numberOfNodes();
//	int64_t m = G.numberOfEdges();
//	file << n << " " << m << " " << (int)weighted << std::endl;

	switch (dim) {
	case 2: {
		file << "*nodedef>name VARCHAR,x DOUBLE,y DOUBLE\n";
		break;
	}
	case 3: {
		file << "nodedef>name VARCHAR,x DOUBLE,y DOUBLE,z DOUBLE\n";
		break;
	}
	default: {
		throw std::runtime_error("Dimension not supported by file format GDF. Skip writing file.");
	}
	}

	Point<float> point; // HM: needs to be float since Graph uses float
	G.forNodes([&](node u) {
		point = G.getCoordinate(u);
		file << u << "," << point.toCsvString() << '\n';
	});

	file << "edgedef>node1 VARCHAR,node2 VARCHAR" << '\n';

	G.forEdges([&](node u, node v) {
		file << u << ',' << v << '\n';
	});
}

void GDFGraphWriter::write(Graph& G, bool weighted, const std::string& path) {
	writeGeneric(G, weighted, path, 2);
}

void GDFGraphWriter::write3D(Graph& G, bool weighted, const std::string& path) {
	writeGeneric(G, weighted, path, 3);
}

} /* namespace NetworKit */
