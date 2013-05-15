/*
 * METISGraphReader.cpp
 *
 *  Created on: 17.01.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include "METISGraphReader.h"

namespace NetworKit {

METISGraphReader::METISGraphReader() {
	// TODO Auto-generated constructor stub

}

METISGraphReader::~METISGraphReader() {
	// TODO Auto-generated destructor stub
}


Graph METISGraphReader::read(std::string path) {

	METISParser parser(path);

	std::tuple<int64_t, int64_t, int64_t> header = parser.getHeader();
	int64_t n = std::get<0>(header);
	int64_t m = std::get<1>(header);
	int64_t weighted = std::get<2>(header);


	Graph G(n);

	std::string graphName = Aux::StringTools::split(Aux::StringTools::split(path, '/').back(), '.').front();

	G.setName(graphName);

	std::cout << "[BEGIN] reading graph G(n=" << n << ", m=" << m << ") from METIS file: " << std::flush;	// progress bar follows

	double p = 0.0; // percentage for progress bar
	node u = 0; // begin with 0

	if (weighted == 0) {
		while (parser.hasNext()) {
			std::vector<node> adjacencies = parser.getNext();
			for (int i=0; i < adjacencies.size(); i++) {
				node v = adjacencies[i] - 1; 	// METIS-indices are 1-based
				assert (v >= 0);
				if (u <= v) { // self-loops are allowed
					G.addEdge(u, v);
				}
			}
			u += 1; // next node
			if ((u % 0x10000) == 0) {
				p = ((double) u / (double) n) * 100;
			}
			std::cout << p << "% " << std::flush;
		}
		std::cout << "[DONE]" << std::endl;
		return G;
	} else {
		while (parser.hasNext()) {
			std::vector<node> adjacencies = parser.getNext();
			for (int i=0; i < adjacencies.size(); i+=2) {
				node v = adjacencies[i] - 1; 	// METIS-indices are 1-based
				assert (v >= 0);
				if (u <= v) { // self-loops are allowed
					G.addEdge(u, v);
					G.setWeight(u, v, adjacencies[i+1]);
					assert(adjacencies[i+1] > 0); // FIXME: store adjacencies and weights separately
				}
			}
			u += 1; // next node
			if ((u % 0x10000) == 0) {
				p = ((double) u / (double) n) * 100;
			}
			std::cout << p << "% " << std::flush;
		}
		std::cout << "[DONE]" << std::endl;
		return G;
	}
}

} /* namespace NetworKit */
