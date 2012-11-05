/*
 * Modularity.cpp
 *
 *  Created on: 15.10.2012
 *      Author: cls
 */

#include "Modularity.h"

namespace EnsembleClustering {

typedef Graph::Node Node;
typedef Graph::Edge Edge;


Modularity::Modularity() {
	// TODO Auto-generated constructor stub

}

Modularity::~Modularity() {
	// TODO Auto-generated destructor stub
}




double Modularity::mod(Clustering clustering) {

	int k; // number of clusters
	int c; // current cluster id

	double mod; // result
	double total; // total weight of all edges
	double cutSum; // sum of cut weights to all other clusters
	double clusterSum; // result of term which is a sum over all clusters
	for (int j = 0; j < k; ++j) {
		// sum over all cuts
		if (j != c) {
			// cutSum += this->cutweight(clustering[c], clustering[j]);
		}
	}
	// this->weight(cluster) * (2 * total - this->weight(cluster)) - 2 * total * cutSum;
	mod = (1 / 4 * total) * clusterSum;
}



} /* namespace EnsembleClustering */
