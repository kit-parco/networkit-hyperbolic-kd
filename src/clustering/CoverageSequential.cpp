/*
 * CoverageSequential.cpp
 *
 *  Created on: 14.02.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include "CoverageSequential.h"

namespace NetworKit {

CoverageSequential::CoverageSequential() {
	// TODO Auto-generated constructor stub

}

CoverageSequential::~CoverageSequential() {
	// TODO Auto-generated destructor stub
}

double CoverageSequential::getQuality(const Clustering& zeta, Graph& G) {

	double cov = 0.0; // term $\frac{\sum_{C \in \zeta} \sum_{ e \in E(C) } \omega(e)}{\sum_{e \in E} \omega(e)}$
	double totalEdgeWeight = G.totalEdgeWeight(); // add edge weight
	DEBUG("total edge weight = " << totalEdgeWeight)


	if (totalEdgeWeight == 0.0) {
		throw std::invalid_argument(
				"Coverage is undefined for graphs without edges (including self-loops).");
	}

	IndexMap<cluster, double> intraEdgeWeight(zeta.upperBound(), 0.0); // cluster -> weight of its internal edges

	// TODO: sum only over intra-cluster edges, not over clusters

	// compute intra-cluster edge weights per cluster
	G.forWeightedEdges(
			[&](node u, node v, edgeweight ew) {
				assert (u < zeta.numberOfNodes());
				assert (v < zeta.numberOfNodes());
				cluster c = zeta[u];
				cluster d = zeta[v];
				if (c == d) {
#ifdef DEBUG
					if ((c >= zeta.upperBound()) || (c < zeta.lowerBound())) {
						ERROR("c=" << c << " = zeta(" << u << ") is larger than upper bound: " << zeta.upperBound());
					}
#endif
					assert ((zeta.lowerBound()) <= c && (c < zeta.upperBound()));
					intraEdgeWeight[c] += ew;
				} // else ignore edge
			});

	double intraEdgeWeightSum = 0.0; //!< term $\sum_{C \in \zeta} \sum_{ e \in E(C) } \omega(e)$
	for (cluster c = zeta.lowerBound(); c < zeta.upperBound(); ++c) {
		intraEdgeWeightSum += intraEdgeWeight[c];
	}
	DEBUG("total intra-cluster edge weight = " << intraEdgeWeightSum);

	cov = intraEdgeWeightSum / totalEdgeWeight;
	DEBUG("coverage = " << cov);

	assert(cov <= 1.0);
	assert(cov >= 0.0);
	return cov;
}

} /* namespace NetworKit */
