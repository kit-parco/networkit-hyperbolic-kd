/*
 * ParallelAgglomerativeClusterer.cpp
 *
 *  Created on: 30.10.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu),
 *      		Henning Meyerhenke (henning.meyerhenke@kit.edu)
 */

#include "ParallelAgglomerativeClusterer.h"
#include "../scoring/ModularityScoring.h"
#include "../matching/PathGrowingMatcher.h"
#include "../coarsening/MatchingContracter.h"
#include "../coarsening/ClusteringProjector.h"

namespace NetworKit {
ParallelAgglomerativeClusterer::ParallelAgglomerativeClusterer(const Graph& G) : CommunityDetectionAlgorithm(G) {};

Partition ParallelAgglomerativeClusterer::run() {
	// copy graph because we make changes due to merges
	Graph Gcopy(G.numberOfNodes(), true); // make weighted copy
	G.forEdges([&](node u, node v, edgeweight w){
		Gcopy.addEdge(u, v, w);
	});

	std::vector<std::vector<node> > mapHierarchy;

	bool repeat = true;
	do {
		// prepare attributes for scoring
		// FIXME: update to new edge attribute system
		//int attrId = Gcopy.addEdgeAttribute_double(0.0);
		int attrId = 0;

		// perform scoring
		TRACE("before scoring graph of size " , Gcopy.numberOfNodes());
		ModularityScoring<double> modScoring(Gcopy);
		modScoring.scoreEdges(attrId);

		// FIXME: so far only sequential
		// compute matching
		PathGrowingMatcher parMatcher;
		Matching M = parMatcher.run(Gcopy);

		// contract graph according to matching, TODO: (and star-like structures)
		MatchingContracter matchingContracter;
		auto GandMap = matchingContracter.run(Gcopy, M);

		// determine if it makes sense to proceed
		count n = Gcopy.numberOfNodes();
		count cn = GandMap.first.numberOfNodes();
		count diff = n - cn;
		repeat = ((diff > 0) &&
				(cn >= MIN_NUM_COMMUNITIES) &&
				((double) diff / (double) n > REL_REPEAT_THRSH)
				); // TODO: last condition: no community becomes too big

		// prepare next iteration if there is one
		if (repeat) {
			Gcopy = GandMap.first;
			mapHierarchy.push_back(GandMap.second);
			TRACE("Repeat agglomeration with graph of size " , Gcopy.numberOfNodes());
		}
	} while (repeat);

	// vertices of coarsest graph are the clusters
	count cn = Gcopy.numberOfNodes();
	Partition zetaCoarse(cn);
	zetaCoarse.allToSingletons();

	// project clustering back to finest graph
	ClusteringProjector projector;
	Partition zeta = projector.projectBackToFinest(zetaCoarse, mapHierarchy,
			G);

	return zeta;
}

std::string ParallelAgglomerativeClusterer::toString() const {
	std::stringstream strm;
	strm << "ParallelAgglomerativeClusterer";
	return strm.str();
}

} /* namespace NetworKit */
