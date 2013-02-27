/*
 * ClusteringProjector.cpp
 *
 *  Created on: 07.01.2013
 *      Author: cls
 */

#include "ClusteringProjector.h"

namespace EnsembleClustering {

ClusteringProjector::ClusteringProjector() {
	// TODO Auto-generated constructor stub

}

ClusteringProjector::~ClusteringProjector() {
	// TODO Auto-generated destructor stub
}

Clustering ClusteringProjector::projectBack(Graph& Gcoarse, Graph& Gfine, NodeMap<node>& fineToCoarse,
		Clustering& zetaCoarse) {

	Clustering zetaFine(Gfine.numberOfNodes());
	// DEBUG
	std::ostringstream oss;	oss << "zeta(" << Gfine.getName() << ")"; zetaFine.setName(oss.str());	//C++??!!
	// DEBUG

	Gfine.forNodes([&](node v) {
		node sv = fineToCoarse[v];
		cluster cv = zetaCoarse.clusterOf(sv);
		zetaFine[v] = cv;
	});

	return zetaFine;
}

Clustering ClusteringProjector::projectBackToFinest(Clustering& zetaCoarse,
		std::vector<NodeMap<node> >& maps, Graph& Gfinest) {
	if (zetaCoarse.numberOfNodes() == Gfinest.numberOfNodes()) {
		return zetaCoarse;
	}

	Clustering zetaFine(Gfinest.numberOfNodes());
	zetaFine.setUpperBound(zetaCoarse.upperBound()); // upper bound for ids in zetaFine must be set to upper bound in zetaCoarse, or modularity assertions fail

	// store temporarily coarsest supernode here
	NodeMap<node> tempMap(Gfinest.numberOfNodes());
	// TODO: add NodeMap.setAll, NodeMap.identity...
	// initialize to identity
	Gfinest.parallelForNodes([&](node v){
		tempMap[v] = v;
	});

	// find coarsest supernode for each node
	for (std::vector<NodeMap<node> >::iterator iter = maps.begin(); iter != maps.end(); ++iter) {
		Gfinest.parallelForNodes([&](node v){
			tempMap[v] = (* iter)[tempMap[v]];
		});
	}


	// set clusters for fine nodes
	Gfinest.parallelForNodes([&](node v) {
		cluster sc = zetaCoarse[tempMap[v]];
		zetaFine[v] = sc;
	});

	return zetaFine;
}

Clustering ClusteringProjector::projectCoarseGraphToFinestClustering(Graph& Gcoarse, Graph& Gfinest, std::vector<NodeMap<node> >& maps) {

	Clustering zeta(Gfinest.numberOfNodes());

	// store temporarily coarsest supernode here
	NodeMap<node> super(Gfinest.numberOfNodes());
	// initialize to identity
	Gfinest.parallelForNodes([&](node v){
		super[v] = v;
	});


	// find coarsest supernode for each node
	for (std::vector<NodeMap<node> >::iterator iter = maps.begin(); iter != maps.end(); ++iter) {
		Gfinest.parallelForNodes([&](node v){
			super[v] = (* iter)[super[v]];
		});
	}

	// assign super node id as cluster id
	Gfinest.parallelForNodes([&](node v) {
		zeta[v] = super[v];
	});

	assert (zeta.numberOfClusters() == Gcoarse.numberOfNodes());

	return zeta;

}

} /* namespace EnsembleClustering */
