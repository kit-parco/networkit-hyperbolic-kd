/*
 * CoarseningGTest.cpp
 *
 *  Created on: 20.12.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef NOGTEST

#include "CoarseningGTest.h"

#include "../../graph/GraphGenerator.h"
#include "../../clustering/ClusteringGenerator.h"
#include "../../coarsening/ClusterContracter.h"
#include "../../coarsening/ClusteringProjector.h"

namespace NetworKit {

TEST_F(CoarseningGTest, testClusterContracter) {
	GraphGenerator graphGen;
	int64_t n = 100;
	Graph G = graphGen.makeErdosRenyiGraph(n, 0.5);

	ClusteringGenerator clusteringGen;
	Partition singleton = clusteringGen.makeSingletonClustering(G);
	DEBUG("singleton(upperBound,numberOfElements,numberOfSubsets): ",singleton.upperBound()," ",singleton.numberOfElements()," ",singleton.numberOfSubsets());


	ClusterContracter contracter;
	auto conSingletonPair = contracter.run(G, singleton);
	Graph Gcon = conSingletonPair.first;

	EXPECT_EQ(G.numberOfNodes(), Gcon.numberOfNodes())
			<< "graph contracted according to singleton clustering should have the same number of nodes as original";
	EXPECT_EQ(G.numberOfEdges(), Gcon.numberOfEdges())
			<< "graph contracted according to singletons clustering should have the same number of nodes as original";

	int k = 2; // number of clusters in random clustering
	Partition random = clusteringGen.makeRandomClustering(G, k);
	DEBUG("random(upperBound,numberOfElements,numberOfSubsets): ",random.upperBound()," ",random.numberOfElements()," ",random.numberOfSubsets());
	auto conRandPair = contracter.run(G, random);
	Graph GconRand = conRandPair.first;

	EXPECT_EQ(k, GconRand.numberOfNodes())
			<< "graph contracted according to random clustering should have the same number of nodes as there are clusters.";

}


TEST_F(CoarseningGTest, testClusteringProjectorWithOneClustering) {
	GraphGenerator graphGen;
	int64_t n = 100;
	Graph G0 = graphGen.makeErdosRenyiGraph(n, 0.5);

	// get 1-clustering of G0
	ClusteringGenerator clusteringGen;
	Partition zeta0 = clusteringGen.makeOneClustering(G0);

	// contract G0 according to 1-clusterings
	ClusterContracter contract;
	auto con = contract.run(G0, zeta0);
	std::vector<NodeMap<node> > maps;
	Graph G1 = con.first;
	maps.push_back(con.second);

	Partition zeta1 = clusteringGen.makeOneClustering(G1);

	ClusteringProjector project;
	Partition zetaBack = project.projectBackToFinest(zeta1, maps, G0);

	EXPECT_TRUE(zeta0.equals(zetaBack, G0)) << "zeta^{1->0} and zeta^{0} should be identical"; //FIXME
}


TEST_F(CoarseningGTest, testClusteringProjectorWithSingletonClustering) {
	GraphGenerator graphGen;
	int64_t n = 100;
	Graph G0 = graphGen.makeErdosRenyiGraph(n, 0.5);

	// get 1-clustering of G0
	ClusteringGenerator clusteringGen;
	Partition zeta0 = clusteringGen.makeSingletonClustering(G0);

	// contract G0 according to 1-clusterings
	ClusterContracter contract;
	auto con = contract.run(G0, zeta0);
	std::vector<NodeMap<node> > maps;
	Graph G1 = con.first;
	maps.push_back(con.second);

	Partition zeta1 = clusteringGen.makeSingletonClustering(G1);

	ClusteringProjector project;
	Partition zetaBack = project.projectBackToFinest(zeta1, maps, G0);

	EXPECT_TRUE(zeta0.equals(zetaBack, G0)) << "zeta^{1->0} and zeta^{0} should be identical"; //FIXME
}



} /* namespace NetworKit */

#endif /*NOGTEST */
