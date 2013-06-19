/*
 * ClusteringAlgoGTest.cpp
 *
 *  Created on: 10.01.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include "ClusteringAlgoGTest.h"

#ifndef NOGTEST

namespace NetworKit {

TEST_F(ClusteringAlgoGTest, testLabelPropagationOnUniformGraph) {
	GraphGenerator graphGenerator;
	int n = 100;
	Graph G = graphGenerator.makeErdosRenyiGraph(n, 0.2);

	LabelPropagation lp;
	Clustering zeta = lp.run(G);

	EXPECT_TRUE(zeta.isProper(G)) << "the resulting partition should be a proper clustering";

	Modularity modularity;
	double mod = modularity.getQuality(zeta, G);
	DEBUG("modularity produced by LabelPropagation: " << mod);
	EXPECT_GE(1.0, mod) << "valid modularity values are in [-0.5, 1]";
	EXPECT_LE(-0.5, mod) << "valid modularity values are in [-0.5, 1]";
}


TEST_F(ClusteringAlgoGTest, testLabelPropagationOnClusteredGraph_ForNumberOfClusters) {
	GraphGenerator graphGenerator;
	int64_t n = 100;
	count k = 3; // number of clusters
	Graph G = graphGenerator.makeClusteredRandomGraph(n, k, 1.0, 0.001);

	LabelPropagation lp;
	Clustering zeta = lp.run(G);

	Modularity modularity;
	double mod = modularity.getQuality(zeta, G);
	DEBUG("modularity produced by LabelPropagation: " << mod);

	EXPECT_TRUE(zeta.isProper(G)) << "the resulting partition should be a proper clustering";
	EXPECT_EQ(k, zeta.numberOfClusters()) << " " << k << " clusters are easy to detect";

}


TEST_F(ClusteringAlgoGTest, testLabelPropagationOnClusteredGraph_ForEquality) {
	int64_t n = 100;

	GraphGenerator graphGen;
	Graph Gtrash = graphGen.makeCompleteGraph(n);

	count k = 3; // number of clusters
	ClusteringGenerator clusteringGen;
	Clustering reference = clusteringGen.makeRandomClustering(Gtrash, 3);
	assert (reference.numberOfClusters() == k);

	Graph G = graphGen.makeClusteredRandomGraph(reference, 1.0, 0.0);	// LabelPropagation is very bad at discerning clusters and needs this large pin/pout difference

	LabelPropagation lp;
	Clustering zeta = lp.run(G);

	Modularity modularity;
	double mod = modularity.getQuality(zeta, G);
	DEBUG("modularity produced by LabelPropagation: " << mod);
	DEBUG("number of clusters produced by LabelPropagation: k=" << zeta.numberOfClusters());

	EXPECT_TRUE(zeta.isProper(G)) << "the resulting partition should be a proper clustering";
	EXPECT_TRUE(zeta.equals(reference, G)) << "LP should detect exactly the reference clustering";

}



TEST_F(ClusteringAlgoGTest, testLabelPropagationOnDisconnectedGraph) {
	GraphGenerator graphGenerator;
	int n = 100;
	int k = 2; // number of clusters
	Graph G = graphGenerator.makeClusteredRandomGraph(n, k, 1.0, 0.0);

	LabelPropagation lp;
	Clustering zeta = lp.run(G);

	Modularity modularity;
	double mod = modularity.getQuality(zeta, G);
	DEBUG("modularity produced by LabelPropagation: " << mod);

	EXPECT_TRUE(zeta.isProper(G)) << "the resulting partition should be a proper clustering";
	EXPECT_EQ(k, zeta.numberOfClusters()) << " " << k << " clusters are easy to detect";

}


TEST_F(ClusteringAlgoGTest, testLabelPropagationOnSingleNodeWithSelfLoop) {
	Graph G(1);
	node v = 0;
	G.setWeight(v, v, 42.0);

	LabelPropagation lp;
	Clustering zeta = lp.run(G);

	EXPECT_TRUE(zeta.isProper(G));
	EXPECT_TRUE(zeta.isSingletonClustering(G));
	EXPECT_TRUE(zeta.isOneClustering(G));

	Modularity modularity;
	double mod = modularity.getQuality(zeta, G);
	DEBUG("modularity produced by LabelPropagation: " << mod);

}



TEST_F(ClusteringAlgoGTest, testLabelPropagationOnManySmallClusters) {
	int64_t n = 1000;
	int k = 100; // number of clusters
	double pin = 1.0;
	double pout = 0.0;

	GraphGenerator graphGen;
	std::pair<Graph, Clustering> G_ref = graphGen.makeClusteredRandomGraphWithReferenceClustering(n, k, pin, pout);


	LabelPropagation lp;
	Clustering zeta = lp.run(G_ref.first);

	Modularity modularity;
	double mod = modularity.getQuality(zeta, G_ref.first);
	DEBUG("modularity produced by LabelPropagation: " << mod);
	DEBUG("number of clusters produced by LabelPropagation: k=" << zeta.numberOfClusters());

	EXPECT_TRUE(zeta.isProper(G_ref.first)) << "the resulting partition should be a proper clustering";
	EXPECT_TRUE(zeta.equals(G_ref.second, G_ref.first)) << "Can LabelPropagation detect the reference clustering?";

}

TEST_F(ClusteringAlgoGTest, testLouvain) {
	count n = 500;
	count k = 25;
	double pin = 0.9;
	double pout = 0.005;
	GraphGenerator graphGen;
	Graph G = graphGen.makeClusteredRandomGraph(n, k, pin, pout);

	Louvain louvain;
	Clustering zeta = louvain.run(G);

	INFO("number of clusters: " << zeta.numberOfClusters());

	Modularity modularity;
	INFO("modularity: " << modularity.getQuality(zeta, G));

}


TEST_F(ClusteringAlgoGTest, testLouvainParallelSimple) {
	count n = 500;
	count k = 25;
	double pin = 0.9;
	double pout = 0.005;
	GraphGenerator graphGen;
	Graph G = graphGen.makeClusteredRandomGraph(n, k, pin, pout);

	Louvain louvain("simple");
	Clustering zeta = louvain.run(G);

	INFO("number of clusters: " << zeta.numberOfClusters());

	Modularity modularity;
	INFO("modularity: " << modularity.getQuality(zeta, G));

}

/*
TEST_F(ClusteringAlgoGTest, testLouvainParallel2Naive) {
	count n = 1000;
	count k = 100;
	double pin = 1.0;
	double pout = 0.005;
	GraphGenerator graphGen;
	Graph G = graphGen.makeClusteredRandomGraph(n, k, pin, pout);

	LouvainParallel louvain;
	Clustering zeta = louvain.run(G);

	INFO("number of clusters: " << zeta.numberOfClusters());

	Modularity modularity;
	INFO("modularity: " << modularity.getQuality(zeta, G));

}
*/


TEST_F(ClusteringAlgoGTest, testLouvainParallelBalanced) {
	count n = 500;
	count k = 25;
	double pin = 0.9;
	double pout = 0.005;
	GraphGenerator graphGen;
	Graph G = graphGen.makeClusteredRandomGraph(n, k, pin, pout);

	Louvain louvain("balanced");
	Clustering zeta = louvain.run(G);

	INFO("number of clusters: " << zeta.numberOfClusters());

	Modularity modularity;
	INFO("modularity: " << modularity.getQuality(zeta, G));

}



TEST_F(ClusteringAlgoGTest, testLouvainIndependent) {
	count n = 500;
	count k = 25;
	double pin = 0.9;
	double pout = 0.005;
	GraphGenerator graphGen;
	Graph G = graphGen.makeClusteredRandomGraph(n, k, pin, pout);

	Louvain louvain("independent");
	Clustering zeta = louvain.run(G);

	INFO("number of clusters: " << zeta.numberOfClusters());

	Modularity modularity;
	INFO("modularity: " << modularity.getQuality(zeta, G));

}


TEST_F(ClusteringAlgoGTest, testCNM) {
	count n = 500;
	count k = 25;
	double pin = 0.9;
	double pout = 0.005;
	GraphGenerator graphGen;
	Graph G = graphGen.makeClusteredRandomGraph(n, k, pin, pout);

	Modularity modularity;
	CNM cnm;

	Clustering clustering = cnm.run(G);
	INFO("CNM number of clusters: " << clustering.numberOfClusters());
	INFO("modularity clustered random graph: " << modularity.getQuality(clustering, G));
}


TEST_F(ClusteringAlgoGTest, testCNMandLouvain) {
	Modularity modularity;
	CNM cnm;
	Louvain louvain;
	METISGraphReader reader;
	Graph jazz = reader.read("input/jazz.graph");
	// this takes much longer than a unit test should
	// Graph blog = reader.read("input/polblogs.graph");



	// *** jazz graph
	// CNM
	Clustering clustering = cnm.run(jazz);
	INFO("CNM number of jazz clusters: " << clustering.numberOfClusters());
	INFO("CNM modularity jazz graph: " << modularity.getQuality(clustering, jazz));

	// Louvain
	clustering = louvain.run(jazz);
	INFO("Louvain number of jazz clusters: " << clustering.numberOfClusters());
	INFO("Louvain modularity jazz graph: " << modularity.getQuality(clustering, jazz));


//	// *** blog graph
//	// CNM
//	clustering = cnm.run(blog);
//	INFO("CNM number of blog clusters: " << clustering.numberOfClusters());
//	INFO("CNM modularity blog graph: " << modularity.getQuality(clustering, jazz));
//
//	// Louvain
//	clustering = louvain.run(blog);
//	INFO("Louvain number of blog clusters: " << clustering.numberOfClusters());
//	INFO("Louvain modularity blog graph: " << modularity.getQuality(clustering, jazz));
}


TEST_F(ClusteringAlgoGTest, testAgglomerativeAndLouvain) {
	Modularity modularity;
	ParallelAgglomerativeClusterer aggl;
	Louvain louvain;
	METISGraphReader reader;
	Graph jazz = reader.read("input/jazz.graph");
	Graph blog = reader.read("input/polblogs.graph");


	// *** jazz graph
	// aggl
	Clustering clustering = aggl.run(jazz);
	INFO("Match-AGGL number of jazz clusters: " << clustering.numberOfClusters());
	INFO("Match-AGGL modularity jazz graph:   " << modularity.getQuality(clustering, jazz));

	// Louvain
	clustering = louvain.run(jazz);
	INFO("Louvain number of jazz clusters: " << clustering.numberOfClusters());
	INFO("Louvain modularity jazz graph:   " << modularity.getQuality(clustering, jazz));


	// *** blog graph
	// CNM
	clustering = aggl.run(blog);
	INFO("Match-AGGL number of blog clusters: " << clustering.numberOfClusters());
	INFO("Match-AGGL modularity blog graph:   " << modularity.getQuality(clustering, jazz));

	// Louvain
	clustering = louvain.run(blog);
	INFO("Louvain number of blog clusters: " << clustering.numberOfClusters());
	INFO("Louvain modularity blog graph:   " << modularity.getQuality(clustering, jazz));
}



} /* namespace NetworKit */

#endif /*NOGTEST */
