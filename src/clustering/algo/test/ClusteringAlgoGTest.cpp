/*
 * ClusteringAlgoGTest.cpp
 *
 *  Created on: 10.01.2013
 *      Author: cls
 */

#include "ClusteringAlgoGTest.h"

namespace EnsembleClustering {

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
	int k = 3; // number of clusters
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

	int k = 3; // number of clusters
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


TEST_F(ClusteringAlgoGTest, testLabelPropagation_OnCliqueGraph_ScaledUp) {
	int64_t n = 1000;

	GraphGenerator graphGen;
	Graph Gtrash(n);

	int k = 100; // number of clusters
	ClusteringGenerator clusteringGen;
	Clustering reference = clusteringGen.makeRandomClustering(Gtrash, k);
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
	G.setWeight(1, 42.0);

	LabelPropagation lp;
	Clustering zeta = lp.run(G);

	EXPECT_TRUE(zeta.isProper(G));
	EXPECT_TRUE(zeta.isSingletonClustering(G));
	EXPECT_TRUE(zeta.isOneClustering(G));

	Modularity modularity;
	double mod = modularity.getQuality(zeta, G);
	DEBUG("modularity produced by LabelPropagation: " << mod);

}


TEST_F(ClusteringAlgoGTest, testLabelPropagationOnUnclearClustering) {
	int64_t n = 100;
	int k = 3; // number of clusters
	double pin = 0.6;
	double pout = 0.2;

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

} /* namespace EnsembleClustering */
