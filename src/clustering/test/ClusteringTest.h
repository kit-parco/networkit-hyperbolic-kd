/*
 * ClusteringTest.h
 *
 *  Created on: 12.12.2012
 *      Author: cls
 */

#ifndef CLUSTERINGTEST_H_
#define CLUSTERINGTEST_H_

#include <gtest/gtest.h>


#include "../../aux/log.h"
#include "../Clustering.h"
#include "../Modularity.h"
#include "../ClusteringGenerator.h"
#include "../../graph/GraphGenerator.h"
#include "../LabelPropagation.h"

namespace EnsembleClustering {

class ClusteringTest: public testing::Test {


};



TEST_F(ClusteringTest, testModularity) {
	GraphGenerator graphGenerator;

	int n = 100;

	DEBUG("testing modularity on clustering of complete graph with " << n << " nodes");


	Graph G = graphGenerator.makeCompleteGraph(n);

	ClusteringGenerator clusteringGenerator;

	Clustering singleton = clusteringGenerator.makeSingletonClustering(G);
	Clustering one = clusteringGenerator.makeOneClustering(G);

	Modularity modularity(G);

	double modSingleton = modularity.getQuality(singleton);
	double modOne = modularity.getQuality(one);

	DEBUG("mod(singleton-clustering) = " << modSingleton);
	DEBUG("mod(1-clustering) = " << modOne);


	EXPECT_EQ(0.0, modOne) << "1-clustering should have modularity of 0.0";
	EXPECT_GE(0.0, modSingleton) << "singleton clustering should have modularity less than 0.0";
	EXPECT_NE(-1* std::numeric_limits<double>::infinity(), modSingleton) << "modularity of singleton clustering became negative infinity - why?";

}


TEST_F(ClusteringTest, testLabelPropagation) {
	GraphGenerator graphGenerator;
	int n = 10;
	Graph G = graphGenerator.makeCompleteGraph(n);

	LabelPropagation lp;
	Clustering zeta = lp.run(G);

	EXPECT_TRUE(zeta.isProper(G)) << "the resulting partition should be a proper clustering";
}

} /* namespace EnsembleClustering */
#endif /* CLUSTERINGTEST_H_ */
