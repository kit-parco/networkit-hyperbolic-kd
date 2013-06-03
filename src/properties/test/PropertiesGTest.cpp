/*
 * PropertiesGTest.cpp
 *
 *  Created on: 03.06.2013
 *      Author: cls
 */

#include "PropertiesGTest.h"

namespace NetworKit {

PropertiesGTest::PropertiesGTest() {
	// TODO Auto-generated constructor stub

}

PropertiesGTest::~PropertiesGTest() {
	// TODO Auto-generated destructor stub
}


TEST_F(PropertiesGTest, testClusteringCoefficient) {

	GraphGenerator gen;
	Graph G = gen.makeErdosRenyiGraph(100, 1.0);

	ClusteringCoefficient clusteringCoefficient;
	double cc = clusteringCoefficient.calculate(G);

	EXPECT_EQ(1.0, cc);

}

TEST_F(PropertiesGTest, testDegreeDistribution) {
	Graph G(3);
	std::vector<count> degreeDist = GraphProperties::degreeDistribution(G);
	G.addEdge(0,1);
	G.addEdge(1,2);
	G.addEdge(2,0);
	degreeDist = GraphProperties::degreeDistribution(G);
	EXPECT_EQ(degreeDist[0], 0);
	EXPECT_EQ(degreeDist[1], 0);
	EXPECT_EQ(degreeDist[2], 3);


}



} /* namespace NetworKit */
