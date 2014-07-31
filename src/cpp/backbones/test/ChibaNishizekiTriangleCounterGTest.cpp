/*
 * ChibaNishizekiCounterGTest.cpp
 *
 *  Created on: 23.05.2014
 *      Author: Gerd Lindner
 */

#ifndef NOGTEST

#include "ChibaNishizekiTriangleCounterGTest.h"

#include "../../backbones/ChibaNishizekiTriangleCounter.h"

namespace NetworKit {

TEST_F(ChibaNishizekiTriangleCounterGTest, testTriangleCountsTrivial) {
	Graph g(5);

	g.addEdge(0,1);
	g.addEdge(0,2);
	g.addEdge(1,2);

	ChibaNishizekiTriangleCounter counter;
	EdgeAttribute counts = counter.getAttribute(g, EdgeAttribute());

	EXPECT_DOUBLE_EQ(1.0, (counts[uEdge(0,1)])) << "wrong triangle count";
	EXPECT_EQ(1.0, (counts[uEdge(0,2)])) << "wrong triangle count";
	EXPECT_EQ(1.0, (counts[uEdge(1,2)])) << "wrong triangle count";
	EXPECT_EQ(0.0, (counts[uEdge(2,3)])) << "wrong triangle count";
}

TEST_F(ChibaNishizekiTriangleCounterGTest, testTriangleCountsSimple) {
	int64_t n = 6;
	Graph g(n);

	g.addEdge(0,1);
	g.addEdge(0,2);
	g.addEdge(1,2);

	g.addEdge(0,4);
	g.addEdge(0,3);
	g.addEdge(3,4);

	g.addEdge(0,5);
	g.addEdge(4,5);

	EXPECT_EQ(8, g.numberOfEdges()) << "wrong edge count";

	ChibaNishizekiTriangleCounter counter;
	EdgeAttribute counts = counter.getAttribute(g, EdgeAttribute());

	EXPECT_EQ(6, g.numberOfNodes()) << "undesired side effect";
	EXPECT_EQ(8, g.numberOfEdges()) << "undesired side effect";

	EXPECT_EQ(8, counts.size()) << "wrong triangle count map size";
	EXPECT_DOUBLE_EQ(1.0, (counts[g.edgeId(0,1)])) << "wrong triangle count";
	EXPECT_DOUBLE_EQ(1.0, (counts[g.edgeId(0,2)])) << "wrong triangle count";
	EXPECT_DOUBLE_EQ(1.0, (counts[g.edgeId(1,2)])) << "wrong triangle count";

	EXPECT_DOUBLE_EQ(1.0, (counts[g.edgeId(0,3)])) << "wrong triangle count";
	EXPECT_DOUBLE_EQ(1.0, (counts[g.edgeId(3,4)])) << "wrong triangle count";

	EXPECT_DOUBLE_EQ(2.0, (counts[g.edgeId(0,4)])) << "wrong triangle count";
	EXPECT_DOUBLE_EQ(1.0, (counts[g.edgeId(0,5)])) << "wrong triangle count";
	EXPECT_DOUBLE_EQ(1.0, (counts[g.edgeId(4,5)])) << "wrong triangle count";

	EXPECT_DOUBLE_EQ(1.0, (counts[g.edgeId(1,0)])) << "wrong triangle count";
	EXPECT_DOUBLE_EQ(1.0, (counts[g.edgeId(2,0)])) << "wrong triangle count";
	EXPECT_DOUBLE_EQ(1.0, (counts[g.edgeId(2,1)])) << "wrong triangle count";

	EXPECT_DOUBLE_EQ(1.0, (counts[g.edgeId(3,0)])) << "wrong triangle count";
	EXPECT_DOUBLE_EQ(1.0, (counts[g.edgeId(4,3)])) << "wrong triangle count";

	EXPECT_DOUBLE_EQ(2.0, (counts[g.edgeId(4,0)])) << "wrong triangle count";
	EXPECT_DOUBLE_EQ(1.0, (counts[g.edgeId(5,0)])) << "wrong triangle count";
	EXPECT_DOUBLE_EQ(1.0, (counts[g.edgeId(5,4)])) << "wrong triangle count";
}

}
/* namespace NetworKit */

#endif /*NOGTEST */
