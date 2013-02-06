/*
 * Graph2GTest.cpp
 *
 *  Created on: 04.02.2013
 *      Author: cls
 */

#include "Graph2GTest.h"

namespace EnsembleClustering {

Graph2GTest::Graph2GTest() {
	// TODO Auto-generated constructor stub

}

Graph2GTest::~Graph2GTest() {
	// TODO Auto-generated destructor stub
}


TEST_F(Graph2GTest, testNumberOfNodes) {
	int64_t n = 10;
	Graph G(n);

	EXPECT_EQ(n, G.numberOfNodes()) << "number of nodes should be " << n << " in the graph";


	// TODO: insert nodes and delete them, test again
}

TEST_F(Graph2GTest, testNumberOfEdges) {
	int64_t n = 10;
	Graph G(n);

	EXPECT_EQ(0, G.numberOfEdges()) << "no edges should be in the graph";

	G.insertEdge(1, 2);
	G.insertEdge(2, 3);
	G.insertEdge(1, 3);

	EXPECT_TRUE(G.hasEdge(1, 2));
	EXPECT_TRUE(G.hasEdge(2, 3));
	EXPECT_TRUE(G.hasEdge(1, 3));
	EXPECT_EQ(3, G.numberOfEdges()) << "3 edges should have been inserted";

	G.removeEdge(3, 1);
	G.removeEdge(2, 1);
	G.removeEdge(3, 2);

	EXPECT_FALSE(G.hasEdge(1, 2));
	EXPECT_FALSE(G.hasEdge(2, 3));
	EXPECT_FALSE(G.hasEdge(1, 3));
	EXPECT_EQ(0, G.numberOfEdges()) << "all edges should have been removed";
}

TEST_F(Graph2GTest, testIsEmpty) {
	Graph G(0);
	EXPECT_EQ(0, G.numberOfNodes()) << "no nodes should be in the graph";
	EXPECT_TRUE(G.isEmpty());

	// TODO: insert node, check if isEmpty is false, delete nodes again, check if true again
}

TEST_F(Graph2GTest, testEdgeInsertionAndRemoval) {

	int64_t n = 10;
	Graph G(n);

	EXPECT_EQ(0, G.numberOfEdges()) << "no edges should be in the graph";

	G.insertEdge(1, 2);
	G.insertEdge(2, 3);
	G.insertEdge(1, 3);

	EXPECT_TRUE(G.hasEdge(1, 2));
	EXPECT_TRUE(G.hasEdge(2, 3));
	EXPECT_TRUE(G.hasEdge(1, 3));
	EXPECT_EQ(3, G.numberOfEdges()) << "3 edges should have been inserted";

	G.removeEdge(3, 1);
	G.removeEdge(2, 1);
	G.removeEdge(3, 2);

	EXPECT_FALSE(G.hasEdge(1, 2));
	EXPECT_FALSE(G.hasEdge(2, 3));
	EXPECT_FALSE(G.hasEdge(1, 3));
	EXPECT_EQ(0, G.numberOfEdges()) << "all edges should have been removed";
}

TEST_F(Graph2GTest, testNodeIteration) {
	int64_t n = 20;
	Graph G(n);

	G.forNodes([&](node v) {
		G.insertEdge(v, (v+1) % n);
	});

	EXPECT_EQ(n, G.numberOfEdges()) << n << " edges should have been inserted";

	G.forNodes([&](node v) {
		G.removeEdge(v, (v+1) % n);
	});

	EXPECT_EQ(0, G.numberOfEdges()) << "all edges should have been removed";
}


TEST_F(Graph2GTest, testParallelNodeIteration) {
	int64_t n = 500;
	int64_t offset = 100;
	Graph G(n);

	G.parallelForNodes([&](node v) {
		G.insertEdge(v, (v+offset) % n);
	});

	EXPECT_EQ(n, G.numberOfEdges()) << n << " edges should have been inserted";

	G.parallelForNodes([&](node v) {
		EXPECT_EQ(2, G.degree(v)) << "degree should be two";
	});

	G.parallelForNodes([&](node v) {
		G.removeEdge(v, (v+offset) % n);
	});

	EXPECT_EQ(0, G.numberOfEdges()) << "all edges should have been removed";
}


TEST_F(Graph2GTest, testEdgeIteration) {
	int64_t n = 500;
	int64_t offset = 100;
	Graph G(n);

	G.forNodes([&](node v) {
		G.insertEdge(v, (v+offset) % n);
	});

	EXPECT_EQ(n, G.numberOfEdges()) << n << " edges should have been inserted";

	G.forEdges([&](node u, node v) {
		G.removeEdge(u, v);
	});

	EXPECT_EQ(0, G.numberOfEdges()) << "all edges should have been removed";
}

TEST_F(Graph2GTest, testParallelEdgeIteration) {
	int64_t n = 500;
	int64_t offset = 100;
	Graph G(n);

	G.parallelForNodes([&](node v) {
		G.insertEdge(v, (v+offset) % n);
	});

	EXPECT_EQ(n, G.numberOfEdges()) << n << " edges should have been inserted";

	G.parallelForEdges([&](node u, node v) {
		G.removeEdge(u, v);
	});

	EXPECT_EQ(0, G.numberOfEdges()) << "all edges should have been removed";
}

TEST_F(Graph2GTest, testNeighborIteration) {
	int64_t n = 50;
	int64_t offset = 10;
	Graph G(n);

	G.forNodes([&](node v) {
		G.insertEdge(v, (v+offset) % n);
		G.insertEdge(v, (v+2*offset) % n);
	});

	G.forNodes([&](node v) {
		G.forNeighborsOf(v, [&](node u) {
			EXPECT_TRUE(G.hasEdge(v, u));
		});
	});
}

TEST_F(Graph2GTest, testParallelSumForNodes) {
	int64_t n = 500;
	int64_t offset = 100;
	Graph G(n);

	G.parallelForNodes([&](node v) {
		G.insertEdge(v, (v+offset) % n);
	});

	EXPECT_EQ(n, G.numberOfEdges()) << n << " edges should have been inserted";

	// sum degrees in parallel
	float sum;
	float summe = G.parallelSumForNodes([&](node v) {
		return (float) G.degree(v);
	}, sum);

	EXPECT_EQ(summe, 2 * G.numberOfEdges()) << " graph volume should be" << (2*G.numberOfEdges());
}

TEST_F(Graph2GTest, testNodePairIteration) {
	int64_t n = 50;
	Graph G(n);

	G.forNodePairs([&](node u, node v) {
		G.insertEdge(u, v);
	});

	EXPECT_EQ(n*(n-1)/2, G.numberOfEdges()) << n*(n-1)/2 << " edges should have been inserted";

	G.forNodePairs([&](node u, node v) {
		G.removeEdge(u, v);
	});

	EXPECT_EQ(0, G.numberOfEdges()) << "all edges should have been removed";
}

TEST_F(Graph2GTest, testConstNodeIteration) {
	int64_t n = 500;
	int64_t offset = 100;
	Graph G2(n);

	G2.forNodes([&](node v) {
		G2.insertEdge(v, (v+offset) % n);
	});

	const Graph G(G2);
	EXPECT_EQ(n, G.numberOfEdges()) << n << " edges should have been inserted";
	G.forNodes([&](node v) {
		EXPECT_EQ(2, G.degree(v)) << "degree should be two";
	});
}


TEST_F(Graph2GTest, testConstParallelNodeIteration) {
	int64_t n = 500;
	int64_t offset = 100;
	Graph G2(n);

	G2.parallelForNodes([&](node v) {
		G2.insertEdge(v, (v+offset) % n);
	});

	const Graph G(G2);
	EXPECT_EQ(n, G.numberOfEdges()) << n << " edges should have been inserted";
	G.parallelForNodes([&](node v) {
		EXPECT_EQ(2, G.degree(v)) << "degree should be two";
	});
}


TEST_F(Graph2GTest, testConstEdgeIteration) {
	int64_t n = 500;
	int64_t offset = 100;
	Graph G(n);

	G.forNodes([&](node v) {
		G.insertEdge(v, (v+offset) % n);
	});

	EXPECT_EQ(n, G.numberOfEdges()) << n << " edges should have been inserted";

	const Graph G2(G);
	G2.forEdges([&](node u, node v) {
		EXPECT_TRUE(G2.hasEdge(v, u));
	});
}

TEST_F(Graph2GTest, testConstParallelEdgeIteration) {
	int64_t n = 500;
	int64_t offset = 100;
	Graph G(n);

	G.parallelForNodes([&](node v) {
		G.insertEdge(v, (v+offset) % n);
	});

	EXPECT_EQ(n, G.numberOfEdges()) << n << " edges should have been inserted";

	const Graph G2(G);
	G2.parallelForEdges([&](node u, node v) {
		EXPECT_TRUE(G2.hasEdge(v, u));
	});
}

TEST_F(Graph2GTest, testConstNeighborIteration) {
	int64_t n = 50;
	int64_t offset = 10;
	Graph G(n);

	G.forNodes([&](node v) {
		G.insertEdge(v, (v+offset) % n);
		G.insertEdge(v, (v+2*offset) % n);
	});

	const Graph G2(G);
	G2.forNodes([&](node v) {
		G2.forNeighborsOf(v, [&](node u) {
			EXPECT_TRUE(G.hasEdge(v, u));
		});
	});
}

TEST_F(Graph2GTest, testConstParallelSumForNodes) {
	int64_t n = 500;
	int64_t offset = 100;
	Graph G(n);

	G.parallelForNodes([&](node v) {
		G.insertEdge(v, (v+offset) % n);
	});

	EXPECT_EQ(n, G.numberOfEdges()) << n << " edges should have been inserted";

	const Graph G2(G);
	// sum degrees in parallel
	float sum;
	float summe = G2.parallelSumForNodes([&](node v) {
		return (float) G.degree(v);
	}, sum);

	EXPECT_EQ(summe, 2 * G2.numberOfEdges()) << " graph volume should be" << (2*G2.numberOfEdges());
}

TEST_F(Graph2GTest, testConstNodePairIteration) {
	int64_t n = 50;
	Graph G(n);

	G.forNodePairs([&](node u, node v) {
		G.insertEdge(u, v);
	});

	EXPECT_EQ(n*(n-1)/2, G.numberOfEdges()) << n*(n-1)/2 << " edges should have been inserted";

	const Graph G2(G);
	G2.forNodePairs([&](node u, node v) {
		EXPECT_TRUE(G2.hasEdge(u, v));
	});
}


TEST_F(Graph2GTest, testAddNode) {
	count n = 10;
	count newN = n + 1;
	Graph G(n);

	node v = G.addNode();
	EXPECT_EQ(v, n) << "new node should have id equal to previous n";
	EXPECT_EQ(newN, G.numberOfNodes()) << "new number of nodes should be " << newN;
	EXPECT_EQ(0, G.degree(v)) << "new node should have degree 0";


	// TODO: insert nodes and delete them, test again
}

} /* namespace EnsembleClustering */
