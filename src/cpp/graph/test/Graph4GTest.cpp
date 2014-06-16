/*
 * BasicGraph.cpp
 *
 *  Created on: 01.06.2014
 *      Author: Klara Reichard (klara.reichard@gmail.com), Marvin Ritter (marvin.ritter@gmail.com)
 */

#ifndef NOGTEST

#include "Graph4GTest.h"

//TODO

	/*
	template<typename L> void forNodesInRandomOrder(L handle) const;
	template<typename L> void balancedParallelForNodes(L handle) const;
	template<typename L> void parallelForNodePairs(L handle) const;
	template<typename L> void forInEdgesOf(node u, L handle) const;
	template<typename L> void forWeightedInEdgesOf(node u, L handle) const;
	template<typename L> double parallelSumForNodes(L handle) const;
	template<typename L> double parallelSumForWeightedEdges(L handle) const;
	template<typename L> void BFSEdgesfrom(node r, L handle) const;
	template<typename L> void DFSfrom(node r, L handle) const;
	template<typename L> void DFSEdgesfrom(node r, L handle) const;
	*/
namespace NetworKit {

INSTANTIATE_TEST_CASE_P(InstantiationName, Graph4GTest, testing::Values(
						std::make_tuple(false, false),
						std::make_tuple(true, false),
						std::make_tuple(false, true),
						std::make_tuple(true, true)));

bool Graph4GTest::isWeightedParameterized() const {
	return std::get<0>(GetParam());
}
bool Graph4GTest::isDirectedParameterized() const {
	return std::get<1>(GetParam());
}

Graph Graph4GTest::createParameterizedGraph(count n) const {
	bool weighted, directed;
	std::tie(weighted, directed) = GetParam();
	Graph G(n, weighted, directed);
	return G;
}

void Graph4GTest::SetUp() {
	/*
	 *    0
	 *   . \
	 *  /   \
	 * /     .
	 * 1 <-- 2
	 * ^ \  .|
	 * |  \/ |
	 * | / \ |
	 * |/   ..
	 * 3 <-- 4
	 *
	 * move you pen from node to node:
	 * 3 -> 1 -> 0 -> 2 -> 1 -> 4 -> 3 -> 2 -> 4
	 */
	n_house = 5;
	m_house = 8;

	Ghouse = createParameterizedGraph(5);
	houseEdgesOut = {
		{0, 2},
		{1, 0},
		{1, 4},
		{2, 1},
		{2, 4},
		{3, 1},
		{3, 2},
		{4, 3}
	};
	Ahouse = {n_house, std::vector<edgeweight>(n_house, 0.0)};
	edgeweight ew = 1.0;
	for (auto& e : houseEdgesOut) {
		node u = e.first;
		node v = e.second;
		Ghouse.addEdge(u, v, ew);
		
		Ahouse[u][v] = ew;
	
		if (!Ghouse.isDirected()) {
			Ahouse[v][u] = ew;
		}
		
		if (Ghouse.isWeighted()) {
			ew += 1.0;
		}

	}
}


/** CONSTRUCTORS **/

TEST_P(Graph4GTest, testCopyConstructor) {
	Graph G = Graph(this->Ghouse, false, false);
	Graph GW = Graph(this->Ghouse, true, false);
	Graph D = Graph(this->Ghouse, false, true);
	Graph DW = Graph(this->Ghouse, true, true);

	ASSERT_FALSE(G.isWeighted());
	ASSERT_FALSE(G.isDirected());
	ASSERT_EQ(this->Ghouse.numberOfNodes(), G.numberOfNodes());
	ASSERT_EQ(this->Ghouse.numberOfEdges(), G.numberOfEdges());

	ASSERT_TRUE(GW.isWeighted());
	ASSERT_FALSE(GW.isDirected());
	ASSERT_EQ(this->Ghouse.numberOfNodes(), GW.numberOfNodes());
	ASSERT_EQ(this->Ghouse.numberOfEdges(), GW.numberOfEdges());

	ASSERT_FALSE(D.isWeighted());
	ASSERT_TRUE(D.isDirected());
	ASSERT_EQ(this->Ghouse.numberOfNodes(), D.numberOfNodes());
	ASSERT_EQ(this->Ghouse.numberOfEdges(), D.numberOfEdges());

	ASSERT_TRUE(DW.isWeighted());
	ASSERT_TRUE(DW.isDirected());
	ASSERT_EQ(this->Ghouse.numberOfNodes(), DW.numberOfNodes());
	ASSERT_EQ(this->Ghouse.numberOfEdges(), DW.numberOfEdges());

	this->Ghouse.forNodes([&](node v) {
		count d = this->Ghouse.degree(v);
		count dUndirected = isDirectedParameterized() ? d + this->Ghouse.degreeIn(v) : d;
		ASSERT_EQ(dUndirected, G.degree(v));
		ASSERT_EQ(dUndirected, GW.degree(v));
		ASSERT_EQ(d, D.degree(v));
		ASSERT_EQ(d, DW.degree(v));
	});

	// if Ghouse was directed we should have an exact copy of it, but if it was undirected
	// we should have edges in both directions
	count m = 0;
	G.forEdges([&](node u, node v) {
		ASSERT_TRUE(G.hasEdge(v, u));
		ASSERT_EQ(defaultEdgeWeight, G.weight(v, u));
		ASSERT_EQ(defaultEdgeWeight, G.weight(u, v));

		auto e = std::make_pair(u, v);
		bool found = std::find(this->houseEdgesOut.begin(), this->houseEdgesOut.end(), e) != this->houseEdgesOut.end();
		if (!found) {
			e = std::make_pair(v, u);
			found = std::find(this->houseEdgesOut.begin(), this->houseEdgesOut.end(), e) != this->houseEdgesOut.end();
		}
		ASSERT_TRUE(found);
		m++;
	});
	ASSERT_EQ(8u, m);	
	
	m = 0;
	GW.forEdges([&](node u, node v) {
		ASSERT_TRUE(GW.hasEdge(v, u));
		ASSERT_EQ(GW.weight(u, v), GW.weight(v, u));

		auto e = std::make_pair(u, v);
		bool found = std::find(this->houseEdgesOut.begin(), this->houseEdgesOut.end(), e) != this->houseEdgesOut.end();
		if (!found) {
			e = std::make_pair(v, u);
			found = std::find(this->houseEdgesOut.begin(), this->houseEdgesOut.end(), e) != this->houseEdgesOut.end();
			ASSERT_EQ(this->Ghouse.weight(v, u), GW.weight(v, u));
		} else {
			ASSERT_EQ(this->Ghouse.weight(u, v), GW.weight(u, v));
		}
		ASSERT_TRUE(found);
		m++;
	});
	ASSERT_EQ(8u, m);

	m = 0;
	D.forEdges([&](node u, node v) {
		ASSERT_EQ(defaultEdgeWeight, D.weight(u, v));

		auto e = std::make_pair(u, v);
		bool found = std::find(this->houseEdgesOut.begin(), this->houseEdgesOut.end(), e) != this->houseEdgesOut.end();
		if (!this->Ghouse.isDirected()) {
			ASSERT_TRUE(D.hasEdge(v, u));

			e = std::make_pair(v, u);
			found = found || (std::find(this->houseEdgesOut.begin(), this->houseEdgesOut.end(), e) != this->houseEdgesOut.end());
		} else {
			ASSERT_FALSE(D.hasEdge(v, u));
		}
		ASSERT_TRUE(found);
		m++;
	});
	count m_expected = isDirectedParameterized() ? 8 : 16;
	ASSERT_EQ(m_expected, m);

	m = 0;
	DW.forEdges([&](node u, node v) {
		ASSERT_EQ(this->Ghouse.weight(u, v), DW.weight(u, v));

		auto e = std::make_pair(u, v);
		bool found = std::find(this->houseEdgesOut.begin(), this->houseEdgesOut.end(), e) != this->houseEdgesOut.end();
		if (!this->Ghouse.isDirected()) {
			ASSERT_TRUE(DW.hasEdge(v, u));
			e = std::make_pair(v, u);
			found = found || (std::find(this->houseEdgesOut.begin(), this->houseEdgesOut.end(), e) != this->houseEdgesOut.end());
		} else {
			ASSERT_FALSE(DW.hasEdge(v, u));
		}
		ASSERT_TRUE(found);
		m++;
	});
	m_expected = isDirectedParameterized() ? 8 : 16;
	ASSERT_EQ(m_expected, m);
}


/** GRAPH INFORMATION **/

TEST_P(Graph4GTest, testGetId) {
	Graph G1 = createParameterizedGraph();
	Graph G2 = createParameterizedGraph(5);

	ASSERT_TRUE(G1.getId() > 0);
	ASSERT_TRUE(G2.getId() > 0);	
	ASSERT_TRUE(G1.getId() < G2.getId());
}

TEST_P(Graph4GTest, testTyp) {
	Graph G = createParameterizedGraph();
	if (isDirectedParameterized()) {
		if (isWeightedParameterized()) {
			ASSERT_EQ("WeightedDirectedGraph", G.typ());
		} else {
			ASSERT_EQ("DirectedGraph", G.typ());
		}
	} else {
		if (isWeightedParameterized()) {
			ASSERT_EQ("WeightedGraph", G.typ());
		} else {
			ASSERT_EQ("Graph", G.typ());
		}
	}
}

TEST_P(Graph4GTest, tesSetName) {
	Graph G1 = createParameterizedGraph(0);
	Graph G2 = createParameterizedGraph(0);
	
	std::string s1 = "Graph 1";
	std::string s2 = "Graph 2";
	G1.setName(s1);
	G2.setName(s2);
	ASSERT_EQ(s1, G1.getName());
	ASSERT_EQ(s2, G2.getName());
}

TEST_P(Graph4GTest, testToString) {
	Graph G1 = createParameterizedGraph(0);
	Graph G2 = createParameterizedGraph(0);

	ASSERT_TRUE(G1.toString() != "");
	ASSERT_TRUE(G2.toString() != "");
}


/** NODE MODIFIERS **/

TEST_P(Graph4GTest, testAddNode) {
	Graph G = createParameterizedGraph();

	ASSERT_FALSE(G.hasNode(0));
	ASSERT_FALSE(G.hasNode(1));
	ASSERT_EQ(0u, G.numberOfNodes());

	G.addNode();
	ASSERT_TRUE(G.hasNode(0));
	ASSERT_FALSE(G.hasNode(1));
	ASSERT_EQ(1u, G.numberOfNodes());

	Graph G2 = createParameterizedGraph(2);
	ASSERT_TRUE(G2.hasNode(0));
	ASSERT_TRUE(G2.hasNode(1));
	ASSERT_FALSE(G2.hasNode(2));
	ASSERT_EQ(2u, G2.numberOfNodes());

	G2.addNode();
	G2.addNode();
	ASSERT_TRUE(G2.hasNode(2));
	ASSERT_TRUE(G2.hasNode(3));
	ASSERT_FALSE(G2.hasNode(4));
	ASSERT_EQ(4u, G2.numberOfNodes());
}

TEST_P(Graph4GTest, testRemoveNode) {
	Graph G = createParameterizedGraph(4);
	G.addEdge(0, 1);

	ASSERT_EQ(4u, G.numberOfNodes());
	ASSERT_TRUE(G.hasNode(0));
	ASSERT_TRUE(G.hasNode(1));
	ASSERT_TRUE(G.hasNode(2));
	ASSERT_TRUE(G.hasNode(3));

	EXPECT_ANY_THROW(G.removeNode(0));
	EXPECT_ANY_THROW(G.removeNode(1));

	G.removeNode(2);

	ASSERT_EQ(3u, G.numberOfNodes());
	ASSERT_TRUE(G.hasNode(0));
	ASSERT_TRUE(G.hasNode(1));
	ASSERT_FALSE(G.hasNode(2));
	ASSERT_TRUE(G.hasNode(3));

	G.removeNode(3);

	ASSERT_EQ(2u, G.numberOfNodes());
	ASSERT_TRUE(G.hasNode(0));
	ASSERT_TRUE(G.hasNode(1));
	ASSERT_FALSE(G.hasNode(2));
	ASSERT_FALSE(G.hasNode(3));
}

TEST_P(Graph4GTest, testHasNode) {
	Graph G = createParameterizedGraph(5);

	ASSERT_TRUE(G.hasNode(0));
	ASSERT_TRUE(G.hasNode(1));
	ASSERT_TRUE(G.hasNode(2));
	ASSERT_TRUE(G.hasNode(3));
	ASSERT_TRUE(G.hasNode(4));
	ASSERT_FALSE(G.hasNode(5));
	ASSERT_FALSE(G.hasNode(6));

	G.removeNode(0);
	G.removeNode(2);
	G.addNode();
	
	ASSERT_FALSE(G.hasNode(0));
	ASSERT_TRUE(G.hasNode(1));
	ASSERT_FALSE(G.hasNode(2));
	ASSERT_TRUE(G.hasNode(3));
	ASSERT_TRUE(G.hasNode(4));
	ASSERT_TRUE(G.hasNode(5));
	ASSERT_FALSE(G.hasNode(6));
}


/** NODE PROPERTIES **/

TEST_P(Graph4GTest, testDegree) {
	if (this->Ghouse.isDirected()) {
		ASSERT_EQ(1u, this->Ghouse.degree(0));
		ASSERT_EQ(2u, this->Ghouse.degree(1));
		ASSERT_EQ(2u, this->Ghouse.degree(2));
		ASSERT_EQ(2u, this->Ghouse.degree(3));
		ASSERT_EQ(1u, this->Ghouse.degree(4));
	} else {
		ASSERT_EQ(2u, this->Ghouse.degree(0));
		ASSERT_EQ(4u, this->Ghouse.degree(1));
		ASSERT_EQ(4u, this->Ghouse.degree(2));
		ASSERT_EQ(3u, this->Ghouse.degree(3));
		ASSERT_EQ(3u, this->Ghouse.degree(4));
	}
}

TEST_P(Graph4GTest, testDegreeIn) {
	if (this->Ghouse.isDirected()) {
		ASSERT_EQ(1u, this->Ghouse.degreeIn(0));
		ASSERT_EQ(2u, this->Ghouse.degreeIn(1));
		ASSERT_EQ(2u, this->Ghouse.degreeIn(2));
		ASSERT_EQ(1u, this->Ghouse.degreeIn(3));
		ASSERT_EQ(2u, this->Ghouse.degreeIn(4));
	} else {
		ASSERT_EQ(2u, this->Ghouse.degreeIn(0));
		ASSERT_EQ(4u, this->Ghouse.degreeIn(1));
		ASSERT_EQ(4u, this->Ghouse.degreeIn(2));
		ASSERT_EQ(3u, this->Ghouse.degreeIn(3));
		ASSERT_EQ(3u, this->Ghouse.degreeIn(4));
	}
}

TEST_P(Graph4GTest, testDegreeOut) {
	if (this->Ghouse.isDirected()) {
		ASSERT_EQ(1u, this->Ghouse.degreeOut(0));
		ASSERT_EQ(2u, this->Ghouse.degreeOut(1));
		ASSERT_EQ(2u, this->Ghouse.degreeOut(2));
		ASSERT_EQ(2u, this->Ghouse.degreeOut(3));
		ASSERT_EQ(1u, this->Ghouse.degreeOut(4));
	} else {
		ASSERT_EQ(2u, this->Ghouse.degreeOut(0));
		ASSERT_EQ(4u, this->Ghouse.degreeOut(1));
		ASSERT_EQ(4u, this->Ghouse.degreeOut(2));
		ASSERT_EQ(3u, this->Ghouse.degreeOut(3));
		ASSERT_EQ(3u, this->Ghouse.degreeOut(4));
	}
}

TEST_P(Graph4GTest, testIsIsolated) {
	ASSERT_FALSE(this->Ghouse.isIsolated(0));
	ASSERT_FALSE(this->Ghouse.isIsolated(1));
	ASSERT_FALSE(this->Ghouse.isIsolated(2));
	ASSERT_FALSE(this->Ghouse.isIsolated(3));
	ASSERT_FALSE(this->Ghouse.isIsolated(4));

	this->Ghouse.addNode();
	ASSERT_TRUE(this->Ghouse.isIsolated(5));

	this->Ghouse.removeEdge(1, 0);
	ASSERT_FALSE(this->Ghouse.isIsolated(0));

	this->Ghouse.removeEdge(0, 2);
	ASSERT_TRUE(this->Ghouse.isIsolated(0));

	this->Ghouse.addEdge(1, 0);
	ASSERT_FALSE(this->Ghouse.isIsolated(0));
}

TEST_P(Graph4GTest, testWeightedDegree) {
	if (this->Ghouse.isWeighted()) {
		if (this->Ghouse.isDirected()) {
			// only sum weight of outgoing edges
			ASSERT_EQ(1.0, this->Ghouse.weightedDegree(0));
			ASSERT_EQ(5.0, this->Ghouse.weightedDegree(1));
			ASSERT_EQ(9.0, this->Ghouse.weightedDegree(2));
			ASSERT_EQ(13.0, this->Ghouse.weightedDegree(3));
			ASSERT_EQ(8.0, this->Ghouse.weightedDegree(4));
		} else {
			ASSERT_EQ(3.0, this->Ghouse.weightedDegree(0));
			ASSERT_EQ(15.0, this->Ghouse.weightedDegree(1));
			ASSERT_EQ(17.0, this->Ghouse.weightedDegree(2));
			ASSERT_EQ(21.0, this->Ghouse.weightedDegree(3));
			ASSERT_EQ(16.0, this->Ghouse.weightedDegree(4));
		}
	} else {
		if (this->Ghouse.isDirected()) {
			// only count outgoing edges
			ASSERT_EQ(1 * defaultEdgeWeight, this->Ghouse.weightedDegree(0));
			ASSERT_EQ(2 * defaultEdgeWeight, this->Ghouse.weightedDegree(1));
			ASSERT_EQ(2 * defaultEdgeWeight, this->Ghouse.weightedDegree(2));
			ASSERT_EQ(2 * defaultEdgeWeight, this->Ghouse.weightedDegree(3));
			ASSERT_EQ(1 * defaultEdgeWeight, this->Ghouse.weightedDegree(4));
		} else {
			ASSERT_EQ(2 * defaultEdgeWeight, this->Ghouse.weightedDegree(0));
			ASSERT_EQ(4 * defaultEdgeWeight, this->Ghouse.weightedDegree(1));
			ASSERT_EQ(4 * defaultEdgeWeight, this->Ghouse.weightedDegree(2));
			ASSERT_EQ(3 * defaultEdgeWeight, this->Ghouse.weightedDegree(3));
			ASSERT_EQ(3 * defaultEdgeWeight, this->Ghouse.weightedDegree(4));
		}
	}
}

TEST_P(Graph4GTest, testVolume) {
	this->Ghouse.addEdge(2, 2, 0.75);

	if (this->Ghouse.isWeighted()) {
		if (this->Ghouse.isDirected()) {
			// only sum weight of outgoing edges
			ASSERT_EQ(1.0, this->Ghouse.volume(0));
			ASSERT_EQ(5.0, this->Ghouse.volume(1));
			ASSERT_EQ(10.5, this->Ghouse.volume(2));
			ASSERT_EQ(13.0, this->Ghouse.volume(3));
			ASSERT_EQ(8.0, this->Ghouse.volume(4));
		} else {
			ASSERT_EQ(3.0, this->Ghouse.volume(0));
			ASSERT_EQ(15.0, this->Ghouse.volume(1));
			ASSERT_EQ(18.5, this->Ghouse.volume(2));
			ASSERT_EQ(21.0, this->Ghouse.volume(3));
			ASSERT_EQ(16.0, this->Ghouse.volume(4));
		}
	} else {
		if (this->Ghouse.isDirected()) {
			// only count outgoing edges
			ASSERT_EQ(1 * defaultEdgeWeight, this->Ghouse.volume(0));
			ASSERT_EQ(2 * defaultEdgeWeight, this->Ghouse.volume(1));
			ASSERT_EQ(4 * defaultEdgeWeight, this->Ghouse.volume(2));
			ASSERT_EQ(2 * defaultEdgeWeight, this->Ghouse.volume(3));
			ASSERT_EQ(1 * defaultEdgeWeight, this->Ghouse.volume(4));
		} else {
			ASSERT_EQ(2 * defaultEdgeWeight, this->Ghouse.volume(0));
			ASSERT_EQ(4 * defaultEdgeWeight, this->Ghouse.volume(1));
			ASSERT_EQ(6 * defaultEdgeWeight, this->Ghouse.volume(2));
			ASSERT_EQ(3 * defaultEdgeWeight, this->Ghouse.volume(3));
			ASSERT_EQ(3 * defaultEdgeWeight, this->Ghouse.volume(4));
		}
	}
}

TEST_P(Graph4GTest, testRandomNode) {
	count n = 4;
	count samples = 100000;
	double maxAbsoluteError = 0.005;

	Graph G = createParameterizedGraph(n);
	std::vector<count> drawCounts(n, 0);
	for (count i = 0; i < samples; i++) {
		node x = G.randomNode();
		drawCounts[x]++;
	}
	for (node v = 0; v < n; v++) {
		double p = drawCounts[v] / (double) samples;
		ASSERT_NEAR(1.0 / n, p, maxAbsoluteError);
	}
}

TEST_P(Graph4GTest, testRandomNeighbor) {
	Graph G = createParameterizedGraph(10);
	G.addEdge(2, 0);
	G.addEdge(2, 1);
	G.addEdge(2, 2);
	G.addEdge(5, 6);

	ASSERT_EQ(none, G.randomNeighbor(3));
	ASSERT_EQ(6u, G.randomNeighbor(5));

	if (G.isDirected()) {
		ASSERT_EQ(none, G.randomNeighbor(1));
	} else {
		ASSERT_EQ(2u, G.randomNeighbor(1));
	}

	count nn = 3;
	count samples = 100000;
	double maxAbsoluteError = 0.005;
	std::vector<count> drawCounts(nn, 0);
	for (count i = 0; i < samples; i++) {
		node x = G.randomNeighbor(2);
		drawCounts[x]++;
	}
	for (node v = 0; v < nn; v++) {
		double p = drawCounts[v] / (double) samples;
		ASSERT_NEAR(1.0 / nn, p, maxAbsoluteError);
	}
}


/** EDGE MODIFIERS **/

TEST_P(Graph4GTest, testAddEdge) {
	Graph G = createParameterizedGraph(3);

	// Graph without edges
	ASSERT_EQ(0u, G.numberOfEdges());
	ASSERT_FALSE(G.hasEdge(0, 2));
	ASSERT_FALSE(G.hasEdge(0, 1));
	ASSERT_FALSE(G.hasEdge(1, 2));
	ASSERT_FALSE(G.hasEdge(2, 2));
	ASSERT_EQ(nullWeight, G.weight(0, 2));
	ASSERT_EQ(nullWeight, G.weight(0, 1));
	ASSERT_EQ(nullWeight, G.weight(1, 2));
	ASSERT_EQ(nullWeight, G.weight(2, 2));

	// Graph with 2 normal edges
	G.addEdge(0, 1, 4.51);
	G.addEdge(1, 2, 2.39);
	ASSERT_EQ(2u, G.numberOfEdges());
	ASSERT_FALSE(G.hasEdge(0, 2)); // was never added
	ASSERT_TRUE(G.hasEdge(0, 1));
	ASSERT_TRUE(G.hasEdge(1, 2));
	ASSERT_FALSE(G.hasEdge(2, 2)); // will be added later

	// check weights
	if (G.isWeighted()) {
		ASSERT_EQ(4.51, G.weight(0, 1));
		ASSERT_EQ(2.39, G.weight(1, 2));
	} else {
		ASSERT_EQ(defaultEdgeWeight, G.weight(0, 1));
		ASSERT_EQ(defaultEdgeWeight, G.weight(1, 2));
	}

	if (G.isDirected()) {
		ASSERT_FALSE(G.hasEdge(1, 0));
		ASSERT_FALSE(G.hasEdge(2, 1));

		// add edge in the other direction
		// note: bidirectional edges are not supported, so both edges have different weights
		G.addEdge(2, 1, 6.23);
		ASSERT_TRUE(G.hasEdge(2, 1));
		if (G.isWeighted()) {
			ASSERT_EQ(2.39, G.weight(1, 2));
			ASSERT_EQ(6.23, G.weight(2, 1));
		} else {
			ASSERT_EQ(defaultEdgeWeight, G.weight(2, 1));
		}
	} else {
		ASSERT_TRUE(G.hasEdge(1, 0));
		ASSERT_TRUE(G.hasEdge(2, 1));
		if (G.isWeighted()) {
			ASSERT_EQ(4.51, G.weight(1, 0));
			ASSERT_EQ(2.39, G.weight(2, 1));
		} else {
			ASSERT_EQ(defaultEdgeWeight, G.weight(1, 0));
			ASSERT_EQ(defaultEdgeWeight, G.weight(2, 1));
		}	
	}

	// add self loop
	G.addEdge(2, 2, 0.72);
	ASSERT_TRUE(G.hasEdge(2, 2));
	if (G.isWeighted()) {
		ASSERT_EQ(0.72, G.weight(2, 2));
	} else {
		ASSERT_EQ(defaultEdgeWeight, G.weight(2, 2));
	}
}

TEST_P(Graph4GTest, testRemoveEdge) {
	Graph G = createParameterizedGraph(3);

	G.addEdge(0, 1);
	G.addEdge(0, 0);

	ASSERT_EQ(2u, G.numberOfEdges());
	ASSERT_TRUE(G.hasEdge(0, 0));
	ASSERT_TRUE(G.hasEdge(0, 1));
	ASSERT_FALSE(G.hasEdge(2, 1));
	G.removeEdge(0, 1);

	ASSERT_EQ(1u, G.numberOfEdges());
	ASSERT_TRUE(G.hasEdge(0, 0));
	ASSERT_FALSE(G.hasEdge(0, 1));
	ASSERT_FALSE(G.hasEdge(2, 1));

	// TODO weights?
}

TEST_P(Graph4GTest, testHasEdge) {
	for (node u = 0; u < this->Ghouse.upperNodeIdBound(); u++) {
		for (node v = 0; v < this->Ghouse.upperNodeIdBound(); v++) {
			auto edge = std::make_pair(u, v);
				auto edgeReverse = std::make_pair(v, u);
			bool hasEdge = std::find(this->houseEdgesOut.begin(), this->houseEdgesOut.end(), edge) != this->houseEdgesOut.end();
			bool hasEdgeReverse = std::find(this->houseEdgesOut.begin(), this->houseEdgesOut.end(), edgeReverse) != this->houseEdgesOut.end();
			if (this->Ghouse.isDirected()) {
				ASSERT_EQ(hasEdge, this->Ghouse.hasEdge(u, v));
			} else {
				ASSERT_EQ(hasEdge || hasEdgeReverse, this->Ghouse.hasEdge(u, v));
			}
		}
	}
}

TEST_P(Graph4GTest, testRandomEdge) {
	// we only test the uniform version

	count n = 4;
	count m = 5;
	count samples = 10000;
	double maxAbsoluteError = 0.01;

	Graph G = createParameterizedGraph(n);
	G.addEdge(0, 1); // 0 * 1 = 0
	G.addEdge(1, 2); // 1 * 2 = 2
	G.addEdge(3, 2); // 3 * 2 = 1 (mod 5)
	G.addEdge(2, 2); // 2 * 2 = 4
	G.addEdge(3, 1); // 3 * 1 = 3
	ASSERT_EQ(m, G.numberOfEdges());

	std::vector<count> drawCounts(m, 0);
	for (count i = 0; i < samples; i++) {
		auto e = G.randomEdge(true);
		count id = (e.first * e.second) % 5;
		drawCounts[id]++;
	}
	for (node id = 0; id < m; id++) {
		double p = drawCounts[id] / (double) samples;
		ASSERT_NEAR(1.0 / m, p, maxAbsoluteError);
	}
}


/** GLOBAL PROPERTIES **/

TEST_P(Graph4GTest, testIsWeighted) {
	ASSERT_EQ(isWeightedParameterized(), this->Ghouse.isWeighted());
}

TEST_P(Graph4GTest, testIsDirected) {
	ASSERT_EQ(isDirectedParameterized(), this->Ghouse.isDirected());
}

TEST_P(Graph4GTest, testIsEmpty) {
	Graph G1 = createParameterizedGraph(0);
	Graph G2 = createParameterizedGraph(2);

	ASSERT_TRUE(G1.isEmpty());
	ASSERT_FALSE(G2.isEmpty());

	node v = G1.addNode();
	G2.removeNode(G2.randomNode());
	ASSERT_FALSE(G1.isEmpty());
	ASSERT_FALSE(G2.isEmpty());

	G1.removeNode(v);
	G2.removeNode(G2.randomNode());
	ASSERT_TRUE(G1.isEmpty());
	ASSERT_TRUE(G2.isEmpty());
}

TEST_P(Graph4GTest, testNumberOfNodes) {
	ASSERT_EQ(this->n_house, this->Ghouse.numberOfNodes());

	Graph G1 = createParameterizedGraph(0);
	ASSERT_EQ(0u, G1.numberOfNodes());
	G1.addNode();
	ASSERT_EQ(1u, G1.numberOfNodes());
	G1.addNode();
	ASSERT_EQ(2u, G1.numberOfNodes());
	G1.removeNode(0);
	ASSERT_EQ(1u, G1.numberOfNodes());
	G1.removeNode(1);
	ASSERT_EQ(0u, G1.numberOfNodes());
}

TEST_P(Graph4GTest, testNumberOfEdges) {
	ASSERT_EQ(this->m_house, this->Ghouse.numberOfEdges());

	Graph G1 = createParameterizedGraph(5);
	ASSERT_EQ(0u, G1.numberOfEdges());
	G1.addEdge(0, 1);
	ASSERT_EQ(1u, G1.numberOfEdges());
	G1.addEdge(1, 2);
	ASSERT_EQ(2u, G1.numberOfEdges());
	G1.removeEdge(0, 1);
	ASSERT_EQ(1u, G1.numberOfEdges());
	G1.removeEdge(1, 2);
	ASSERT_EQ(0u, G1.numberOfEdges());
}

TEST_P(Graph4GTest, testNumberOfSelfLoops) {
	Graph G = createParameterizedGraph(3);
	G.addEdge(0, 1);
	ASSERT_EQ(0u, G.numberOfSelfLoops());
	G.addEdge(0, 0);
	ASSERT_EQ(1u, G.numberOfSelfLoops());
	G.addEdge(1, 1);
	G.addEdge(1, 2);
	ASSERT_EQ(2u, G.numberOfSelfLoops());
	G.removeEdge(0, 0);
	ASSERT_EQ(1u, G.numberOfSelfLoops());
}

TEST_P(Graph4GTest, testUpperNodeIdBound) {
	ASSERT_EQ(5u, this->Ghouse.upperNodeIdBound());

	Graph G1 = createParameterizedGraph(0);
	ASSERT_EQ(0u, G1.upperNodeIdBound());
	G1.addNode();
	ASSERT_EQ(1u, G1.upperNodeIdBound());
	G1.addNode();
	ASSERT_EQ(2u, G1.upperNodeIdBound());
	G1.removeNode(1);
	ASSERT_EQ(2u, G1.upperNodeIdBound());
	G1.addNode();
	ASSERT_EQ(3u, G1.upperNodeIdBound());
}


/** DYNAMICS **/

TEST_P(Graph4GTest, testTime) {
	ASSERT_EQ(0u, this->Ghouse.time());
	this->Ghouse.timeStep();
	ASSERT_EQ(1u, this->Ghouse.time());
	this->Ghouse.timeStep();
	this->Ghouse.timeStep();
	ASSERT_EQ(3u, this->Ghouse.time());
}


/** EDGE ATTRIBUTES **/

TEST_P(Graph4GTest, testWeight) {
	this->Ghouse.forNodes([&](node u) {
		this->Ghouse.forNodes([&](node v) {
			ASSERT_EQ(this->Ahouse[u][v], this->Ghouse.weight(u, v));
		});
	});
}

TEST_P(Graph4GTest, testSetWeight) {
	Graph G = createParameterizedGraph(10);
	G.addEdge(0, 1);
	G.addEdge(1, 2);

	if (isWeightedParameterized()) {
		// edges should get weight defaultWeight on creation and setWeight should overwrite this
		G.setWeight(1, 2, 2.718);
		EXPECT_EQ(defaultEdgeWeight, G.weight(0, 1));
		EXPECT_EQ(2.718, G.weight(1, 2));
		if (isDirectedParameterized()) {
			EXPECT_EQ(nullWeight, G.weight(1, 0));
			EXPECT_EQ(nullWeight, G.weight(2, 1));
		} else {
			// undirected graph is symmetric
			EXPECT_EQ(defaultEdgeWeight, G.weight(1, 0));
			EXPECT_EQ(2.718, G.weight(2, 1));
		}

		// setting an edge weight should create the edge if it doesn't exists
		ASSERT_FALSE(G.hasEdge(5, 6));
		G.setWeight(5, 6, 56.0);
		ASSERT_EQ(56.0, G.weight(5, 6));
		ASSERT_EQ(isDirectedParameterized() ? nullWeight : 56.0, G.weight(6, 5));
		ASSERT_TRUE(G.hasEdge(5, 6));

		// directed graphs are not symmetric, undirected are
		G.setWeight(2, 1, 5.243);
		if (isDirectedParameterized()) {
			EXPECT_EQ(2.718, G.weight(1, 2));
			EXPECT_EQ(5.243, G.weight(2, 1));
		} else {
			EXPECT_EQ(5.243, G.weight(1, 2));
			EXPECT_EQ(5.243, G.weight(2, 1));
		}
		
		// self-loop
		G.addEdge(4, 4, 2.5);
		ASSERT_EQ(2.5, G.weight(4, 4));
		G.setWeight(4, 4, 3.14);
		ASSERT_EQ(3.14, G.weight(4, 4));
	} else {
		EXPECT_ANY_THROW(G.setWeight(0, 1, 1.5));
	}
}

TEST_P(Graph4GTest, increaseWeight) {
	Graph G = createParameterizedGraph(5);
	G.addEdge(0, 1);
	G.addEdge(1, 2);
	G.addEdge(3, 4, 3.14);

	if (G.isWeighted()) {
		G.increaseWeight(1, 2, 0.5);
		G.increaseWeight(3, 4, - 0.5);

		ASSERT_EQ(defaultEdgeWeight, G.weight(0, 1));
		ASSERT_EQ(defaultEdgeWeight + 0.5, G.weight(1, 2));
		ASSERT_EQ(3.14 - 0.5, G.weight(3, 4));

		if (G.isDirected()) {	
			// reverse edges do net exist => weight should be nullWeight
			ASSERT_EQ(nullWeight, G.weight(1, 0));
			ASSERT_EQ(nullWeight, G.weight(2, 1));
			ASSERT_EQ(nullWeight, G.weight(4, 3));
		} else {
			ASSERT_EQ(defaultEdgeWeight, G.weight(1, 0));
			ASSERT_EQ(defaultEdgeWeight + 0.5, G.weight(2, 1));
			ASSERT_EQ(3.14 - 0.5, G.weight(3, 4));
		}
	} else {
		EXPECT_ANY_THROW(G.increaseWeight(1, 2, 0.3));
		EXPECT_ANY_THROW(G.increaseWeight(2, 3, 0.3)); // edge does not exists
	}	
}

TEST_P(Graph4GTest, testEdgeAttributes) {
	count n = 5;
	Graph G(n);

	int attrId = G.addEdgeAttribute_double(0.0);

	G.forNodePairs([&](node u, node v){
		G.addEdge(u, v);
	});

	G.forEdges([&](node u, node v){
		EXPECT_EQ(0.0, G.attribute_double(u, v, attrId));
	});

	G.forEdges([&](node u, node v){
		G.setAttribute_double(u, v, attrId, 42.0);
	});

	G.forEdges([&](node u, node v){
		EXPECT_EQ(42.0, G.attribute_double(u, v, attrId));
	});

	node v = G.addNode();
	G.addEdge(v, 0);

	EXPECT_EQ(0.0, G.attribute_double(v, 0, attrId));

}
// done!
// int addEdgeAttribute_double(double defaultValue);

// double attribute_double(node u, node v, int attrId) const;

// void setAttribute_double(node u, node v, int attrId, double attr);


/** SUMS **/

TEST_P(Graph4GTest, testTotalEdgeWeight) {
	Graph G1 = createParameterizedGraph(5);
	Graph G2 = createParameterizedGraph(5);
	G2.addEdge(0, 1, 3.14);

	if (this->Ghouse.isWeighted()) {
		ASSERT_EQ(0.0, G1.totalEdgeWeight());
		ASSERT_EQ(3.14, G2.totalEdgeWeight());
		ASSERT_EQ(36.0, this->Ghouse.totalEdgeWeight());
	} else {
		ASSERT_EQ(0 * defaultEdgeWeight, G1.totalEdgeWeight());
		ASSERT_EQ(1 * defaultEdgeWeight, G2.totalEdgeWeight());
		ASSERT_EQ(8 * defaultEdgeWeight, this->Ghouse.totalEdgeWeight());
	}
}


/** Collections **/

TEST_P(Graph4GTest, testNodes) {
	Graph G = createParameterizedGraph(3);
	G.addNode();
	G.removeNode(2);
	G.addNode();
	auto nodes = G.nodes();

	ASSERT_EQ(4u, nodes.size());
	ASSERT_EQ(0u, nodes[0]);
	ASSERT_EQ(1u, nodes[1]);
	ASSERT_EQ(3u, nodes[2]);
	ASSERT_EQ(4u, nodes[3]);
}

TEST_P(Graph4GTest, testNeighbors) {
	Graph G = this->Ghouse;
	auto neighbors = G.neighbors(1);	
	ASSERT_EQ(0, neighbors[0]);
	ASSERT_EQ(4, neighbors[1]);



}

TEST_P(Graph4GTest, testEdges){

	Graph G = this->Ghouse;
	auto edges = G.edges();
	
	if(G.isDirected()){
	
		ASSERT_EQ(std::make_pair(0ul,2ul), edges[0]);
		ASSERT_EQ(std::make_pair(1ul,0ul), edges[1]);
		ASSERT_EQ(std::make_pair(1ul,4ul), edges[2]);
		ASSERT_EQ(std::make_pair(2ul,1ul), edges[3]);
		ASSERT_EQ(std::make_pair(2ul,4ul), edges[4]);
		ASSERT_EQ(std::make_pair(3ul,1ul), edges[5]);
		ASSERT_EQ(std::make_pair(3ul,2ul), edges[6]);
		ASSERT_EQ(std::make_pair(4ul,3ul), edges[7]);
	}else{
		ASSERT_EQ(std::make_pair(1ul,0ul), edges[0]);
		ASSERT_EQ(std::make_pair(2ul,0ul), edges[1]);
		ASSERT_EQ(std::make_pair(2ul,1ul), edges[2]);
		ASSERT_EQ(std::make_pair(3ul,1ul), edges[3]);
		ASSERT_EQ(std::make_pair(3ul,2ul), edges[4]);
		ASSERT_EQ(std::make_pair(4ul,1ul), edges[5]);
		ASSERT_EQ(std::make_pair(4ul,2ul), edges[6]);
		ASSERT_EQ(std::make_pair(4ul,3ul), edges[7]);

	}





}

/** NODE ITERATORS **/

TEST_P(Graph4GTest, testForNodes) {
	Graph G = createParameterizedGraph(3);
	std::vector<bool> visited(4, false);
	G.forNodes([&](node v) {
		ASSERT_FALSE(visited[v]);
		if (v == 2) {
			G.addNode();
		}
		visited[v] = true;
	});
	for (bool b : visited) {
		ASSERT_TRUE(b);
	}
}

//TODO
TEST_P(Graph4GTest, testParallelForNodes) {

	Graph G = this->Ghouse;
	std::vector<node> visited;
	std::vector<bool> visit(G.numberOfNodes(), false);
	G.parallelForNodes([&](node u){
		
		#pragma omp critical
		visited.push_back(u);
	});
	
	for(int i = 0; i < G.numberOfNodes(); i++)
	{
		if(!visit[visited[i]]){
			visit[visited[i]] = true;	
		}else{
			visit[visited[i]] = false;
		}
	}
	for( bool b: visit){
		ASSERT_TRUE(b);
	}

}

TEST_P(Graph4GTest, forNodesWhile) {

	count n = 100;
	Graph G(n);

	count nc = 0;
	G.forNodesWhile([&](){ return (nc < n); }, [&](node u){
		nc += 1;
	});

	ASSERT_EQ(n, nc);
}



// template<typename L> void forNodesInRandomOrder(L handle) const;

// template<typename L> void balancedParallelForNodes(L handle) const;

TEST_P(Graph4GTest, testForNodePairs) {
	count n = 10;
	count m = n * (n - 1) / 2;
	Graph G = createParameterizedGraph(n);

	// add all edges
	G.forNodePairs([&](node u, node v) {
		ASSERT_FALSE(G.hasEdge(u, v));
		G.addEdge(u, v);
		ASSERT_TRUE(G.hasEdge(u, v));
	});

	EXPECT_EQ(m, G.numberOfEdges());

	// remove all edges
	G.forNodePairs([&](node u, node v) {
		ASSERT_TRUE(G.hasEdge(u, v));
		G.removeEdge(u, v);
		ASSERT_FALSE(G.hasEdge(u, v));
	});

	EXPECT_EQ(0u, G.numberOfEdges());
}

// template<typename L> void parallelForNodePairs(L handle) const;


/** EDGE ITERATORS **/

TEST_P(Graph4GTest, testForEdges){

	Graph G = this->Ghouse;
	auto edges = G.edges();
	std::vector<std::pair<node, node> > edgescopy;
	G.forEdges([&](node u, node v){
		edgescopy.push_back(std::make_pair(u,v));
	});

	for(int i=0; i< G.numberOfEdges(); i++)
	{
		ASSERT_EQ(edges[i], edgescopy[i]);
	}
}

TEST_P(Graph4GTest, testForWeightedEdges){
	
	count n = 4;
	Graph G(n);
	G.forNodePairs([&](node u, node v){
		G.addEdge(u, v, 1.0);
	});

	edgeweight weightSum = 0.0;
	G.forWeightedEdges([&](node u, node v, node w){
		weightSum += w;
	});
	
	if(!G.isDirected()){
		ASSERT_EQ(6.0, weightSum) << "sum of edge weights should be 6 in undirected case";
	}else{
		ASSERT_EQ(12.0, weightSum) << "sum of edge weights should be 12 in directed case";
	}
}

TEST_P(Graph4GTest, testParallelForWeightedEdges){
	
	count n = 4;
	Graph G(n);
	G.forNodePairs([&](node u, node v){
		G.addEdge(u, v, 1.0);
	});

	edgeweight weightSum = 0.0;
	G.parallelForWeightedEdges([&](node u, node v, node w){
		#pragma omp atomic
		weightSum += w;
	});
	
	if(!G.isDirected()){
		ASSERT_EQ(6.0, weightSum) << "sum of edge weights should be 6 in undirected case";
	}else{
		ASSERT_EQ(12.0, weightSum) << "sum of edge weights should be 12 in directed case";
	}
}

TEST_P(Graph4GTest, testParallelForEdges){

	count n = 4;
	Graph G(n);
	G.forNodePairs([&](node u, node v){
		G.addEdge(u, v);
	});

	edgeweight weightSum = 0.0;
	G.parallelForEdges([&](node u, node v){
		#pragma omp atomic
		weightSum += 1;
	});
	
	if(!G.isDirected()){
		ASSERT_EQ(6.0, weightSum) << "sum of edge weights should be 6 in undirected case";
	}else{
		ASSERT_EQ(12.0, weightSum) << "sum of edge weights should be 12 in directed case";
	}

}

// template<typename L> void forEdgesWithAttribute_double(int attrId, L handle) const;

/** NEIGHBORHOOD ITERATORS **/


TEST_P(Graph4GTest, testForNeighborsOf){

	Graph G = this->Ghouse;
	std::vector<node> visited;
	G.forNeighborsOf(3, [&](node u){
		
		visited.push_back(u);
	});
	ASSERT_EQ(visited[0], 1);
	ASSERT_EQ(visited[1], 2);
	
}

TEST_P(Graph4GTest, testForWeightedNeighborsOf){

	Graph G = this->Ghouse;
	std::vector<std::pair<node, edgeweight> > visited;
	G.forWeightedNeighborsOf(3, [&](node u, edgeweight ew){
		
		visited.push_back(std::make_pair(u, ew));
	});
	if(G.isWeighted()&&G.isDirected()){
		ASSERT_EQ(visited[0].first, 1);
		ASSERT_EQ(visited[1].first, 2);
		ASSERT_EQ(visited[0].second, 6);
		ASSERT_EQ(visited[1].second, 7);
	} if(G.isWeighted()&& !G.isDirected()){
		ASSERT_EQ(visited[0].first, 1);
		ASSERT_EQ(visited[1].first, 2);
		ASSERT_EQ(visited[2].first, 4);
		ASSERT_EQ(visited[0].second, 6);
		ASSERT_EQ(visited[1].second, 7);
	}if(!G.isWeighted()&&G.isDirected()){
		ASSERT_EQ(visited[0].first, 1);
		ASSERT_EQ(visited[1].first, 2);
		ASSERT_EQ(visited[0].second, defaultEdgeWeight);
		ASSERT_EQ(visited[1].second, defaultEdgeWeight);
	}if(!G.isWeighted() && !G.isDirected()){
		ASSERT_EQ(visited[0].first, 1);
		ASSERT_EQ(visited[1].first, 2);
		ASSERT_EQ(visited[2].first, 4);
		ASSERT_EQ(visited[0].second, defaultEdgeWeight);
		ASSERT_EQ(visited[1].second, defaultEdgeWeight);
		
	}
	
}


TEST_P(Graph4GTest, testForEdgesOf) {
	count m = 0;
	std::vector<int> visited(this->m_house, 0);

	this->Ghouse.forNodes([&](node u) {
		this->Ghouse.forEdgesOf(u, [&](node v, node w) {
			// edges should be v to w, so if we iterate over edges from u, u should be equal v
			EXPECT_EQ(u, v);
			
			auto e = std::make_pair(v, w);
			auto it = std::find(this->houseEdgesOut.begin(), this->houseEdgesOut.end(), e);
			if (!isDirectedParameterized() && it == this->houseEdgesOut.end()) {
				auto e2 = std::make_pair(w, v);
				it = std::find(this->houseEdgesOut.begin(), this->houseEdgesOut.end(), e2);
			}

			EXPECT_TRUE(it != this->houseEdgesOut.end());
		
			// find index in edge array
			int i = std::distance(this->houseEdgesOut.begin(), it);
			if (isDirectedParameterized()) {
				// make sure edge was not visited before (would be visited twice)
				EXPECT_EQ(0, visited[i]);
			}
			
			// mark edge as visited
			visited[i]++;
			m++;
		});
	});

	if (isDirectedParameterized()) {
		// we iterated over all outgoing edges once
		EXPECT_EQ(this->m_house, m);
		for (auto c : visited) {
			EXPECT_EQ(1, c);
		}
	} else {
		// we iterated over all edges in both directions
		EXPECT_EQ(2 * this->m_house, m);
		for (auto c : visited) {
			EXPECT_EQ(2, c);
		}
	}
}

TEST_P(Graph4GTest, testForWeightedEdgesOf) {

	count m = 0;
	std::vector<int> visited(this->m_house, 0);
	double sumOfWeights = 0;

	this->Ghouse.forNodes([&](node u) {
		this->Ghouse.forWeightedEdgesOf(u, [&](node v, node w, edgeweight ew) {
			// edges should be v to w, so if we iterate over edges from u, u should be equal v
			EXPECT_EQ(u, v);
			sumOfWeights+= ew;
			auto e = std::make_pair(v, w);
			auto it = std::find(this->houseEdgesOut.begin(), this->houseEdgesOut.end(), e);
			if (!isDirectedParameterized() && it == this->houseEdgesOut.end()) {
				auto e2 = std::make_pair(w, v);
				it = std::find(this->houseEdgesOut.begin(), this->houseEdgesOut.end(), e2);
			}

			EXPECT_TRUE(it != this->houseEdgesOut.end());
		
			// find index in edge array
			int i = std::distance(this->houseEdgesOut.begin(), it);
			if (isDirectedParameterized()) {
				// make sure edge was not visited before (would be visited twice)
				EXPECT_EQ(0, visited[i]);
			}
			
			// mark edge as visited
			visited[i]++;
			m++;
		});
	});

	if (isDirectedParameterized() && !isWeightedParameterized()) {
		// we iterated over all outgoing edges once
		EXPECT_EQ(this->m_house, m);
		EXPECT_EQ(sumOfWeights, m);
		for (auto c : visited) {
			EXPECT_EQ(1, c);
		}
	} else if (isWeightedParameterized() && !isDirectedParameterized()) {
		// we iterated over all edges in both directions
		EXPECT_EQ(2 * this->m_house, m);
		EXPECT_EQ(sumOfWeights, 72);
		for (auto c : visited) {
			EXPECT_EQ(2, c);
		}
	} else if (isWeightedParameterized() && isDirectedParameterized()) {

		EXPECT_EQ(sumOfWeights, 36);
		EXPECT_EQ(this->m_house, m);
		for (auto c : visited) {
			EXPECT_EQ(1, c);
		}
	} else if(!isWeightedParameterized() && !isDirectedParameterized()) {

		EXPECT_EQ(sumOfWeights, m);
		EXPECT_EQ(2 * this->m_house, m);
		for (auto c : visited) {
			EXPECT_EQ(2, c);
		}
	}
}

TEST_P(Graph4GTest, testForInNeighborsOf) {
	std::vector<int> visited(this->n_house, 0);
	this->Ghouse.forInNeighborsOf(3, [&](node v){
		
		visited[v] = 1;
	});
	if( isDirectedParameterized()){

		EXPECT_EQ(visited[2], 0);
		EXPECT_EQ(visited[4], 1);
		EXPECT_EQ(visited[1], 0);
	}else{

		EXPECT_EQ(visited[2], 1);
		EXPECT_EQ(visited[4], 1);
		EXPECT_EQ(visited[1], 1);
	}

}


TEST_P(Graph4GTest, testForWeightedInNeighborsOf){

	std::vector<int> visited(this->n_house, 0);
	this->Ghouse.forWeightedInNeighborsOf(3,[&](node v, edgeweight ew){
		
		visited[v] = ew;
	});
	if( isDirectedParameterized()&& !isWeightedParameterized()){

		EXPECT_EQ(visited[2], 0);
		EXPECT_EQ(visited[4], defaultEdgeWeight);
		EXPECT_EQ(visited[1], 0);

	}if(!isDirectedParameterized()&& !isWeightedParameterized()){

		EXPECT_EQ(visited[2], defaultEdgeWeight);
		EXPECT_EQ(visited[4], defaultEdgeWeight);
		EXPECT_EQ(visited[1], defaultEdgeWeight);

	}if(!isDirectedParameterized()&& isWeightedParameterized()){
	
		EXPECT_EQ(visited[2], 7);
		EXPECT_EQ(visited[4], 8);
		EXPECT_EQ(visited[1], 6);
	
	}if( isDirectedParameterized()&& isWeightedParameterized()){
	
		EXPECT_EQ(visited[2], 0);
		EXPECT_EQ(visited[4], 8);
		EXPECT_EQ(visited[1], 0);
	
	}

}

TEST_P(Graph4GTest, forInEdgesOf) {

	std::vector<int> visited(this->n_house, 0);
	this->Ghouse.forInEdgesOf(3, [&](node u, node v){
		
		visited[v] = 1;
	});
	if( isDirectedParameterized()){

		EXPECT_EQ(visited[2], 0);
		EXPECT_EQ(visited[4], 1);
		EXPECT_EQ(visited[1], 0);
	}else{

		EXPECT_EQ(visited[2], 1);
		EXPECT_EQ(visited[4], 1);
		EXPECT_EQ(visited[1], 1);
	}



}

TEST_P(Graph4GTest, testForWeightedInEdgesOf){

	std::vector<int> visited(this->n_house, 0);
	this->Ghouse.forWeightedInEdgesOf(3,[&](node u, node v, edgeweight ew){
		
		visited[v] = ew;
	});
	if( isDirectedParameterized()&& !isWeightedParameterized()){

		EXPECT_EQ(visited[2], 0);
		EXPECT_EQ(visited[4], defaultEdgeWeight);
		EXPECT_EQ(visited[1], 0);

	}if(!isDirectedParameterized()&& !isWeightedParameterized()){

		EXPECT_EQ(visited[2], defaultEdgeWeight);
		EXPECT_EQ(visited[4], defaultEdgeWeight);
		EXPECT_EQ(visited[1], defaultEdgeWeight);

	}if(!isDirectedParameterized()&& isWeightedParameterized()){
	
		EXPECT_EQ(visited[2], 7);
		EXPECT_EQ(visited[4], 8);
		EXPECT_EQ(visited[1], 6);
	
	}if( isDirectedParameterized()&& isWeightedParameterized()){
	
		EXPECT_EQ(visited[2], 0);
		EXPECT_EQ(visited[4], 8);
		EXPECT_EQ(visited[1], 0);
	
	}

}


/** REDUCTION ITERATORS **/

// template<typename L> double parallelSumForNodes(L handle) const;

// template<typename L> double parallelSumForWeightedEdges(L handle) const;


/** GRAPH SEARCHES **/

TEST_P(Graph4GTest, testBFSfrom) {
		std::vector<count> visitedOrder(5, none);
		index i = 0;
		this->Ghouse.BFSfrom(3, [&](node v) {
			EXPECT_EQ(none, visitedOrder[v]); // visit every node once
			visitedOrder[v] = i++;
		});
		// have we visited all nodes
		for (count l : visitedOrder) {
			EXPECT_TRUE(l != none);
		}

		if(isDirectedParameterized()){
			// root on level 0
			EXPECT_EQ(0u, visitedOrder[3]);

			// level 1
			EXPECT_TRUE( (visitedOrder[1] == 1) ^ (visitedOrder[1] == 2) );
			EXPECT_TRUE( (visitedOrder[2] == 1) ^ (visitedOrder[2] == 2) );

			// level 2
			EXPECT_TRUE( (visitedOrder[0] == 3) ^ (visitedOrder[0] == 4) );
			EXPECT_TRUE( (visitedOrder[4] == 3) ^ (visitedOrder[4] == 4) );
		} else {
			EXPECT_EQ(0u, visitedOrder[3]); 
			EXPECT_TRUE( (visitedOrder[1] == 1) ^ (visitedOrder[1] == 2) ^ (visitedOrder[1] == 3 ));
			EXPECT_TRUE( (visitedOrder[2] == 1) ^ (visitedOrder[2] == 2) ^ (visitedOrder[2] == 3));
			EXPECT_TRUE( (visitedOrder[4] == 1) ^ (visitedOrder[4] == 2) ^ (visitedOrder[4] == 3 ));
			EXPECT_TRUE( (visitedOrder[0] == 4) );

		
		}
}

TEST_P(Graph4GTest, testDFSfrom) {
	if (isDirectedParameterized()) {
		std::vector<count> visitedOrder(5, none);
		index i = 0;
		this->Ghouse.DFSfrom(3, [&](node v) {
			EXPECT_EQ(none, visitedOrder[v]); // visit every node once
			visitedOrder[v] = i++;
		});

		// have we visited all nodes
		for (count l : visitedOrder) {
			EXPECT_TRUE(l != none);
		}

		// root on level 0
		EXPECT_EQ(0u, visitedOrder[3]);

		// level 1
		EXPECT_TRUE( (visitedOrder[1] == 1) ^ (visitedOrder[2] == 1) );

		// level 2
		EXPECT_TRUE( (visitedOrder[0] == 2) ^ (visitedOrder[1] == 2) ^ (visitedOrder[4] == 2) );

		// level 3
		EXPECT_TRUE( (visitedOrder[2] == 3) ^ (visitedOrder[0] == 3) ^ (visitedOrder[4] == 3) ^ (visitedOrder[1] == 3) );

		// level 4
		EXPECT_TRUE( (visitedOrder[2] == 4) ^ (visitedOrder[4] == 4) ^ (visitedOrder[0] == 4) );
	} else {

		std::vector<count> visitedOrder(5, none);
		Graph G(5);
		G.addEdge(0,1);
		G.addEdge(0,2);
		G.addEdge(2,3);
		G.addEdge(3,4);
		index i = 0;
		G.DFSfrom(0, [&](node v) {
			visitedOrder[v] = i++;
		});
		for (count l : visitedOrder) {
			EXPECT_TRUE(l != none);
		}
		EXPECT_EQ(0u, visitedOrder[0]);

		EXPECT_TRUE((visitedOrder[1] == 1) ^ (visitedOrder[2] == 1) );

		EXPECT_TRUE((visitedOrder[2] == 2) ^ (visitedOrder[3] == 2) );

		EXPECT_TRUE((visitedOrder[3] == 3) ^ (visitedOrder[4] == 3) );

		EXPECT_TRUE((visitedOrder[4] == 4) ^ (visitedOrder[1] == 4) );
		
		
	}
}

} /* namespace NetworKit */

#endif /*NOGTEST */
