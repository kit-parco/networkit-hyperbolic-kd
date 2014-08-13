/*
Dy * GeneratorsTest.cpp
 *
 *  Created on: 09.04.2013
 *      Author: cls
 */

#ifndef NOGTEST

#include <numeric>

#include "GeneratorsGTest.h"

#include "../DynamicPathGenerator.h"
#include "../DynamicForestFireGenerator.h"
#include "../DynamicDorogovtsevMendesGenerator.h"
#include "../DorogovtsevMendesGenerator.h"
#include "../WattsStrogatzGenerator.h"
#include "../RegularRingLatticeGenerator.h"
#include "../../properties/ClusteringCoefficient.h"
#include "../../community/PLM.h"
#include "../../community/Modularity.h"
#include "../StochasticBlockmodel.h"


namespace NetworKit {

GeneratorsGTest::GeneratorsGTest() {

}

GeneratorsGTest::~GeneratorsGTest() {

}


TEST_F(GeneratorsGTest, testDynamicBarabasiAlbertGeneratorSingleStep) {
	count k = 2; // number of edges added per node
	DynamicGraphSource* gen = new DynamicBarabasiAlbertGenerator(k);
	GraphEventProxy* Gproxy = gen->newGraph();
	Graph* G = Gproxy->G;

	gen->initializeGraph();

	count nPre = G->numberOfNodes();
	count mPre = G->numberOfEdges();
	EXPECT_EQ(k, nPre) << "graph should have been initialized to k nodes";
	EXPECT_EQ(k - 1, mPre) << "graph should have been initialized to a path of k nodes which means k-1 edges";

	// perform single preferential attachment step
	gen->generate();

	count nPost = G->numberOfNodes();
	count mPost = G->numberOfEdges();
	EXPECT_EQ(nPre + 1, nPost) << "one more node should have been added";
	EXPECT_EQ(mPre + k, mPost) << "k edges should have been added";

	delete G;
}

TEST_F(GeneratorsGTest, testDynamicBarabasiAlbertGenerator) {


	DynamicGraphSource* gen = new DynamicBarabasiAlbertGenerator(2);

	GraphEventProxy* Gproxy = gen->newGraph();
	Graph* G = Gproxy->G;

	gen->initializeGraph();

	EXPECT_EQ(2u, G->numberOfNodes()) << "initially the generator creates two connected nodes";
	EXPECT_EQ(1u, G->numberOfEdges()) << "initially the generator creates two connected nodes";

	count n = 100;

	gen->generateWhile([&]() {
				return ( G->numberOfNodes() < n );
			});

	EXPECT_EQ(n, G->numberOfNodes());
	DEBUG("m = " , G->numberOfEdges());

	// resume generator

	gen->generateWhile([&]() {
		return (G->numberOfNodes() < 2 * n);
	});
	EXPECT_EQ(2 * n, G->numberOfNodes());
}


TEST_F(GeneratorsGTest, viewDynamicBarabasiAlbertGenerator) {
	DynamicGraphSource* gen = new DynamicBarabasiAlbertGenerator(2);
	GraphEventProxy* Gproxy = gen->newGraph();
	Graph* G = Gproxy->G;
	gen->initializeGraph();
	count n = 42;
	gen->generateWhile([&]() {
				return ( G->numberOfNodes() < n );
			});
	METISGraphWriter writer;
	writer.write(*G, "output/BATest.graph");

	delete G;
}


TEST_F(GeneratorsGTest, testStaticPubWebGenerator) {
	count n = 1800;
	count numCluster = 24;
	count maxNumNeighbors = 36;
	float rad = 0.075;

	PubWebGenerator gen(n, numCluster, rad, maxNumNeighbors);
	Graph G = gen.generate();
	EXPECT_EQ(n, G.numberOfNodes()) << "number of generated nodes";

	// check degree
	G.forNodes([&](node v) {
		EXPECT_LE(G.degree(v), maxNumNeighbors) << "maximum degree";
	});

	// 1-clustering
	ClusteringGenerator clusterGen;
	Partition oneClustering = clusterGen.makeOneClustering(G);
	EXPECT_EQ(G.numberOfNodes(),oneClustering.numberOfElements());

	// output to EPS file
	PostscriptWriter psWriter(true);
	psWriter.write(G, oneClustering, "output/pubweb.eps");

	// clustering
	PLM clusterAlgo;
	Partition clustering = clusterAlgo.run(G);
	EXPECT_EQ(G.numberOfNodes(),clustering.numberOfElements());
	psWriter.write(G, clustering, "output/pubweb-clustered-PLM.eps");

	Modularity mod;
	double modVal = mod.getQuality(clustering, G);
	EXPECT_GE(modVal, 0.2) << "modularity of clustering";
	DEBUG("Modularity of clustering: " , modVal);
	DEBUG("Total edge weight: " , G.totalEdgeWeight());
}


TEST_F(GeneratorsGTest, testDynamicPubWebGenerator) {
//	count nSteps = 100;
//	count n = 1200;
	count nSteps = 15;
	count n = 300;
	count numCluster = 30;
	count maxNumNeighbors = 40;
	float rad = 0.08;

	DynamicPubWebGenerator dynGen(n, numCluster, rad, maxNumNeighbors, false);
	Graph G = dynGen.getGraph();
	GraphUpdater gu(G);
	std::vector<GraphEvent> stream;

	// static clustering algorithm for better visual output
	PostscriptWriter psWriter(true);
	psWriter.write(G, "output/pubweb-0000.eps");

	for (index i = 1; i <= nSteps; ++i) {
		stream = dynGen.generate(1);
		DEBUG("updating graph");
		gu.update(stream);
		G.initCoordinates();

		DEBUG("updated graph, new (n, m) = (" , G.numberOfNodes() , ", " , G.numberOfEdges() , ")");
		edgeweight tew = G.totalEdgeWeight();
		DEBUG("1/2 graph volume: ", tew);
		EXPECT_GT(tew, 0);

		// update coordinates
		std::map<node, Point<float> > newCoordinates = dynGen.getNewCoordinates();
		for (std::map<node, Point<float> >::iterator iter = newCoordinates.begin();
				iter != newCoordinates.end(); ++iter) {
			node v = iter->first;
			Point<float> p = iter->second;
			G.setCoordinate(v, p);
		}

		// output for visual inspection
		char path[23];
		sprintf(path, "output/pubweb-%04llu.eps", static_cast<unsigned long long>(i));
		TRACE("path: " , path);
		psWriter.write(G, path);
	}
}


TEST_F(GeneratorsGTest, testBarabasiAlbertGenerator) {
	count k = 3;
	count nMax = 100;
	count n0 = 3;

	BarabasiAlbertGenerator BarabasiAlbert(k, nMax, n0);
	Graph G(0);
	EXPECT_TRUE(G.isEmpty());

	G = BarabasiAlbert.generate();
	EXPECT_FALSE(G.isEmpty());

	EXPECT_EQ(nMax, G.numberOfNodes());
	EXPECT_EQ( ((n0-1) + ((nMax - n0) * k)), G.numberOfEdges());



}

TEST_F(GeneratorsGTest, generatetBarabasiAlbertGeneratorGraph) {
		count k = 3;
		count nMax = 1000;
		count n0 = 3;

		BarabasiAlbertGenerator BarabasiAlbert(k, nMax, n0);

		Graph G = BarabasiAlbert.generate();
		GraphIO io;
		io.writeAdjacencyList(G, "output/"
				"BarabasiGraph.txt");
}


TEST_F(GeneratorsGTest, testDynamicPathGenerator) {
	DynamicPathGenerator gen;
	auto stream = gen.generate(42);
	for (auto ev : stream) {
		TRACE(ev.toString());
	}
}


TEST_F(GeneratorsGTest, testErdosRenyiGenerator) {
	count n = 2000;
	double p = 1.5 * (log(n) / (double) n);

	ErdosRenyiGenerator generator(n, p);
	Graph G = generator.generate();
	EXPECT_EQ(n, G.numberOfNodes());
	EXPECT_FALSE(G.isEmpty());
	EXPECT_TRUE(G.consistencyCheck());

	count nPairs = (n * (n-1)) / 2;
	count nEdges = G.numberOfEdges();
	EXPECT_GE(nEdges, 0.75 * p * nPairs);
	EXPECT_LE(nEdges, 1.25 * p * nPairs);

	DEBUG("Number of edges with probability " , p , " (actual/expected): " , nEdges , " / " , (nPairs * p));
}

TEST_F(GeneratorsGTest, testRmatGenerator) {
	count scale = 9;
	count n = (1 << scale);
	count edgeFactor = 12;
	double a = 0.51;
	double b = 0.12;
	double c = 0.12;
	double d = 0.2;

	RmatGenerator rmat(scale, edgeFactor, a, b, c, d);
	Graph G = rmat.generate();

	EXPECT_EQ(G.numberOfNodes(), n);
	EXPECT_LE(G.numberOfEdges(), n * edgeFactor);

	ClusteringCoefficient cc;
	double ccex = cc.exactGlobal(G);
	EXPECT_LE(ccex, 0.4);

	PLM clusterer(true);
	Partition zeta = clusterer.run(G);
	Modularity mod;
	double modVal = mod.getQuality(zeta, G);
	INFO("Modularity of R-MAT graph clustering: ", modVal);
	EXPECT_GE(modVal, 0.0);
}


TEST_F(GeneratorsGTest, testChungLuGenerator) {
	count n = 400;
	count maxDegree = n / 8;
	std::vector<unsigned int> sequence(n); // TODO: revert to count when cython issue fixed
	count expVolume = 0;
	count actualVolume = 0;

	// fill sequence with random values (this is not power-law, of course!)
	for (index i = 0; i < n; ++i) {
		sequence[i] = rand() % maxDegree;
		expVolume += sequence[i];
	}

	ChungLuGenerator gen(sequence);
	Graph G = gen.generate();

	EXPECT_EQ(n, G.numberOfNodes());
	G.forNodes([&](node v) {
		actualVolume += G.degree(v);
	});

	INFO("expected volume: ", expVolume, ", actual volume: ", actualVolume);
}


TEST_F(GeneratorsGTest, testHavelHakimiGeneratorOnRandomSequence) {
	count n = 400;
	count maxDegree = n / 10;
	std::vector<unsigned int> sequence(n); // TODO: revert to count when cython issue fixed
//	std::vector<count> sequence = {5, 4, 4, 3, 2, 2, 2, 2, 2, 2};
	bool realizable = false;

	do {
		// fill sequence with random values (this is not power-law, of course!)
		for (index i = 0; i < n; ++i) {
			sequence[i] = rand() % maxDegree;
		}

		// check if sequence is realizable
		bool skipTest = false;
		HavelHakimiGenerator hhgen(sequence, skipTest);
		realizable = hhgen.getRealizable();

		if (realizable) {
			Graph G = hhgen.generate();
			count volume = std::accumulate(sequence.begin(), sequence.end(), 0);
			EXPECT_EQ(volume, 2 * G.numberOfEdges());
		}
	} while (! realizable);
}

TEST_F(GeneratorsGTest, testHavelHakimiGeneratorOnRealSequence) {
	METISGraphReader reader;
	std::vector<std::string> graphs = {"input/jazz.graph",
			"input/lesmis.graph"}; //, "input/PGPgiantcompo.graph", "input/coAuthorsDBLP.graph"};

	for (auto path : graphs) {
		Graph G = reader.read(path);
		count n = G.numberOfNodes();
		std::vector<unsigned int> sequence = GraphProperties::degreeSequence(G); // TODO: revert to count when cython issue fixed

		bool skipTest = false;
		HavelHakimiGenerator hhgen(sequence, skipTest);
		Graph G2 = hhgen.generate();

		count volume = std::accumulate(sequence.begin(), sequence.end(), 0);
		EXPECT_EQ(volume, 2 * G2.numberOfEdges());

		if (volume < 50000) {
			std::vector<unsigned int> testSequence = GraphProperties::degreeSequence(G2);
			std::sort(testSequence.begin(), testSequence.end(), std::greater<unsigned int>());
			std::sort(sequence.begin(), sequence.end(), std::greater<unsigned int>());

			for (index i = 0; i < n; ++i) {
				EXPECT_EQ(sequence[i], testSequence[i]);
			}
		}
	}
}


TEST_F(GeneratorsGTest, testDynamicForestFireGenerator) {
	Graph G1(0);
	GraphUpdater gu1(G1);
	std::vector<GraphEvent> stream;
	DynamicForestFireGenerator ffg1(0.0, false);
	stream = ffg1.generate(10);
	gu1.update(stream);
	EXPECT_EQ(11, G1.numberOfNodes());
	G1.forNodes([&](node u) {
		int count = 0;
		G1.forNeighborsOf(u, [&](node v) {
			if (v < u) {
				count += 1;
			}
		});
		if (u == 0) {
			EXPECT_EQ(0, count);
		} else {
			EXPECT_EQ(1, count);
		}
	});

	Graph G2(0);
	GraphUpdater gu2(G2);
	DynamicForestFireGenerator ffg2(1.0, true, 1.0);
	stream = ffg2.generate(10);
	gu2.update(stream);
	EXPECT_EQ(11, G2.numberOfNodes());
	G2.forNodePairs([&](node u, node v) {
		if (v < u) {
			EXPECT_TRUE(G2.hasEdge(u,v));
		}
	});
}

TEST_F(GeneratorsGTest, testRegularRingLatticeGenerator) {
	int n0 = 10;
	int neighbors = 2;
	auto testRingLattice = [&](Graph G) {
		EXPECT_EQ(n0, G.numberOfNodes());
		EXPECT_EQ(n0*neighbors, G.numberOfEdges());
		G.forNodePairs([&](node u, node v) {
			int diff = std::abs((int) u- (int) v);
			if (u != v && (diff <= neighbors || diff >= n0 - neighbors)) {
				EXPECT_TRUE(G.hasEdge(u,v));
			} else {
				EXPECT_FALSE(G.hasEdge(u,v));
			}
		});
	};

	RegularRingLatticeGenerator rrlg = RegularRingLatticeGenerator(n0, neighbors);
	testRingLattice(rrlg.generate());
}

TEST_F(GeneratorsGTest, testWattsStrogatzGenerator) {
	int n0 = 10;
	int neighbors = 2;
	auto testRingLattice = [&](Graph G) {
		G.forNodePairs([&](node u, node v) {
			int diff = std::abs((int) u- (int) v);
			if (u != v && (diff <= neighbors || diff >= n0 - neighbors)) {
				EXPECT_TRUE(G.hasEdge(u,v));
			} else {
				EXPECT_FALSE(G.hasEdge(u,v));
			}
		});
	};

	WattsStrogatzGenerator wsg1 = WattsStrogatzGenerator(n0, neighbors, 0.0);
	testRingLattice(wsg1.generate());

	WattsStrogatzGenerator wsg2 = WattsStrogatzGenerator(n0, neighbors, 0.3);
	Graph G = wsg2.generate();
	EXPECT_EQ(n0, G.numberOfNodes());
	EXPECT_EQ(n0*neighbors, G.numberOfEdges());
}

TEST_F(GeneratorsGTest, testDorogovtsevMendesGenerator) {
	int n0 = 20;
	DorogovtsevMendesGenerator dmg = DorogovtsevMendesGenerator(n0);
	Graph G = dmg.generate();

	EXPECT_EQ(n0, G.numberOfNodes());
	EXPECT_EQ(2*n0-3, G.numberOfEdges());
	G.forNodes([&](node u) {
		int count = 0;
		G.forNeighborsOf(u, [&](node v) {
			if (v < u) {
				count += 1;
			}
		});
		if (u <= 2) {
			EXPECT_EQ(u, count);
		} else {
			EXPECT_EQ(2, count);
		}
	});
}

TEST_F(GeneratorsGTest, testDynamicDorogovtsevMendesGenerator) {
	count n0 = 20;
	DynamicDorogovtsevMendesGenerator ddmg = DynamicDorogovtsevMendesGenerator();
	Graph G(0);
	GraphUpdater gu(G);
	std::vector<GraphEvent> stream;
	stream = ddmg.generate(n0 - 3);
	gu.update(stream);

	EXPECT_EQ(n0, G.numberOfNodes());
	EXPECT_EQ(2*n0-3, G.numberOfEdges());
	G.forNodes([&](node u) {
		int count = 0;
		G.forNeighborsOf(u, [&](node v) {
			if (v < u) {
				count += 1;
			}
		});
		if (u <= 2) {
			EXPECT_EQ(u, count);
		} else {
			EXPECT_EQ(2, count);
		}
	});
}



TEST_F(GeneratorsGTest, testStochasticBlockmodel) {
	count n = 10;
	count nBlocks = 2;
	std::vector<index> membership = {0, 0, 0, 0, 0, 1, 1, 1, 1, 1};
	std::vector<std::vector<double> > affinity = {{1.0, 0.0}, {0.0, 1.0}};
	StochasticBlockmodel sbm(n, nBlocks, membership, affinity);
	Graph G = sbm.generate();

	EXPECT_EQ(n, G.numberOfNodes());
	EXPECT_EQ(20, G.numberOfEdges());
}

} /* namespace NetworKit */

#endif /*NOGTEST */
