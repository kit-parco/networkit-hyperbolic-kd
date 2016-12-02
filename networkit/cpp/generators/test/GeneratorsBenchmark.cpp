/*
 * GeneratorsBenchmark.cpp
 *
 *  Created on: May 29, 2013
 *      Author: forigem
 */

#ifndef NOGTEST

#include <omp.h>
#include <random>
#include <functional>
#include <cstdlib>
#include <string>

#include "GeneratorsBenchmark.h"
#include "../../auxiliary/Log.h"
#include "../../auxiliary/Parallel.h"
#include "../../auxiliary/Parallelism.h"

#include "../HyperbolicGenerator.h"
#include "../DynamicHyperbolicGenerator.h"
#include "../BarabasiAlbertGenerator.h"
#include "../ChungLuGenerator.h"
#include "../../graph/GraphBuilder.h"
#include "../../io/HyperbolicGraphReader.h"


namespace NetworKit {


TEST_F(GeneratorsBenchmark, benchmarkGraphBuilder) {
	// parameters for Erdös-Renyi
	count n = 25000;
	double p = 0.001;
	count m_expected = p * n * (n + 1) / 2;

	Graph G;
	GraphBuilder builder;

	// prepare a random generator for each possible thread
	int maxThreads = omp_get_max_threads();
	std::vector< std::function<double()> > randomPerThread;
	std::random_device device;
	std::uniform_int_distribution<uint64_t> intDist;
	for (int tid = 0; tid < maxThreads; tid++) {
		auto seed = intDist(device);
		std::mt19937_64 gen(seed);
		std::uniform_real_distribution<double> dist{0.0, std::nexttoward(1.0, 2.0)};
		auto rdn = std::bind(dist, gen);
		randomPerThread.push_back(rdn);
	}

	count m_actual;
	uint64_t t1, t2;

	// half parallel way
	m_actual = 0;
	t1 = timeOnce([&]() {
		builder = GraphBuilder(n);
		builder.parallelForNodePairs([&](node u, node v) {
			int tid = omp_get_thread_num();
			double rdn = randomPerThread[tid]();
			if (rdn <= p) {
				builder.addHalfEdge(u, v);
			}
		});
	});
	t2 = timeOnce([&]() {
		G = builder.toGraph(true);
	});
	m_actual = G.numberOfEdges();
	EXPECT_NEAR(m_actual / (double) m_expected, 1.0, 0.1);
	std::cout << "parallelForNodePairs + toGraphSequentiel:\t\t" << t1 << " + " << t2 << " = " << (t1 + t2) << " ms\n";

	// fully parallel way
	m_actual = 0;
	t1 = timeOnce([&]() {
		builder = GraphBuilder(n);
		builder.parallelForNodePairs([&](node u, node v) {
			int tid = omp_get_thread_num();
			double rdn = randomPerThread[tid]();
			if (rdn <= p) {
				builder.addHalfEdge(u, v);
			}
		});
	});
	t2 = timeOnce([&]() {
		G = builder.toGraph(true, false);
	});
	m_actual = G.numberOfEdges();
	EXPECT_NEAR(m_actual / (double) m_expected, 1.0, 0.1);
	std::cout << "parallelForNodePairs + toGraphParallel:\t\t" << t1 << " + " << t2 << " = " << (t1 + t2) << " ms\n";

	// old way
	t1 = timeOnce([&]() {
		G = Graph(n);
		G.forNodePairs([&](node u, node v) {
			if (randomPerThread[0]() <= p) {
				G.addEdge(u, v);
			}
		});
	});
	m_actual = G.numberOfEdges();
	EXPECT_NEAR(m_actual / (double) m_expected, 1.0, 0.1);
	std::cout << "forNodePairs + Graph.addEdge:\t\t\t\t" << t1 << " ms\n";
}

TEST_F(GeneratorsBenchmark, benchmarkBarabasiAlbertGenerator) {
	count k = 2;
	count nMax = 100000;
	count n0 = 2;

	BarabasiAlbertGenerator BarabasiAlbert(k, nMax, n0, false);
	Graph G(0);
	EXPECT_TRUE(G.isEmpty());

	G = BarabasiAlbert.generate();
	EXPECT_FALSE(G.isEmpty());

	EXPECT_EQ(nMax, G.numberOfNodes());
	EXPECT_EQ( ((n0-1) + ((nMax - n0) * k)), G.numberOfEdges());
}

TEST_F(GeneratorsBenchmark, benchBarabasiAlbertGeneratorBatagelj) {
	for (index i = 0; i < 10; ++i) {
		Aux::Random::setSeed(i, false);
		count n = Aux::Random::integer(100, 10000);
		count k = n / Aux::Random::integer(5, 20);
		BarabasiAlbertGenerator gen(k, n, 0);
		auto G = gen.generate();
		//EXPECT_TRUE(G.checkConsistency());
		//INFO(G.toString());
	}
}

TEST_F(GeneratorsBenchmark, benchBarabasiAlbertGenerator2) {
	for (index i = 0; i < 10; ++i) {
		Aux::Random::setSeed(i, false);
		count n = Aux::Random::integer(100, 10000);
		count k = n / Aux::Random::integer(5, 20);
		BarabasiAlbertGenerator gen(k, n, k, false);
		auto G = gen.generate();
		//EXPECT_TRUE(G.checkConsistency());
		//INFO(G.toString());
	}
}

TEST_F(GeneratorsBenchmark, benchmarkHyperbolicGenerator) {
	count n = 100000;
	HyperbolicGenerator gen(n);
	Graph G = gen.generate();
	EXPECT_EQ(G.numberOfNodes(), n);
}

TEST_F(GeneratorsBenchmark, benchmarkHyperbolicGeneratorWithSortedNodes) {
	count n = 100000;
	double s = 1.0;
	double alpha = 1.0;

	vector<double> angles(n);
	vector<double> radii(n);
	double R = s*HyperbolicSpace::hyperbolicAreaToRadius(n);
	double r = HyperbolicSpace::hyperbolicRadiusToEuclidean(R);
	//sample points randomly

	HyperbolicSpace::fillPoints(angles, radii, R, alpha);
	vector<index> permutation(n);

	index p = 0;
	std::generate(permutation.begin(), permutation.end(), [&p](){return p++;});

	//can probably be parallelized easily, but doesn't bring much benefit
	Aux::Parallel::sort(permutation.begin(), permutation.end(), [&angles,&radii](index i, index j){return angles[i] < angles[j] || (angles[i] == angles[j] && radii[i] < radii[j]);});

	vector<double> anglecopy(n);
	vector<double> radiicopy(n);

	#pragma omp parallel for
	for (index j = 0; j < n; j++) {
		anglecopy[j] = angles[permutation[j]];
		radiicopy[j] = radii[permutation[j]];
	}

	Graph G = HyperbolicGenerator().generate(anglecopy, radiicopy, r);
	EXPECT_EQ(G.numberOfNodes(), n);
}

TEST_F(GeneratorsBenchmark, benchmarkDynamicHyperbolicGeneratorOnNodeMovement) {
	const count runs = 100;
	const double fractionStep = 0.01;
	const count stepCount = 100;
	const count n = 1000000;
	const double k = 6;
	const double exp = 3;
	const double moveDistance = 0.1;

	for (index i = 0; i < stepCount; i++) {
		double moveFraction = fractionStep * i;
		for (index j = 0; j < runs; j++) {
			DynamicHyperbolicGenerator dyngen(n, k, exp, moveFraction, moveDistance);
			dyngen.generate(1);
		}
	}
}

TEST_F(GeneratorsBenchmark, benchmarkExternalEmbedderCall) {
	const count iterations = 1000;
	const count minN = 1 << 13;
	const count maxN = 1 << 23;

	const double alpha = 0.75;
	const double T = 0.1;
	const double C = -1;

	for (count n = minN; n <= maxN; n *= 2) {
		std::cout << "Started Test Case with " << n << " nodes." << std::endl;
		std::string commandstring = std::string("/home/moritzl/Gadgets/hyperbolic-embedder/embedder") + std::string(" --generate test-") + std::to_string(n)
				+ std::string(" --n ") + std::to_string(n) + std::string(" --C ") + std::to_string(C) + std::string(" --T ") + std::to_string(T) + std::string(" --alpha ") + std::to_string(alpha);

		if (std::system(NULL)) {
			int returnValue = system(commandstring.c_str());
			if (returnValue != 0) {
				DEBUG("Return value was ", returnValue);
			}
		} else {
			throw std::runtime_error("No system access!");
		}

		std::string filenamePrefix = "test-"+std::to_string(n);
		Graph G;
		vector<double> angles;
		vector<double> radii;
		HyperbolicGraphReader::readGraph(filenamePrefix, G, angles, radii);
		EXPECT_EQ(n, angles.size());
		EXPECT_EQ(n, radii.size());
		EXPECT_EQ(n, G.numberOfNodes());

		const double R = 2*log(n)+C;
		Quadtree<index, false> quad(R, true, alpha, 20, 0.999);
		for (index i = 0; i < n; i++) {
			quad.addContent(i, angles[i], radii[i]);
		}

		quad.trim();

		double beta = 1/T;
		assert(beta == beta);
		auto edgeProb = [beta, R](double distance) -> double {return 1 / (exp(beta*(distance-R)/2)+1);};

		double maxcdf = cosh(alpha*R);
		std::uniform_real_distribution<double> phidist{0, 2*M_PI};
		std::uniform_real_distribution<double> rdist{1, maxcdf};

		int seed = std::time(NULL);
		Aux::Random::setSeed(seed, false);
		std::cout << "Used seed " << seed << " for dynamic benchmarks." << std::endl;

		//now measure the dynamic part
		Aux::Timer timer;
		timer.start();

		for (index i = 0; i < iterations; i++) {
			index toMove = Aux::Random::integer(n);

			//remove old position and nodes
			quad.removeContent(toMove, angles[toMove], radii[toMove]);
			for (index neighbor : G.neighbors(toMove)) {
				G.removeEdge(toMove, neighbor);
			}

			//get new position
			angles[toMove] = phidist(Aux::Random::getURNG());
			double random = rdist(Aux::Random::getURNG());
			radii[toMove] = (acosh(random)/alpha);
			if (radii[toMove] == R) radii[toMove] = std::nextafter(radii[toMove], 0);

			//get new edges
			vector<index> newNeighbors;
			quad.getElementsProbabilistically({angles[toMove], radii[toMove]}, edgeProb, newNeighbors);
			for (index u : newNeighbors) {
				G.addEdge(u, toMove);
			}
		}
		timer.stop();
		std::cout << iterations << " iterations took " << timer.elapsedMilliseconds() << " milliseconds." << std::endl;
	}
}

TEST_F(GeneratorsBenchmark, benchmarkSingleNodeMovements) {
	const count iterations = 1000;
	const count minN = 1 << 13;
	const count maxN = 1 << 23;

	const double alpha = 0.75;
	const double T = 0.1;
	const double C = -1;

	for (count n = minN; n <= maxN; n *= 2) {
		const double R = 2*log(n)+C;
		INFO("Started Test Case with ", n, " nodes.");

		//setup coordinates and quadtree
		Aux::Timer graphTimer;
		graphTimer.start();
		vector<double> angles(n);
		vector<double> radii(n);
		HyperbolicSpace::fillPoints(angles, radii, R, alpha);
		INFO("Sampled point positions");
		Quadtree<index, false> quad(R, true, alpha, 20, 0.1);
		for (index i = 0; i < n; i++) {
			quad.addContent(i, angles[i], radii[i]);
		}

		quad.trim();
		//now define lambda
		double beta = 1/T;
		assert(beta == beta);
		auto edgeProb = [beta, R](double distance) -> double {return 1 / (exp(beta*(distance-R)/2)+1);};

		//get Graph
		GraphBuilder result(n, false, false);//no direct swap with probabilistic graphs
		count totalCandidates = 0;
		#pragma omp parallel for
		for (index i = 0; i < n; i++) {
			vector<index> near;
			totalCandidates += quad.getElementsProbabilistically({angles[i], radii[i]}, edgeProb, near);
			for (index j : near) {
				if (j >= n) ERROR("Node ", j, " prospective neighbour of ", i, " does not actually exist. Oops.");
				if (j > i) {
					result.addHalfEdge(i, j);
				}
			}
		}

		Graph G = result.toGraph(true, true);
		graphTimer.stop();
		INFO("Generated ", G.numberOfEdges(), " edges in ", graphTimer.elapsedMilliseconds(), " milliseconds, k:", 2*G.numberOfEdges() / (double)n);

		double maxcdf = cosh(alpha*R);
		std::uniform_real_distribution<double> phidist{0, 2*M_PI};
		std::uniform_real_distribution<double> rdist{1, maxcdf};

		//now measure the dynamic part
		Aux::Timer timer;
		timer.start();

		for (index i = 0; i < iterations; i++) {
			index toMove = Aux::Random::integer(n);

			//remove old position and nodes
			quad.removeContent(toMove, angles[toMove], radii[toMove]);
			for (index neighbor : G.neighbors(toMove)) {
				G.removeEdge(toMove, neighbor);
			}

			//get new position
			angles[toMove] = phidist(Aux::Random::getURNG());
			double random = rdist(Aux::Random::getURNG());
			radii[toMove] = (acosh(random)/alpha);
			if (radii[toMove] == R) radii[toMove] = std::nextafter(radii[toMove], 0);

			//get new edges
			vector<index> newNeighbors;
			quad.getElementsProbabilistically({angles[toMove], radii[toMove]}, edgeProb, newNeighbors);
			for (index u : newNeighbors) {
				G.addEdge(u, toMove);
			}
		}
		timer.stop();
		INFO(iterations, " iterations took ", timer.elapsedMilliseconds(), " milliseconds.");
	}
}

TEST_F(GeneratorsBenchmark, benchmarkSingleNodeMovementsCold) {
	const count iterations = 1000;
	const count minN = 1 << 13;
	const count maxN = 1 << 23;

	const double alpha = 0.75;
	const double C = -1;

	for (count n = minN; n <= maxN; n *= 2) {
		const double R = 2*log(n)+C;
		INFO("Started Test Case with ", n, " nodes.");

		//setup coordinates and quadtree
		Aux::Timer graphTimer;
		graphTimer.start();
		vector<double> angles(n);
		vector<double> radii(n);
		HyperbolicSpace::fillPoints(angles, radii, R, alpha);
		INFO("Sampled point positions");
		Quadtree<index, false> quad(R, false, alpha);
		for (index i = 0; i < n; i++) {
			quad.addContent(i, angles[i], radii[i]);
		}
		INFO("Filled Quadtree.");

		quad.trim();

		//get Graph
		GraphBuilder result(n, false, false);
		#pragma omp parallel for
		for (index i = 0; i < n; i++) {
			vector<index> near;

			quad.getElementsInHyperbolicCircle({angles[i], radii[i]}, R, near);
			for (index j : near) {
				if (j >= n) ERROR("Node ", j, " prospective neighbour of ", i, " does not actually exist. Oops.");
				if (j > i) {
					result.addHalfEdge(i, j);
				}
			}
		}

		Graph G = result.toGraph(true, true);
		graphTimer.stop();
		INFO("Generated ", G.numberOfEdges(), " edges in ", graphTimer.elapsedMilliseconds(), " milliseconds, k:", 2*G.numberOfEdges() / (double)n);

		double maxcdf = cosh(alpha*R);
		std::uniform_real_distribution<double> phidist{0, 2*M_PI};
		std::uniform_real_distribution<double> rdist{1, maxcdf};

		//now measure the dynamic part
		Aux::Timer timer;
		timer.start();

		for (index i = 0; i < iterations; i++) {
			index toMove = Aux::Random::integer(n);

			//remove old position and nodes
			quad.removeContent(toMove, angles[toMove], radii[toMove]);
			for (index neighbor : G.neighbors(toMove)) {
				G.removeEdge(toMove, neighbor);
			}

			//get new position
			angles[toMove] = phidist(Aux::Random::getURNG());
			double random = rdist(Aux::Random::getURNG());
			radii[toMove] = (acosh(random)/alpha);
			if (radii[toMove] == R) radii[toMove] = std::nextafter(radii[toMove], 0);

			//get new edges
			vector<index> newNeighbors;
			quad.getElementsInHyperbolicCircle({angles[i], radii[i]}, R, newNeighbors);
			for (index u : newNeighbors) {
				G.addEdge(u, toMove);
			}
		}
		timer.stop();
		INFO(iterations, " iterations took ", timer.elapsedMilliseconds(), " milliseconds.");
	}

}

TEST_F(GeneratorsBenchmark, benchmarkParallelQuadtreeConstruction) {
	count n = 33554432;
	Quadtree<index> quad(n,1.0);
	EXPECT_EQ(quad.size(), n);
}

TEST_F(GeneratorsBenchmark, benchmarkSequentialQuadtreeConstruction) {
	count n = 33554432;
	count capacity = 1000;
	double s =1;
	double alpha = 1;
	double R = s*HyperbolicSpace::hyperbolicAreaToRadius(n);
	vector<double> angles(n);
	vector<double> radii(n);
	HyperbolicSpace::fillPoints(angles, radii, s, alpha);

	Quadtree<index> quad(HyperbolicSpace::hyperbolicRadiusToEuclidean(R),false,alpha,capacity);

	for (index i = 0; i < n; i++) {
		quad.addContent(i, angles[i], radii[i]);
	}
	EXPECT_EQ(quad.size(), n);
}

TEST_F(GeneratorsBenchmark, benchmarkHyperbolicGeneratorMechanicGraphs) {
	count n = 10000;
	double k = 6;
	count m = n*k/2;
	HyperbolicGenerator gen(n, k, 3, 0.14);
	gen.setLeafCapacity(10);
	Graph G = gen.generate();
	EXPECT_NEAR(G.numberOfEdges(), m, m/10);
}

TEST_F(GeneratorsBenchmark, benchmarkChungLuGenerator) {
	count n = 100000;
    int maxDegree = 100;
	std::vector<count> vec;
	/* Creates a random weight list */
	for (int i = 0; i < n; i++){
	int grad = Aux::Random::integer(1, maxDegree);
		vec.push_back(grad);
	}
	ChungLuGenerator generator(vec);
	Graph G = generator.generate();
	EXPECT_EQ(G.numberOfNodes(), n);
}

} /* namespace NetworKit */

#endif /*NOGTEST */
