/*
 * GeneratorsBenchmark.cpp
 *
 *  Created on: May 29, 2013
 *      Author: forigem
 */

#ifndef NOGTEST

#include "GeneratorsBenchmark.h"

namespace NetworKit {

GeneratorsBenchmark::GeneratorsBenchmark() {
	// TODO Auto-generated constructor stub

}

GeneratorsBenchmark::~GeneratorsBenchmark() {
	// TODO Auto-generated destructor stub
}

TEST_F(GeneratorsBenchmark, benchmarkBarabasiAlbertGenerator) {
	count k = 2;
	count nMax = 100000;
	count n0 = 2;

	BarabasiAlbertGenerator BarabasiAlbert(k, nMax, n0);
	Graph G(0);
	EXPECT_TRUE(G.isEmpty());

	G = BarabasiAlbert.generate();
	EXPECT_FALSE(G.isEmpty());

	EXPECT_EQ(nMax, G.numberOfNodes());
	EXPECT_EQ( ((n0-1) + ((nMax - n0) * k)), G.numberOfEdges());

}

TEST_F(GeneratorsBenchmark, benchmarkHyperbolicGenerator) {
	count n = 100000;
	HyperbolicGenerator gen(n,1);
	Graph G = gen.generate();
	EXPECT_EQ(G.numberOfNodes(), n);
}

} /* namespace NetworKit */

#endif /*NOGTEST */
