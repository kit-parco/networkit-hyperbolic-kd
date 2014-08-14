/*
 * GraphBuilderGTest.h
 *
 *  Created on: 14.08.2014
 *      Author: Marvin Ritter (marvin.ritter@gmail.com)
 */

#ifndef NOGTEST

#ifndef GRAPHBUILDERGTEST_H_
#define GRAPHBUILDERGTEST_H_

#include <tuple>
#include <gtest/gtest.h>

#include "../Graph.h"
#include "../GraphBuilder.h"

namespace NetworKit {

class GraphBuilderGTest: public testing::TestWithParam< std::tuple<bool, bool, bool> > {
public:
	virtual void SetUp();

protected:
	GraphBuilder bHouse;
	Graph Ghouse;
	std::vector< std::pair<node, node> > houseEdgesOut;
	std::vector< std::vector<edgeweight> > Ahouse;
	count n_house;
	count m_house;

	bool isGraph() const { return !isWeighted() && !isDirected(); }
	bool isWeightedGraph() const { return isWeighted() && !isDirected(); }
	bool isDirectedGraph() const { return !isWeighted() && isDirected(); }
	bool isWeightedDirectedGraph() const { return isWeighted() && isDirected(); }

	bool isWeighted() const;
	bool isDirected() const;
	bool useParallel() const;
	GraphBuilder createGraphBuilder(count n = 0) const;
};

} /* namespace NetworKit */

#endif /* GRAPHBUILDERGTEST_H_ */

#endif /* NOGTEST */
