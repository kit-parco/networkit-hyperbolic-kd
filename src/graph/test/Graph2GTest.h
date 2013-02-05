/*
 * Graph2GTest.h
 *
 *  Created on: 04.02.2013
 *      Author: cls
 */

#ifndef GRAPH2GTEST_H_
#define GRAPH2GTEST_H_

#include <gtest/gtest.h>

#include "../Graph.h"

namespace EnsembleClustering {

class Graph2GTest: public testing::Test {
public:
	Graph2GTest();
	virtual ~Graph2GTest();
};

} /* namespace EnsembleClustering */
#endif /* GRAPH2GTEST_H_ */
