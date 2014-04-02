/*
 * AdjacencyMatrixGTest.h
 *
 *  Created on: 02.04.2014
 *      Author: Michael
 */

#ifndef ADJACENCYMATRIXGTEST_H_
#define ADJACENCYMATRIXGTEST_H_

#include "gtest/gtest.h"
#include "../AdjacencyMatrix.h"
#include "../../graph/Graph.h"
#include "../../io/METISGraphReader.h"

class AdjacencyMatrixGTest : public testing::Test {
public:
	AdjacencyMatrixGTest();
	virtual ~AdjacencyMatrixGTest();
};

#endif /* ADJACENCYMATRIXGTEST_H_ */
