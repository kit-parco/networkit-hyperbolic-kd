/*
 * GraphTest.h
 *
 *  Created on: 06.12.2012
 *      Author: cls
 */

#ifndef GRAPHTEST_H_
#define GRAPHTEST_H_

#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>

#include "../Graph.h"
#include "../GraphGenerator.h"
#include "../../aux/log.h"

namespace EnsembleClustering {

class GraphTest : public CPPUNIT_NS::TestFixture {

	CPPUNIT_TEST_SUITE(GraphTest);
	CPPUNIT_TEST(testCppUnit);
	CPPUNIT_TEST(testIteration);
	CPPUNIT_TEST(testForallEdges);
	CPPUNIT_TEST_SUITE_END();

protected:

	GraphGenerator gen;
	Graph randomGraph;


public:

	GraphTest();

	virtual ~GraphTest();

	void setUp();

	void tearDown();

	/*** Tests **/

	void testIteration();

	void testCppUnit();

	void testForallEdges();
};

} /* namespace EnsembleClustering */
#endif /* GRAPHTEST_H_ */
