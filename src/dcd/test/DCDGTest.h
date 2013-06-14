/*
 * DCDGTest.h
 *
 *  Created on: 10.04.2013
 *      Author: cls
 */

#ifndef NOGTEST

#ifndef DCDGTEST_H_
#define DCDGTEST_H_

#include <gtest/gtest.h>

#include "../DynamicLabelPropagation.h"
#include "../../community/LabelPropagation.h"
#include "../../community/Louvain.h"
#include "../../generators/DynamicBarabasiAlbertGenerator.h"

namespace NetworKit {

class DCDGTest: public testing::Test {
public:
	DCDGTest();
	virtual ~DCDGTest();
};

} /* namespace NetworKit */
#endif /* DCDGTEST_H_ */

#endif /*NOGTEST */
