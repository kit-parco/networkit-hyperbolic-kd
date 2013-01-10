/*
 * ClusteringTest.h
 *
 *  Created on: 12.12.2012
 *      Author: cls
 */

#ifndef CLUSTERINGGTEST_H_
#define CLUSTERINGGTEST_H_

#include <gtest/gtest.h>


#include "../../aux/log.h"
#include "../base/Clustering.h"
#include "../base/Modularity.h"
#include "../base/ClusteringGenerator.h"
#include "../../graph/GraphGenerator.h"
#include "../algo/LabelPropagation.h"

namespace EnsembleClustering {

class ClusteringGTest: public testing::Test {


};



} /* namespace EnsembleClustering */
#endif /* CLUSTERINGTEST_H_ */
