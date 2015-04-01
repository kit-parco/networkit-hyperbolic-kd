/*
 * LinkPredictionGTest.cpp
 *
 *  Created on: 07.12.2014
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#ifndef NOGTEST

#include "LinkPredictionGTest.h"
#include "../../io/METISGraphReader.h"
#include "../KatzIndex.h"
#include "../CommonNeighborsIndex.h"
#include "../JaccardIndex.h"
#include "../ROCMetric.h"
#include "../PrecisionRecallMetric.h"
#include "../TrainingGraphGenerator.h"
#include "../KFoldCrossValidator.h"
#include "../UnconnectedNodesFinder.h"
#include "../UDegreeIndex.h"
#include "../VDegreeIndex.h"

namespace NetworKit {

LinkPredictionGTest::LinkPredictionGTest() : G(7) {
}

void LinkPredictionGTest::SetUp() {
  G.addEdge(0, 1);
  G.addEdge(0, 3);
  G.addEdge(1, 2);
  G.addEdge(1, 4);
  G.addEdge(2, 3);
  G.addEdge(2, 4);
  G.addEdge(2, 5);
  G.addEdge(3, 4);
  G.addEdge(3, 5);
  G.addEdge(4, 5);
}

TEST_F(LinkPredictionGTest, testTrainingGraphGenerator) {
  Graph trainingGraph = TrainingGraphGenerator::byPercentage(G, 0.7);
  EXPECT_EQ(7, trainingGraph.numberOfEdges());
}

TEST_F(LinkPredictionGTest, testUDegreeIndexRun) {
  UDegreeIndex uDegreeIndex(G);
  EXPECT_EQ(2, uDegreeIndex.run(0, 3));
  EXPECT_EQ(3, uDegreeIndex.run(1, 0));
  EXPECT_EQ(4, uDegreeIndex.run(2, 3));
  EXPECT_EQ(4, uDegreeIndex.run(3, 4));
  EXPECT_EQ(4, uDegreeIndex.run(4, 5));
  EXPECT_EQ(3, uDegreeIndex.run(5, 4));
}

TEST_F(LinkPredictionGTest, testVDegreeIndexRun) {
  VDegreeIndex vDegreeIndex(G);
  EXPECT_EQ(2, vDegreeIndex.run(3, 0));
  EXPECT_EQ(3, vDegreeIndex.run(0, 1));
  EXPECT_EQ(4, vDegreeIndex.run(3, 2));
  EXPECT_EQ(4, vDegreeIndex.run(4, 3));
  EXPECT_EQ(4, vDegreeIndex.run(5, 4));
  EXPECT_EQ(3, vDegreeIndex.run(4, 5));
}


/*TEST_F(LinkPredictionGTest, testCommonNeighborsRunOn) {
  METISGraphReader graphReader;
  Graph newG = graphReader.read("input/caidaRouterLevel.graph");
  RandomEdgePartitioner partitioner(newG);
  std::pair<Graph, Graph> graphPartitions = partitioner.partitionByPercentage(0.1);

  CommonNeighborsIndex cn(graphPartitions.first);

  UnconnectedNodesFinder unf(graphPartitions.first);
  std::vector<std::pair<node, node>> missingLinks = unf.findAll(2);
  INFO("Found ", missingLinks.size(), " missing links with distance 2");
  std::vector<LinkPredictor::node_dyad_score_pair> scoresParallel = cn.runOnParallel(missingLinks);
  //std::vector<LinkPredictor::node_dyad_score_pair> scores = cn.runOn(missingLinks);

  INFO("Size = ", scoresParallel.size());

  //for (index i = 0; i < scores.size(); ++i) {
    //EXPECT_EQ(scoresParallel[i].second, scores[i].second);
    //INFO("entries[", i, "] = ((", scores[i].first.first, ", ", scores[i].first.second, "), ", scores[i].second, ")");
  //}
}*/

/*TEST_F(LinkPredictionGTest, testKFoldCrossValidator) {
  KatzIndex katzIndex;
  METISGraphReader graphReader;
  Graph jazz = graphReader.read("input/jazz.graph");
  ROC roc;

  KFoldCrossValidator validator(jazz, &katzIndex, &roc);
  double average = validator.crossValidate(10);
  EXPECT_NEAR(average, 0.78, 0.03);
}*/

/*TEST_F(LinkPredictionGTest, testMissingLinkFinder) {
  METISGraphReader graphReader;
  Graph jazz = graphReader.read("input/PGPgiantcompo.graph");
  UnconnectedNodesFinder unf(jazz);
  std::vector<std::pair<node, node>> missingEdges = unf.findAll(2);
  INFO(missingEdges.size());
}*/

/*TEST_F(LinkPredictionGTest, testPrecisionRecall) {
  METISGraphReader graphReader;
  Graph newG = graphReader.read("input/hep-th.graph");
  RandomEdgePartitioner partitioner(newG);
  std::pair<Graph, Graph> graphPartitions = partitioner.partitionByPercentage(0.3);

  for (std::pair<node, node> e : graphPartitions.second.edges()) {
    INFO("Removed edge (", e.first, ", ", e.second, ").");
  }

  UnconnectedNodesFinder unf(graphPartitions.first);
  std::vector<std::pair<node, node>> nodePairs = unf.findAll(2);
  INFO("nodePairs.size() = ", nodePairs.size());

  CommonNeighborsIndex cn(graphPartitions.first);
  std::vector<LinkPredictor::node_dyad_score_pair> scores = cn.runOnParallel(nodePairs);
  LinkPredictor::sortByScore(scores);
  for (index i = 0; i < 50; ++i) {
    INFO("entries[", i, "] = ((", scores[i].first.first, ", ", scores[i].first.second, "), ", scores[i].second, ")");
  }
  PrecisionRecallMetric pr(graphPartitions.second, scores);
  pr.generatePoints();
  std::pair<std::vector<double>, std::vector<double>> points = pr.getPoints();
  INFO("Size = ", points.first.size());
  INFO("Last point = (", points.first.back(), ", ", points.second.back(), ").");
  for (index i = 0; i < points.first.size(); ++i) {
    //INFO("Point[", i, "] = (", points.first[i], ", ", points.second[i], ").");
  }
}*/

/*TEST_F(LinkPredictionGTest, testUnconnectedNodesFinder) {
  METISGraphReader graphReader;
  Graph newG = graphReader.read("input/caidaRouterLevel.graph");
  UnconnectedNodesFinder unf(newG);
  std::vector<std::pair<node, node>> missingLinks = unf.findAll(2);
  INFO("Size = ", missingLinks.size());
  //for (std::pair<node, node> p : missingLinks) {
  //  INFO("(", p.first, ", ", p.second, ")");
  //}
}*/

/*TEST_F(LinkPredictionGTest, testReceiverOperatingCharacteristic) {
  METISGraphReader graphReader;
  Graph newG = graphReader.read("input/PGPgiantcompo.graph");
  RandomEdgePartitioner partitioner(G);
  std::pair<Graph, Graph> graphPartitions = partitioner.partitionByPercentage(0.3);

  UnconnectedNodesFinder unf(graphPartitions.first);
  std::vector<std::pair<node, node>> nodePairs = unf.findAll(2);
  INFO("nodePairs.size() = ", nodePairs.size());

  KatzIndex katz(graphPartitions.first, 2, 1);
  std::vector<LinkPredictor::node_dyad_score_pair> scores = katz.runOnParallel(nodePairs);
  for (index i = 0; i < scores.size(); ++i) {
    INFO("entries[", i, "] = ((", scores[i].first.first, ", ", scores[i].first.second, "), ", scores[i].second, ")");
  }
  ROCMetric roc(graphPartitions.second, scores);
  roc.generatePoints();
  std::pair<std::vector<double>, std::vector<double>> points = roc.getPoints();
  INFO("Size = ", points.first.size());
  for (index i = 0; i < points.first.size(); ++i) {
    INFO("Point[", i, "] = (", points.first[i], ", ", points.second[i], ").");
  }
}*/


/*
TEST_F(LinkPredictionGTest, testCommonNeighborsRun) {
  CommonNeighborsIndex cni(G);
  EXPECT_EQ(3.0, cni.run(1, 3));
  EXPECT_EQ(2.0, cni.run(3, 5));
  EXPECT_EQ(0.0, cni.run(4, 6));
}

TEST_F(LinkPredictionGTest, testCommonNeighborsRunAll) {
  CommonNeighborsIndex cni(G);
  ScoreCollection scores = cni.runAll();
  EXPECT_EQ(3.0, scores.getScore(1, 3));
  EXPECT_EQ(2.0, scores.getScore(3, 5));
  EXPECT_EQ(0.0, scores.getScore(4, 6));
}*/
/*
TEST_F(LinkPredictionGTest, testPreferentialAttachmentRun) {
  PreferentialAttachmentIndex pai(G);
  EXPECT_EQ(12.0, pai.run(1, 3));
  EXPECT_EQ(16.0, pai.run(2, 3));
  EXPECT_EQ(0.0, pai.run(4, 6));
}

TEST_F(LinkPredictionGTest, testPreferentialAttachmentRunAll) {
  PreferentialAttachmentIndex pai(G);
  ScoreCollection scores = pai.runAll();
  EXPECT_EQ(12.0, scores.getScore(1, 3));
  EXPECT_EQ(16.0, scores.getScore(2, 3));
  EXPECT_EQ(0.0, scores.getScore(4, 6));
}

TEST_F(LinkPredictionGTest, testJaccardCoefficientRun) {
  JaccardCoefficientIndex jci(G);
  EXPECT_DOUBLE_EQ(0.75, jci.run(1, 3));
  EXPECT_DOUBLE_EQ((double) 1 / 3, jci.run(2, 3));
  EXPECT_EQ(0.0, jci.run(4, 6));
}

TEST_F(LinkPredictionGTest, testJaccardCoefficientRunAll) {
  JaccardCoefficientIndex jci(G);
  ScoreCollection scores = jci.runAll();
  EXPECT_DOUBLE_EQ(0.75, scores.getScore(1, 3));
  EXPECT_DOUBLE_EQ((double) 1 / 3, scores.getScore(2, 3));
  EXPECT_EQ(0.0, scores.getScore(4, 6));
}
*/
/*TEST_F(LinkPredictionGTest, testNonexistentResultsException) {
  PreferentialAttachmentIndex pai(G);
  EXPECT_THROW(pai.getResult(1, 2), std::runtime_error);
}*/

// TODO: Write a test to make sure the score of a node with itself
// is [...] (maybe 0, MISSING_SCORE or a new INVALID_SCORE)
// Maybe it should throw an exception?

} // namespace NetworKit

#endif /* NOGTEST */