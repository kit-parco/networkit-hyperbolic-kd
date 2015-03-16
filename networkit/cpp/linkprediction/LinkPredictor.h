/*
 * LinkPredictor.h
 *
 *  Created on: 28.02.2015
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#ifndef LINKPREDICTOR_H_
#define LINKPREDICTOR_H_

#include <memory>

#include "../graph/Graph.h"

namespace NetworKit {

/**
 * @ingroup linkprediction
 *
 * Abstract base class for link predictors.
 */
class LinkPredictor {
public:
  // Declare typedef in advance for later use
  typedef std::pair<std::pair<node, node>, double> node_dyad_score_pair;

private:
  struct SecondGreater {
    bool operator()(const node_dyad_score_pair& a, const node_dyad_score_pair& b) const {
      return (a.second < b.second) || (a.second == b.second && a.first > b.first);
    }
  };

protected:
  const Graph& G; //!< Graph to operate on

public:
  /**
   * Constructs a new LinkPredictor instance for the graph @a G.
   *
   * @param G The graph to use
   */
  explicit LinkPredictor(const Graph& G);

  /**
   * Default destructor.
   */
  virtual ~LinkPredictor() = default;

  /**
   * Predicts the likelihood of an edge between the given nodes.
   *
   * @param u node in graph
   * @param v node in graph
   * @return a prediction-score indicating the likelihood of an
   * edge between the given nodes
   */
  virtual double run(node u, node v) = 0;

  virtual std::vector<node_dyad_score_pair> runAll(count limit = 0);
};

} // namespace NetworKit

#endif /* LINKPREDICTOR_H_ */