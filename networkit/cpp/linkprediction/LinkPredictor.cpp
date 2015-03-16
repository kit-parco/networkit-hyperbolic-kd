/*
 * LinkPredictor.cpp
 *
 *  Created on: 28.02.2015
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#include "LinkPredictor.h"

namespace NetworKit {

LinkPredictor::LinkPredictor(const Graph& G) : G(G) {
}

std::vector<LinkPredictor::node_dyad_score_pair> LinkPredictor::runAll(count limit) {
  std::priority_queue<node_dyad_score_pair, std::vector<node_dyad_score_pair>, SecondGreater> pairQueue;
  std::vector<node_dyad_score_pair> result;
  G.forNodes([&](node u) {
    G.forNodes([&](node v) {
      if (u < v && !G.hasEdge(u, v)) {
        double score = run(u, v);
        if (limit == 0 || pairQueue.size() < limit) {
          pairQueue.push(std::make_pair(std::make_pair(u, v), score));
        } else if (score > pairQueue.top().second) {
          pairQueue.pop();
          pairQueue.push(std::make_pair(std::make_pair(u, v), score));
        }
      }
    });
  });
  count numEntries = pairQueue.size();
  for (index i = 0; i < numEntries; ++i) {
    result.push_back(pairQueue.top());
    pairQueue.pop();
  }
  return result;
}

} // namespace NetworKit