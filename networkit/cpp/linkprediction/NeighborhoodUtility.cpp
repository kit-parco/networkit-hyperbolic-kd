/*
 * NeighborhoodUtility.cpp
 *
 *  Created on: 06.04.2015
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#include "NeighborhoodUtility.h"

namespace NetworKit {

std::pair<std::vector<node>, std::vector<node>> NeighborhoodUtility::prepareNeighbors(const Graph& G, node u, node v) {
  std::vector<node> uNeighbors = G.neighbors(u);
  std::vector<node> vNeighbors = G.neighbors(v);
  // We have no guarantee that the neighbor-vectors are sorted so we have to
  // sort them in order for set-functions to work properly.
  std::sort(uNeighbors.begin(), uNeighbors.end());
  std::sort(vNeighbors.begin(), vNeighbors.end());
  return std::make_pair(uNeighbors, vNeighbors);
}

std::vector<node> NeighborhoodUtility::getNeighborsUnion(const Graph& G, node u, node v) {
  if (!G.hasNode(u) || !G.hasNode(v)) {
    throw std::invalid_argument("Invalid node provided.");
  }
  std::pair<std::vector<node>, std::vector<node>> neighborhoods = prepareNeighbors(G, u, v);
  std::vector<node> neighborsUnion;
  std::set_union(neighborhoods.first.begin(), neighborhoods.first.end(), neighborhoods.second.begin(),
    neighborhoods.second.end(), std::back_inserter(neighborsUnion));
  return neighborsUnion;
}

std::vector<node> NeighborhoodUtility::getCommonNeighbors(const Graph& G, node u, node v) {
  if (!G.hasNode(u) || !G.hasNode(v)) {
    throw std::invalid_argument("Invalid node provided.");
  }
  std::pair<std::vector<node>, std::vector<node>> neighborhoods = prepareNeighbors(G, u, v);
  std::vector<node> commonNeighbors;
  std::set_intersection(neighborhoods.first.begin(), neighborhoods.first.end(), neighborhoods.second.begin(),
    neighborhoods.second.end(), std::back_inserter(commonNeighbors));
  return commonNeighbors;
}

} // namespace NetworKit
