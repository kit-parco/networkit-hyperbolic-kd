/*
 * Dijkstra.cpp
 *
 *  Created on: Jul 23, 2013
 *      Author: Henning, Christian Staudt
 */

#include "Dijkstra.h"

#include <algorithm>

namespace NetworKit {

Dijkstra::Dijkstra(const Graph& G, node source, bool storePaths, bool storeStack) : SSSP(G, source, storePaths, storeStack) {
}

Dijkstra::Dijkstra(const Graph& G, const std::vector<node>& sources) : SSSP(G, sources) {
}




void Dijkstra::run(node t) {

  DEBUG("initializing Dijkstra data structures");
  // init distances
  edgeweight infDist = std::numeric_limits<edgeweight>::max();
  distances.clear();
  distances.resize(G.upperNodeIdBound(), infDist);
  if (storePaths) {
    previous.clear();
    previous.resize(G.upperNodeIdBound());
    npaths.clear();
    npaths.resize(G.upperNodeIdBound(), 0);
    for (node source: sources) {
      npaths[source] = 1;
    }
  }

  if (storeStack) {
    std::stack<node> empty;
    std::swap(stack, empty);
  }

  // priority queue with distance-node pairs
  Aux::PrioQueue<edgeweight, node> pq(distances);

  for (node source: sources) {
    distances[source] = 0;
  }

  auto relax([&](node u, node v, edgeweight w) {
    if (distances[v] > distances[u] + w) {
      distances[v] = distances[u] + w;
      if (storePaths) {
        previous[v] = {u}; // new predecessor on shortest path
        npaths[v] = npaths[u];
      }
      pq.decreaseKey(distances[v], v);
    } else if (storePaths && (distances[v] == distances[u] + w)) {
      previous[v].push_back(u); 	// additional predecessor
      npaths[v] += npaths[u]; 	// all the shortest paths to u are also shortest paths to v now
    }
  });

  bool breakWhenFound = (t != none);
  DEBUG("traversing graph");
  while (pq.size() > 0) {
//		DEBUG("pq size: ", pq.size());

    node current = pq.extractMin().second;
//		DEBUG("pq size: ", pq.size());
//		TRACE("current node in Dijkstra: " , current);
    if (breakWhenFound && t == current) {
      break;
    }

    if (storeStack) {
      stack.push(current);
    }

    G.forEdgesOf(current, relax);
  }

}

void Dijkstra::runUntil(node t) {

  DEBUG("initializing Dijkstra data structures");
  // init distances
  edgeweight infDist = std::numeric_limits<edgeweight>::max();
  distances.clear();
  distances.resize(G.upperNodeIdBound(), infDist);
  if (storePaths) {
    previous.clear();
    previous.resize(G.upperNodeIdBound());
    npaths.clear();
    npaths.resize(G.upperNodeIdBound(), 0);
    for (node source : sources) {
      npaths[source] = 1;
    }
  }

  if (storeStack) {
    std::stack<node> empty;
    std::swap(stack, empty);
  }

  // priority queue with distance-node pairs
  Aux::PrioQueue<edgeweight, node> pq(distances);

  for (node source : sources) {
    distances[source] = 0;
  }


  auto relax([&](node u, node v, edgeweight w) {
    if (distances[v] > distances[u] + w) {
      distances[v] = distances[u] + w;
      if (storePaths) {
        previous[v] = {u}; // new predecessor on shortest path
        npaths[v] = npaths[u];
      }
      pq.decreaseKey(distances[v], v);
    } else if (storePaths && (distances[v] == distances[u] + w)) {
      previous[v].push_back(u); 	// additional predecessor
      npaths[v] += npaths[u]; 	// all the shortest paths to u are also shortest paths to v now
    }
  });


  DEBUG("traversing graph");
  while (pq.size() > 0) {
//		DEBUG("pq size: ", pq.size());

    node current = pq.extractMin().second;
//		DEBUG("pq size: ", pq.size());
//		TRACE("current node in Dijkstra: " , current);

    if (storeStack) {
      stack.push(current);
    }

    if (current == t) {
      break;
    }

    G.forEdgesOf(current, relax);
  }

}


} /* namespace NetworKit */
