/*
 * GraphEventProxy.h
 *
 *  Created on: 03.04.2013
 *      Author: cls
 */

#ifndef GRAPHEVENTPROXY_H_
#define GRAPHEVENTPROXY_H_

#include "../graph/Graph.h"
#include "GraphEventHandler.h"

namespace NetworKit {

/**
 * This class enables the observer pattern for dynamic graphs: It has the same modifiers as a Graph object.
 * When these modifiers are called, they are also called on the underlying graphs. Also, all registered
 * observers (type GraphEventHandler) are notified.
 */
class GraphEventProxy {

protected:

	std::vector<GraphEventHandler*> observers;


public:

	Graph* G;

	GraphEventProxy(Graph& G);

	virtual ~GraphEventProxy();

	void registerObserver(GraphEventHandler& observer);

	node addNode();

	void removeNode(node u);

	void addEdge(node u, node v);

	void removeEdge(node u, node v);

	void setWeight(node u, node v, edgeweight w);
};

} /* namespace NetworKit */
#endif /* GRAPHEVENTPROXY_H_ */
