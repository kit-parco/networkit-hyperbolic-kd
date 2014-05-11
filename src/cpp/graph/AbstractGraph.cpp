/*
 * AbstractGraph.cpp
 *
 *  Created on: 10.05.2014
 *      Author: Klara Reichard (klara.reichard@gmail.com), Marvin Ritter (marvin.ritter@gmail.com)
 */

#include <stdexcept>

#include "AbstractGraph.h"
#include "../auxiliary/Random.h"

namespace NetworKit {

AbstractGraph::AbstractGraph(count n, bool weighted):
	n(n),
	m(0),
	z(n),
	t(0),
	weighted(weighted),
	exists(z, true) {

	// set name from global id
	std::stringstream sstm;
	sstm << "G#" << graphId;
	this->name = sstm.str();
}

void AbstractGraph::setName(std::string name) {
	this->name = name;
}

std::string AbstractGraph::getName() const {
	return name;
}

std::string AbstractGraph::toString() const {
	std::stringstream strm;
	strm << "Graph(name=" << this->getName() << ", n=" << this->numberOfNodes()
			<< ", m=" << this->numberOfEdges() << ")";
	return strm.str();
}


/** NODE MODIFIERS **/

node AbstractGraph::addNode() {
	node v = this->z;	// node gets maximum id
	this->z++;	// increment node range
	this->n++;	// increment node count

	// update per node data structures
	this->exists.push_back(true);

	return v;
}

node AbstractGraph::addNode(float x, float y) {
	node v = addNode();
	std::vector<float> coords = {x, y};
	coordinates.addCoordinates(coords);
	return v;
}

void AbstractGraph::removeNode(node v) {
	assert (v < z);
	assert (this->exists[v]);

	if (degree(v) > 0) {
		throw std::runtime_error("nodes must have degree 0 before they can be removed");
	}

	this->exists[v] = false;
	this->n--;
}

bool AbstractGraph::hasNode(node v) const {
	return (v < z) && this->exists[v];	// exists array determines whether node is present
}

node AbstractGraph::randomNode() const {
	assert (this->numberOfNodes() > 0);
	node u;
	do {
		u = Aux::Random::integer(z);
	} while (!this->hasNode(u));
	return u;
}


/** GLOBAL PROPERTIES **/

bool AbstractGraph::isWeighted() const {
	return weighted;
}

bool AbstractGraph::isEmpty() const {
	return (n == 0);
}

count AbstractGraph::numberOfNodes() const {
	return n;
}

count AbstractGraph::numberOfEdges() const {
	return m;
}

count AbstractGraph::numberOfSelfLoops() const {
	count c = 0;
	this->forEdges([&](node u, node v) {
		if (u == v) {
			c += 1;
		}
	});
	return c;
}

index AbstractGraph::upperNodeIdBound() const {
	return z;
}

/** DYNAMICS **/

void AbstractGraph::timeStep() {
	this->t++;
}

count AbstractGraph::time() {
	return this->t;
}


/** COORDINATES **/

void AbstractGraph::setCoordinate(node v, Point<float> value) {
	coordinates.setCoordinate(v, value);
}

Point<float>& AbstractGraph::getCoordinate(node v) {
	return coordinates.getCoordinate(v);
}

float AbstractGraph::minCoordinate(count dim) {
	return coordinates.minCoordinate(dim);
}

float AbstractGraph::maxCoordinate(count dim) {
	return coordinates.maxCoordinate(dim);
}

void AbstractGraph::initCoordinates() {
	coordinates.init(z);
}

} /* namespace NetworKit */
