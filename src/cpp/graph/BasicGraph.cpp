/*
 * BasicGraph.cpp
 *
 *  Created on: 01.06.2014
 *      Author: Klara Reichard (klara.reichard@gmail.com), Marvin Ritter (marvin.ritter@gmail.com)
 */

#include <stdexcept>
#include <sstream>

#include "BasicGraph.h"
#include "../auxiliary/Random.h"

namespace NetworKit {

namespace graph_impl {

template<Weighted w, Directed d>
BasicGraph<w, d>::BasicGraph(count n) :
	DData(n),
	n(n),
	m(0),
	z(n),
	t(0),
	exists(n, true) {
	
	// set name from global id
	static count nextGraphId = 1;
	std::stringstream sstm;
	sstm << "G#" << nextGraphId++;
	this->name = sstm.str();
}


/** GRAPH INFORMATION **/

template<>
std::string BasicGraph<Weighted::unweighted, Directed::undirected>::typ() const { return "Graph"; }

template<>
std::string BasicGraph<Weighted::weighted, Directed::undirected>::typ() const { return "WeightedGraph"; }

template<>
std::string BasicGraph<Weighted::unweighted, Directed::directed>::typ() const { return "DirectedGraph"; }

template<>
std::string BasicGraph<Weighted::weighted, Directed::directed>::typ() const { return "WeightedDirectedGraph"; }

template<Weighted w, Directed d>
count BasicGraph<w, d>::getMemoryUsage() const {
	count mem = WData::getMemoryUsage() + DData::getMemoryUsage();

	mem += exists.capacity() / 8;

	for (auto& map : edgeMaps_double) {
		for (auto& a : map) {
			mem += sizeof(double) * a.capacity();
		}
	}

	return mem;
}

template<Weighted w, Directed d>
void BasicGraph<w, d>::shrinkToFit() {
	WData::shrinkToFit();
	DData::shrinkToFit();

	exists.shrink_to_fit();

	edgeMaps_double.shrink_to_fit();
	for (auto& map : edgeMaps_double) {
		map.shrink_to_fit();
		for (auto& a : map) {
			a.shrink_to_fit();
		}
	}

}

template<Weighted w, Directed d>
std::string BasicGraph<w, d>::toString() const {
	std::stringstream strm;
	strm << typ() << "(name=" << getName() << ", n=" << numberOfNodes()
			<< ", m=" << numberOfEdges() << ")";
	return strm.str();
}


/** NODE MODIFIERS **/

template<Weighted w, Directed d>
node BasicGraph<w, d>::addNode() {
	node v = z;	// node gets maximum id
	z++;	// increment node range
	n++;	// increment node count

	// update per node data structures
	exists.push_back(true);

	// update per node data structures
	WData::addNode();
	DData::addNode();

	// update edge attribute data structures
	for (int attrId = 0; attrId < this->edgeMaps_double.size(); attrId++) {
		std::vector<double> attrVector;
		this->edgeMaps_double[attrId].push_back(attrVector);
	}

	return v;
}

template<Weighted w, Directed d>
node BasicGraph<w, d>::addNode(float x, float y) {
	node v = addNode();
	std::vector<float> coords = {x, y};
	coordinates.addCoordinates(coords);
	return v;
}

template<Weighted w, Directed d>
void BasicGraph<w, d>::removeNode(node v) {
	assert (v < z);
	assert (exists[v]);

	if (!isIsolated(v)) {
		throw std::runtime_error("nodes must have degree 0 before they can be removed");
	}

	exists[v] = false;
	n--;
}

/** NODE PROPERTIES **/

template<Weighted w>
bool isIsolated_impl(const BasicGraph<w, Directed::undirected>& G, node v) {
	return G.deg[v] == 0;
}

template<Weighted w>
bool isIsolated_impl(const BasicGraph<w, Directed::directed>& G, node v) {
	return G.inDeg[v] == 0 && G.outDeg[v];
}

template<Weighted w, Directed d>
edgeweight BasicGraph<w, d>::weightedDegree(node v) const {
	if (w) {
		edgeweight sum = 0.0;
		forWeightedNeighborsOf(v, [&](node u, edgeweight ew) {
			sum += ew;
		});
		return sum;
	}
	return defaultEdgeWeight * degree(v);
}

template<Weighted w, Directed d>
node BasicGraph<w, d>::randomNode() const {
	assert (numberOfNodes() > 0);

	node v;
	do {
		v = Aux::Random::integer(z);
	} while (!exists[v]);

	return v;
}


/** EDGE MODIFIERS **/

template<Weighted w>
void addEdge_impl(BasicGraph<w, Directed::undirected>& G, node u, node v, edgeweight ew) {
	assert (u >= 0);
	assert (u < G.z);
	assert (G.exists[u]);
	assert (v >= 0);
	assert (v < G.z);
	assert (G.exists[v]);

	if (u == v) { // self-loop case
		G.adja[u].push_back(u);
		G.deg[u] += 1;
		if (w == Weighted::weighted) {
			G.edgeWeights[u].push_back(ew);
		}
		for (index attrId = 0; attrId < G.edgeMaps_double.size(); ++attrId) {
			double defaultAttr = G.edgeAttrDefaults_double[attrId];
			G.edgeMaps_double[attrId][u].push_back(defaultAttr);
		}
	} else {
		// set adjacency
		G.adja[u].push_back(v);
		G.adja[v].push_back(u);
		// increment degree counters
		G.deg[u]++;
		G.deg[v]++;
		// set edge weight
		if (w == Weighted::weighted) {
			G.edgeWeights[u].push_back(ew);
			G.edgeWeights[v].push_back(ew);	
		}
		// loop over all attributes, setting default attr
		for (index attrId = 0; attrId < G.edgeMaps_double.size(); ++attrId) {
			double defaultAttr = G.edgeAttrDefaults_double[attrId];
			G.edgeMaps_double[attrId][u].push_back(defaultAttr);
			G.edgeMaps_double[attrId][v].push_back(defaultAttr);
		}
	}

	G.m++; // increasing the number of edges
}

template<Weighted w>
void addEdge_impl(BasicGraph<w, Directed::directed>& G, node u, node v, edgeweight ew) {
	assert (u >= 0);
	assert (u < G.z);
	assert (G.exists[u]);
	assert (v >= 0);
	assert (v < G.z);
	assert (G.exists[v]);

	G.m++; // increasing the number of edges
	G.inDeg[u]++;
	G.outDeg[v]++;

	G.outEdges[u].push_back(v);
	if (w == Weighted::weighted) {
		G.edgeWeights[u].push_back(ew);
	}
	G.inEdges[v].push_back(u);

	// loop over all attributes, setting default attr
	for (index attrId = 0; attrId < G.edgeMaps_double.size(); ++attrId) {
		double defaultAttr = G.edgeAttrDefaults_double[attrId];
		G.edgeMaps_double[attrId][u].push_back(defaultAttr);
		G.edgeMaps_double[attrId][v].push_back(defaultAttr);
	}

}

template<Weighted w>
void removeEdge_impl(const BasicGraph<w, Directed::undirected>& G, node u, node v) {
	index ui = G.indexInEdgeArray(v, u);
	index vi = G.indexInEdgeArray(u, v);

	if (vi == none) {
		std::stringstream strm;
		strm << "edge (" << u << "," << v << ") does not exist";
		throw std::runtime_error(strm.str());
	} else {
		G.m--; // decreasing the number of edges
		G.adja[u][vi] = none;
		G.adja[v][ui] = none;
		// decrement degree counters
		G.deg[u]--;
		if (u != v) { // self-loops are counted only once
			G.deg[v]--;
		}
		if (w == Weighted::weighted) {
			// remove edge weight
			G.weights[u][vi] = G.nullWeight;
			G.weights[v][ui] = G.nullWeight;
		}

		// dose not make a lot of sense do remove attributes
	}
}

template<Weighted w>
void removeEdge_impl(const BasicGraph<w, Directed::directed>& G, node u, node v) {
	// remove adjacency, for u it is on outgoing edge, at v an incoming
	index ui = G.indexInInEdgeArray(v, u);
	index vi = G.indexInOutEdgeArray(u, v);

	if (vi == none) {
		std::stringstream strm;
		strm << "edge (" << u << "," << v << ") does not exist";
		throw std::runtime_error(strm.str());
	} else {
		G.m--; // decreasing the number of edges
		G.outEdges[u][vi] = none;
		G.inEdges[v][ui] = none;
		// decrement degree counters
		G.deg[u].out--;
		G.deg[v].in--;
		if (w == Weighted::weighted) {
			// remove edge weight
			G.edgeWeights[u][vi] = G.nullWeight;
			G.edgeWeights[v][ui] = G.nullWeight;
		}

		// dose not make a lot of sense do remove attributes
	}
}

template<Weighted w, Directed d>
bool BasicGraph<w, d>::hasEdge(node u, node v) const {
	return find(u, v) != none;
}

template<Weighted w, Directed d>
node BasicGraph<w, d>::mergeEdge(node u, node v, bool discardSelfLoop) {
	// ToDo implement
	throw std::runtime_error("TODO");
}


/** GLOBAL PROPERTIES **/

template<Weighted w, Directed d>
count BasicGraph<w, d>::numberOfSelfLoops() const {
	count c = 0;
	forEdges([&](node u, node v) {
		if (u == v) {
			c += 1;
		}
	});
	return c;
}


/** COORDINATES **/

template<Weighted w, Directed d>
void BasicGraph<w, d>::setCoordinate(node v, Point<float> value) {
	// ToDo implement
	throw std::runtime_error("TODO");
}

template<Weighted w, Directed d>
Point<float>& BasicGraph<w, d>::getCoordinate(node v) {
	// ToDo implement
	throw std::runtime_error("TODO");
}

template<Weighted w, Directed d>
float BasicGraph<w, d>::minCoordinate(count dim) {
	// ToDo implement
	throw std::runtime_error("TODO");
}

template<Weighted w, Directed d>
float BasicGraph<w, d>::maxCoordinate(count dim) {
	// ToDo implement
	throw std::runtime_error("TODO");
}

template<Weighted w, Directed d>
void BasicGraph<w, d>::initCoordinates() {
	// ToDo implement
	throw std::runtime_error("TODO");
}


/** EDGE ATTRIBUTES **/

template<Directed d>
edgeweight weight_impl(const BasicGraph<Weighted::unweighted, d>& G, node u, node v) {
	return G.hasEdge(u, v) ? defaultEdgeWeight : nullWeight;
}

template<Directed d>
edgeweight weight_impl(const BasicGraph<Weighted::weighted, d>& G, node u, node v) {
	index vi = G.find(u, v);
	if (vi != none) {
		return G.edgeWeights[u][vi];
	} else {
		return nullWeight;
	}
}

template<>
void BasicGraph<Weighted::unweighted, Directed::undirected>::setWeight(node u, node v, edgeweight ew) = delete;

template<>
void BasicGraph<Weighted::weighted, Directed::undirected>::setWeight(node u, node v, edgeweight ew) {
	if (u == v) {
		// self-loop case
		index ui = find(u, u);
		if (ui != none) {
			edgeWeights[u][ui] = ew;
		} else {
			addEdge(u, u, ew);
		}
	} else {
		index vi = find(u, v);
		index ui = find(v, u);
		if ((vi != none) && (ui != none)) {
			edgeWeights[u][vi] = ew;
			edgeWeights[v][ui] = ew;
		} else {
			addEdge(u, v, ew);
		}
	}
}

template<>
void BasicGraph<Weighted::unweighted, Directed::directed>::setWeight(node u, node v, edgeweight ew) = delete;

template<>
void BasicGraph<Weighted::weighted, Directed::directed>::setWeight(node u, node v, edgeweight ew) {
	index vi = find(u, v);
	index ui = find(v, u); // ToDo use findIn
	if ((vi != none) && (ui != none)) {
		edgeWeights[u][vi] = ew;
		edgeWeights[v][ui] = ew;
	} else {
		addEdge(u, v, ew);
	}
}

template<Weighted w, Directed d>
void BasicGraph<w, d>::increaseWeight(node u, node v, edgeweight ew) {
	setWeight(u, v, weight(u, v) + ew);
}

template<Weighted w, Directed d>
int BasicGraph<w, d>::addEdgeAttribute_double(double defaultValue) {
	// ToDo implement
	throw std::runtime_error("TODO");
}

template<Weighted w, Directed d>
double BasicGraph<w, d>::attribute_double(node u, node v, int attrId) const {
	// ToDo implement
	throw std::runtime_error("TODO");
}

template<Weighted w, Directed d>
void BasicGraph<w, d>::setAttribute_double(node u, node v, int attrId, double attr) {
	// ToDo implement
	throw std::runtime_error("TODO");
}


/** SUMS **/

template<Weighted w, Directed d>
edgeweight BasicGraph<w, d>::totalEdgeWeight() const {
	// ToDo implement
	throw std::runtime_error("TODO");
}


/** Collections **/

template<Weighted w, Directed d>
std::vector<node> BasicGraph<w, d>::nodes() const {
	// ToDo implement
	throw std::runtime_error("TODO");
}

template<Weighted w, Directed d>
std::vector<std::pair<node, node> > BasicGraph<w, d>::edges() const {
	// ToDo implement
	throw std::runtime_error("TODO");
}


template<Weighted w, Directed d>
std::vector<node> BasicGraph<w, d>::neighbors(node u) const {
	// ToDo implement
	throw std::runtime_error("TODO");
}


/** NODE ITERATORS **/

template<Weighted w, Directed d>
template<typename L>
void BasicGraph<w, d>::forNodes(L handle) const {
	for (node v = 0; v < z; ++v) {
		if (exists[v]) {
			handle(v);
		}
	}
}

template<Weighted w, Directed d>
template<typename L>
void BasicGraph<w, d>::parallelForNodes(L handle) const {
	#pragma omp parallel for
	for (node v = 0; v < z; ++v) {
		if (exists[v]) {
			handle(v);
		}
	}
}

template<Weighted w, Directed d>
template<typename C, typename L>
void BasicGraph<w, d>::forNodesWhile(C condition, L handle) const {
	for (node v = 0; v < z; ++v) {
		if (exists[v]) {
			if (!condition()) {
				break;
			}
			handle(v);

		}
	}
}

template<Weighted w, Directed d>
template<typename C, typename L>
void BasicGraph<w, d>::forNodes(C condition, L handle) const {
	for (node v = 0; v < z; ++v) {
		if (exists[v]) {
			if (condition()) {
				handle(v);
			}
		}
	}
}

template<Weighted w, Directed d>
template<typename L>
void BasicGraph<w, d>::forNodesInRandomOrder(L handle) const {
	std::vector<node> randVec = nodes();
	random_shuffle(randVec.begin(), randVec.end());
	for (node v : randVec) {
		handle(v);
	}
}

template<Weighted w, Directed d>
template<typename L>
void BasicGraph<w, d>::balancedParallelForNodes(L handle) const {
	#pragma omp parallel for schedule(guided) // TODO: define min block size (and test it!)
	for (node v = 0; v < z; ++v) {
		if (exists[v]) {
			handle(v);
		}
	}
}

template<Weighted w, Directed d>
template<typename L>
void BasicGraph<w, d>::forNodePairs(L handle) const {
	for (node u = 0; u < z; ++u) {
		if (exists[u]) {
			for (node v = u + 1; v < z; ++v) {
				if (exists[v]) {
					handle(u, v);
				}
			}
		}
	}
}

template<Weighted w, Directed d>
template<typename L>
void BasicGraph<w, d>::parallelForNodePairs(L handle) const {
	#pragma omp parallel for
	for (node u = 0; u < z; ++u) {
		if (exists[u]) {
			for (node v = u + 1; v < z; ++v) {
				if (exists[v]) {
					handle(u, v);
				}
			}
		}
	}
}

	
/** EDGE ITERATORS **/

template<Weighted w, typename L>
void forEdges_impl(const BasicGraph<w, Directed::undirected>& G, L handle) {
	for (node u = 0; u < G.z; ++u) {
		for (index i = 0; i < G.adja[u].size(); i++) {
			node v = G.adja[u][i];
			if (v != none) {
				handle(u, v);
			}
		}
	}
}

template<Weighted w, typename L>
void forEdges_impl(const BasicGraph<w, Directed::directed>& G, L handle) {
	for (node u = 0; u < G.z; ++u) {
		for (index i = 0; i < G.outEdges.size(); i++) {
			node v = G.outEdges[i];
			if (v != none) {
				handle(u, v);
			}
		}
	}
}

template<Weighted w, typename L>
void parallelForEdges_impl(const BasicGraph<w, Directed::undirected>& G, L handle) {
	#pragma omp parallel for
	for (node u = 0; u < G.z; ++u) {
		for (index i = 0; i < G.adja[u].size(); i++) {
			node v = G.adja[u][i];
			if (v != none) {
				handle(u, v);
			}
		}
	}
}

template<Weighted w, typename L>
void parallelForEdges_impl(const BasicGraph<w, Directed::directed>& G, L handle) {
	#pragma omp parallel for
	for (node u = 0; u < G.z; ++u) {
		for (index i = 0; i < G.outEdges.size(); i++) {
			node v = G.outEdges[i];
			if (v != none) {
				handle(u, v);
			}
		}
	}
}

template<Weighted w, Directed d>
template<typename L>
void BasicGraph<w, d>::forWeightedEdges(L handle) const {
	// ToDo implement
	throw std::runtime_error("TODO");
}

template<Weighted w, Directed d>
template<typename L>
void BasicGraph<w, d>::parallelForWeightedEdges(L handle) const {
	// ToDo implement
	throw std::runtime_error("TODO");
}

template<Weighted w, Directed d>
template<typename L>
void BasicGraph<w, d>::forEdgesWithAttribute_double(int attrId, L handle) const {
	// ToDo implement
	throw std::runtime_error("TODO");
}


/** NEIGHBORHOOD ITERATORS **/

template<Weighted w, Directed d>
template<typename L>
void BasicGraph<w, d>::forNeighborsOf(node u, L handle) const {
	// ToDo implement
	throw std::runtime_error("TODO");
}

template<Weighted w, Directed d>
template<typename L>
void BasicGraph<w, d>::forWeightedNeighborsOf(node u, L handle) const {
	// ToDo implement
	throw std::runtime_error("TODO");
}

template<Weighted w, Directed d>
template<typename L>
void BasicGraph<w, d>::forEdgesOf(node u, L handle) const {
	// ToDo implement
	throw std::runtime_error("TODO");
}

template<Weighted w, Directed d>
template<typename L>
void BasicGraph<w, d>::forWeightedEdgesOf(node u, L handle) const {
	// ToDo implement
	throw std::runtime_error("TODO");
}


/** REDUCTION ITERATORS **/

template<Weighted w, Directed d>
template<typename L>
double BasicGraph<w, d>::parallelSumForNodes(L handle) const {
	double sum = 0.0;
	#pragma omp parallel for reduction(+:sum)
	for (node v = 0; v < z; ++v) {
		if (exists[v]) {
			sum += handle(v);
		}
	}
	return sum;
}

template<Weighted w, typename L>
double parallelSumForNodes_impl(const BasicGraph<w, Directed::undirected>& G, L handle) {
	double sum = 0.0;
	#pragma omp parallel for reduction(+:sum)
	for (node u = 0; u < G.z; u++) {
		for (index i = 0; i < G.adja[u].size(); i++) {
			node v = G.adja[u][i];
			edgeweight ew = (w == Weighted::weighted) ? G.edgeWeights[u][i] : defaultEdgeWeight;
			if (v != none) {
				sum += handle(u, v, ew);
			}
		}
	}
	return sum;
}

template<Weighted w, typename L>
double parallelSumForNodes_impl(const BasicGraph<w, Directed::directed>& G, L handle) {
	double sum = 0.0;
	#pragma omp parallel for reduction(+:sum)
	for (node u = 0; u < G.z; u++) {
		for (index i = 0; i < G.outEdges[u].size(); i++) {
			node v = G.outEdges[u][i];
			edgeweight ew = (w == Weighted::weighted) ? G.edgeWeights[u][i] : defaultEdgeWeight;
			if (v != none) {
				sum += handle(u, v, ew);
			}
		}
	}
	return sum;
}

/** GRAPH SEARCHES **/

template<Weighted w, Directed d>
template<typename L>
void BasicGraph<w, d>::BFSfrom(node r, L handle) const {
	// ToDo implement
	throw std::runtime_error("TODO");
}

template<Weighted w, Directed d>
template<typename L>
void BasicGraph<w, d>::BFSEdgesfrom(node r, L handle) const {
	// ToDo implement
	throw std::runtime_error("TODO");
}

template<Weighted w, Directed d>
template<typename L>
void BasicGraph<w, d>::DFSfrom(node r, L handle) const {
	// ToDo implement
	throw std::runtime_error("TODO");
}

template<Weighted w, Directed d>
template<typename L>
void BasicGraph<w, d>::DFSEdgesfrom(node r, L handle) const {
	// ToDo implement
	throw std::runtime_error("TODO");
}

} /* namespace graph_impl */

} /* namespace NetworKit */
