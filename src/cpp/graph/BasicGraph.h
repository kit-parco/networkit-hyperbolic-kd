/*
 * BasicGraph.h
 *
 *  Created on: 20.05.2014
 *      Author: Klara Reichard (klara.reichard@gmail.com), Marvin Ritter (marvin.ritter@gmail.com)
 */

#ifndef BASICGRAPH_H
#define BASICGRAPH_H

#include <algorithm>
#include <limits>
#include <vector>
#include <type_traits>
#include <utility>

// #include "../Globals.h"
// #include "../viz/Point.h"

namespace NetworKit {

/** Typedefs **/

typedef uint64_t index; // more expressive name for an index into an array
typedef uint64_t count; // more expressive name for an integer quantity
typedef index node; // node indices are 0-based

constexpr index none = std::numeric_limits<index>::max(); // value for not existing nodes/edges

typedef double edgeweight; // edge weight type
constexpr edgeweight defaultEdgeWeight = 1.0;
constexpr edgeweight nullWeight = 0.0;

typedef std::function<void(node)> FNode;
typedef std::function<void(node, node)> FNodePair;
typedef std::function<void(node, node)> FEdge;
typedef std::function<void(node, node, double)> FEdgeWithWeight;
typedef std::function<bool()> FCondition;

enum class Weighted {
	weighted,
	unweighted
};

enum class Directed {
	directed,
	undirected
};

// hide implementation details in extra namespace
namespace graph_impl {

// data structures for special graph classes, the BasicGraph class will inherit this private as needed
// TODO: comment variables
struct UnweightedData {
};
struct WeightedData {
	std::vector< std::vector<double> > weights;
};

struct UndirectedData {
	UndirectedData(count n = 0) :
		deg(n, 0),
		adja(n)
	{}
	std::vector<count> deg;
	std::vector< std::vector<node> > adja;
};
struct DirectedData {
	DirectedData(count n = 0) :
		inDeg(n, 0),
		outDeg(n, 0),
		inEdges(n),
		outEdges(n)
	{}
	std::vector<count> inDeg;
	std::vector<count> outDeg;
	std::vector< std::vector<node> > inEdges;
	std::vector< std::vector<node> > outEdges;
};

// choose between types
template<bool b, class T1, class T2>
using either = typename std::conditional<b, T1, T2>::type;

// class for all graphs in NetworKit, some methods have special implementation depending on the template parameters
template<Weighted weighted, Directed directed>
class BasicGraph :
	private either<weighted == Weighted::weighted, WeightedData, UnweightedData>,
	private either<directed == Directed::directed, DirectedData, UndirectedData>
{
public:
	BasicGraph(count n = 0);

	count degree(node v) const { return degree_impl(*this, v); }

	template<Weighted w>
	friend count degree_impl(const BasicGraph<w, directed>& G, node v);
	
	count numberOfNodes() const { return n; }
	count numberOfEdges() const { return m; }

	bool isWeighted() const { return weighted == Weighted::weighted; }
	bool isDirected() const { return directed == Directed::directed; }

	int typ() const;

	inline void forNodes(FNode f) const {
		for (node v = 0; v < z; ++v) {
			if (exists[v]) {
				f(v);
			}
		}
	}

	/*EDGE-ITERATORS*/

	void add_Edge(node u, node v, edgeweight weight) { return add_Edge_impl(*this, u, v, weight);}
	
	template<Weighted w>
	friend void add_Edge_impl(const BasicGraph<w, directed> &G, node u, node v, edgeweight weight);

	index find(node u, node v) { return find_impl(*this, u, v);}

	template<Weighted w>
	friend index find_impl(const BasicGraph<w, directed> &G, node u, node v);

	void remove_Edge(node u, node v, edgeweight weight) { return add_Edge_impl(*this, u, v, weight);}
	
	template<Weighted w>
	friend void remove_Edge_impl(const BasicGraph<w, directed> &G, node u, node v, edgeweight weight);

	template<Weighted w, Directed d>
	friend bool hasEdge(const BasicGraph<w,d>, node u, node v);

	

private:

	using DData = either<directed == Directed::directed, DirectedData, UndirectedData>;
	using WData = either<weighted == Weighted::weighted, WeightedData, UnweightedData>;

	// scalars
	count n; //!< current number of nodes
	count m; //!< current number of edges
	node z; //!< current upper bound of node ids
	count t; //!< current time step

	

	// per node data
	std::vector<bool> exists; //!< exists[v] is true if node v has not been removed from the graph

	std::string name;

	std::vector<std::vector<std::vector<double> > > edgeMaps_double; // contains edge maps (u, v) -> double
	std::vector<double> edgeAttrDefaults_double; // stores default value for edgeMaps_double[i] at index i
};

} // namespace graph_impl

using graph_impl::BasicGraph;
// using Graph = BasicGraph<Weighted::unweighted, Directed::undirected>;
// using WeightedGraph = BasicGraph<Weighted::weighted, Directed::undirected>;
// using DirectedGraph = BasicGraph<Weighted::unweighted, Directed::directed>;
// using WeightedDirectedGraph = BasicGraph<Weighted::weighted, Directed::directed>;

template<Weighted w>
using IUndirectedGraph = BasicGraph<w, Directed::undirected>;

template<Weighted w>
using IDirectedGraph = BasicGraph<w, Directed::directed>;

} // namespace NetworKit


// algorithm examples
namespace NetworKit {

// algorithm for all graph classes
template<Weighted w, Directed d>
count degreeSum(const BasicGraph<w, d>& G);

// algorithm for all undirected graph classes
template<Weighted w>
count undirected_algo(const IUndirectedGraph<w>& G);

template<Weighted w>
count directed_algo(const IDirectedGraph<w>& G);

} // namespace NetworKit

#endif /* BASICGRAPH_H */
