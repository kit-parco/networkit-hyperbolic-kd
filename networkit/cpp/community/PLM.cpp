/*
 * MLPLM.cpp
 *
 *  Created on: 20.11.2013
 *      Author: cls
 */

#include "PLM.h"
#include <omp.h>
#include "../coarsening/ParallelPartitionCoarsening.h"
#include "../coarsening/ClusterContractor.h"
#include "../coarsening/ClusteringProjector.h"
#include "../auxiliary/Log.h"

#include <sstream>

namespace NetworKit {

PLM::PLM(const Graph& G, bool refine, double gamma, std::string par, count maxIter, bool parallelCoarsening, bool turbo) : CommunityDetectionAlgorithm(G), parallelism(par), refine(refine), gamma(gamma), maxIter(maxIter), parallelCoarsening(parallelCoarsening), turbo(turbo) {

}

void PLM::run() {
	DEBUG("calling run method on " , G.toString());

	count z = G.upperNodeIdBound();


	// init communities to singletons
	Partition zeta(z);
	zeta.allToSingletons();
	index o = zeta.upperBound();

	// init graph-dependent temporaries
	std::vector<double> volNode(z, 0.0);
	// $\omega(E)$
	edgeweight total = G.totalEdgeWeight();
	DEBUG("total edge weight: " , total);
	edgeweight divisor = (2 * total * total); // needed in modularity calculation

	G.parallelForNodes([&](node u) { // calculate and store volume of each node
		volNode[u] += G.weightedDegree(u);
		volNode[u] += G.weight(u, u); // consider self-loop twice
		// TRACE("init volNode[" , u , "] to " , volNode[u]);
	});

	// init community-dependent temporaries
	std::vector<double> volCommunity(o, 0.0);
	zeta.parallelForEntries([&](node u, index C) { 	// set volume for all communities
		if (C != none)
			volCommunity[C] = volNode[u];
	});

	// first move phase
	bool moved = false; // indicates whether any node has been moved in the last pass
	bool change = false; // indicates whether the communities have changed at all

	// stores the affinity for each neighboring community (index), one vector per thread
	std::vector<std::vector<edgeweight> > turboAffinity;
	// stores the list of neighboring communities, one vector per thread
	std::vector<std::vector<index> > neigh_comm;


	if (turbo) {
		if (this->parallelism != "none") { // initialize arrays for all threads only when actually needed
			turboAffinity.resize(omp_get_max_threads());
			neigh_comm.resize(omp_get_max_threads());
			for (auto &it : turboAffinity) {
				// resize to maximum community id
				it.resize(zeta.upperBound());
			}
			for (auto &it : neigh_comm) {
				// resize to maximum size so no further memory allocation is needed
				it.resize(G.upperNodeIdBound());
			}
		} else { // initialize array only for first thread
			turboAffinity.emplace_back(zeta.upperBound());
			neigh_comm.emplace_back(G.upperNodeIdBound());
		}
	}

	// try to improve modularity by moving a node to neighboring clusters
	auto tryMove = [&](node u) {
		// TRACE("trying to move node " , u);
		index tid = omp_get_thread_num();

		// collect edge weight to neighbor clusters
		std::map<index, edgeweight> affinity;
		// only for turbo mode, stores the number of neighbors
		index numNeighbors = 0;

		if (turbo) {
			G.forNeighborsOf(u, [&](node v) {
				turboAffinity[tid][zeta[v]] = -1; // set all to -1 so we can see when we get to it the first time
			});
			turboAffinity[tid][zeta[u]] = 0;
			G.forNeighborsOf(u, [&](node v, edgeweight weight) {
				if (u != v) {
					index C = zeta[v];
					if (turboAffinity[tid][C] == -1) {
						// found the neighbor for the first time, initialize to 0 and add to list of neighboring communities
						turboAffinity[tid][C] = 0;
						neigh_comm[tid][numNeighbors++] = C;
					}
					turboAffinity[tid][C] += weight;
				}
			});
		} else {
			G.forNeighborsOf(u, [&](node v, edgeweight weight) {
				if (u != v) {
					index C = zeta[v];
					affinity[C] += weight;
				}
			});
		}


		// sub-functions

		// $\vol(C \ {x})$ - volume of cluster C excluding node x
		auto volCommunityMinusNode = [&](index C, node x) {
			double volC = 0.0;
			double volN = 0.0;
			volC = volCommunity[C];
			if (zeta[x] == C) {
				volN = volNode[x];
				return volC - volN;
			} else {
				return volC;
			}
		};

		// // $\omega(u | C \ u)$
		// auto omegaCut = [&](node u, index C) {
		// 	return affinity[C];
		// };

		auto modGain = [&](node u, index C, index D, edgeweight affinityC, edgeweight affinityD) {
			double volN = 0.0;
			volN = volNode[u];
			double delta = (affinityD - affinityC) / total + this->gamma * ((volCommunityMinusNode(C, u) - volCommunityMinusNode(D, u)) * volN) / divisor;
			//TRACE("(" , affinity[D] , " - " , affinity[C] , ") / " , total , " + " , this->gamma , " * ((" , volCommunityMinusNode(C, u) , " - " , volCommunityMinusNode(D, u) , ") *" , volN , ") / 2 * " , (total * total));
			return delta;
		};

		index best = none;
		index C = none;
		index D = none;
		double deltaBest = -1;

		C = zeta[u];

		if (turbo) {
			edgeweight affinityC = turboAffinity[tid][C];

			for (index i = 0; i < numNeighbors; ++i) {
				D = neigh_comm[tid][i];

				if (D != C) { // consider only nodes in other clusters (and implicitly only nodes other than u)
					double delta = modGain(u, C, D, affinityC, turboAffinity[tid][D]);

					// TRACE("mod gain: " , delta);
					if (delta > deltaBest) {
						deltaBest = delta;
						best = D;
					}
				}
			}
		} else {
			edgeweight affinityC = affinity[C];

//			TRACE("Processing neighborhood of node " , u , ", which is in cluster " , C);
			for (auto it : affinity) {
				D = it.first;
				if (D != C) { // consider only nodes in other clusters (and implicitly only nodes other than u)
					double delta = modGain(u, C, D, affinityC, it.second);
					// TRACE("mod gain: " , delta);
					if (delta > deltaBest) {
						deltaBest = delta;
						best = D;
					}
				}
			}
		}

		// TRACE("deltaBest=" , deltaBest);
		if (deltaBest > 0) { // if modularity improvement possible
			assert (best != C && best != none);// do not "move" to original cluster

			zeta[u] = best; // move to best cluster
			// TRACE("node " , u , " moved");

			// mod update
			double volN = 0.0;
			volN = volNode[u];
			// update the volume of the two clusters
			#pragma omp atomic update
			volCommunity[C] -= volN;
			#pragma omp atomic update
			volCommunity[best] += volN;

			moved = true; // change to clustering has been made

		} else {
			// TRACE("node " , u , " not moved");
		}
	};

	// performs node moves
	auto movePhase = [&](){
		count iter = 0;
		do {
			moved = false;
			// apply node movement according to parallelization strategy
			if (this->parallelism == "none") {
				G.forNodes(tryMove);
			} else if (this->parallelism == "simple") {
				G.parallelForNodes(tryMove);
			} else if (this->parallelism == "balanced") {
				G.balancedParallelForNodes(tryMove);
			} else {
				ERROR("unknown parallelization strategy: " , this->parallelism);
				throw std::runtime_error("unknown parallelization strategy");
			}
			if (moved) change = true;

			if (iter == maxIter) {
				WARN("move phase aborted after ", maxIter, " iterations");
			}
			iter += 1;
		} while (moved && (iter <= maxIter));
		DEBUG("iterations in move phase: ", iter);
	};

	// first move phase
	movePhase();

	if (change) {
		INFO("nodes moved, so begin coarsening and recursive call");
		std::pair<Graph, std::vector<node>> coarsened = coarsen(G, zeta, parallelCoarsening);	// coarsen graph according to communitites
		PLM onCoarsened(coarsened.first);
		onCoarsened.run();
		Partition zetaCoarse = onCoarsened.getPartition();

		INFO("coarse graph has ", coarsened.first.numberOfEdges(), " edges");
		zeta = prolong(coarsened.first, zetaCoarse, G, coarsened.second); // unpack communities in coarse graph onto fine graph
		// refinement phase
		if (refine) {
			DEBUG("refinement phase");
			// reinit community-dependent temporaries
			o = zeta.upperBound();
			volCommunity.clear();
			volCommunity.resize(o, 0.0);
			zeta.parallelForEntries([&](node u, index C) { 	// set volume for all communities
				if (C != none) {
					edgeweight volN = volNode[u];
					#pragma omp atomic update
					volCommunity[C] += volN;
				}
			});

			// second move phase
			movePhase();
		}
	}
	result = std::move(zeta);
	hasRun = true;
}

std::string NetworKit::PLM::toString() const {
	std::string refined;
	if (refine) {
		refined = "refinement";
	} else {
		refined = "";
	}
	std::string parCoarsening;
	if (parallelCoarsening) {
		parCoarsening = "parallel coarsening";
	} else {
		parCoarsening = "";
	}

	std::stringstream stream;
	stream << "PLM(" << parallelism << "," << refined << "," << parCoarsening << ")";

	return stream.str();
}

std::pair<Graph, std::vector<node> > PLM::coarsen(const Graph& G, const Partition& zeta, bool parallel) {
	if (parallel) {
		ParallelPartitionCoarsening parCoarsening(true);
		return parCoarsening.run(G, zeta);
	} else {
		ClusterContractor seqCoarsening;
		return seqCoarsening.run(G, zeta);
	}


}

Partition PLM::prolong(const Graph& Gcoarse, const Partition& zetaCoarse, const Graph& Gfine, std::vector<node> nodeToMetaNode) {
	Partition zetaFine(Gfine.upperNodeIdBound());
	zetaFine.setUpperBound(zetaCoarse.upperBound());

	Gfine.forNodes([&](node v) {
		node mv = nodeToMetaNode[v];
		index cv = zetaCoarse[mv];
		zetaFine[v] = cv;
	});


	return zetaFine;
}

} /* namespace NetworKit */
