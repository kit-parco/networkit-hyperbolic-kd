/*
 * TSelectiveSCAN.h
 *
 *  Created on: 24.06.2013
 *      Author: cls
 */

#ifndef TSELECTIVESCAN_H_
#define TSELECTIVESCAN_H_

#include "SelectiveCommunityDetector.h"

namespace NetworKit {

template <class Distance> class TSelectiveSCAN: public SelectiveCommunityDetector {
public:

	TSelectiveSCAN(Graph& G, Parameters& param, double epsilon=0.5, double mu=2);

	virtual ~TSelectiveSCAN();

	virtual std::unordered_map<node, std::unordered_set<node> > run(std::unordered_set<node> seeds);

protected:


	void expandCore(node core, node label, std::unordered_set<node>* community,
					std::unordered_map<node, node>* nodesState, std::unordered_set<node>* candidates);

	std::pair<bool,std::unordered_set<node>> isCore(node u);


	Parameters param; //!< parameter class needed to supply custom arguments to Distance while having the same constructor signature across Distance implementations
	double epsilon; //!< threshold value for maximal node distances for nodes to be similar
	double mu;		//!< threshold for number of similar neighbors to define a core
	Distance distMeasure; //!< template distance measure
};




template<class Distance>
inline TSelectiveSCAN<Distance>::TSelectiveSCAN(Graph& G, Parameters& param, double epsilon, double mu) : SelectiveCommunityDetector(G), param(param), epsilon(epsilon), mu(mu), distMeasure(G) {
}


template<class Distance>
inline TSelectiveSCAN<Distance>::~TSelectiveSCAN() {
}



template<class Distance>
inline std::unordered_map<node, std::unordered_set<node> > TSelectiveSCAN<Distance>::run(std::unordered_set<node> seeds) {

	// initialize distance measure
	distMeasure.initialize(param);

	std::unordered_map<node, node> nodesState;
	std::unordered_map<node, std::unordered_set<node>> communities;
	std::unordered_set<node> community;

	G.forNodes([&](node u) {
		nodesState.insert(std::pair<node,int>(u, -1));
	});

	for (node u : seeds) {
		std::pair<bool, std::unordered_set<node>> isCore = this->isCore(u);
		if ((nodesState.find(u))->second == -1 && isCore.first) {
			expandCore(u, u, &community, &nodesState, &isCore.second);
		} else if (((nodesState.find(u))->second != -1)
				&& ((nodesState.find(u))->second != -2)) {
			community = communities.find((nodesState.find(u))->second)->second;
		} else {
			bool clustered = false;
			nodesState.find(u)->second = -2;
			std::pair<bool, std::unordered_set<node>> isCore = this->isCore(u);
			std::pair<bool, std::unordered_set<node>> candidates;

			while (!isCore.second.empty() && !clustered) {
				node core = *isCore.second.begin();
				candidates = this->isCore(core);
				if (candidates.first) {
					expandCore(core, u, &community, &nodesState,
							&candidates.second);
					clustered = true;
				} else {
					isCore.second.erase(isCore.second.begin());
				}
			}
		}
		communities.insert( { u, community });
		community.clear();
	}
	return communities;
}


template<class Distance>
inline void TSelectiveSCAN<Distance>::expandCore(node core,
		node label, std::unordered_set<node>* community,
		std::unordered_map<node, node>* nodesState,
		std::unordered_set<node>* candidates) {
	std::pair<bool,std::unordered_set<node>> isCore;
	community->insert(core);
	nodesState->find(core)->second = label;

	while (!candidates->empty()) {
		node v = *(candidates->begin());
		nodesState->find(v)->second = label;
		community->insert(v);
		isCore = this->isCore(v);
		if (isCore.first) {
			for (node x : isCore.second) {
				int tmp = nodesState->find(x)->second;
				if (tmp < 0) {
					nodesState->find(x)->second = label;
					community->insert(x);
				}
				if (tmp == -1) {
					candidates->insert(x);
				}
			}
		}
		candidates->erase(candidates->find(v));
		// commit test
	}
}

template<class Distance>
inline std::pair<bool, std::unordered_set<node> > TSelectiveSCAN<Distance>::isCore(node u) {
	bool core = false;
	std::unordered_set<node> similarNeighbors;
	int count = 0;
	G.forNeighborsOf(u, [&](node v){
		if (this->distMeasure.distance(u, v) <= this->epsilon) {
			count++;
			similarNeighbors.insert(v);
		}
	});
	if (count >= this->mu) {
		core = true;
	}
	return std::pair<bool,std::unordered_set<node>>(core, similarNeighbors);
}

} /* namespace NetworKit */



#endif /* TSELECTIVESCAN_H_ */
