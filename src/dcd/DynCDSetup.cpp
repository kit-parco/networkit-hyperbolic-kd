/*
 * DynCDSetup.cpp
 *
 *  Created on: 16.06.2013
 *      Author: cls
 */

#include "DynCDSetup.h"

namespace NetworKit {

DynCDSetup::DynCDSetup(DynamicGraphSource& dynGen, std::vector<DynamicCommunityDetector*>& dynDetectors, count tMax, count deltaT) :
		gen(&dynGen),
		detectors(dynDetectors),
		tMax(tMax),
		deltaT(deltaT),
		staticAlgo(NULL) {
	if (deltaT >= tMax) {
		ERROR("deltaT >= tMax");
		throw std::runtime_error("deltaT must be smaller than tMax");
	} else {
		INFO("deltaT " << deltaT);
		INFO("tMax " << tMax);

	}

	// create graph and proxy instances
	Gproxy = this->gen->newGraph();
	G = Gproxy->G;

	DEBUG("setting up " << detectors.size() << " community detection algorithms");
	// register algorithms
	for (DynamicCommunityDetector* dynCD : detectors) {
		Gproxy->registerObserver(dynCD); 	// register community detection algorithm as observer
		dynCD->setGraph(*G);				// point the community detection algorithm to the graph
	}

}

DynCDSetup::~DynCDSetup() {
	// TODO Auto-generated destructor stub
}

void DynCDSetup::run() {
	Aux::Timer runtime;
	runtime.start();

	// helpers
	Modularity modularity;
	DynamicNMIDistance NMID;
	SampledRandMeasure sampledRand(500);

	// initialize graph
	gen->initializeGraph();

	// store the resulting clusterings in results

	for (count i = 0; i < this->detectors.size(); ++i) {
		std::vector<Clustering> dynZeta;
		results.push_back(dynZeta);
	}


	auto runIt = [&]() {
		// inspect the current graph
		INFO("current graph : " << G->toString());

		// run the dynamic community detectors
		for (index detectorIndex = 0; detectorIndex < this->detectors.size(); ++detectorIndex) {
			DynamicCommunityDetector* dynCD = this->detectors.at(detectorIndex);

			INFO("running dynamic community detector " << dynCD->toString());
			results.at(detectorIndex).push_back(dynCD->run());

			// evaluations which need the current graph
			if (checkNumCom) {

				INFO("calculating number of clusters");
				INFO("[RESULT] number of communities: " << results.at(detectorIndex).back().numberOfClusters());
			}

			// modularity
			if (checkMod) {
				double mod = modularity.getQuality(results.at(detectorIndex).back(), *G);
				INFO("[RESULT] communities have modularity: " << mod);
			}

			// NMID
			if (checkNMID && (results[detectorIndex].size() >= 2)) {
				double nmid = NMID.getDissimilarity(*G, results.at(detectorIndex).at(results.at(detectorIndex).size() - 2), results.at(detectorIndex).back());
				INFO("[RESULT] NMID for communities at t=" << G->time() << " vs t=" << (G->time() - deltaT) << ": " << nmid);
			}

			// continuity by sampling
			if (checkSampledRand  && (results[detectorIndex].size() >= 2)) {
				double cont = sampledRand.getDissimilarity(*G, results.at(detectorIndex).at(results.at(detectorIndex).size() - 2), results.at(detectorIndex).back());
				INFO("[RESULT] continuity for " << detectors.at(detectorIndex)->toString() << " at t=" << G->time() << " vs t=" << (G->time() - deltaT) << ": " << cont);

			}
		} // end for multiple detectors

		// optionally also run a static community detector
		if (staticAlgo != NULL) {
			staticClusterings.push_back(staticAlgo->run(*G));

			if (checkSampledRand  && (staticClusterings.size() >= 2)) { // if static algo has been set,
				double contStatic = sampledRand.getDissimilarity(*G, staticClusterings.at(staticClusterings.size() - 2), staticClusterings.back());
				INFO("[RESULT] continuity for " << staticAlgo->toString() <<  " at t=" << G->time() << " vs t=" << (G->time() - deltaT) << ": " << contStatic);
			}
		}



	};

	// for all community detectors, perform run
	bool sourceEnd = false;
	while (G->time() < tMax) {
		try {
			 // try generating the next batch of events
			gen->generateTimeSteps(G->time() + deltaT);

			INFO("=========================== current time step: " << G->time() << " of " << tMax << " ========================================");

		} catch (std::logic_error& e) {
			INFO("exception caught: " << e.what());
			sourceEnd = true;
		}

		// run algorithms and evaluation
		runIt();
		// break from the loop if source cannot produce more events
		if (sourceEnd) {
			INFO("breaking from loop because of end of source");
			break;
		}

	} // end while

	runtime.stop();
	INFO("setup runtime: " << runtime.elapsedTag());


	for (DynamicCommunityDetector* dynCD : this->detectors){
		INFO("timer history for algorithm " << dynCD->toString() << ": " << Aux::vectorToString(dynCD->getTimerHistory()));
	}

}

Graph* DynCDSetup::getGraph() {
	return this->G;
}

void DynCDSetup::setStatic(Clusterer* staticAlgo) {
	this->staticAlgo = staticAlgo;
}

Graph DynCDSetup::getGraphCopy() {
	return *(this->G);
}

void DynCDSetup::checkModularity() {
	this->checkMod = true;
}

void DynCDSetup::checkNumberOfCommunities() {
	this->checkNumCom = true;
}

void DynCDSetup::checkNMIDistance() {
	this->checkNMID = true;
}

void DynCDSetup::checkContinuity() {
	this->checkSampledRand = true;
}

} /* namespace NetworKit */
