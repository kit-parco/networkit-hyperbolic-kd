/*
 * DynamicCommunityDetection.cpp
 *
 *  Created on: 10.01.2014
 *      Author: cls
 */

#include "DynamicCommunityDetection.h"

#include "DynPLP.h"
#include "DynPLM.h"
#include "StaticAdapter.h"
#include "../dynamics/DGSStreamParser.h"
#include "../dynamics/GraphUpdater.h"
#include "../clustering/Modularity.h"
#include "../auxiliary/Timer.h"
#include "../community/PLP.h"
#include "../community/PLM2.h"
#include "../clustering/SampledRandMeasure.h"


namespace NetworKit {

DynamicCommunityDetection::DynamicCommunityDetection(std::string inputPath, std::string algoName, std::string updateStrategy, count interval,
	std::vector<std::string> recordSettings) :
	 inputPath(inputPath), algoName(algoName), updateStrategy(updateStrategy), interval(interval), recordSettings(recordSettings) {
}

void DynamicCommunityDetection::run() {

	// check if property with the given name should be recorded
	auto record = [&](std::string key) {
		return std::find(recordSettings.begin(), recordSettings.end(), key) != recordSettings.end();
	};

	DEBUG("reading full event stream from file");
	DGSStreamParser parser(inputPath, true, 1);
	std::vector<GraphEvent> stream = parser.getStream();
	DEBUG("setting up graph");
	Graph G;
	GraphUpdater gu(G);

	DEBUG("setting up algorithms");
	DynCommunityDetector* algo = NULL;
	if (algoName == "DynPLP") {
		algo = new DynPLP();
	} else if (algoName == "DynPLM") {
		algo = new DynPLM(updateStrategy);
	} else if (algoName == "PLP") {
		algo = new StaticAdapter(new PLP());
	} else if (algoName == "PLM") {
		algo = new StaticAdapter(new PLM2());
	}

	DEBUG("attaching graph");
	algo->attachGraph(G);

	Aux::Timer timer;


	// slicing stream by timesteps
	auto i = stream.begin();
	auto j = stream.begin();

	count run = 0;

	// advance the iterators forward to the next <interval> time steps
	auto advance = [&](){
		count steps = 0;
		i = j;
		for (; (steps < interval) && (j < stream.end()); ++j) {
			if (j->type == GraphEvent::TIME_STEP) {
				steps++;
			}
		}
	};

	DEBUG("advancing dynamic graph");
	while (j != stream.end()) {
		run++; // next run
		advance(); // advance iterators
		std::vector<GraphEvent> slice(i, j);
		DEBUG("updating graph");
		gu.update(slice);
		DEBUG("number of nodes: " << G.numberOfNodes());


		timer.start();
		//
		algo->update(slice);
		//
		timer.stop();
		updateTime.push_back(timer.elapsedMilliseconds());

		timer.start();
		//
		Clustering zeta = algo->detect();
		//
		timer.stop();
		detectTime.push_back(timer.elapsedMilliseconds());

		DEBUG("upper community id bound: " << zeta.upperBound());
		// DEBUG("zeta at time " << G.time() << ": " << Aux::vectorToString(zeta.getVector()));

		if (record("quality")) {
			Modularity mod;
			quality.push_back(mod.getQuality(zeta, G));
			INFO("modularity at time " << G.time() << ": " << quality.back());
		}

		if (record("communitySize")) {
			std::vector<count> sizes = zeta.clusterSizes();
			// double avgSize = std::accumulate(sizes.begin(), sizes.end(), 0) / (double) sizes.size();
			communitySizes.push_back(sizes);
		}

		if (record("continuity") && (run >= 2)) {
			SampledRandMeasure sampledRand(1000);
			double cont = sampledRand.getDissimilarity(G, zeta, previous);
			continuity.push_back(cont);
		}

		previous = zeta; // save last solution for next run

	}


	// record graph size
	size = gu.getSizeTimeline();

}

std::vector<count> DynamicCommunityDetection::getUpdateTimeline() {
	return updateTime;
}

std::vector<count> DynamicCommunityDetection::getDetectTimeline() {
	return detectTime;
}

std::vector<double> DynamicCommunityDetection::getQualityTimeline() {
	return quality;
}

std::vector<double> DynamicCommunityDetection::getContinuityTimeline() {
	return continuity;
}

std::vector<std::pair<count, count> > DynamicCommunityDetection::getGraphSizeTimeline() {
	return size;
}

} /* namespace NetworKit */
