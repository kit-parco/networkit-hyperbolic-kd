/*
 * EnsemblePreprocessing.cpp
 *
 *  Created on: 26.02.2013
 *      Author: cls
 */

#include "EnsemblePreprocessing.h"

namespace EnsembleClustering {

EnsemblePreprocessing::EnsemblePreprocessing() : Clusterer() {
	this->finalClusterer = NULL;
}

EnsemblePreprocessing::~EnsemblePreprocessing() {
	// TODO Auto-generated destructor stub
}

void EnsemblePreprocessing::addBaseClusterer(Clusterer& base) {
	this->baseClusterers.push_back(&base);
}

void EnsemblePreprocessing::setFinalClusterer(Clusterer& final) {
	this->finalClusterer = &final;
}

void EnsemblePreprocessing::setOverlapper(Overlapper& overlap) {
	this->overlap = &overlap;
}

Clustering EnsemblePreprocessing::run(Graph& G) {
	INFO("STARTING EnsemblePreprocessing on G=" << G.toString());

	// fixed sub-algorithms
	ClusterContracter contracter;
	ClusteringProjector projector;

	// data
	std::vector<Clustering> baseClusterings(baseClusterers.size(), Clustering(0)); // collection of base clusterings - fill with empty clustering

	// run base clusterers in parallel
	#pragma omp parallel for
	for (int b = 0; b < baseClusterers.size(); b += 1) {
		baseClusterings.at(b) = baseClusterers.at(b)->run(G);
	}

	// ANALYSIS
	if (CALC_DISSIMILARITY) {
		JaccardMeasure dm;
		double dissimilaritySum = 0.0;
		for (int b = 0; b < baseClusterings.size(); b += 1) {
			for (int c = b + 1; c < baseClusterings.size(); c += 1) {
				double d = dm.getDissimilarity(G, baseClusterings.at(b), baseClusterings.at(c));
				dissimilaritySum += d;
			}
		}
		double avgDissimilarity = dissimilaritySum / (baseClusterings.size() * (baseClusterings.size() - 1) / 2.0);
		std::cout << "[INFO] avg. base clustering dissimilarity: " << avgDissimilarity << std::endl;
	}
	//

	// create core clustering
	Clustering core = this->overlap->run(G, baseClusterings);
	// contract graph according to core clustering
	std::pair<Graph, NodeMap<node> > contraction = contracter.run(G, core);
	Graph Gcore = contraction.first;
	NodeMap<node> fineToCoarse = contraction.second;
	// send contracted graph to final clusterer
	Clustering finalCoarse = this->finalClusterer->run(Gcore);

	// project clustering of contracted graph back to original graph
	Clustering final = projector.projectBack(Gcore, G, fineToCoarse, finalCoarse);
	// return clustering
	return final;
}

std::string EnsemblePreprocessing::toString() const {
	std::stringstream strm;
	strm << "TODO: string representation";
	strm << "EnsemblePreprocessing(" << "base=" << this->baseClusterers.front()->toString() << ",ensemble=" << this->baseClusterers.size() << ",final=" << this->finalClusterer->toString() << ")";
	return strm.str();
}

} /* namespace EnsembleClustering */
