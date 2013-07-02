/*
 * BalancedLabelPropagation.cpp
 *
 *  Created on: Jun 19, 2013
 *      Author: Henning
 */

#include "BalancedLabelPropagation.h"

namespace NetworKit {

BalancedLabelPropagation::BalancedLabelPropagation(double myexponent): exponent(myexponent) {

}

BalancedLabelPropagation::~BalancedLabelPropagation() {

}

Clustering BalancedLabelPropagation::run(Graph& graph, count numParts) {
	ClusteringGenerator gen;
	// FIXME: change to region growing from random
	Clustering partition = gen.makeRandomClustering(graph, numParts);
	partition = rerun(graph, numParts, partition);
	return partition;
}


Clustering& BalancedLabelPropagation::rerun(Graph& graph, count numParts, Clustering& partition) {
	float avg = ceil((float) graph.numberOfNodes() / (float) numParts);
	std::vector<float> adjustByFactor(numParts);

	EdgeCut edgeCut;

	std::vector<count> clusterSizes = partition.clusterSizes();
	for (index p = 0; p < clusterSizes.size(); ++p) {
		DEBUG("size of cluster " << p << ": " << clusterSizes[p]);
	}

	count numIters = 15;
	if (exponent >= 4.0) {
		numIters = 3;
	}
	for (index i = 0; i < numIters; ++i) { // FIXME: different termination criterion
		// read cluster sizes and compute scale values
		clusterSizes = partition.clusterSizes();
		for (index p = 0; p < numParts; ++p) {
			adjustByFactor[p] = avg / clusterSizes[p];
		}

		// perform LP
		graph.forNodes([&](node v) {

			std::map<cluster, double> labelWeights; // neighborLabelCounts maps label -> frequency in the neighbors
			cluster heaviest = none;

			// weigh the labels in the neighborhood of v
			graph.forWeightedNeighborsOf(v, [&](node w, edgeweight weight) {
				cluster lw = partition[w];
				labelWeights[lw] += weight * adjustByFactor[lw];
				heaviest = lw; // init heaviest by some neighbor's cluster
			});

//			heaviest = std::max_element(labelWeights.begin(),
//							labelWeights.end(),
//							[](const std::pair<cluster, edgeweight>& p1, const std::pair<cluster, edgeweight>& p2) {
//								return p1.second < p2.second;})->first;

			if (labelWeights.size() > 1) {

				std::vector<std::pair<cluster, double> > neighborhood;
				for (std::map<cluster, double>::iterator iter = labelWeights.begin();
						iter != labelWeights.end(); ++iter) {
					neighborhood.push_back(std::make_pair(iter->first, iter->second));
				}

				double sum = 0.0;
				for (std::vector<std::pair<cluster, double>>::iterator iter = neighborhood.begin();
						iter != neighborhood.end(); ++iter) {
					sum += pow(iter->second, exponent); // favor large int-/ext-degrees
				}

				// normalize to sum 1 and prefix sum
				for (std::vector<std::pair<cluster, double>>::iterator iter = neighborhood.begin();
						iter != neighborhood.end(); ++iter) {
					iter->second = pow(iter->second, exponent) / sum;
				}

				double prefixSum = 0.0;
				double temp = 0.0;
				for (std::vector<std::pair<cluster, double>>::iterator iter = neighborhood.begin();
						iter != neighborhood.end(); ++iter) {
					temp = iter->second;
					iter->second += prefixSum;
					prefixSum += temp;
				}

				// get heaviest label
				Aux::RandomProbability randGen;
				double prob = randGen.generate();

				for (std::vector<std::pair<cluster, double>>::iterator iter = neighborhood.begin();
						iter != neighborhood.end(); ++iter) {
					if (iter->second >= prob) {
						heaviest = iter->first;
						break;
					}
				}
			}

			if (partition[v] != heaviest) { // UPDATE
				partition[v] = heaviest;
			}


//			DEBUG("partition[" << v << "]: " << partition[v]);
		});

		DEBUG("cut in iter " << i << ": " << edgeCut.getQuality(partition, graph));
	}

	clusterSizes = partition.clusterSizes();
	for (index p = 0; p < clusterSizes.size(); ++p) {
		DEBUG("after: size of cluster " << p << ": " << clusterSizes[p]);
	}

	return partition;
}

Clustering& BalancedLabelPropagation::postsmooth(Graph& graph, count numBlocks,
		Clustering& partition) {
	double tmpExp = this->exponent;
	this->setExponent(8.0); // quasi-deterministic
	partition = this->rerun(graph, numBlocks, partition);
	this->exponent = tmpExp;
	return partition;
}

} /* namespace NetworKit */
