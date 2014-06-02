/*
 * NormalizedLaplacianMatrix.cpp
 *
 *  Created on: 20.03.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#include "NormalizedLaplacianMatrix.h"

namespace NetworKit {

NormalizedLaplacianMatrix::NormalizedLaplacianMatrix() : Matrix() {
}

NormalizedLaplacianMatrix::NormalizedLaplacianMatrix(const Graph &graph) : Matrix(graph.numberOfNodes()) {
	graph.forNodes([&](const index i){
		double weightedDegree = graph.weightedDegree(i);

		graph.forWeightedNeighborsOf(i, [&](const index j, double weight){
			if (i != j) {
				double weightedNeighborDegree = graph.weightedDegree(j);
				setValue(i, j, -weight/sqrt(weightedDegree * weightedNeighborDegree));
			}
		});

		if (weightedDegree != 0.0) {
			if (graph.isWeighted()) {
				setValue(i, i, 1-(graph.weight(i, i)) / weightedDegree);
			} else {
				setValue(i, i, 1);
			}
		}
	});
}

NormalizedLaplacianMatrix::~NormalizedLaplacianMatrix() {
}

} /* namespace NetworKit */
