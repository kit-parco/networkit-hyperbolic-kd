/*
 * DynamicHyperbolicGenerator.h
 *
 *  Created on: 29.07.2014
 *      Author: moritzl
 */

#ifndef DYNAMICHYPERBOLICGENERATOR_H_
#define DYNAMICHYPERBOLICGENERATOR_H_

#include <map>

#include "DynamicGraphGenerator.h"
#include "Quadtree/Quadtree.h"


namespace NetworKit {

class DynamicHyperbolicGenerator: public NetworKit::DynamicGraphGenerator  {
	friend class GeneratorsGTest;
public:
	/**
	 * Initialize a dynamic hyperbolic generator and generate initial node positions
	 * Node movement and change of neighborhood disks happens in discrete time steps
	 *
	 * @param n number of nodes
	 * @param initialFactor initial value of thresholdFactor
	 * @param alpha dispersion parameter, stays fixed
	 * @param stretch geometric stretch factor s, stays fixed
	 * @param moveEachStep fraction of nodes which are moved at each time step, should be non-negative
	 * @param factorGrowth increment added to the value of thresholdFactor at each step, should be non-negative
	 * @param moveDistance base value for the node movements
	 */

	DynamicHyperbolicGenerator(count n = 10, double avgDegree=6, double exp=3, double moveEachStep = 0, double moveDistance = 0);

	/**
	 * Initialize a dynamic hyperbolic generator with given initial node positions in polar coordinates
	 * Node movement and change of neighborhood disks happens in discrete time steps
	 *
	 * @param n number of nodes
	 * @param angles angular coordinates of initial positions
	 * @param radii radial coordinates of initial positions
	 * @param stretch geometric stretch factor s, stays fixed
	 * @param initialFactor initial value of thresholdFactor
	 * @param moveEachStep fraction of nodes which are moved at each time step, should be non-negative
	 * @param factorGrowth increment added to the value of thresholdFactor at each step, should be non-negative
	 * @param moveDistance base value for the node movements
	 */
	DynamicHyperbolicGenerator(std::vector<double> &angles, std::vector<double> &radii,  double R, double alpha, double moveEachStep = 0, double moveDistance = 0);

	/**
	 * Default constructor
	 */
	DynamicHyperbolicGenerator();

	/**
	 * Run the dynamic, changes state of generator
	 *
	 * @param nSteps number of time steps to iterate over
	 * @return changed edges
	 */
	std::vector<GraphEvent> generate(count nSteps) override;

	/**
	 * Get the graph corresponding to the current state of the generator. Does not change the generator
	 *
	 * @return graph at the current state
	 */
	Graph getGraph() const;

	/**
	 * Get coordinates in native representation
	 * @return vector of 2D-Points in Cartesian coordinates
	 */
	std::vector<Point<float> > getCoordinates() const;

private:
	/**
	 * Generate initial node positions and fill the quadtree with them
	 */
	void initializePoints();

	/**
	 * Generate initial movement vectors for all points
	 */
	void initializeMovement();

	/**
	 * Generate initial movement vectors for all points
	 */
	void recomputeBands();

	vector<index> getNeighborsInBands(index i, bool bothDirections=true);

	/**
	 * Execute node movement part of time step
	 *
	 * @param result vector to store GraphEvents in
	 */
	void getEventsFromNodeMovement(vector<GraphEvent> &result);

	/**
     * Execute node movement part of time step
	 *
	 * @param result vector to store GraphEvents in
	 */
	void moveNode(index node);

	count nodeCount;
	double alpha;
	double moveEachStep;
	double moveDistance;
	vector<double> angles;
	vector<double> radii;
	vector<double> angularMovement;
	vector<double> radialMovement;
	vector<double> bandRadii;
	vector<vector<Point2D<double>>> bands;
	vector<vector<double> > bandAngles;
	double R;
	bool initialized;
};

} /* namespace NetworKit */
#endif /* DYNAMICHYPERBOLICGENERATOR_H_ */

