/*
 * StrongConnectedComponents.h
 *
 *  Created on: 01.06.2014
 *      Author: Klara Reichard (klara.reichard@gmail.com), Marvin Ritter (marvin.ritter@gmail.com)
 */

#ifndef STRONGCONNECTEDCOMPONENTS_H_
#define STRONGCONNECTEDCOMPONENTS_H_

#include "../graph/Graph.h"
#include "../structures/Partition.h"

namespace NetworKit {

/**
 * Determines the strongly connected components of an directed graph.
 */
class StrongConnectedComponents {
public:
	StrongConnectedComponents(const Graph& G);
	virtual ~StrongConnectedComponents();

	/**
	 * This method determines the connected components for the graph g.
	 */
	void run();

	/**
	 * This method returns the number of connected components.
	 */
	count numberOfComponents();

	/**
	 * This method returns the the component in which node query is situated.
	 *
	 * @param[in]	query	the node whose component is asked for
	 */
	count componentOfNode(node u);


	/** 
	 * Return a Partition that represents the components
	 */
	Partition getPartition();


private:
	const Graph& G;
	Partition component;
};

}


#endif /* STRONGSrongCONNECTEDCOMPONENTS_H_ */
