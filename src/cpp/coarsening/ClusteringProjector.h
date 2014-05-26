/*
 * ClusteringProjector.h
 *
 *  Created on: 07.01.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef CLUSTERINGPROJECTOR_H_
#define CLUSTERINGPROJECTOR_H_

#include "../structures/Partition.h"

namespace NetworKit {

class ClusteringProjector {

public:

	ClusteringProjector();

	virtual ~ClusteringProjector();

	/**
	 * DEPRECATED
	 *
	 * Given
	 * 		@param[in]	Gcoarse
	 * 		@param[in] 	Gfine
	 * 		@param[in]	fineToCoarse
	 * 		@param[in]	zetaCoarse	a clustering of the coarse graph
	 *
	 * 	, project the clustering back to the fine graph to create a clustering of the fine graph.
	 * 		@param[out] 			a clustering of the fine graph
	 **/
	//virtual Partition projectBack(Graph& Gcoarse, Graph& Gfine, std::vector<node>& fineToCoarse,	Partition& zetaCoarse);


	/**
	 * Given
	 * 		@param[in]	Gcoarse
	 * 		@param[in] 	Gfine
	 * 		@param[in]	fineToCoarse
	 * 		@param[in]	zetaCoarse	a clustering of the coarse graph
	 *
	 * 	, project the clustering back to the fine graph to create a clustering of the fine graph.
	 * 		@param[out] 			a clustering of the fine graph
	 **/
	virtual Partition projectBack(Graph& Gcoarse, Graph& Gfine, std::vector<node>& fineToCoarse,
			Partition& zetaCoarse);



	/**
	 * Project a clustering \zeta^{i} of the coarse graph G^{i}�back to
	 * the finest graph G^{0}, using the hierarchy of fine->coarse maps
	 */
	virtual Partition projectBackToFinest(Partition& zetaCoarse,
			std::vector<std::vector<node> >& maps, Graph& Gfinest);


	/**
	 * Assuming that the coarse graph resulted from contracting and represents a clustering of the finest graph
	 *
	 * @param[in]	Gcoarse		coarse graph
	 * @param[in]	Gfinest		finest graph
	 * @param[in]	maps		hierarchy of maps M^{i->i+1} mapping nodes in finer graph to supernodes in coarser graph
	 */
	virtual Partition projectCoarseGraphToFinestClustering(Graph& Gcoarse, Graph& Gfinest, std::vector<std::vector<node> >& maps);

};

} /* namespace NetworKit */
#endif /* CLUSTERINGPROJECTOR_H_ */
