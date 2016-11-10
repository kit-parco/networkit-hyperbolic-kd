/*
 * QuadNode.h
 *
 *  Created on: 21.05.2014
 *      Author: Moritz v. Looz (moritz.looz-corswarem@kit.edu)
 */

#ifndef QUADNODE_H_
#define QUADNODE_H_

#include <vector>
#include <algorithm>
#include <functional>
#include <assert.h>

#include "SpatialCell.h"
#include "../../auxiliary/Log.h"
#include "../../auxiliary/Parallel.h"
#include "../../geometric/HyperbolicSpace.h"

using std::vector;
using std::min;
using std::max;
using std::cos;

namespace NetworKit {

template <class T, bool poincare = true>
class QuadNode : public NetworKit::SpatialCell<T> {
	friend class QuadTreeGTest;
private:
	/**
	double leftAngle;
	double minR;
	double rightAngle;
	double maxR;
	Point2D<double> a,b,c,d;
	unsigned capacity;
	*/
	static const long unsigned sanityNodeLimit = 10E15; //just assuming, for debug purposes, that this algorithm never runs on machines with more than 4 Petabyte RAM
	/**
	count subTreeSize;
	std::vector<T> content;
	std::vector<Point2D<double> > positions;
	std::vector<double> angles;
	std::vector<double> radii;
	bool isLeaf;
		*/
	bool splitTheoretical;
	double alpha;
	double balance;
	index ID;
	//double lowerBoundR;

public:
	//std::vector<QuadNode> children;

	QuadNode() = default;

	/**
	 * Construct a QuadNode for polar coordinates.
	 *
	 *
	 * @param leftAngle Minimal angular coordinate of region, in radians from 0 to 2\pi
	 * @param minR Minimal radial coordinate of region, between 0 and 1
	 * @param rightAngle Maximal angular coordinate of region, in radians from 0 to 2\pi
	 * @param maxR Maximal radial coordinate of region, between 0 and 1
	 * @param capacity Number of points a leaf cell can store before splitting
	 * @param splitTheoretical Whether to split in a theoretically optimal way or in a way to decrease measured running times
	 * @param alpha dispersion Parameter of the point distribution. Only has an effect if theoretical split is true
	 * @param balance
	 *
	 */
	QuadNode(const Point<double> &minCoords, const Point<double> &maxCoords, count capacity=1000, bool splitTheoretical = false, double alpha = 1, double balance = 0.5) {
		if (minCoords.getDimensions() != 2) throw std::runtime_error("Quadnode currently only supported for 2D");
		if (maxCoords.getDimensions() != 2) throw std::runtime_error("Quadnode currently only supported for 2D");
		if (poincare && maxCoords[1] >= 1) throw std::runtime_error("The Poincare disk has a radius of 1, cannot create quadtree larger than that!");
		this->isLeaf = true;
		this->subTreeSize = 0;
		this->minCoords = minCoords;
		this->maxCoords = maxCoords;
		this->capacity = capacity;
		this->alpha = alpha;
		this->splitTheoretical = splitTheoretical;
		this->balance = balance;
		this->ID = 0;
		if (balance <= 0 || balance >= 1) throw std::runtime_error("Quadtree balance parameter must be between 0 and 1.");

	}

	void split() {
		assert(this->isLeaf);
		//heavy lifting: split up!
		double leftAngle = this->minCoords[0];
		double rightAngle = this->maxCoords[0];
		double minR = this->minCoords[1];
		double maxR = this->maxCoords[1];
		double middleAngle = (rightAngle - leftAngle) / 2 + leftAngle;
		/**
		 * we want to make sure the space is evenly divided to obtain a balanced tree
		 * Simply halving the radius will cause a larger space for the outer Quadnode, resulting in an unbalanced tree
		 */

		double middleR;

		if (poincare) {
			if (splitTheoretical) {
				double hyperbolicOuter = HyperbolicSpace::EuclideanRadiusToHyperbolic(maxR);
				double hyperbolicInner = HyperbolicSpace::EuclideanRadiusToHyperbolic(minR);
				double hyperbolicMiddle = acosh((1-balance)*cosh(alpha*hyperbolicOuter) + balance*cosh(alpha*hyperbolicInner))/alpha;
				middleR = HyperbolicSpace::hyperbolicRadiusToEuclidean(hyperbolicMiddle);
			} else {
				double nom = maxR - minR;
				double denom = pow((1-maxR*maxR)/(1-minR*minR), 0.5)+1;
				middleR = nom/denom + minR;
			}
		} else {
			middleR = acosh((1-balance)*cosh(alpha*maxR) + balance*cosh(alpha*minR))/alpha;
		}

		//one could also use the median here. Results in worse asymptotical complexity, but maybe better runtime?

		assert(middleR < maxR);
		assert(middleR > minR);

		std::shared_ptr<QuadNode<index,poincare> > southwest(new QuadNode<index,poincare>({leftAngle, minR}, {middleAngle, middleR}, this->capacity, splitTheoretical, alpha, balance));
		std::shared_ptr<QuadNode<index,poincare> > southeast(new QuadNode<index,poincare>({middleAngle, minR}, {rightAngle, middleR}, this->capacity, splitTheoretical, alpha, balance));
		std::shared_ptr<QuadNode<index,poincare> > northwest(new QuadNode<index,poincare>({leftAngle, middleR}, {middleAngle, maxR}, this->capacity, splitTheoretical, alpha, balance));
		std::shared_ptr<QuadNode<index,poincare> > northeast(new QuadNode<index,poincare>({middleAngle, middleR}, {rightAngle, maxR}, this->capacity, splitTheoretical, alpha, balance));
		this->children = {southwest, southeast, northwest, northeast};
		this->isLeaf = false;
	}

	std::pair<double, double> distances(const Point<double> &query) const override {
		return this->hyperbolicPolarDistances(query, poincare);
	}

	double distance(const Point<double> &query, index k) const override {
		double result;
		Point<double> pos = this->positions[k];
		if (poincare) {
			result = HyperbolicSpace::poincareMetric(HyperbolicSpace::polarToCartesian(pos[0], pos[1]), HyperbolicSpace::polarToCartesian(query[0], query[1]));
		} else {
			result = HyperbolicSpace::nativeDistance(pos[0], pos[1], query[0], query[1]);
		}
		return result;
	}

	double getLeftAngle() const {
		return this->minCoords[0];
	}

	double getRightAngle() const {
		return this->maxCoords[0];
	}

	double getMinR() const {
		return this->minCoords[1];
	}

	double getMaxR() const {
		return this->maxCoords[1];
	}

	index getID() const {
		return ID;
	}

	/**
	index indexSubtree(index nextID) {
		index result = nextID;
		assert(this->children.size() == 4 || this->children.size() == 0);
		for (int i = 0; i < this->children.size(); i++) {
			result = this->children[i]->indexSubtree(result);
		}
		this->ID = result;
		return result+1;
	}

	index getCellID(const Point<double> query) const {
		if (!this->responsible(query)) return -1;
		if (this->isLeaf) return getID();
		else {
			for (int i = 0; i < 4; i++) {
				index childresult = this->children[i]->getCellID(query);
				if (childresult >= 0) return childresult;
			}
			assert(false); //if responsible
			return -1;
		}
	}

	index getMaxIDInSubtree() const {
		if (this->isLeaf) return getID();
		else {
			index result = -1;
			for (int i = 0; i < 4; i++) {
				result = std::max(this->children[i]->getMaxIDInSubtree(), result);
			}
			return std::max(result, getID());
		}
	}

	count reindex(count offset) {
		if (this->isLeaf)
		{
			#pragma omp task
			{
				index p = offset;
				std::generate(this->content.begin(), this->content.end(), [&p](){return p++;});
			}
			offset += this->size();
		} else {
			for (int i = 0; i < 4; i++) {
				offset = this->children[i]->reindex(offset);
			}
		}
		return offset;
	}
	*/

};
}

#endif /* QUADNODE_H_ */
