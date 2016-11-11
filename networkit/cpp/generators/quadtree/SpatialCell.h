/*
 * SpatialCell.h
 *
 *  Created on: 10.11.2016
 *      Author: moritzl
 */

#ifndef SPATIALCELL_H_
#define SPATIALCELL_H_

#include <vector>
#include <functional>
#include <memory>

#include "../../geometric/HyperbolicSpace.h"
#include "../../auxiliary/Log.h"

#include "../../viz/Point.h"

using std::min;
using std::max;
using std::cos;

namespace NetworKit {

template <class T>
class SpatialCell {
	friend class QuadTreeGTest;
public:
	SpatialCell() = default;
	virtual ~SpatialCell() = default;
	SpatialCell(const Point<double> &minCoords, const Point<double> &maxCoords, count capacity=1000) : minCoords(minCoords), maxCoords(maxCoords), capacity(capacity) {
		isLeaf = true;
		subTreeSize = 0;
		ID = 0;
	}

	virtual void split() = 0;
	virtual std::pair<double, double> distances(const Point<double> &query) const = 0;
	virtual double distance(const Point<double> &query, index k) const = 0;

	void getCoordinates(vector<double> &anglesContainer, vector<double> &radiiContainer) const {
		assert(minCoords.getDimensions() == 2);
		if (this->isLeaf) {
			for (Point<double> pos : positions) {
				anglesContainer.push_back(pos[0]);
				radiiContainer.push_back(pos[1]);
			}
		}
		else {
			assert(this->content.size() == 0);

			for (index i = 0; i < children.size(); i++) {
				this->children[i]->getCoordinates(anglesContainer, radiiContainer);
			}
		}
	}

	virtual void coarsen() {
		assert(this->height() == 2);
		assert(content.size() == 0);
		assert(positions.size() == 0);

		vector<T> allContent;
		vector<Point<double> > allPositions;
		for (index i = 0; i < this->children.size(); i++) {
			allContent.insert(allContent.end(), children[i]->content.begin(), children[i]->content.end());
			allPositions.insert(allPositions.end(), children[i]->positions.begin(), children[i]->positions.end());
		}
		assert(this->subTreeSize == allContent.size());
		assert(this->subTreeSize == allPositions.size());

		this->children.clear();
		this->content.swap(allContent);
		this->positions.swap(allPositions);
		this->isLeaf = true;
	}

	/**
	 * Remove content at polar coordinates (angle, R). May cause coarsening of the quadtree
	 *
	 * @param input Content to be removed
	 * @param angle Angular coordinate
	 * @param R Radial coordinate
	 *
	 * @return True if content was found and removed, false otherwise
	 */
	bool removeContent(T input, const Point<double> &pos) {
		if (!this->responsible(pos)) return false;
		if (this->isLeaf) {
			index i = 0;
			for (; i < this->content.size(); i++) {
				if (this->content[i] == input) break;
			}
			if (i < this->content.size()) {
				//remove element
				this->content.erase(this->content.begin()+i);
				this->positions.erase(this->positions.begin()+i);
				return true;
			} else {
				return false;
			}
		}
		else {
			bool removed = false;
			bool allLeaves = true;
			assert(this->children.size() > 0);
			for (index i = 0; i < children.size(); i++) {
				if (!children[i]->isLeaf) allLeaves = false;
				if (children[i]->removeContent(input, pos)) {
					assert(!removed);
					removed = true;
				}
			}
			if (removed) this->subTreeSize--;
			//coarsen?
			if (removed && allLeaves && this->size() < this->coarsenLimit) {
				this->coarsen();
			}

			return removed;
		}
	}

	void recount() {
		this->subTreeSize = 0;
		for (index i = 0; i < children.size(); i++) {
			this->children[i]->recount();
			this->subTreeSize += this->children[i]->size();
		}
	}

	virtual bool outOfReach(Point<double> query, double radius) const {
		return distances(query).first > radius;
	}

	virtual void getElementsInCircle(Point<double> center, double radius, vector<T> &result) const {
		if (outOfReach(center, radius)) {
			return;
		}

		if (this->isLeaf) {
			const count cSize = content.size();

			for (index i = 0; i < cSize; i++) {
				if (distance(center, i) < radius) {
					result.push_back((content[i]));
				}
			}
		} else {
			for (index i = 0; i < children.size(); i++) {
				children[i]->getElementsInCircle(center, radius, result);
			}
		}
	}

	virtual void addContent(T input, Point<double> coords) {
		assert(content.size() == positions.size());
		assert(this->responsible(coords));
		if (isLeaf) {
			if (content.size() + 1 < capacity) {
				content.push_back(input);
				positions.push_back(coords);
			} else {
				split();

				for (index i = 0; i < content.size(); i++) {
					this->addContent(content[i], positions[i]);
				}
				assert(subTreeSize == content.size());//we have added everything twice
				subTreeSize = content.size();
				content.clear();
				positions.clear();
				this->addContent(input, coords);
			}
		}
		else {
			assert(children.size() > 0);
			bool foundResponsibleChild = false;
			for (index i = 0; i < children.size(); i++) {
				if (children[i]->responsible(coords)) {
					foundResponsibleChild = true;
					children[i]->addContent(input, coords);
					break;
				}
			}
			assert(foundResponsibleChild);
			subTreeSize++;
		}
		TRACE("Element added.");
	}

	virtual bool responsible(const Point<double> &point) const {
		const index d = minCoords.getDimensions();
		assert(point.getDimensions() == d);
		for (index i = 0; i < d; i++) {
			if (point[i] < minCoords[i] || point[i] >= maxCoords[i]) return false;
		}
		return true;
	}

	virtual count size() const {
		return isLeaf ? content.size() : subTreeSize;
	}

	virtual void trim() {
		content.shrink_to_fit();
		positions.shrink_to_fit();

		for (index i = 0; i < children.size(); i++) {
			children[i]->trim();
		}
	}

	virtual count getElementsProbabilistically(const Point<double> &query, std::function<double(double)> prob, std::vector<T> &result) const {
		auto distancePair = distances(query);
		double probUB = prob(distancePair.first);
		double probLB = prob(distancePair.second);
		assert(probLB <= probUB);
		if (probUB > 0.5) probUB = 1;//if we are going to take every second element anyway, no use in calculating expensive jumps
		if (probUB == 0) return 0;
		//TODO: return whole if probLB == 1
		double probdenom = std::log(1-probUB);
		if (probdenom == 0) {
			DEBUG(probUB, " not zero, but too small too process. Ignoring.");
			return 0;
		}
		TRACE("probUB: ", probUB, ", probdenom: ", probdenom);

		count expectedNeighbours = probUB*size();
		count candidatesTested = 0;

		if (isLeaf) {
			const count lsize = content.size();
			TRACE("Leaf of size ", lsize);
			for (index i = 0; i < lsize; i++) {
				//jump!
				if (probUB < 1) {
					double random = Aux::Random::real();
					double delta = std::log(random) / probdenom;
					assert(delta == delta);
					assert(delta >= 0);
					i += delta;
					if (i >= lsize) break;
					TRACE("Jumped with delta ", delta, " arrived at ", i);
				}

				//see where we've arrived
				candidatesTested++;
				double dist = distance(query, i);
				assert(dist >= distancePair.first);

				double q = prob(dist);
				q = q / probUB; //since the candidate was selected by the jumping process, we have to adjust the probabilities
				assert(q <= 1);
				assert(q >= 0);

				//accept?
				double acc = Aux::Random::real();
				if (acc < q) {
					TRACE("Accepted node ", i, " with probability ", q, ".");
					result.push_back(content[i]);
				}
			}
		}	else {
			if (expectedNeighbours < 1) {//select candidates directly instead of calling recursively
				TRACE("probUB = ", probUB,  ", switching to direct candidate selection.");
				assert(probUB < 1);
				const count stsize = size();
				for (index i = 0; i < stsize; i++) {
					double delta = std::log(Aux::Random::real()) / probdenom;
					assert(delta >= 0);
					i += delta;
					TRACE("Jumped with delta ", delta, " arrived at ", i, ". Calling maybeGetKthElement.");
					if (i < size()) maybeGetKthElement(probUB, query, prob, i, result);//this could be optimized. As of now, the offset is subtracted separately for each point
					else break;
					candidatesTested++;
				}
			} else {//carry on as normal
				for (index i = 0; i < children.size(); i++) {
					TRACE("Recursively calling child ", i);
					candidatesTested += children[i]->getElementsProbabilistically(query, prob, result);
				}
			}
		}
		//DEBUG("Expected at most ", expectedNeighbours, " neighbours, got ", result.size() - offset);
		return candidatesTested;
	}

	virtual void maybeGetKthElement(double upperBound, Point<double> query, std::function<double(double)> prob, index k, std::vector<T> &circleDenizens) const {
			TRACE("Maybe get element ", k, " with upper Bound ", upperBound);
			assert(k < size());
			if (isLeaf) {
				double dist = distance(query, k);

				double acceptance = prob(dist)/upperBound;
				TRACE("Is leaf, accept with ", acceptance);
				if (Aux::Random::real() < acceptance) circleDenizens.push_back(content[k]);
			} else {
				TRACE("Call recursively.");
				index offset = 0;
				for (index i = 0; i < children.size(); i++) {
					count childsize = children[i]->size();
					if (k - offset < childsize) {
						children[i]->maybeGetKthElement(upperBound, query, prob, k - offset, circleDenizens);
						break;
					}
					offset += childsize;
				}
			}
		}

	//actually, move this thing to hyperbolic space
	std::pair<double, double> hyperbolicPolarDistances(Point<double> query, bool poincare) const {
		double phi = query[0];
		double r = query[1];

		double leftAngle = this->minCoords[0];
		double minR = this->minCoords[1];
		double rightAngle = this->maxCoords[0];
		double maxR = this->maxCoords[1];

		assert(leftAngle >= 0);
		assert(rightAngle <= 2*M_PI);
		assert(leftAngle < rightAngle);

		assert(minR >= 0);
		assert(maxR > minR);
		assert(phi >= 0);
		assert(phi < 2*M_PI);
		assert(r > 0);

		double minRHyper, maxRHyper, r_h;
		if (poincare) {
			minRHyper=HyperbolicSpace::EuclideanRadiusToHyperbolic(minR);
			maxRHyper=HyperbolicSpace::EuclideanRadiusToHyperbolic(maxR);
			r_h = HyperbolicSpace::EuclideanRadiusToHyperbolic(r);
		} else {
			minRHyper=minR;
			maxRHyper=maxR;
			r_h = r;
		}

		double coshr = cosh(r_h);
		double sinhr = sinh(r_h);
		double coshMinR = cosh(minRHyper);
		double coshMaxR = cosh(maxRHyper);
		double sinhMinR = sinh(minRHyper);
		double sinhMaxR = sinh(maxRHyper);
		double cosDiffLeft = cos(phi - leftAngle);
		double cosDiffRight = cos(phi - rightAngle);

		/**
		 * If the query point is not within the quadnode, the distance minimum is on the border.
		 * Need to check whether extremum is between corners:
		 */

		double coshMinDistance, coshMaxDistance;

		//Left border
		double lowerLeftDistance = coshMinR*coshr-sinhMinR*sinhr*cosDiffLeft;
		double upperLeftDistance = coshMaxR*coshr-sinhMaxR*sinhr*cosDiffLeft;
		if (this->responsible(query)) coshMinDistance = 1; //strictly speaking, this is wrong
		else coshMinDistance = min(lowerLeftDistance, upperLeftDistance);

		coshMaxDistance = max(lowerLeftDistance, upperLeftDistance);
		//double a = cosh(r_h);
		double b, extremum;

		if (phi != leftAngle) {
			b = sinhr*cosDiffLeft;
			extremum = log((coshr+b)/(coshr-b))/2;
			if (extremum < maxRHyper && extremum >= minRHyper) {
				double extremeDistance = cosh(extremum)*coshr-sinh(extremum)*sinhr*cosDiffLeft;
				coshMinDistance = min(coshMinDistance, extremeDistance);
				coshMaxDistance = max(coshMaxDistance, extremeDistance);
			}
		} else {
			coshMinDistance = 1;
		}

		/**
		 * cosh is a function from [0,\infty) to [1, \infty)
		 * Variables thus need
		 */
		assert(coshMaxDistance >= 1);
		assert(coshMinDistance >= 1);

		//Right border
		double lowerRightDistance = coshMinR*coshr-sinhMinR*sinhr*cosDiffRight;
		double upperRightDistance = coshMaxR*coshr-sinhMaxR*sinhr*cosDiffRight;
		coshMinDistance = min(coshMinDistance, lowerRightDistance);
		coshMinDistance = min(coshMinDistance, upperRightDistance);
		coshMaxDistance = max(coshMaxDistance, lowerRightDistance);
		coshMaxDistance = max(coshMaxDistance, upperRightDistance);

		assert(coshMaxDistance >= 1);
		assert(coshMinDistance >= 1);

		if (phi != rightAngle) {
			b = sinhr*cosDiffRight;
			extremum = log((coshr+b)/(coshr-b))/2;
			if (extremum < maxRHyper && extremum >= minRHyper) {
				double extremeDistance = cosh(extremum)*coshr-sinh(extremum)*sinhr*cosDiffRight;
				coshMinDistance = min(coshMinDistance, extremeDistance);
				coshMaxDistance = max(coshMaxDistance, extremeDistance);
			}
		} else {
			coshMinDistance = 1;
		}

		assert(coshMaxDistance >= 1);
		assert(coshMinDistance >= 1);

		//upper and lower borders
		if (phi >= leftAngle && phi < rightAngle) {
			double lower = cosh(abs(r_h-minRHyper));
			double upper = cosh(abs(r_h-maxRHyper));
			coshMinDistance = min(coshMinDistance, lower);
			coshMinDistance = min(coshMinDistance, upper);
			coshMaxDistance = max(coshMaxDistance, upper);
			coshMaxDistance = max(coshMaxDistance, lower);
		}

		assert(coshMaxDistance >= 1);
		assert(coshMinDistance >= 1);

		//again with mirrored phi
		double mirrorphi;
		if (phi >= M_PI) mirrorphi = phi - M_PI;
		else mirrorphi = phi + M_PI;
		if (mirrorphi >= leftAngle && mirrorphi < rightAngle) {
			double lower = coshMinR*coshr+sinhMinR*sinhr;
			double upper = coshMaxR*coshr+sinhMaxR*sinhr;
			coshMinDistance = min(coshMinDistance, lower);
			coshMinDistance = min(coshMinDistance, upper);
			coshMaxDistance = max(coshMaxDistance, upper);
			coshMaxDistance = max(coshMaxDistance, lower);
		}

		assert(coshMaxDistance >= 1);
		assert(coshMinDistance >= 1);

		double minDistance, maxDistance;
		minDistance = acosh(coshMinDistance);
		maxDistance = acosh(coshMaxDistance);
		assert(maxDistance >= 0);
		assert(minDistance >= 0);
		return std::pair<double, double>(minDistance, maxDistance);
	}

	/**
	 * Get all Elements in this QuadNode or a descendant of it
	 *
	 * @return vector of content type T
	 */
	std::vector<T> getElements() const {
		if (isLeaf) {
			return content;
		} else {
			assert(content.size() == 0);
			assert(positions.size() == 0);

			std::vector<T> result;
			for (index i = 0; i < children.size(); i++) {
				std::vector<T> subresult = children[i]->getElements();
				result.insert(result.end(), subresult.begin(), subresult.end());
			}
			return result;
		}
	}

	count height() const {
		count result = 1;//if leaf node, the children loop will not execute
		for (auto child : children) result = std::max(result, child->height()+1);
		return result;
	}

	/**
	 * Leaf cells in the subtree hanging from this QuadNode
	 */
	count countLeaves() const {
		if (isLeaf) return 1;
		count result = 0;
		for (index i = 0; i < children.size(); i++) {
			result += children[i]->countLeaves();
		}
		return result;
	}

	index getID() const {
			return ID;
		}


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

protected:
	Point<double> minCoords;
	Point<double> maxCoords;
	std::vector<T> content;
	std::vector<Point<double> > positions;
	std::vector<std::shared_ptr<SpatialCell<T> > > children;
	bool isLeaf;
	count capacity;
	index subTreeSize;
	index ID;

private:
	static const unsigned coarsenLimit = 4;
};

} /* namespace NetworKit */
#endif /* SPATIALCELL_H_ */
