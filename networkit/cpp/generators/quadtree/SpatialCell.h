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

#include "../../viz/Point.h"

namespace NetworKit {

template <class T>
class SpatialCell {
public:
	SpatialCell() = default;
	virtual ~SpatialCell() = default;
	SpatialCell(const Point<double> &minCoords, const Point<double> &maxCoords, count capacity=1000) : minCoords(minCoords), maxCoords(maxCoords), capacity(capacity) {
		isLeaf = true;
		subTreeSize = 0;
	}

	virtual void split() = 0;

	virtual std::pair<double, double> distances(const Point<double> &query) const = 0;

	virtual double distance(const Point<double> &query, index k) const = 0;

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

	virtual count getElementsProbabilistically(const Point<double> &query, std::function<double(double)> prob, vector<T> &result) const {
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

protected:
	Point<double> minCoords;
	Point<double> maxCoords;
	std::vector<T> content;
	std::vector<Point<double> > positions;
	std::vector<std::shared_ptr<SpatialCell<T> > > children;
	bool isLeaf;
	count capacity;
	index subTreeSize;
};

} /* namespace NetworKit */
#endif /* SPATIALCELL_H_ */
