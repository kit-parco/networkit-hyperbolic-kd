/*
 * SpatialTree.h
 *
 *  Created on: 10.11.2016
 *      Author: moritzl
 */

#ifndef SPATIALTREE_H_
#define SPATIALTREE_H_



namespace NetworKit {

template <class T>
class SpatialTree {
public:
	SpatialTree() = default;
	virtual ~SpatialTree() = default;

	//SpatialTree(const Point<double> &minCoords, const Point<double> &maxCoords, count capacity=1000) = 0;

	void addContent(T content, const Point<double> &coords) {
		root->addContent(content, coords);
	}

	void getElementsInCircle(const Point<double> query, const double radius, vector<T> &circleDenizens) const {
		root->getElementsInCircle(query, radius, circleDenizens);
	}

	count getElementsProbabilistically(Point<double> query, std::function<double(double)> prob, vector<T> &circleDenizens) {
		return root->getElementsProbabilistically(query, prob, circleDenizens);
	}

	count size() const {
		return root->size();
	}

	count height() const {
		return root->height();
	}

	void trim() {
		root->trim();
	}

protected:
	std::shared_ptr<SpatialCell<T>> root;
};

} /* namespace NetworKit */
#endif /* SPATIALTREE_H_ */
