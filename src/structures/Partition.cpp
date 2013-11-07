/*
 * Partition.cpp
 *
 *  Created on: 03.10.2013
 *      Author: cls
 */

#include "Partition.h"

namespace NetworKit {

Partition::Partition(index z) : z(z), omega(0), data(z+1, none) {

}

void Partition::addToSubset(index s, index e) {
	assert (data[e] == none);	// guarantee that element was unasigned
	assert (s <= omega); // do not create new subset ids
	data[e] = s;
}

void Partition::moveToSubset(index s, index e) {
	assert (s <= omega); // do not create new subset ids
	data[e] = s;
}

void Partition::toSingleton(index e) {
	data[e] = newSubsetId();
}

void Partition::allToSingletons() {
	for (index e = 0; e <= this->z; ++e) {
		toSingleton(e);
	}
}

void Partition::mergeSubets(index s, index t) {
	assert (s <= omega);
	assert (t <= omega);
	index m = newSubsetId(); // new id for merged set
	for (index e = 0; e <= this->z; ++e) {
		if (data[e] == t || data[e] == t) {
			data[e] = m;
		}
	}
}

bool Partition::isOnePartition(const std::set<index>& elements) {
	index one = data[0];	// first subset id should be equal to all others
	for (index e = 0; e <= this->z; ++e) {
		if (data[e] != one) {
			return false;
		}
	}
	return true;
}

bool Partition::isSingletonPartition(const std::set<index>& elements) const {
	return (numberOfElements() == numberOfSubsets());
}

count Partition::numberOfSubsets() const {
	// TODO:
}

index Partition::upperBound() const {
	return omega + 1;
}

index Partition::lowerBound() const {
	return 0;
}

void Partition::compact() {
	// TODO:
}

bool Partition::contains(index e) const {
	return (e <= z) && (data[e] != none);	// e is in the element index range and the entry is not empty
}

bool Partition::inSameSubset(index e1, index e2) const {
	assert (data[e1] != none);
	assert (data[e2] != none);
	return (data[e1] == data[e2]);
}

bool Partition::equals(const Partition& other,
		const std::set<index>& elements) const {
}

std::vector<count> Partition::subsetSizes() const {
	// TODO:
}

std::map<index, count> Partition::subsetSizeMap() const {
	// TODO:
}

std::set<index> Partition::getMembers(const index s) const {
	assert (s <= omega);
	std::set<index> subset;
	for (index e = 0; e <= this->z; ++e) {
		if (data[e] == s) {
			subset.insert(e);
		}
	}
	return subset;
}

count Partition::numberOfElements() const {
	return z;	// z is the maximum element id
}


} /* namespace NetworKit */
