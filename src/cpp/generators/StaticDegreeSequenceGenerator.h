/*
 * StaticDegreeSequenceGenerator.h
 *
 *  Created on: 24.02.2014
 *      Author: Henning
 */

#ifndef STATICDEGREESEQUENCEGENERATOR_H_
#define STATICDEGREESEQUENCEGENERATOR_H_

#include "StaticGraphGenerator.h"

namespace NetworKit {

const short NO = 0;
const short YES = 1;
const short UNKNOWN = 2;

class StaticDegreeSequenceGenerator: public NetworKit::StaticGraphGenerator {
protected:
	std::vector<unsigned int> seq;
	short realizable;


public:
	StaticDegreeSequenceGenerator(const std::vector<unsigned int>& sequence);

	/**
	 * Erdoes-Gallai test if degree sequence seq is realizable.
	 */
	virtual bool isRealizable();

	virtual bool getRealizable() const;


	virtual Graph generate() = 0;
};

} /* namespace NetworKit */
#endif /* STATICDEGREESEQUENCEGENERATOR_H_ */
