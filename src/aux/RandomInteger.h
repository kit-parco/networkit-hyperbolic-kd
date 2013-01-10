/*
 * RandomInteger.h
 *
 *  Created on: 10.01.2013
 *      Author: cls
 */

#ifndef RANDOMINTEGER_H_
#define RANDOMINTEGER_H_

#include <random>
#include <ctime>

namespace EnsembleClustering {

class RandomInteger {

protected:

	std::default_random_engine randomEngine;
	std::uniform_int_distribution<> distribution;

public:

	RandomInteger(int64_t lower, int64_t upper);

	virtual ~RandomInteger();

	virtual int64_t generate();
};

} /* namespace EnsembleClustering */
#endif /* RANDOMINTEGER_H_ */
