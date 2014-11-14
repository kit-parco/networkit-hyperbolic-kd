/*
 * Smoother.h
 *
 *  Created on: 31.10.2014
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#ifndef SMOOTHER_H_
#define SMOOTHER_H_

#include "../algebraic/Matrix.h"
#include "../algebraic/Vector.h"

#include <limits>

namespace NetworKit {

class Smoother {
public:
	Smoother() {}
	virtual ~Smoother(){}

	virtual Vector relax(const Matrix &A, const Vector &b, const Vector &initialGuess, const count maxIterations = std::numeric_limits<count>::max()) const = 0;
};

} /* namespace NetworKit */

#endif /* SMOOTHER_H_ */
