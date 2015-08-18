#include "Algorithm.h"
#include "../auxiliary/SignalHandling.h"
#include "../auxiliary/Log.h"
#include <exception>

namespace NetworKit {
	Algorithm::Algorithm() : hasRun(false) {

	}

	bool Algorithm::hasFinished() const {
		return hasRun;
	}

	std::string Algorithm::toString() const {
		return "Algorithm base class";
	}

	bool Algorithm::isParallel() const {
		throw std::runtime_error("TODO: Implement in subclass");
		return false;
	}

} /* NetworKit */
