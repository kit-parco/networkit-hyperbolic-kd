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
		throw std::runtime_error("TODO: implement in subclass and return string representation");
	}

} /* NetworKit */
