/*
 * Timer.cpp
 *
 *  Created on: 14.01.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include "Timer.h"

namespace Aux {

Timer::Timer() : running(false) {
}

std::chrono::steady_clock::time_point Timer::start() {
	this->started = std::chrono::steady_clock::now();
	running = true;
	return this->started;
}

std::chrono::steady_clock::time_point Timer::stop() {
	this->stopped = std::chrono::steady_clock::now();
	running = false;
	return this->stopped;
}

std::chrono::duration<uint64_t, std::milli> Timer::elapsed() {
	if (running) {
		return std::chrono::duration_cast<std::chrono::duration<uint64_t, std::milli>>(std::chrono::steady_clock::now() - this->started);
	}

	std::chrono::duration<uint64_t, std::milli> elapsed = std::chrono::duration_cast<std::chrono::duration<uint64_t, std::milli>>(this->stopped - this->started);
	return elapsed;
}

std::chrono::steady_clock::time_point Timer::startTime() {
	return this->started;
}

std::chrono::steady_clock::time_point Timer::stopTime() {
	return this->stopped;
}

uint64_t Timer::elapsedMilliseconds() {
	return this->elapsed().count();
}

uint64_t Timer::elapsedMicroseconds() {
	return std::chrono::duration_cast<std::chrono::duration<uint64_t, std::micro>>(this->stopped - this->started).count();
}

uint64_t Timer::elapsedNanoseconds() {
	return std::chrono::duration_cast<std::chrono::duration<uint64_t, std::nano>>(this->stopped - this->started).count();
}

std::string Timer::elapsedTag() {
	std::stringstream s;
	s << "(" << this->elapsed().count() << " ms) ";
	return s.str();
}


} /* namespace Aux */

