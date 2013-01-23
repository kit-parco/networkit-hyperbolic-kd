/*
 * METISParser.cpp
 *
 *  Created on: 27.11.2012
 *      Author: cls
 */

#include "METISParser.h"


namespace EnsembleClustering {


/**
 * Extract a vector of indices from a line in the file.
 *
 * @param[in]	line		line from input file containing node indices
 *
 * @param[out]	indices		node indices extracted from line
 */
static std::vector<node> parseLine(std::string line) {

	std::stringstream stream(line);
	std::string token;
	char delim = ' ';
	std::vector<node> adjacencies;

	// split string and push adjacent nodes
	while (std::getline(stream, token, delim)) {
		node v = atoi(token.c_str());
		adjacencies.push_back(v);
	}

	TRACE("line parsed");
	return adjacencies;
}


METISParser::METISParser(std::string path) : graphFile(path) {
	if (!(this->graphFile)) {
		ERROR("invalid graph file: " << path);
		throw std::runtime_error("invalid graph file");
	}
}



METISParser::~METISParser() {
	// TODO Auto-generated destructor stub
}




std::pair<int, int> METISParser::getHeader() {

	// handle header line
	int n;  // number of nodes
	int m;	// number of edges

	std::string line = "";
	assert (this->graphFile);
	if (std::getline(this->graphFile, line)) {
		std::vector<node> tokens = parseLine(line);
		n = tokens[0];
		m = tokens[1];

		TRACE("n = " << n << " m = " << m );

		return std::make_pair(n, m);
	} else {
		ERROR("getline not successful");
	}

}


bool METISParser::hasNext() {
	// if graph file has lines left, return true
	return this->graphFile.good();
}


std::vector<node> METISParser::getNext() {

	std::string line;
	bool comment = false;
	do {
		comment = false;
		std::getline(this->graphFile, line);
		TRACE("reading line: " << line);
		// check for comment line starting with '%'
		if (line[0] == '%') {
			comment = true;
			TRACE("comment line found");
		} else {
			return parseLine(line);
		}

	} while (comment);
}




} /* namespace EnsembleClustering */
