/*
 * HyperbolicGraphReader.cpp
 *
 *  Created on: 29.11.2016
 *      Author: moritzl
 */

#include "HyperbolicGraphReader.h"
#include "EdgeListReader.h"
#include "LineFileReader.h"

namespace NetworKit {
	void HyperbolicGraphReader::readGraph(std::string filenamePrefix, Graph &G, std::vector<double> &angleOutput, std::vector<double> &radiiOutput) {

		/**
		 * read coordinates and parameters
		 */
		std::vector<std::string> coordinates =  LineFileReader::read(filenamePrefix+"-coordinates.txt");
		if (coordinates.size() == 0) {
			throw std::runtime_error("File "+filenamePrefix+"-coordinates.txt"+" not found or empty.");
		}

		while (std::string("").compare(coordinates[coordinates.size() -1]) == 0) {
			coordinates.pop_back();
		}

		//get parameters

		std::stringstream stream(coordinates[1]);
		std::string item;
		const char delim = '\t';

		std::getline(stream, item, delim);
		const count n = std::stoul(item);

		std::getline(stream, item, delim);
		const double R = std::stod(item);

		if (coordinates.size() != n+2) {
			std::string lastElement = coordinates[coordinates.size() -1];
			throw std::runtime_error("Coordinate vector has " + std::to_string(coordinates.size()) + " elements, should have " + std::to_string(n+2) + ". Last element is " + lastElement);
		}

		for (index i = 2; i < coordinates.size(); i++) {
			std::stringstream ss(coordinates[i]);

			std::getline(ss, item, delim);
			index node = std::stoul(item);
			if (node != i-2) {
				throw std::runtime_error("Node should be " + std::to_string(i-2) + ", is " + std::to_string(node));
			}

			std::getline(ss, item, delim);
			double radius = std::stod(item);
			if (radius < 0 || radius > R) {
				throw std::runtime_error("Read radius " + std::to_string(radius) + ". (R is " + std::to_string(R) + ")");
			}
			radiiOutput.push_back(radius);

			std::getline(ss, item, delim);
			double angle = (2*M_PI) * std::stod(item) / 360;
			if (angle < 0 || angle > 2*M_PI) {
				throw std::runtime_error("Read illegal angle: " + std::to_string(angle));
			}
			angleOutput.push_back(angle);
		}

		/**
		 * read links
		 */

		EdgeListReader reader('\t', 0);
		G = reader.read(filenamePrefix+"-links.txt");
		if (G.numberOfNodes() > n) {
			throw std::runtime_error("Graph has " + std::to_string(G.numberOfNodes()) + " nodes, but only " + std::to_string(n) + " coordinates were given.");
		}
		while (G.numberOfNodes() < n) {
			INFO("Adding nodes for additional coordinates.");
			G.addNode();
		}

	}
} /* namespace NetworKit */
