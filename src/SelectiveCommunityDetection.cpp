//============================================================================
// Name        : SelectiveCommunityDetection.cpp
// Author      : Christian Staudt (christian.staudt@kit.edu),
//				 Henning Meyerhenke (henning.meyerhenke@kit.edu)
// Version     :
// Copyright   : � 2013, Christian Staudt, Henning Meyerhenke
// Description : Framework for engineering selective community detection algorithms
//============================================================================

// includes


// includes
#include <iostream>
#include <fstream>
#include <utility>
#include <cfenv>	// floating point exceptions
#include <stdexcept>

// log4cxx
#ifndef NOLOG4CXX
#include "log4cxx/logger.h"
#include "log4cxx/basicconfigurator.h"
#endif

// GoogleTest
#ifndef NOGTEST
#include <gtest/gtest.h>
#endif

// OpenMP
#include <omp.h>


// EnsembleClustering
#include "Globals.h"
#include "ext/optionparser.h"
#include "auxiliary/Log.h"
#include "auxiliary/Timer.h"
#include "auxiliary/Functions.h"
#include "auxiliary/StringTools.h"
#include "graph/Graph.h"



using namespace NetworKit;






#ifndef NOLOGGING
#ifndef NOLOG4CXX
/**
 * Call this first to configure logging output.
 */
void configureLogging(const std::string& loglevel = "INFO") {
	// configure logging
	log4cxx::BasicConfigurator::configure();
	if (loglevel == "TRACE") {
		log4cxx::Logger::getRootLogger()->setLevel(log4cxx::Level::getTrace());
	} else if (loglevel == "DEBUG") {
		log4cxx::Logger::getRootLogger()->setLevel(log4cxx::Level::getDebug());
	} else if (loglevel == "INFO") {
		log4cxx::Logger::getRootLogger()->setLevel(log4cxx::Level::getInfo());
	} else if (loglevel == "WARN") {
		log4cxx::Logger::getRootLogger()->setLevel(log4cxx::Level::getWarn());
	} else if (loglevel == "ERROR") {
		log4cxx::Logger::getRootLogger()->setLevel(log4cxx::Level::getError());
	} else {
		ERROR("unknown loglevel: " << loglevel);
		exit(1);
	}
}
#endif
#endif


/**
 *  Set the number of threads available to the program.
 */
void setNumberOfThreads(int nThreads) {
#ifdef _OPENMP
		omp_set_num_threads(nThreads);
#else
		WARN("Thread option ignored since OpenMP is deactivated.");
#endif
}



// *** Option Parser Configuration ***//

class Arg: public OptionParser::Arg {

static OptionParser::ArgStatus Required(const OptionParser::Option& option, bool msg)
{
  if (option.arg != 0)
    return OptionParser::ARG_OK;

  if (msg) {
	  std::cout << "Option '" << option << "' requires an argument" << std::endl;
  }
  return OptionParser::ARG_ILLEGAL;
}

};

// TODO: clean up obsolete parameters
enum  optionIndex { UNKNOWN, HELP, LOGLEVEL, THREADS, TESTS, ALGORITHM, RUNS, SAVE_GRAPH, PROGRESS, SUMMARY, SCALETHREADS, UPDATE_THRESHOLD};
const OptionParser::Descriptor usage[] =
{
 {UNKNOWN, 0,"" , ""    ,OptionParser::Arg::None, "USAGE: EnsembleClustering [options]\n\n"
                                            "Options:" },
 {HELP,    0,"h" , "help",OptionParser::Arg::None, "  --help  \t Print usage and exit." },
 {LOGLEVEL,    0, "" , "loglevel", OptionParser::Arg::Required, "  --loglevel=<LEVEL>  \t set the log level" },
 {THREADS,    0, "" , "threads", OptionParser::Arg::Required, "  --threads=<NUM>  \t set the maximum number of threads" },
 {TESTS, 0, "t", "tests", OptionParser::Arg::None, "  --tests \t Run unit tests"},
 {ALGORITHM, 0, "", "algorithm", OptionParser::Arg::Required, "  --algorithm=<ALGORITHM>:<PARAMS> \t select clustering algorithm"},
 {RUNS, 0, "", "runs", OptionParser::Arg::Required, "  --runs=<NUMBER> \t set number of clusterer runs"},
 {SAVE_GRAPH, 0, "", "saveGraph", OptionParser::Arg::Required, "  --saveGraph=<PATH> \t write the graph to a file"},
 {PROGRESS, 0, "", "progress", OptionParser::Arg::None, "  --progress \t print progress bar"},
 {SUMMARY, 0, "", "summary", OptionParser::Arg::Required, "  --summary=<PATH> \t append summary as a .csv line to this file"},
 {SCALETHREADS, 0, "", "scaleThreads", OptionParser::Arg::Required, "  --scaleThreads=<MAXTHREADS> \t scale number of threads by factor 2 until maximum is reached"},
 {UNKNOWN, 0,"" ,  ""   ,OptionParser::Arg::None, "\nExamples:\n"
                                            " TODO" },
 {0,0,0,0,0,0}
};


// MAIN FUNCTIONS






int main(int argc, char **argv) {
	std::cout << "=== NetworKit - Selective Community Detection ===" << std::endl;

	// ENABLE FLOATING POINT EXCEPTIONS (needs GNU extension, apparently only available on Linux)
#ifdef _GNU_SOURCE
	// feenableexcept(FE_ALL_EXCEPT);
#endif


	/// PARSE OPTIONS

	argc-=(argc>0); argv+=(argc>0); // skip program name argv[0] if present
	OptionParser::Stats  stats(usage, argc, argv);
	OptionParser::Option options[stats.options_max], buffer[stats.buffer_max];
	OptionParser::Parser parse(usage, argc, argv, options, buffer);

	if (parse.error())
	 return 1;

	if (options[HELP] || argc == 0) {
	 OptionParser::printUsage(std::cout, usage);
	 return 0;
	}

	for (OptionParser::Option* opt = options[UNKNOWN]; opt; opt = opt->next())
	 std::cout << "Unknown option: " << opt->name << "\n";

	for (int i = 0; i < parse.nonOptionsCount(); ++i)
	 std::cout << "Non-option #" << i << ": " << parse.nonOption(i) << "\n";




	// CONFIGURE LOGGING


#ifndef NOLOGGING
#ifndef NOLOG4CXX
	if (options[LOGLEVEL]) {
		configureLogging(options[LOGLEVEL].arg);
	} else {
		configureLogging();	// with default level
	}
#endif
#endif


	// CONFIGURE PARALLELISM

#ifdef _OPENMP
	omp_set_nested(1); // enable nested parallelism
#endif

	if (options[THREADS]) {
		// set number of threads
		int nThreads = std::atoi(options[THREADS].arg);
		setNumberOfThreads(nThreads);
	}

	// CONFIGURE OUTPUT
	if (options[PROGRESS]) {
		PRINT_PROGRESS = true;
	} else {
		PRINT_PROGRESS = false;
	}


	// TODO: how to do summary?
	if (options[SUMMARY]) {
		// check if file exists by trying to read from it
		std::ifstream summary(options[SUMMARY].arg);
		if (summary) {
			// file exists
		} else {
			// create it and write CSV header
			std::ofstream summary(options[SUMMARY].arg);
			summary << "threads;algo;revision;graph;running;eps;clusters;mod" << std::endl; // TODO: update header
		}
	}



	// RUN UNIT TESTS

	if (options[TESTS]) {

#ifndef NOGTEST
	   ::testing::InitGoogleTest(&argc, argv);
	   INFO("=== starting unit tests ===");
	   return RUN_ALL_TESTS();
#else
	   throw std::runtime_error("unit tests are excluded from build by the NOGTEST preprocessor directive");
#endif
	}


	// SET NUMBER OF CLUSTERER RUNS
	int runs = 1;
	if (options[RUNS]) {
		runs = std::atoi(options[RUNS].arg);
	}



	// RUN PROGRAM
	// TODO: get source

	// allow for scripted thread scaling
	if (options[SCALETHREADS]) {
		// perform scaling
		int maxThreads = std::atoi(options[SCALETHREADS].arg);
		for (int nThreads = 1; nThreads <= maxThreads; nThreads *= 2) {
			setNumberOfThreads(nThreads);
			// allow for multiple runs
			for (int run = 0; run < runs; run++) {
				// TODO: run
				// TODO: inspect(G, clustering, options);
			}

		}
	} else {
		// allow for multiple runs
		for (int run = 0; run < runs; run++) {
			// TODO: run
			// TODO: inspect
		}
	}

	std::cout << "[EXIT] terminated normally" << std::endl;
	return 0;
}





