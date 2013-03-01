//============================================================================
// Name        : EnsembleClustering.cpp
// Author      : Christian Staudt
// Version     :
// Copyright   : � 2012, Christian Staudt
// Description : Combining parallel clustering and ensemble learning
//============================================================================

// includes
#include <iostream>
#include <utility>
#include <cfenv>	// floating point exceptions

// log4cxx
#ifndef NOLOG4CXX
#include "log4cxx/logger.h"
#include "log4cxx/basicconfigurator.h"
#endif

// GoogleTest
#include <gtest/gtest.h>

// OpenMP
#include <omp.h>


// EnsembleClustering
#include "Globals.h"
#include "aux/optionparser.h"
#include "aux/Log.h"
#include "aux/Timer.h"
#include "aux/Functions.h"
#include "aux/StringTools.h"
#include "graph/Graph.h"
#include "graph/GraphGenerator.h"
#include "ensemble/EnsembleClusterer.h"
#include "ensemble/EnsemblePreprocessing.h"
#include "clustering/algo/LabelPropagation.h"
#include "clustering/algo/ParallelAgglomerativeClusterer.h"
#include "clustering/algo/RandomClusterer.h"
#include "clustering/algo/Louvain.h"
#include "clustering/algo/LouvainParallel.h"
#include "clustering/base/Clustering.h"
#include "clustering/base/Modularity.h"
#include "clustering/base/Coverage.h"
#include "clustering/base/ClusteringGenerator.h"
#include "io/METISGraphReader.h"
#include "io/METISGraphWriter.h"
#include "io/ClusteringWriter.h"


// revision
static const std::string REVISION = "r007";


using namespace EnsembleClustering;



// graph functions

Graph generateRandomGraph(count n, double p) {

	std::cout << "[BEGIN] generating random graph..." << std::flush;
	Aux::Timer running;
	running.start();
	GraphGenerator graphGen;
	Graph G = graphGen.makeRandomGraph(n, p);
	running.stop();
	std::cout << "[DONE] (" << running.elapsed().count() << " ms)" << std::endl;
	return G;
}

Graph generateClusteredRandomGraph(count n, count k, double pin, double pout) {

	// prepare generated clustered random graph (planted partition)
	std::cout << "[BEGIN] generating random clustering..." << std::flush;
	Graph emptyG(n);	// clustering generator needs a dummy graph
	ClusteringGenerator clusteringGen;
	Clustering planted = clusteringGen.makeRandomClustering(emptyG, k);
	std::cout << "[DONE]" << std::endl;

	std::cout << "[BEGIN] generating clustered random graph..." << std::flush;
	GraphGenerator graphGen;

	Aux::Timer running;
	running.start();
	//
	Graph G = graphGen.makeClusteredRandomGraph(planted, pin, pout);
	//
	running.stop();
	std::cout << "[DONE] (" << running.elapsed().count() << " ms)" << std::endl;

	return G;
}

Graph generatePreferentialAttachmentGraph(count n, count a) {

	std::cout << "[BEGIN] generating preferential attachment graph..." << std::flush;
	GraphGenerator graphGen;

	Aux::Timer running;
	running.start();
	//
	Graph G = graphGen.makePreferentialAttachmentGraph(n, a);
	//
	running.stop();
	std::cout << "[DONE] (" << running.elapsed().count() << " ms)" << std::endl;

	return G;
}


/**
 * Read a graph from a file.
 */
Graph readGraph(const std::string& graphPath) {

	// READ GRAPH

	METISGraphReader reader;	// TODO: add support for multiple graph file formats

	// TIMING
	Aux::Timer readTimer;
	readTimer.start();
	//
	std::cout << "[BEGIN] reading file: " << graphPath << std::endl;

	Graph G = reader.read(graphPath);
	//
	readTimer.stop();
	std::cout << "[DONE] read graph file " << readTimer.elapsedTag() << std::endl;
	// TIMING

	return G;

}




/**
 * Split a string at delimiter and return vector of parts.
 *
 */
std::vector<std::string> splitAString(const std::string& s, char delim = ' ') {
	std::stringstream stream(s);
	std::string token;
	std::vector<std::string> tokens;

	// split string and push adjacent nodes
	while (std::getline(stream, token, delim)) {
		tokens.push_back(token);
	}

	return tokens;
}


/**
 * Set LabelPropagation parameter: randOrder
 */
void configureRandOrder(const std::string& isOrderRandomized) {
	if (isOrderRandomized == "no") {
		RAND_ORDER = false;
	}
	else if (isOrderRandomized == "yes") {
		RAND_ORDER = true;
	} else {
		// TODO: print warning?
	}
}


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


enum  optionIndex { UNKNOWN, HELP, LOGLEVEL, THREADS, TESTS, GRAPH, GENERATE, ALGORITHM, RUNS,
	ENSEMBLE, SOLO, NOREC, NORM_VOTES, SCALESTRENGTH,
	WRITEGRAPH, SAVE_CLUSTERING, PROGRESS, SUMMARY, RANDORDER, UPDATE_THRESHOLD, OVERLAP, DISSIMILARITY};
const OptionParser::Descriptor usage[] =
{
 {UNKNOWN, 0,"" , ""    ,OptionParser::Arg::None, "USAGE: EnsembleClustering [options]\n\n"
                                            "Options:" },
 {HELP,    0,"h" , "help",OptionParser::Arg::None, "  --help  \t Print usage and exit." },
 {LOGLEVEL,    0, "" , "loglevel", OptionParser::Arg::Required, "  --loglevel  \t set the log level" },
 {THREADS,    0, "" , "threads", OptionParser::Arg::Required, "  --threads=<NUM>  \t set the maximum number of threads" },
 {TESTS, 0, "t", "tests", OptionParser::Arg::None, "  --tests \t Run unit tests"},
 {GRAPH, 0, "g", "graph", OptionParser::Arg::Required, "  --graph \t Run ensemble clusterer on graph"},
 {GENERATE, 0, "", "generate", OptionParser::Arg::Required, "  --generate \t Run ensemble clusterer on generated graph with planted partition"},
 {ALGORITHM, 0, "", "algorithm", OptionParser::Arg::Required, "  --algorithm=<NAME>:<PARAMS> \t select clustering algorithm"},
 {RUNS, 0, "", "runs", OptionParser::Arg::Required, "  --runs=<NUMBER> \t set number of clusterer runs"},
 {ENSEMBLE, 0, "", "ensemble", OptionParser::Arg::Required, "  --ensemble=<b>*<BASE>+<FINAL> \t <b>: number of base clusterers in the ensemble, <BASE>: base clusterer name, <FINAL>: final clusterer name"},
 {SOLO, 0, "", "solo", OptionParser::Arg::Required, "  --solo=<Algorithm> \t run only a single base algorithm"},
 {NOREC, 0, "", "noRecursion", OptionParser::Arg::None, "  --noRecursion \t run only on the finest graph, even if ensemble is used"},
 {NORM_VOTES, 0, "", "normalizeVotes", OptionParser::Arg::None, "  --normalizeVotes \t normalize votes in label propagation by weighted degree"},
 {SCALESTRENGTH, 0, "", "scaleStrength", OptionParser::Arg::Required, "  --scaleStrength=<value in [0,1]> \t scale cluster strengths"},
 {WRITEGRAPH, 0, "", "writeGraph", OptionParser::Arg::Required, "  --writegraph=<PATH> \t write the graph to a file"},
 {SAVE_CLUSTERING, 0, "", "saveClustering", OptionParser::Arg::Required, "  --saveClustering=<PATH> \t save the clustering to a file"},
 {PROGRESS, 0, "", "progress", OptionParser::Arg::None, "  --progress \t print progress bar"},
 {SUMMARY, 0, "", "summary", OptionParser::Arg::Required, "  --summary=<PATH> \t append summary as a .csv line to this file"},
 {RANDORDER, 0, "", "randOrder", OptionParser::Arg::Required, "  --randOrder=<yes,no> \t don't randomize vertex processing order"},
 {UPDATE_THRESHOLD, 0, "", "updateThreshold", OptionParser::Arg::Required, "  --updateThreshold=<N> or --updateThreshold=auto \t number of updated nodes below which label propagation terminates - auto determines this automatically from the size of the graph"},
 {OVERLAP, 0, "", "overlap", OptionParser::Arg::Required, "  --overlap=<Algorithm> set overlap algorithm which combines the base clusterings"},
 {DISSIMILARITY, 0, "", "dissimilarity", OptionParser::Arg::None, "  --dissimilarity \t calculate and print base clustering dissimilarity for ensemble (expensive!)"},
 {UNKNOWN, 0,"" ,  ""   ,OptionParser::Arg::None, "\nExamples:\n"
                                            " TODO" },
 {0,0,0,0,0,0}
};


// MAIN FUNCTIONS


Graph getGraph(OptionParser::Option* options) {

	// generated graph
	if (options[GENERATE]) {
		std::string genOption = options[GENERATE].arg;
		std::cout << "\t --generate=" << options[GENERATE].arg << std::endl;

		// get graph generation model
		std::vector<std::string> genArgs = splitAString(genOption, ':');
		std::vector<std::string> genNumArgs = splitAString(genArgs.at(1), ',');
		assert (genArgs.size() == 2);
		std::string model = genArgs.at(0);
		if (model == "RG") { // random graph (Erdos-Renyi)
			assert (genNumArgs.size() == 2);
			count n = std::atoi(genNumArgs.at(0).c_str());
			double p = std::atoi(genNumArgs.at(1).c_str());
			return generateRandomGraph(n, p);
		} else if (model == "CR") {	// clustered random graph
			assert (genNumArgs.size() == 4);
			count n = std::atoi(genNumArgs.at(0).c_str());
			count k = std::atoi(genNumArgs.at(1).c_str());
			double pin = std::atof(genNumArgs.at(2).c_str());
			double pout = std::atof(genNumArgs.at(3).c_str());
			return generateClusteredRandomGraph(n, k, pin, pout);
		} else if (model == "PA") {	// preferential attachment (Barabasi-Albert)
			assert (genNumArgs.size() == 2);
			count n = std::atoi(genNumArgs.at(0).c_str());
			count a = std::atoi(genNumArgs.at(1).c_str());
			return generatePreferentialAttachmentGraph(n, a);
		} else {
			std::cout << "[ERROR] unknown graph generation model: " << model << " [EXIT]" << std::endl;
			exit(1);
		}


	} else if (options[GRAPH]) { // graph from file
		std::string graphPath = options[GRAPH].arg;
		std::cout << "\t --graph=" << graphPath << std::endl;

		Graph G = readGraph(graphPath);
		return G;
	} else {
		Graph G(0);	// return empty graph
		G.setName("NONE");
		return G;
	}
}


Clustering startClusterer(Graph& G, OptionParser::Option* options) {

	// if getGraph returns empty graph, abort
	if (G.isEmpty() && (G.getName() == "NONE")) {
		std::cout << "[ERROR]�no graph instance provided." << std::endl;
		std::cout << "[EXIT]" << std::endl;
		exit(1);
	}


	if (options[WRITEGRAPH]) {
		// write input graph to file
		std::cout << "\t --writegraph=" << options[WRITEGRAPH].arg << std::endl;

		std::string path = options[WRITEGRAPH].arg;
		METISGraphWriter writer;
		std::cout << "[BEGIN] writing graph to file: " << path << std::endl;
		writer.write(G, path);
		std::cout << "[DONE]" << std::endl;

	}


	// determine update threshold / abort criterion for LabelPropagation
	count updateThreshold = 0;
	if (options[UPDATE_THRESHOLD]) {
		std::string updateThresholdArg = options[UPDATE_THRESHOLD].arg;
		if (updateThresholdArg == "auto") {
			updateThreshold = (count) (G.numberOfNodes() / 1e3);
		} else {
			updateThreshold = std::atoi(updateThresholdArg.c_str());
		}
	} else {
		// default = "auto"
		updateThreshold = (count) (G.numberOfNodes() / 1e3);
	}

	// prepare clusterer run

	Clusterer* algo = NULL; // the clusterer

//	std::pair<Clustering, Graph> result = std::make_pair(Clustering(0), G); // this will be returned
	Aux::Timer running; // measures running time of clusterer

	if (options[ALGORITHM]) {
		std::string algoArg = options[ALGORITHM].arg;
		std::string algoName = Aux::StringTools::split(algoArg, ':').front();
		std::string algoParams = Aux::StringTools::split(algoArg, ':').back();

		if (algoName == "LabelPropagation") {
			algo = new LabelPropagation(updateThreshold);
		} else if (algoName == "Agglomerative") {
			algo = new ParallelAgglomerativeClusterer();
		} else if (algoName == "RandomClusterer") {
			algo = new RandomClusterer();
		} else if (algoName == "Louvain") {
			// algoParams is parallelization strategy
			if (algoParams.empty()) {
				algo = new Louvain();
			} else {
				algo = new Louvain(algoParams);
			}
		} else if (algoName == "LouvainParallel") {
			algo = new LouvainParallel();
		} else if (algoName == "EnsembleMultiLevel") {
			// TODO:
		} else if (algoName == "EnsemblePreprocessing") {
			EnsemblePreprocessing* ensemblePre = new EnsemblePreprocessing();
			// parse params
			std::string ensembleFrontArg = Aux::StringTools::split(algoParams, '+').front();
			std::string finalClustererArg = Aux::StringTools::split(algoParams, '+').back();
			std::string ensembleSizeArg = Aux::StringTools::split(ensembleFrontArg, '*').front();
			std::string baseClustererArg = Aux::StringTools::split(ensembleFrontArg, '*').back();

			int ensembleSize = std::atoi(ensembleSizeArg.c_str());
			// 1. add base clusterers
			for (int i = 0; i < ensembleSize; i += 1) {
				Clusterer* base = NULL;
				if (baseClustererArg == "LabelPropagation") {
					base = new LabelPropagation(updateThreshold);
				} else if (baseClustererArg == "Agglomerative") {
					base = new ParallelAgglomerativeClusterer();
				} else {
					std::cout << "[ERROR]�unknown base clusterer: " << baseClustererArg << std::endl;
					exit(1);
				}
				ensemblePre->addBaseClusterer(*base);
			}
			// 2. overlap algorithm
			Overlapper* overlap = NULL;
			if (options[OVERLAP]) {
				std::string overlapArg = options[OVERLAP].arg;
				if (overlapArg == "Hashing") {
					overlap = new HashingOverlapper();
				} else if (overlapArg == "RegionGrowing") {
					overlap = new RegionGrowingOverlapper();
				} else {
					std::cout << "[ERROR] unknown overlap algorithm: " << overlapArg << std::endl;
					exit(1);
				}
			} else {
				// default
				overlap = new RegionGrowingOverlapper();
			}
			ensemblePre->setOverlapper(*overlap);
			// 3. Final Clusterer
			Clusterer* final = NULL;
			if (finalClustererArg == "LabelPropagation") {
				final = new LabelPropagation();
			} else if (finalClustererArg == "Agglomerative") {
				final = new ParallelAgglomerativeClusterer();
			} else if (finalClustererArg == "Louvain") {
				final = new Louvain();
			} else if (finalClustererArg == "LouvainParallel") {
				final = new LouvainParallel();
			} else {
				std::cout << "[ERROR] unknown final clusterer: " << finalClustererArg << std::endl;
				exit(1);
			}

			ensemblePre->setFinalClusterer(*final);

			algo = ensemblePre;

		} else {
			std::cout << "[ERROR] unknown algorithm: " << algoName << std::endl;
			std::cout << "[EXIT]" << std::endl;
			exit(1);
		}

		// START CLUSTERER
			std::cout << "[BEGIN] clusterer: " << algo->toString() << std::endl;
			running.start();
			Clustering resultClustering = algo->run(G);
			running.stop();
			//
			std::cout << "[DONE] " << algo->toString() << " ran: \t" << running.elapsedTag() << std::endl;


			// result = std::make_pair(resultClustering, G);

			// print speed info
			double eps = (G.numberOfEdges() / ((double) running.elapsed().count() / 1000.0));	// edges per second
			std::cout << "\t # edges per second:\t" << eps << std::endl;


		 	if (options[SUMMARY]) {
		 		char sep = ';';
		 		std::ofstream summary(options[SUMMARY].arg, std::ios::app); // open summary file to append to
		 		// APPEND number of threads available
		 		#ifdef _OPENMP
		 		summary << omp_get_max_threads() << sep;
		 		#else
		 		summary << "NoOpenMP" << sep;
		 		#endif
		 		summary << algo->toString() << sep;				// APPEND algorithm description
		 		summary << REVISION << sep;						// APPEND revision number
		 		summary << G.getName() << sep;					// APPEND graph description
		 		summary << running.elapsed().count() << sep;	// APPEND running time
		 		summary << eps << sep;							// APPEND edges per second
		 	}


		//	return result;	// return empty clustering
		 	return resultClustering;

	}

	// TODO: REMOVE OLD PARAMETERS!
	if (options[SOLO]) {
		std::cout << "[WARNING]�DEPRECATED: use --algorithm=<NAME>:<PARAMS> instead!" << std::endl;

		std::cout << "\t --solo=" << options[SOLO].arg << std::endl;
		// RUN ONLY SINGLE BASE ALGORITHM



		// get specified base algorithm
		std::string algoName = options[SOLO].arg;
		if (algoName == "LabelPropagation") {

			LabelPropagation* lp = new LabelPropagation(updateThreshold);

			algo = lp;
		} else if (algoName == "Agglomerative") {
			ParallelAgglomerativeClusterer* agglo = new ParallelAgglomerativeClusterer();
			algo = agglo;
		} else if (algoName == "RandomClusterer") {
			algo = new RandomClusterer();
		} else if (algoName == "Louvain") {
			algo = new Louvain();
		}else {
			std::cout << "[ERROR] unknown base algorithm: " << algoName << std::endl;
			std::cout << "[EXIT]" << std::endl;
			exit(1);
		}


	} else if (options[ENSEMBLE]) {
		// RUN ENSEMBLE

		std::cout << "[WARNING]�DEPRECATED: use --algorithm=<NAME>:<PARAMS> instead!" << std::endl;

		std::string ensembleArg = options[ENSEMBLE].arg;

		std::string ensembleFrontArg = Aux::StringTools::split(ensembleArg, '+').front();
		std::string finalClustererArg = Aux::StringTools::split(ensembleArg, '+').back();
		std::string ensembleSizeArg = Aux::StringTools::split(ensembleFrontArg, '*').front();
		std::string baseClustererArg = Aux::StringTools::split(ensembleFrontArg, '*').back();

		int ensembleSize = std::atoi(ensembleSizeArg.c_str());

		EnsembleClusterer* ensemble = new EnsembleClusterer();

		// CONFIGURE ENSEMBLE CLUSTERER

		// 1. Quality Measure
		QualityMeasure* qm = new Modularity();
		ensemble->setQualityMeasure(*qm);

		// 2. Base Clusterers
		for (int i = 0; i < ensembleSize; i += 1) {
			Clusterer* base = NULL;
			if (baseClustererArg == "LabelPropagation") {
				base = new LabelPropagation(updateThreshold);
			} else if (baseClustererArg == "Agglomerative") {
				base = new ParallelAgglomerativeClusterer();
			} else {
				std::cout << "[ERROR]�unknown base clusterer: " << baseClustererArg << std::endl;
				exit(1);
			}
			ensemble->addBaseClusterer(*base);
		}

		// 3. Overlap
		Overlapper* overlap = NULL;

		// select overlap algorithm
		if (options[OVERLAP]) {
			std::string overlapArg = options[OVERLAP].arg;
			if (overlapArg == "Hashing") {
				overlap = new HashingOverlapper();
			} else if (overlapArg == "RegionGrowing") {
				overlap = new RegionGrowingOverlapper();
			} else {
				std::cout << "[ERROR] unknown overlap algorithm: " << overlapArg << std::endl;
				exit(1);
			}
		} else {
			// default
			overlap = new RegionGrowingOverlapper();
		}
		ensemble->setOverlapper(*overlap);


		// 4. Final Clusterer
		Clusterer* final = NULL;
		if (finalClustererArg == "LabelPropagation") {
			final = new LabelPropagation();
		} else if (finalClustererArg == "Agglomerative") {
			final = new ParallelAgglomerativeClusterer();
		} else {
			std::cout << "[ERROR] unknown final clusterer: " << finalClustererArg << std::endl;
			exit(1);
		}

		ensemble->setFinalClusterer(*final);

		algo = ensemble;
	}
	else {
		// if no algorithm specified, return empty clustering
		std::cout << "[INFO] no algorithm specified - returning empty clustering" << std::endl;
		return Clustering(0);
	}

	// recurse on ensemble
	if (options[NOREC]) {
		NO_RECURSION = true;
	}


	// START CLUSTERER
	// start solo base algorithm
	std::cout << "[BEGIN] clusterer: " << algo->toString() << std::endl;
	running.start();
	Clustering resultClustering = algo->run(G);
	running.stop();
	//
	std::cout << "[DONE] " << algo->toString() << " ran: \t" << running.elapsedTag() << std::endl;


	// result = std::make_pair(resultClustering, G);

	// print speed info
	double eps = (G.numberOfEdges() / ((double) running.elapsed().count() / 1000.0));	// edges per second
	std::cout << "\t # edges per second:\t" << eps << std::endl;


 	if (options[SUMMARY]) {
 		char sep = ';';
 		std::ofstream summary(options[SUMMARY].arg, std::ios::app); // open summary file to append to
 		summary << algo->toString() << sep;				// APPEND algorithm description
 		summary << REVISION << sep;						// APPEND revision number
 		summary << G.getName() << sep;					// APPEND graph description
 		summary << running.elapsed().count() << sep;	// APPEND running time
 		summary << eps << sep;							// APPEND edges per second
 	}


//	return result;	// return empty clustering
 	return resultClustering;
}

/**
 * Examine the returned clustering and print stats.
 */
bool inspect(Graph& G, Clustering& clustering, OptionParser::Option* options) {

	std::cout << "[INFO] Graph: " << G.toString() << std::endl;

	if (clustering.numberOfEntries() == 0) {
		std::cout << "[EXIT] no clusterer specified" << std::endl;
		return true; // no inspection
	}
	std::cout << "[INFO] inspecting result clustering " << std::endl;

	count k = clustering.numberOfClusters();

	Aux::Timer running;
	running.start();
	Modularity modularity;
	double mod = modularity.getQuality(clustering, G);
	running.stop();
	std::cout << "[DONE] calculating modularity " << running.elapsedTag() << std::endl;

	if ((mod > 1.0) || (mod < -0.5)) {
		std::cout << "[ERROR] modularity calculation went wrong: " << mod << " is not in range [-0.5, 1.0]";
	}


 	std::cout << "\t # clusters:\t" << k << std::endl;
 	std::cout << "\t modularity:\t" << mod << std::endl;


 	if (options[SUMMARY]) {
 		char sep = ';';
 		std::ofstream summary(options[SUMMARY].arg, std::ios::app); // open summary file to append to
 		summary << k << sep;
 		summary << mod << std::endl;  		// since this is the last column, append endl instead of separator
 	}

 	if (options[SAVE_CLUSTERING]) {
 		std::string clusteringFilePath = options[SAVE_CLUSTERING].arg;
 		ClusteringWriter writer;
 		std::cout << "[BEGIN]�saving clustering to file" << std::endl;
 		writer.write(clustering, clusteringFilePath);
 		std::cout << "[DONE] saved clustering to file: " << clusteringFilePath << std::endl;
 	}

 	return true;
}




int main(int argc, char **argv) {
	std::cout << "*** EnsembleClustering: combining parallel clustering algorithms with an ensemble learning strategy *** " << std::endl;

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

#ifndef NOLOG4CXX
	if (options[LOGLEVEL]) {
		configureLogging(options[LOGLEVEL].arg);
	} else {
		configureLogging();	// with default level
	}
#endif


	// CONFIGURE PARALLELISM

#ifdef _OPENMP
	omp_set_nested(1); // enable nested parallelism
#endif

	if (options[THREADS]) {
		// set number of threads
		int nThreads = std::atoi(options[THREADS].arg);
#ifdef _OPENMP
		omp_set_num_threads(nThreads);
#else
		WARN("Thread option ignored since OpenMP is deactivated.");
#endif
	}

	// CONFIGURE OUTPUT
	if (options[PROGRESS]) {
		PRINT_PROGRESS = true;
	} else {
		PRINT_PROGRESS = false;
	}


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

	// CONFIGURE GLOBAL FLAGS

	// CONFIGURE RANDOM ORDER

	if (options[RANDORDER]) {
		configureRandOrder(options[RANDORDER].arg);
	}


	// CALCULATE BASE CLUSTERING DISSIMILARITY?
	if (options[DISSIMILARITY]) {
		CALC_DISSIMILARITY = true;
	}


	// CONFIGURE NORMALIZED VOTES

	if (options[NORM_VOTES]) {
		NORMALIZE_VOTES = true;
	}

	if (options[SCALESTRENGTH]) {
		SCALE_STRENGTH = options[SCALE_STRENGTH].arg;
	}


	// RUN UNIT TESTS

	if (options[TESTS]) {

	   ::testing::InitGoogleTest(&argc, argv);
	   INFO("=== starting unit tests ===");
	   return RUN_ALL_TESTS();
	}


	// SET NUMBER OF CLUSTERER RUNS
	int runs = 1;
	if (options[RUNS]) {
		runs = std::atoi(options[RUNS].arg);
	}



	// RUN PROGRAM
	Graph G = getGraph(options);

	// allow for multiple runs
	for (int run = 0; run < runs; run++) {
		Clustering clustering = startClusterer(G, options);
		inspect(G, clustering, options);
	}

	std::cout << "[EXIT] terminated normally" << std::endl;
	return 0;
}
