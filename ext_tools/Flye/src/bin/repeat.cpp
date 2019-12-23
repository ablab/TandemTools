//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <iostream>
#include <signal.h>
#include <stdlib.h>
#include <unistd.h>
#include <cmath>

#include "../sequence/sequence_container.h"
#include "../common/config.h"
#include "../common/logger.h"
#include "../common/utils.h"
#include "../common/memory_info.h"

#include "../repeat_graph/repeat_graph.h"
#include "../repeat_graph/multiplicity_inferer.h"
#include "../repeat_graph/graph_processing.h"
#include "../repeat_graph/repeat_resolver.h"
#include "../repeat_graph/output_generator.h"

#include <getopt.h>

bool parseArgs(int argc, char** argv, std::string& readsFasta, 
			   std::string& outFolder, std::string& logFile, 
			   std::string& inAssembly, int& kmerSize,
			   int& minOverlap, bool& debug, size_t& numThreads, 
			   std::string& configPath, bool& unevenCov)
{
	auto printUsage = [argv]()
	{
		std::cerr << "Usage: flye-repeat "
				  << " --disjointigs path --reads path --out-dir path --config path\n"
				  << "\t\t[--log path] [--treads num] [--kmer size] [--meta]\n"
				  << "\t\t[--min-ovlp size] [--debug] [-h]\n\n"
				  << "Required arguments:\n"
				  << "  --disjointigs path\tpath to disjointigs file\n"
				  << "  --reads path\tcomma-separated list of read files\n"
				  << "  --out-dir path\tpath to output directory\n"
				  << "  --config path\tpath to the config file\n\n"
				  << "Optional arguments:\n"
				  << "  --kmer size\tk-mer size [default = 15] \n"
				  << "  --min-ovlp size\tminimum overlap between reads "
				  << "[default = 5000] \n"
				  << "  --debug \t\tenable debug output "
				  << "[default = false] \n"
				  << "  --meta \t\tenable uneven coverage (metagenome) mode "
				  << "[default = false] \n"
				  << "  --log log_file\toutput log to file "
				  << "[default = not set] \n"
				  << "  --threads num_threads\tnumber of parallel threads "
				  << "[default = 1] \n";
	};
	
	int optionIndex = 0;
	static option longOptions[] =
	{
		{"disjointigs", required_argument, 0, 0},
		{"reads", required_argument, 0, 0},
		{"out-dir", required_argument, 0, 0},
		{"config", required_argument, 0, 0},
		{"log", required_argument, 0, 0},
		{"threads", required_argument, 0, 0},
		{"kmer", required_argument, 0, 0},
		{"min-ovlp", required_argument, 0, 0},
		{"meta", no_argument, 0, 0},
		{"debug", no_argument, 0, 0},
		{0, 0, 0, 0}
	};

	int opt = 0;
	while ((opt = getopt_long(argc, argv, "h", longOptions, &optionIndex)) != -1)
	{
		switch(opt)
		{
		case 0:
			if (!strcmp(longOptions[optionIndex].name, "kmer"))
				kmerSize = atoi(optarg);
			else if (!strcmp(longOptions[optionIndex].name, "threads"))
				numThreads = atoi(optarg);
			else if (!strcmp(longOptions[optionIndex].name, "min-ovlp"))
				minOverlap = atoi(optarg);
			else if (!strcmp(longOptions[optionIndex].name, "log"))
				logFile = optarg;
			else if (!strcmp(longOptions[optionIndex].name, "debug"))
				debug = true;
			else if (!strcmp(longOptions[optionIndex].name, "meta"))
				unevenCov = true;
			else if (!strcmp(longOptions[optionIndex].name, "reads"))
				readsFasta = optarg;
			else if (!strcmp(longOptions[optionIndex].name, "out-dir"))
				outFolder = optarg;
			else if (!strcmp(longOptions[optionIndex].name, "disjointigs"))
				inAssembly = optarg;
			else if (!strcmp(longOptions[optionIndex].name, "config"))
				configPath = optarg;
			break;

		case 'h':
			printUsage();
			exit(0);
		}
	}
	if (readsFasta.empty() || outFolder.empty() || 
		inAssembly.empty() || configPath.empty())
	{
		printUsage();
		return false;
	}

	return true;
}

int main(int argc, char** argv)
{
	#ifdef NDEBUG
	signal(SIGSEGV, segfaultHandler);
	std::set_terminate(exceptionHandler);
	#endif

	bool debugging = false;
	size_t numThreads = 1;
	int kmerSize = 15;
	int minOverlap = 5000;
	bool unevenCov = false;
	std::string readsFasta;
	std::string inAssembly;
	std::string outFolder;
	std::string logFile;
	std::string configPath;
	if (!parseArgs(argc, argv, readsFasta, outFolder, logFile, inAssembly,
				   kmerSize, minOverlap, debugging, 
				   numThreads, configPath, unevenCov))  return 1;
	
	Logger::get().setDebugging(debugging);
	if (!logFile.empty()) Logger::get().setOutputFile(logFile);
	Logger::get().debug() << "Build date: " << __DATE__ << " " << __TIME__;
	std::ios::sync_with_stdio(false);
	
	Logger::get().debug() << "Total RAM: " 
		<< getMemorySize() / 1024 / 1024 / 1024 << " Gb";
	Logger::get().debug() << "Available RAM: " 
		<< getFreeMemorySize() / 1024 / 1024 / 1024 << " Gb";
	Logger::get().debug() << "Total CPUs: " << std::thread::hardware_concurrency();

	Config::load(configPath);
	Parameters::get().numThreads = numThreads;
	Parameters::get().kmerSize = kmerSize;
	Parameters::get().minimumOverlap = minOverlap;
	Parameters::get().unevenCoverage = unevenCov;
	Logger::get().debug() << "Running with k-mer size: " << 
		Parameters::get().kmerSize; 
	Logger::get().debug() << "Selected minimum overlap " << minOverlap;
	Logger::get().debug() << "Metagenome mode: " << "NY"[unevenCov];

	Logger::get().info() << "Reading sequences";
	SequenceContainer seqAssembly; 
	SequenceContainer seqReads;
	std::vector<std::string> readsList = splitString(readsFasta, ',');
	try
	{
		seqAssembly.loadFromFile(inAssembly);
		for (auto& readsFile : readsList)
		{
			seqReads.loadFromFile(readsFile);
		}
	}
	catch (SequenceContainer::ParseException& e)
	{
		Logger::get().error() << e.what();
		return 1;
	}
	seqReads.buildPositionIndex();
	seqAssembly.buildPositionIndex();

	SequenceContainer edgeSequences;
	RepeatGraph rg(seqAssembly, &edgeSequences);
	GraphProcessor proc(rg, seqAssembly);
	ReadAligner aligner(rg, seqReads);
	OutputGenerator outGen(rg, aligner, seqReads);

	Logger::get().info() << "Building repeat graph";
	rg.build();
	//outGen.outputDot(proc.getEdgesPaths(), outFolder + "/graph_raw.gv");

	Logger::get().info() << "Aligning reads to the graph";
	aligner.alignReads();

	MultiplicityInferer multInf(rg, aligner, seqAssembly, seqReads);
	multInf.estimateCoverage();

	multInf.removeUnsupportedConnections();
	multInf.maskUnsupportedEdges();
	multInf.collapseHeterozygousLoops(/*remove alternatives*/ false);
	multInf.collapseHeterozygousBulges(/*remove alternatives*/ false);
	
	//aligner.storeAlignments(outFolder + "/read_alignment_before_rr");

	Logger::get().info() << "Resolving repeats";
	RepeatResolver resolver(rg, seqAssembly, seqReads, aligner, multInf);
	resolver.findRepeats();
	outGen.outputDot(proc.getEdgesPaths(), outFolder + "/graph_before_rr.gv");
	//outGen.outputGfa(proc.getEdgesPaths(), outFolder + "/graph_before_rr.gfa");
	outGen.outputFasta(proc.getEdgesPaths(), outFolder + "/graph_before_rr.fasta");
	resolver.resolveRepeats();

	//clean graph again after repeat resolution
	multInf.removeUnsupportedEdges();
	multInf.collapseHeterozygousLoops(/*remove alternatives*/ true);
	multInf.collapseHeterozygousBulges(/*remove alternatives*/ true);
	resolver.findRepeats();
	
	resolver.finalizeGraph();

	outGen.outputDot(proc.getEdgesPaths(), outFolder + "/graph_after_rr.gv");
	rg.storeGraph(outFolder + "/repeat_graph_dump");
	aligner.storeAlignments(outFolder + "/read_alignment_dump");
	SequenceContainer::writeFasta(edgeSequences.iterSeqs(), 
								  outFolder + "/repeat_graph_edges.fasta",
								  /*only pos strand*/ true);

	Logger::get().debug() << "Peak RAM usage: " 
		<< getPeakRSS() / 1024 / 1024 / 1024 << " Gb";

	return 0;
}
