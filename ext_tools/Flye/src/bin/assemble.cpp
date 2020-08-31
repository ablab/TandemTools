//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <iostream>
#include <signal.h>
#include <stdlib.h>
#include <unistd.h>
#include <cmath>
#include <execinfo.h>

#include "../sequence/vertex_index.h"
#include "../sequence/sequence_container.h"
#include "../sequence/overlap.h"
#include "../sequence/consensus_generator.h"
#include "../common/config.h"
#include "../assemble/extender.h"
#include "../assemble/parameters_estimator.h"
#include "../common/logger.h"
#include "../common/utils.h"
#include "../common/memory_info.h"

#include <getopt.h>

bool parseArgs(int argc, char** argv, std::string& readsFasta,
               std::string& inAssembly, std::string& outFile, std::string& kmersList,
               std::string& outAssembly, std::string& logFile, size_t& genomeSize,
			   int& kmerSize, float& maxDiff, bool& debug, size_t& numThreads, int& minOverlap,
			   std::string& configPath, int& minReadLength, size_t& minKmers, bool& unevenCov, bool& polishSam)
{
	auto printUsage = [argv]()
	{
		std::cerr << "Usage: flye-assemble "
				  << " --reads path --out-asm path --genome-size size --config path\n"
				  << "\t\t[--min-read length] [--log path] [--treads num]\n"
				  << "\t\t[--kmer size] [--meta] [--min-ovlp size] [--debug] [-h]\n\n"
				  << "Required arguments:\n"
				  << "  --reads path\tcomma-separated list of read files\n"
				  << "  --out-asm path\tpath to output file\n"
				  << "  --genome-size size\tgenome size in bytes\n"
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
		{"reads", required_argument, 0, 0},
		{"asm", required_argument, 0, 0},
		{"kmers", required_argument, 0, 0},
		{"min-kmers", required_argument, 0, 0},
		{"max-diff", required_argument, 0, 0},
        {"out-file", required_argument, 0, 0},
		{"out-asm", required_argument, 0, 0},
		{"genome-size", required_argument, 0, 0},
		{"config", required_argument, 0, 0},
		{"min-read", required_argument, 0, 0},
		{"log", required_argument, 0, 0},
		{"threads", required_argument, 0, 0},
		{"kmer", required_argument, 0, 0},
		{"min-ovlp", required_argument, 0, 0},
		{"meta", no_argument, 0, 0},
		{"debug", no_argument, 0, 0},
		{"polish", no_argument, 0, 0},
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
			else if (!strcmp(longOptions[optionIndex].name, "max-diff"))
				maxDiff = atof(optarg);
			else if (!strcmp(longOptions[optionIndex].name, "min-read"))
				minReadLength = atoi(optarg);
			else if (!strcmp(longOptions[optionIndex].name, "min-kmers"))
				minKmers = atoi(optarg);
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
			else if (!strcmp(longOptions[optionIndex].name, "polish"))
				polishSam = true;
			else if (!strcmp(longOptions[optionIndex].name, "reads"))
				readsFasta = optarg;
			else if (!strcmp(longOptions[optionIndex].name, "asm"))
				inAssembly = optarg;
			else if (!strcmp(longOptions[optionIndex].name, "out-file"))
				outFile = optarg;
			else if (!strcmp(longOptions[optionIndex].name, "kmers"))
				kmersList = optarg;
			else if (!strcmp(longOptions[optionIndex].name, "out-asm"))
				outAssembly = optarg;
			else if (!strcmp(longOptions[optionIndex].name, "genome-size"))
				genomeSize = atoll(optarg);
			else if (!strcmp(longOptions[optionIndex].name, "config"))
				configPath = optarg;
			break;

		case 'h':
			printUsage();
			exit(0);
		}
	}
	if (readsFasta.empty() || outAssembly.empty() || 
		genomeSize == 0 || configPath.empty())
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

	int kmerSize = 15;
	int minReadLength = 0;
	size_t minKmers = 0;
	size_t genomeSize = 0;
	int minOverlap = 5000;
	float maxDiff = 0.15;
	bool debugging = false;
	bool unevenCov = false;
	bool polishSam = false;
	size_t numThreads = 1;
	std::string readsFasta;
	std::string inAssembly;
	std::string outFile;
	std::string kmersList;
	std::string outAssembly;
	std::string logFile;
	std::string configPath;

	if (!parseArgs(argc, argv, readsFasta, inAssembly, outFile, kmersList, outAssembly, logFile, genomeSize,
				   kmerSize, maxDiff, debugging, numThreads, minOverlap, configPath,
				   minReadLength, minKmers, unevenCov, polishSam)) return 1;

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
	Parameters::get().polishSam = polishSam;

	Parameters::get().maxDiff = maxDiff;
	Parameters::get().minKmers = minKmers;
	Logger::get().debug() << "Running with k-mer size: " << 
		Parameters::get().kmerSize; 
	Logger::get().debug() << "Running with minimum overlap " << minOverlap;
	Logger::get().debug() << "Metagenome mode: " << "NY"[unevenCov];

	SequenceContainer readsContainer;
	SequenceContainer assemblyContainer;
	try
	{
	    assemblyContainer.loadFromFile(inAssembly, minReadLength);
	}
	catch (SequenceContainer::ParseException& e)
	{
		Logger::get().error() << e.what();
		return 1;
	}
	assemblyContainer.buildPositionIndex();

	std::vector<std::string> readsList = splitString(readsFasta, ',');
	Logger::get().info() << "Reading sequences";
	try
	{
		for (auto& readsFile : readsList)
		{
			readsContainer.loadFromFile(readsFile, minReadLength);
		}
	}
	catch (SequenceContainer::ParseException& e)
	{
		Logger::get().error() << e.what();
		return 1;
	}
	readsContainer.buildPositionIndex();
	VertexIndex vertexIndex(assemblyContainer,
							(int)Config::get("assemble_kmer_sample"));
	vertexIndex.outputProgress(true);

	int64_t sumLength = 0;
	for (auto& seq : readsContainer.iterSeqs())
	{
		sumLength += seq.sequence.length();
	}
	int coverage = sumLength / 2 / genomeSize;
	Logger::get().debug() << "Expected read coverage: " << coverage;

	Logger::get().info() << "Generating solid k-mer index";
	if (!Parameters::get().unevenCoverage)
	{
		size_t hardThreshold = std::min(5, std::max(2, 
				coverage / (int)Config::get("hard_min_coverage_rate")));
		vertexIndex.countKmers(1, genomeSize);
	}
	else
	{
		vertexIndex.countKmers(/*hard threshold*/ 1, genomeSize);
	}
	//ParametersEstimator estimator(readsContainer, vertexIndex, genomeSize);
	//estimator.estimateMinKmerCount();
	//int minKmerCov = estimator.minKmerCount();
	vertexIndex.setRepeatCutoff(1);
	if (!Parameters::get().unevenCoverage)
	{
		vertexIndex.buildIndex(1, kmersList);
	}
	else
	{
		static const float SELECT_RATE = 0.25;
		static const int TANDEM_FREQ = 10;
		vertexIndex.buildIndexUnevenCoverage(/*min coverage*/ 2, SELECT_RATE, 
											 TANDEM_FREQ);
	}

	Logger::get().debug() << "Peak RAM usage: " 
		<< getPeakRSS() / 1024 / 1024 / 1024 << " Gb";

	std::vector<int> kmerPositions;

	int64_t assemblyLength = 0;
	for (auto& seq : assemblyContainer.iterSeqs())
	{
		assemblyLength += seq.sequence.length();
	}
	kmerPositions.assign(assemblyLength, 0);

	for (auto& seq : assemblyContainer.iterSeqs())
	{
	    size_t numKmers = 0;
	    for (auto kmerPos : IterKmers(seq.sequence))
		{
		    if (vertexIndex.kmerFreq(kmerPos.kmer)) {
                //Logger::get().debug() << "pos " << kmerPos.position << " freq " << vertexIndex.kmerFreq(kmerPos.kmer);
                //for (const auto& extReadPos : vertexIndex.iterKmerPos(itKmer.second))
                numKmers++;
                kmerPositions[kmerPos.position] = numKmers;
		    }
		}
	}

	//int maxOverlapsNum = !Parameters::get().unevenCoverage ? 5 * coverage : 0;
	OverlapDetector ovlp(assemblyContainer, vertexIndex,
						 (int)Config::get("maximum_jump"), 
						 Parameters::get().minimumOverlap,
						 (int)Config::get("maximum_overhang"),
						 /*no max overlaps*/ 0, 
						 /*store alignment*/ true,
						 /*only max*/ true,
						 /*no div threshold*/ 1.0f,
						 /* bad end adjustment*/ 0.0f,
						 /* nucl alignent*/ true);
	OverlapContainer readOverlaps(ovlp, readsContainer);
	readOverlaps.estimateOverlaperParameters(kmerPositions, outFile);
	readOverlaps.setRelativeDivergenceThreshold(
		(float)Config::get("assemble_ovlp_relative_divergence"));

	Extender extender(readsContainer, readOverlaps);
	extender.assembleDisjointigs();
	vertexIndex.clear();

	ConsensusGenerator consGen;
	auto disjointigsFasta = 
		consGen.generateConsensuses(extender.getDisjointigPaths());
	SequenceContainer::writeFasta(disjointigsFasta, outAssembly);

	Logger::get().debug() << "Peak RAM usage: " 
		<< getPeakRSS() / 1024 / 1024 / 1024 << " Gb";

	return 0;
}
