//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <ctime>
#include <iostream>
#include <getopt.h>
#include <cstring>

#include "../polishing/bubble_processor.h"


bool parseArgs(int argc, char** argv, std::string& bubblesFile, 
			   std::string& scoringMatrix, std::string& hopoMatrix,
			   std::string& outConsensus, std::string& outVerbose,
			   int& numThreads, bool& quiet)
{
	auto printUsage = [argv]()
	{
		std::cerr << "Usage: flye-polish "
				  << " --bubbles path --subs-mat path --hopo-mat size --out path\n"
				  << "\t\t[--treads num] [--quiet] [--debug] [-h]\n\n"
				  << "Required arguments:\n"
				  << "  --bubbles path\tpath to bubbles file\n"
				  << "  --subs-mat path\tpath to substitution matrix\n"
				  << "  --hopo-mat size\tpath to homopolymer matrix\n"
				  << "  --out path\tpath to output file\n\n"
				  << "Optional arguments:\n"
				  << "  --quiet \t\tno terminal output "
				  << "[default = false] \n"
				  << "  --debug \t\textra debug output "
				  << "[default = false] \n"
				  << "  --threads num_threads\tnumber of parallel threads "
				  << "[default = 1] \n";
	};
	
	int optionIndex = 0;
	static option longOptions[] =
	{
		{"bubbles", required_argument, 0, 0},
		{"subs-mat", required_argument, 0, 0},
		{"hopo-mat", required_argument, 0, 0},
		{"out", required_argument, 0, 0},
		{"threads", required_argument, 0, 0},
		{"debug", no_argument, 0, 0},
		{"quiet", no_argument, 0, 0},
		{0, 0, 0, 0}
	};

	int opt = 0;
	while ((opt = getopt_long(argc, argv, "h", longOptions, &optionIndex)) != -1)
	{
		switch(opt)
		{
		case 0:
			if (!strcmp(longOptions[optionIndex].name, "threads"))
				numThreads = atoi(optarg);
			else if (!strcmp(longOptions[optionIndex].name, "debug"))
				outVerbose = true;
			else if (!strcmp(longOptions[optionIndex].name, "quiet"))
				quiet = true;
			else if (!strcmp(longOptions[optionIndex].name, "bubbles"))
				bubblesFile = optarg;
			else if (!strcmp(longOptions[optionIndex].name, "subs-mat"))
				scoringMatrix = optarg;
			else if (!strcmp(longOptions[optionIndex].name, "hopo-mat"))
				hopoMatrix = optarg;
			else if (!strcmp(longOptions[optionIndex].name, "out"))
				outConsensus = optarg;
			break;

		case 'h':
			printUsage();
			exit(0);
		}
	}
	if (bubblesFile.empty() || scoringMatrix.empty() || 
		hopoMatrix.empty() || outConsensus.empty())
	{
		printUsage();
		return false;
	}

	return true;
}

int main(int argc, char* argv[]) 
{
	std::string bubblesFile;
	std::string scoringMatrix;
	std::string hopoMatrix;
	std::string outConsensus;
	std::string outVerbose;
	int  numThreads = 0;
	bool quiet = false;
	if (!parseArgs(argc, argv, bubblesFile, scoringMatrix, 
				   hopoMatrix, outConsensus, outVerbose, numThreads,
				   quiet))
		return 1;

	BubbleProcessor bp(scoringMatrix, hopoMatrix, !quiet);
	if (!outVerbose.empty())
		bp.enableVerboseOutput(outVerbose);
	bp.polishAll(bubblesFile, outConsensus, numThreads); 

	return 0;
}
