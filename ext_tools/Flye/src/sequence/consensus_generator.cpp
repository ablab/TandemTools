//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <cassert>
#include <algorithm>
#include <fstream>

#include "consensus_generator.h"

#include "../common/config.h"
#include "../common/logger.h"
#include "../common/parallel.h"
#include "../common/matrix.h"


namespace
{
	//banded glocal alignment
	void pairwiseAlignment(const std::string& seqOne, const std::string& seqTwo,
						   std::string& outOne, std::string& outTwo, 
						   int bandWidth)
	{
		static const int32_t MATCH = 5;
		static const int32_t SUBST = -5;
		static const int32_t INDEL = -3;

		static const int32_t matchScore[] = {SUBST, MATCH};
		static const int32_t indelScore[] = {INDEL, 0};

		const size_t numRows = seqOne.length() + 1;
		const size_t numCols = 2 * bandWidth + 1;

		Matrix<char> backtrackMatrix(numRows, numCols);
		std::vector<int32_t> scoreRowOne(numCols, 0);
		std::vector<int32_t> scoreRowTwo(numCols, 0);


		for (size_t i = 0; i < numRows; ++i) 
		{
			size_t j = std::max(0, bandWidth - (int)i);
			backtrackMatrix.at(i, j) = 1;
		}
		for (size_t i = 0; i < numCols; ++i) 
		{
			backtrackMatrix.at(0, i) = 0;
		}

		//filling DP matrices
		for (size_t i = 1; i < numRows; ++i)
		{
			int leftOverhang = bandWidth - (int)i + 1;
			int rightOverhand = (int)i + bandWidth - (int)seqTwo.length();
			size_t colLeft = std::max(0, leftOverhang);
			size_t colRight = std::min((int)numCols, (int)numCols - rightOverhand);

			for (int j = colLeft; j < (int)colRight; ++j)
			{
				size_t twoCoord = j + i - bandWidth;
				int32_t cross = scoreRowOne[j] + 
								matchScore[seqOne[i - 1] == seqTwo[twoCoord - 1]];
				char maxStep = 2;
				int32_t maxScore = cross;

				if (j < (int)numCols - 1) //up
				{
					int32_t up = scoreRowOne[j + 1] +
								 indelScore[twoCoord == seqTwo.length()];
					if (up > maxScore)
					{
						maxStep = 1;
						maxScore = up;
					}
				}

				if (j > 0) //left
				{
					int32_t left = scoreRowTwo[j - 1] + 
								   indelScore[i == seqOne.length()];
					if (left > maxScore)
					{
						maxStep = 0;
						maxScore = left;
					}
				}

				scoreRowTwo[j] = maxScore;
				backtrackMatrix.at(i, j) = maxStep;
			}
			scoreRowOne.swap(scoreRowTwo);
		}

		//backtrack
		outOne.reserve(seqOne.length() * 3 / 2);
		outTwo.reserve(seqTwo.length() * 3 / 2);

		int i = numRows - 1;
		int j = bandWidth - (int)numRows + (int)seqTwo.length() + 1;
		while (i != 0 || j != bandWidth) 
		{
			size_t twoCoord = j + i - bandWidth;
			if(backtrackMatrix.at(i, j) == 1) //up
			{
				outOne += seqOne[i - 1];
				outTwo += '-';
				i -= 1;
				j += 1;
			}
			else if (backtrackMatrix.at(i, j) == 0) //left
			{
				outOne += '-';
				outTwo += seqTwo[twoCoord - 1];
				j -= 1;
			}
			else	//cross
			{
				outOne += seqOne[i - 1];
				outTwo += seqTwo[twoCoord - 1];
				i -= 1;
			}
		}
		std::reverse(outOne.begin(), outOne.end());
		std::reverse(outTwo.begin(), outTwo.end());
	}
}


std::vector<FastaRecord> 
	ConsensusGenerator::generateConsensuses(const std::vector<ContigPath>& contigs, 
											bool verbose)
{
	if (verbose) Logger::get().info() << "Generating sequence";
	std::vector<std::vector<AlignmentInfo>> allAlignments;
	std::vector<FastaRecord> consensuses;

	auto alnMap = this->generateAlignments(contigs, verbose);
	//then, generate contig sequences
	for (size_t i = 0; i < contigs.size(); ++i)
	{
		if (contigs[i].sequences.empty()) continue;
		if (contigs[i].sequences.size() == 1)
		{
			FastaRecord rec(contigs[i].sequences.front(), contigs[i].name, 
							FastaRecord::ID_NONE);
			consensuses.push_back(rec);
		}
		else
		{
			consensuses.push_back(this->generateLinear(contigs[i], alnMap));
		}
	}
	return consensuses;
}

FastaRecord ConsensusGenerator::generateLinear(const ContigPath& path, 
											   const AlignmentsMap& alnMap)
{
	std::vector<FastaRecord> contigParts;

	auto prevSwitch = std::make_pair(0, 0);
	std::string contigSequence;
	for (size_t i = 0; i < path.sequences.size(); ++i)
	{
		auto& sequence = path.sequences[i];
		int32_t leftCut = prevSwitch.second;
		int32_t rightCut = sequence.length();
		if (i != path.sequences.size() - 1)
		{
			auto curSwitch = 
				this->getSwitchPositions(alnMap.at(&path.overlaps[i]),
										 prevSwitch.second);
			rightCut = curSwitch.first;
			prevSwitch = curSwitch;
		}

		if (rightCut - leftCut > 0)	//shoudn't happen, but just in case
		{
			contigSequence += sequence.substr(leftCut, rightCut - leftCut).str();
		}
	}
	int32_t cutLen = contigSequence.length() - (path.trimLeft + path.trimRight);
	if (cutLen > 0)
	{
		contigSequence = contigSequence.substr(path.trimLeft, cutLen);
	}
	return FastaRecord(DnaSequence(contigSequence), path.name, 
					   FastaRecord::ID_NONE);
}


ConsensusGenerator::AlignmentsMap 
	ConsensusGenerator::generateAlignments(const std::vector<ContigPath>& contigs,
										   bool verbose)
{
	typedef std::pair<const ContigPath*, size_t> AlnTask;

	AlignmentsMap alnMap;
	std::mutex mapMutex;
	std::function<void(const AlnTask&)> alnFunc =
	[this, &alnMap, &mapMutex](const AlnTask& task)
	{
		const ContigPath* path = task.first;
		size_t i = task.second;

		int32_t leftStart = path->overlaps[i].curBegin;
		int32_t leftLen = path->overlaps[i].curRange();
		std::string leftSeq = path->sequences[i].substr(leftStart, leftLen).str();

		int32_t rightStart = path->overlaps[i].extBegin;
		int32_t rightLen = path->overlaps[i].extRange();
		std::string rightSeq = path->sequences[i + 1].substr(rightStart, 
															 rightLen).str();
		
		const int bandWidth = abs((int)leftSeq.length() - 
								  (int)rightSeq.length()) + 
								  		Config::get("maximum_jump");
		/*if (abs((int)leftSeq.length() - (int)rightSeq.length()) >
			std::min((int)leftSeq.length(), (int)rightSeq.length()))
		{
			Logger::get().warning() << "Aligning sequence that are too "
				<< " different - something is terribly wrong!";
		}*/

		std::string alignedLeft;
		std::string alignedRight;
		pairwiseAlignment(leftSeq, rightSeq, alignedLeft, 
						  alignedRight, bandWidth);

		{
			std::lock_guard<std::mutex> lock(mapMutex);
			alnMap[&path->overlaps[i]] = {alignedLeft, alignedRight, 
						 				  leftStart, rightStart};
		}
	};

	std::vector<AlnTask> tasks;
	for (auto& path : contigs)
	{
		for (size_t i = 0; i < path.sequences.size() - 1; ++i)
		{
			tasks.emplace_back(&path, i);
		}
	}
	processInParallel(tasks, alnFunc, Parameters::get().numThreads, verbose);

	return alnMap;
}


std::pair<int32_t, int32_t> 
ConsensusGenerator::getSwitchPositions(const AlignmentInfo& aln,
									int32_t prevSwitch)
{
	int leftPos = aln.startOne;
	int rightPos = aln.startTwo;
	int matchRun = 0;
	for (size_t i = 0; i < aln.alnOne.length(); ++i)
	{
		if (aln.alnOne[i] != '-') ++leftPos;
		if (aln.alnTwo[i] != '-') ++rightPos;

		if (aln.alnOne[i] == aln.alnTwo[i] &&
			leftPos > prevSwitch + Config::get("maximum_jump"))
		{
			++matchRun;
		}
		else
		{
			matchRun = 0;
		}
		if (matchRun == (int)Parameters::get().kmerSize)
		{
			return {leftPos, rightPos};
		}
	}

	//Logger::get().debug() << "No jump found!";
	prevSwitch = std::max(prevSwitch + 1, aln.startOne);
	return {prevSwitch, aln.startTwo};
}
