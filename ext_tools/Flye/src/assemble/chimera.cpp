//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <vector>
#include <iostream>
#include <unordered_set>
#include <unordered_map>
#include <iomanip>
#include <cmath>

#include "../common/config.h"
#include "../common/logger.h"
#include "chimera.h"



bool ChimeraDetector::isChimeric(FastaRecord::Id readId)
{
	if (!_chimeras.contains(readId))
	{
		const auto& ovlps = _ovlpContainer.lazySeqOverlaps(readId);
		bool result = this->testReadByCoverage(readId, ovlps) ||
					  _ovlpContainer.hasSelfOverlaps(readId);
		_chimeras.insert(readId, result);
		_chimeras.insert(readId.rc(), result);
	}
	return _chimeras.find(readId);
}

bool ChimeraDetector::isChimeric(FastaRecord::Id readId,
								 const std::vector<OverlapRange>& readOvlps)
{
	const int JUMP = Config::get("maximum_jump");
	if (!_chimeras.contains(readId))
	{
		bool result = this->testReadByCoverage(readId, readOvlps);
		for (const auto& ovlp : readOvlps)
		{
			if (ovlp.curId == ovlp.extId.rc()) 
			{
				int32_t projEnd = ovlp.extLen - ovlp.extEnd - 1;
				if (abs(ovlp.curEnd - projEnd) < JUMP)
				{
					result = true;
				}
			}
		}
		_chimeras.insert(readId, result);
		_chimeras.insert(readId.rc(), result);
	}
	return _chimeras.find(readId);
}

void ChimeraDetector::estimateGlobalCoverage()
{
	Logger::get().debug() << "Estimating overlap coverage";

	int numSamples = std::min(1000, (int)_seqContainer.iterSeqs().size());
	int sampleRate = (int)_seqContainer.iterSeqs().size() / numSamples;
	//int minCoverage = _inputCoverage / 
	//				(int)Config::get("max_coverage_drop_rate") + 1;
	//int maxCoverage = _inputCoverage * 
	//				(int)Config::get("max_coverage_drop_rate");
	int flankSize = 0;

	std::unordered_map<int32_t, int32_t> readHist;
	std::vector<int32_t> covList;
	
	//std::ofstream fout("../cov_hist.txt");

	int64_t sum = 0;
	int64_t num = 0;
	for (const auto& seq : _seqContainer.iterSeqs())
	{
		if (rand() % sampleRate) continue;
		const auto& overlaps = _ovlpContainer.lazySeqOverlaps(seq.id);
		auto coverage = this->getReadCoverage(seq.id, overlaps);
		bool nonZero = false;
		for (auto c : coverage) nonZero |= (c != 0);
		if (!nonZero) continue;

		for (size_t i = flankSize; i < coverage.size() - flankSize; ++i)
		{
			{
				++readHist[coverage[i]];
				sum += coverage[i];
				++num;
				covList.push_back(coverage[i]);
			}
		}
	}

	if (readHist.empty())
	{
		Logger::get().warning() << "No overlaps found!";
		_overlapCoverage = 0;
	}
	else
	{
		_overlapCoverage = median(covList);
	}

	Logger::get().info() << "Overlap-based coverage: " << _overlapCoverage;
}

std::vector<int32_t> 
	ChimeraDetector::getReadCoverage(FastaRecord::Id readId,
									 const std::vector<OverlapRange>& readOverlaps)
{
	static const int WINDOW = Config::get("chimera_window");
	//const int FLANK = (int)Config::get("maximum_overhang") / WINDOW;
	const int FLANK = 1;

	std::vector<int> coverage;
	int numWindows = std::ceil((float)_seqContainer.seqLen(readId) / WINDOW) + 1;
	if (numWindows - 2 * FLANK <= 0) return {0};

	coverage.assign(numWindows - 2 * FLANK, 0);
	for (const auto& ovlp : readOverlaps)
	{
		if (ovlp.curId == ovlp.extId.rc() ||
			ovlp.curId == ovlp.extId) continue;

		//skip 2 first/last windows of overlap to be more robust to
		//possible coorinate shifts
		for (int pos = ovlp.curBegin / WINDOW + FLANK; 		
			 pos <= ovlp.curEnd / WINDOW - FLANK; ++pos)
		{
			if (pos - FLANK >= 0 && 					//shouldn't happen, but just in case..
				pos - FLANK < (int)coverage.size())
			{
				//assert(pos - FLANK >= 0 && pos - FLANK < (int)coverage.size());
				++coverage[pos - FLANK];
			}
		}
	}

	return coverage;
}


bool ChimeraDetector::testReadByCoverage(FastaRecord::Id readId,
										 const std::vector<OverlapRange>& readOvlps)
{
	const float MAX_DROP_RATE = Config::get("max_coverage_drop_rate");

	auto coverage = this->getReadCoverage(readId, readOvlps);
	if (coverage.empty()) return false;

	int32_t maxCov = 0;
	int64_t sumCov = 0;
	for (auto cov : coverage)
	{
		maxCov = std::max(maxCov, cov);
		sumCov += cov;
	}
	//int32_t meanCoverage = sumCov / coverage.size();
	int32_t medianCoverage = median(coverage);
	if (sumCov == 0) return false;	//no overlaps found, but it's not chimeric either

	int threshold = 0;	
	if (!Parameters::get().unevenCoverage)
	{
		threshold = std::max(1L, std::lround((float)_overlapCoverage / 
											 MAX_DROP_RATE));
	}
	else
	{
		/*threshold = std::round((float)std::min(_overlapCoverage, maxCov) /
							   MAX_DROP_RATE);*/
		/*threshold = 1;*/
		threshold = std::max(1L, std::lround(medianCoverage / MAX_DROP_RATE));
	}

	int32_t goodStart = 0;
	int32_t goodEnd = coverage.size() - 1;
	//if (!Parameters::get().unevenCoverage)
	//{
	const int MAX_FLANK = (int)Config::get("maximum_overhang") / 
							Config::get("chimera_window");
	goodStart = MAX_FLANK;
	goodEnd = coverage.size() - MAX_FLANK - 1;
	//}
	/*else
	{
		for (goodEnd = 0; goodStart < (int)coverage.size(); ++goodStart)
		{
			if (coverage[goodStart] >= threshold) break;
			if (goodStart > 0 && coverage[goodStart] < coverage[goodStart - 1]) break;
		}
		for (goodEnd = coverage.size() - 1; goodEnd >= 0; --goodEnd)
		{
			if (coverage[goodEnd] >= threshold) break;
			if (goodEnd > 0 && coverage[goodEnd] > coverage[goodEnd - 1]) break;
		}
	}*/
	
	bool lowCoverage = false;
	if (goodEnd <= goodStart) lowCoverage = true;
	for (int32_t i = goodStart; i <= goodEnd; ++i)
	{
		if (coverage[i] < threshold)
		{
			lowCoverage = true;
			break;
		}
	}

	/*static std::ofstream fout("chimera_dump.txt");
	static std::mutex logLock;
	logLock.lock();
	std::string covStr;
	for (int32_t i = 0; i < (int)coverage.size() - 1; ++i)
	{
		if (i == goodStart) covStr += " { ";
		covStr += std::to_string(coverage[i]) + " ";
		if (i == goodEnd) covStr += " } ";
	}
	fout << _seqContainer.seqName(readId) << " " 
		<< _seqContainer.seqLen(readId) << std::endl;
	fout << covStr << std::endl;
	fout << "max: " << maxCov << " mean: " << meanCoverage << " median: " << medianCoverage
		<< " threshold: " << threshold << " chim:" << lowCoverage << std::endl << std::endl;
	logLock.unlock();*/

	return lowCoverage;
}
