//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <set>
#include <iostream>
#include <cassert>
#include <algorithm>
#include <thread>
#include <ctime>
#include <queue>
#include <cmath>
#include <chrono>
#include <ctime>
#include <cstring>
#include <iomanip>
#include <numeric>

#define HAVE_KALLOC
#include "kalloc.h"
#include "ksw2.h"
#undef HAVE_KALLOC

#include "overlap.h"
#include "../common/config.h"
#include "../common/utils.h"
#include "../common/parallel.h"
#include "../common/disjoint_set.h"
#include "../common/bfcontainer.h"


namespace
{
	using namespace std::chrono;
	struct ThreadMemPool
	{
		ThreadMemPool():
			prevCleanup(system_clock::now() + seconds(rand() % 60))
		{
		   memPool = km_init();
		}
		~ThreadMemPool()
		{
		   km_destroy(memPool);
		}
		void cleanIter()
		{
			if ((system_clock::now() - prevCleanup) > seconds(60))
			{
				km_destroy(memPool);
               	memPool = km_init();
				prevCleanup = system_clock::now();
			}
		}

		time_point<system_clock> prevCleanup;
		void* memPool;
    };

    std::string kswAlign(const DnaSequence& trgSeq, std::string trgName, size_t trgBegin, size_t trgLen,
                   const DnaSequence& qrySeq, std::string qryName, size_t qryBegin, size_t qryLen,
				   int matchScore, int misScore, int gapOpen, int gapExtend,
				   bool showAlignment)
	{
		static const int32_t MAX_JUMP = 1000;
		const int KMER_SIZE = Parameters::get().kmerSize;
		const int IS_FLYE_POLISH = Parameters::get().polishSam;

		thread_local ThreadMemPool buf;
		thread_local std::vector<uint8_t> trgByte;
		thread_local std::vector<uint8_t> qryByte;
		buf.cleanIter();
		trgByte.assign(trgLen, 0);
		qryByte.assign(qryLen, 0);

		for (size_t i = 0; i < trgLen; ++i)
		{
			trgByte[i] = trgSeq.atRaw(i + trgBegin);
		}
		for (size_t i = 0; i < qryLen; ++i)
		{
			qryByte[i] = qrySeq.atRaw(i + qryBegin);
		}

		int seqDiff = abs((int)trgByte.size() - (int)qryByte.size());
		int bandWidth = seqDiff + MAX_JUMP;

		//substitution matrix
		int8_t a = matchScore;
		int8_t b = misScore < 0 ? misScore : -misScore; // a > 0 and b < 0
		int8_t subsMat[] = {a, b, b, b, 0, 
						  	b, a, b, b, 0, 
						  	b, b, a, b, 0, 
						  	b, b, b, a, 0, 
						  	0, 0, 0, 0, 0};

		ksw_extz_t ez;
		memset(&ez, 0, sizeof(ksw_extz_t));
		const int NUM_NUCL = 5;
		const int Z_DROP = -1;
		const int FLAG = KSW_EZ_APPROX_MAX | KSW_EZ_APPROX_DROP;
		const int END_BONUS = 0;
		//ksw_extf2_sse(0, qseq.size(), &qseq[0], tseq.size(), &tseq[0], matchScore,
		//		 	  misScore, gapOpen, bandWidth, Z_DROP, &ez);
		ksw_extz2_sse(buf.memPool, qryByte.size(), &qryByte[0], 
					  trgByte.size(), &trgByte[0], NUM_NUCL,
				 	  subsMat, gapOpen, gapExtend, bandWidth, Z_DROP, 
					  END_BONUS, FLAG, &ez);

		//decode CIGAR
		std::string strQ;
		std::string strT;
		std::string alnQry;
		std::string alnTrg;
		if (true)
		{
			for (auto x : qryByte) strQ += "ACGT"[x];
			for (auto x : trgByte) strT += "ACGT"[x];
		}

		size_t posQry = 0;
		size_t posTrg = 0;
        int matches = 0;
		int alnLength = 0;
		int condensedLength = 0;
		std::string cigar_str;
        cigar_str += std::to_string(qryBegin);
        cigar_str += "S";

		for (size_t i = 0; i < (size_t)ez.n_cigar; ++i)
		{
			int size = ez.cigar[i] >> 4;
			char op = "MID"[ez.cigar[i] & 0xf];
			alnLength += size;
            cigar_str += std::to_string(size);
            cigar_str += op;

        	if (op == 'M')
			{
				for (size_t i = 0; i < (size_t)size; ++i)
				{
					if (trgByte[posTrg + i] == qryByte[posQry + i]) ++matches;
				}
				if (showAlignment)
				{
					alnQry += strQ.substr(posQry, size);
					alnTrg += strT.substr(posTrg, size);
				}
				posQry += size;
				posTrg += size;
				condensedLength += size;
			}
            else if (op == 'I')
			{
				if (showAlignment)
				{
					alnQry += strQ.substr(posQry, size);
					if (IS_FLYE_POLISH) alnTrg += std::string(size, '-');
				}
                posQry += size;
				condensedLength += std::min(size, KMER_SIZE);
			}
            else //D
			{
				if (showAlignment)
				{
				    if (IS_FLYE_POLISH) alnQry += std::string(size, '-');
					alnTrg += strT.substr(posTrg, size);
				}
				posTrg += size;
				condensedLength += std::min(size, KMER_SIZE);
			}
		}
        cigar_str += std::to_string(qrySeq.length() - qryBegin - qryLen);
        cigar_str += "S";

        float errRate = 1 - (float(matches) / condensedLength);
        int errors = qryLen - matches;

        std::stringstream ss;
        if (errRate<=0.2) {
            std::string readFlag = qryName.front() == '+' ? "\t0\t" : "\t16\t";
            if (IS_FLYE_POLISH) {
            ss << qryName.substr(1) << readFlag << trgName.substr(1) << "\t" << (trgBegin+1) <<"\t60\t" << cigar_str << "\t*\t0\t0\t" << alnTrg << "\t" << alnQry << "\t*\tNM:i:" << errors << "\n";
            }
            else {
            ss << qryName.substr(1) << readFlag << trgName.substr(1) << "\t" << (trgBegin+1) <<"\t60\t" << cigar_str << "\t*\t0\t0\t" << qrySeq.str() << "\t*\tNM:i:" << errors << "\n";
            }

        }
        //float errRate = 1 - float(matches) / alnLength;
        //float errRate = 1 - float(matches) / std::max(tseq.size(), qseq.size());
		kfree(buf.memPool, ez.cigar);

		return ss.str();
	}
}

//Check if it is a proper overlap
bool OverlapDetector::overlapTest(const OverlapRange& ovlp,
								  bool& outSuggestChimeric) const
{
	if (ovlp.curRange() < _minOverlap || 
		ovlp.extRange() < _minOverlap) 
	{
		return false;
	}

	//filter overlaps that way to divergent in length.
	//theoretically, they should not pass sequence divergence filter,
	//but just in case
	static const float OVLP_DIV = 0.5;
	float lengthDiff = abs(ovlp.curRange() - ovlp.extRange());
	if (lengthDiff > OVLP_DIV * std::min(ovlp.curRange(), ovlp.extRange()))
	{
		return false;
	}

	//check if it's "almost trivial" match with intersecting sequence
	if (ovlp.curId == ovlp.extId)
	{
		int32_t intersect = std::min(ovlp.curEnd, ovlp.extEnd) - 
			   				std::max(ovlp.curBegin, ovlp.extBegin);
		if (intersect > ovlp.curRange() / 2) return false;
	}

	//check "strand skipping" PacBio pattern
	if (ovlp.curId == ovlp.extId.rc())
	{
		int32_t intersect = std::min(ovlp.curEnd, ovlp.extLen - ovlp.extBegin) - 
			   				std::max(ovlp.curBegin, ovlp.extLen - ovlp.extEnd);

		if (intersect > -_maxJump) outSuggestChimeric = true;
		if (intersect > ovlp.curRange() / 2) return false;
	}

	if (_checkOverhang)
	{
		if (std::min(ovlp.curBegin, ovlp.extBegin) > 
			_maxOverhang) 
		{
			return false;
		}
		if (std::min(ovlp.curLen - ovlp.curEnd, ovlp.extLen - ovlp.extEnd) > 
			_maxOverhang)
		{
			return false;
		}
	}

	return true;
}


namespace
{

	template <class T>
	void shrinkAndClear(std::vector<T>& vec, float rate)
	{
		if (vec.empty()) return;
		//if ((float)vec.capacity() / vec.size() < rate) return;
		//Logger::get().debug() << "Shrink: " << vec.capacity() << " " << vec.size();

		size_t newCapacity = std::roundf((float)vec.capacity() / rate);
		vec = std::vector<T>();
		vec.reserve(newCapacity);
	}
}

//This implementation was inspired by Heng Li's minimap2 paper
//might be used in parallel
std::vector<OverlapRange> 
OverlapDetector::getSeqOverlaps(const FastaRecord& fastaRec,
								bool& outSuggestChimeric,
								OvlpDivStats& divStats) const
{
	//static std::ofstream fout("../kmers.txt");
	
	//const int MAX_LOOK_BACK = 50;
	const int kmerSize = Parameters::get().kmerSize;
	//const float minKmerSruvivalRate = std::exp(-_maxDivergence * kmerSize);
	const float minKmerSruvivalRate = 0.001;
	const float LG_GAP = 10;
	const float SM_GAP = 0;
	const float SM_DIFF = 10;
	const float MAX_DIFF = Parameters::get().maxDiff;  // max difference between distances in a read and in an assembly
	const float KMER_BONUS = 10; // bonus for each added non-overlapping k-mer
    int maxOverlaps = 0;

	outSuggestChimeric = false;
	int32_t curLen = fastaRec.sequence.length();
	std::vector<int32_t> curFilteredPos;

	//cache memory-intensive containers as
	//many parallel memory allocations slow us down significantly
	//thread_local std::vector<KmerMatch> vecMatches;
	thread_local std::vector<KmerMatch> matchesList;
	thread_local std::vector<int32_t> scoreTable;
	thread_local std::vector<int32_t> backtrackTable;

	static ChunkPool<KmerMatch> sharedChunkPool;	//shared accoress threads
	BFContainer<KmerMatch> vecMatches(sharedChunkPool);

	//speed benchmarks
	thread_local float timeMemory = 0;
	thread_local float timeKmerIndexFirst = 0;
	thread_local float timeKmerIndexSecond = 0;
	thread_local float timeDp = 0;
	auto timeStart = std::chrono::system_clock::now();

	/*static std::mutex reportMutex;
	thread_local int threadId = rand() % 1000; 
	thread_local int numTicks = 0;
	++numTicks;
	thread_local auto prevReport = std::chrono::system_clock::now();
	if ((std::chrono::system_clock::now() - prevReport) > 
		std::chrono::seconds(60))
	{
		std::lock_guard<std::mutex> lock(reportMutex);
		prevReport = std::chrono::system_clock::now();
		Logger::get().debug() << ">Perf " << threadId
			<< " mem:" << timeMemory
			<< " kmFst:" << timeKmerIndexFirst 
			<< " kmSnd:" << timeKmerIndexSecond 
			<< " dp:" << timeDp << " ticks: " << numTicks; 
		Logger::get().debug() << ">Mem  " << threadId 
			<< " chunks:" << sharedChunkPool.numberChunks();

		timeMemory = 0;
		timeKmerIndexFirst = 0;
		timeKmerIndexSecond = 0;
		timeDp = 0;
		numTicks = 0;
	}*/
	timeStart = std::chrono::system_clock::now();

	//although once in a while shrink allocated memory size
	//thread_local auto prevCleanup = 
	//	std::chrono::system_clock::now() + std::chrono::seconds(rand() % 60);
	thread_local int prevCleanup = 0;
	if (++prevCleanup > 50)
	{
		prevCleanup = 0;
		shrinkAndClear(matchesList, 2);
		shrinkAndClear(scoreTable, 2);
		shrinkAndClear(backtrackTable, 2);
	}
	timeMemory += std::chrono::duration_cast<std::chrono::duration<float>>
						(std::chrono::system_clock::now() - timeStart).count();
	timeStart = std::chrono::system_clock::now();

	cuckoohash_map<Kmer, char> _goodKmers;
	cuckoohash_map<Kmer, char> _badKmers;
	for (const auto& curKmerPos : IterKmers(fastaRec.sequence))
	{
		if (!_vertexIndex.kmerFreq(curKmerPos.kmer)) continue;

        if (_goodKmers.contains(curKmerPos.kmer) && !_badKmers.contains(curKmerPos.kmer))
        {
            _badKmers.insert(curKmerPos.kmer, true);
        }
        else
        {
            _goodKmers.insert(curKmerPos.kmer, true);
        }
	}

	for (const auto& curKmerPos : IterKmers(fastaRec.sequence))
	{
		if (!_goodKmers.contains(curKmerPos.kmer)) continue;

		//FastaRecord::Id prevSeqId = FastaRecord::ID_NONE;
        Kmer k = curKmerPos.kmer;
        k.standardForm();
		for (const auto& extReadPos : _vertexIndex.iterKmerPos(k))
		{
			//no trivial matches
			if ((extReadPos.readId == fastaRec.id &&
				extReadPos.position == curKmerPos.position)) continue;

            int32_t freq = _vertexIndex.kmerFreq(curKmerPos.kmer);
            int bonus = freq > 1 ? std::log(51-freq) : 10.0;
                bonus = std::max(bonus, 1); //std::log(52-freq)
			vecMatches.emplace_back(curKmerPos.position, 
									extReadPos.position,
									bonus,
									extReadPos.readId);
		}
	}
	timeKmerIndexFirst += std::chrono::duration_cast<std::chrono::duration<float>>
							(std::chrono::system_clock::now() - timeStart).count();
	timeStart = std::chrono::system_clock::now();

	std::sort(vecMatches.begin(), vecMatches.end(),
			  [](const KmerMatch& k1, const KmerMatch& k2)
			  {return k1.extId != k2.extId ? k1.extId < k2.extId : 
			  								 k1.curPos < k2.curPos;});

	timeKmerIndexSecond += std::chrono::duration_cast<std::chrono::duration<float>>
								(std::chrono::system_clock::now() - timeStart).count();
	timeStart = std::chrono::system_clock::now();

	const int STAT_WND = 10000;
	std::vector<OverlapRange> divStatWindows(curLen / STAT_WND + 1);

	std::vector<OverlapRange> detectedOverlaps;
	size_t extRangeBegin = 0;
	size_t extRangeEnd = 0;
	while(extRangeEnd < vecMatches.size())
	{
		if (maxOverlaps != 0 &&
			detectedOverlaps.size() >= (size_t)maxOverlaps) break;

		extRangeBegin = extRangeEnd;
		size_t uniqueMatches = 0;
		int32_t prevPos = 0;
		while (extRangeEnd < vecMatches.size() &&
			   vecMatches[extRangeBegin].extId == 
			   vecMatches[extRangeEnd].extId)
		{
			if (vecMatches[extRangeEnd].curPos != prevPos)
			{
				++uniqueMatches;
				prevPos = vecMatches[extRangeEnd].curPos;
			}
			++extRangeEnd;
		}
		if (uniqueMatches < minKmerSruvivalRate * _minOverlap) continue;

		matchesList.assign(vecMatches.begin() + extRangeBegin,
						   vecMatches.begin() + extRangeEnd);
		assert(matchesList.size() > 0 && 
			   matchesList.size() < (size_t)std::numeric_limits<int32_t>::max());

		FastaRecord::Id extId = matchesList.front().extId;
		int32_t extLen = _seqContainer.seqLen(extId);

		//pre-filtering
		int32_t minCur = matchesList.front().curPos;
		int32_t maxCur = matchesList.back().curPos;
		int32_t minExt = std::numeric_limits<int32_t>::max();
		int32_t maxExt = std::numeric_limits<int32_t>::min();
		for (const auto& match : matchesList)
		{
			minExt = std::min(minExt, match.extPos);
			maxExt = std::max(maxExt, match.extPos);
		}
		if (maxCur - minCur < _minOverlap || 
			maxExt - minExt < _minOverlap) continue;
		if (_checkOverhang)
		{
			if (std::min(minCur, minExt) > _maxOverhang) continue;
			if (std::min(curLen - maxCur, 
						 extLen - maxExt) > _maxOverhang) continue;
		}
		//++uniqueCandidates;

		//chain matiching positions with DP
		scoreTable.assign(matchesList.size(), 0);
		backtrackTable.assign(matchesList.size(), -1);

		bool extSorted = extLen > curLen;
		if (extSorted)
		{
			std::sort(matchesList.begin(), matchesList.end(),
					  [](const KmerMatch& k1, const KmerMatch& k2)
					  {return k1.extPos < k2.extPos;});
		}

		for (int32_t i = 1; i < (int32_t)scoreTable.size(); ++i)
		{
			int32_t maxScore = 0;
			int32_t maxId = 0;
			int32_t curNext = matchesList[i].curPos;
			int32_t extNext = matchesList[i].extPos;
			//int32_t noImprovement = 0;

			for (int32_t j = i - 1; j >= 0; --j)
			{
				int32_t curPrev = matchesList[j].curPos;
				int32_t extPrev = matchesList[j].extPos;
				if (0 < curNext - curPrev && curNext - curPrev < _maxJump &&
					0 < extNext - extPrev && extNext - extPrev < _maxJump)
				{
					int32_t matchScore = std::min(abs(curNext - curPrev), abs(extNext - extPrev));
					int32_t distDiff = abs(abs(curNext - curPrev) - abs(extNext - extPrev));

					if (distDiff< 20000) { //  MAX_DIFF*matchScore) {
                        int32_t diffPenalty = distDiff > SM_DIFF ? (distDiff*100)/matchScore : std::min(1, distDiff);
                        int32_t maxBonus = matchesList[i].bonus; //matchScore >= kmerSize ? matchesList[i].bonus : 1;

                        int32_t nextScore = std::max(scoreTable[j]-20, scoreTable[j] + maxBonus - diffPenalty);
                        if (nextScore > maxScore && nextScore - scoreTable[j] > -100)
                        {
                            maxScore = nextScore;
                            maxId = j;
                            //noImprovement = 0;

                            //if (jumpDiv == 0 && curNext - curPrev < kmerSize) break;
                        }
                        /*else
                        {
                            if (++noImprovement > MAX_LOOK_BACK) break;
                        }*/
					}
				}
				if (extSorted && extNext - extPrev > _maxJump) break;
				if (!extSorted && curNext - curPrev > _maxJump) break;
			}

			scoreTable[i] = std::max(maxScore, kmerSize);
			if (maxScore > kmerSize)
			{
				backtrackTable[i] = maxId;
			}
		}

		//backtracking
		std::vector<OverlapRange> extOverlaps;
		std::vector<int32_t> shifts;
		std::vector<KmerMatch> kmerMatches;

		for (int32_t chainStart = backtrackTable.size() - 1;
			 chainStart > 0; --chainStart)
		{
			if (backtrackTable[chainStart] == -1) continue;

			int32_t chainMaxScore = scoreTable[chainStart];
			int32_t lastMatch = chainStart;
			int32_t firstMatch = 0;

			int32_t chainLength = 0;
			shifts.clear();
			kmerMatches.clear();
			//int32_t totalMatch = kmerSize;
			//int32_t totalGap = 0;
			
			int32_t pos = chainStart;
			while (pos != -1)
			{
				//found a new maximum, shorten the chain end
				if (scoreTable[pos] > chainMaxScore)
				{
					chainMaxScore = scoreTable[pos];
					lastMatch = pos;
					chainLength = 0;
					shifts.clear();
					kmerMatches.clear();
				}

				firstMatch = pos;
				shifts.push_back(matchesList[pos].curPos - 
								 matchesList[pos].extPos);
				++chainLength;

				/*int32_t prevPos = backtrackTable[pos];
				if (prevPos != -1)
				{
					int32_t curNext = matchesList[pos].curPos;
					int32_t extNext = matchesList[pos].extPos;
					int32_t curPrev = matchesList[prevPos].curPos;
					int32_t extPrev = matchesList[prevPos].extPos;
					int32_t matchScore = 
							std::min(std::min(curNext - curPrev, extNext - extPrev),
											  kmerSize);
					int32_t jumpDiv = abs((curNext - curPrev) - 
										  (extNext - extPrev));
					int32_t gapCost = (jumpDiv > 50) ? 2 * jumpDiv : jumpDiv;

					totalMatch += matchScore;
					totalGap += gapCost;
				}*/
				if (_keepAlignment)
				{
                    kmerMatches.emplace_back(matchesList[pos].curPos,
                                             matchesList[pos].extPos,
                                             scoreTable[pos],
                                             matchesList[pos].extId);
				}

				assert(pos >= 0 && pos < (int32_t)backtrackTable.size());
				int32_t newPos = backtrackTable[pos];
				//backtrackTable[pos] = -1;
				pos = newPos;
			}

			//Logger::get().debug() << chainStart - firstMatch << " " << lastMatch - firstMatch;

			OverlapRange ovlp(fastaRec.id, matchesList.front().extId,
							  matchesList[firstMatch].curPos, 
							  matchesList[firstMatch].extPos,
							  curLen, extLen);
			ovlp.curEnd = matchesList[lastMatch].curPos + kmerSize - 1;
			ovlp.extEnd = matchesList[lastMatch].extPos + kmerSize - 1;
			ovlp.score = scoreTable[lastMatch] - scoreTable[firstMatch] + kmerSize - 1;

			if (this->overlapTest(ovlp, outSuggestChimeric))
			{
				if (_keepAlignment)
				{
					//kmerMatches.emplace_back(ovlp.curBegin, ovlp.extBegin);
					std::reverse(kmerMatches.begin(), kmerMatches.end());
					//kmerMatches.emplace_back(ovlp.curEnd, ovlp.extEnd);
					ovlp.kmerMatches = kmerMatches;
				}
				ovlp.leftShift = median(shifts);
				ovlp.rightShift = extLen - curLen + ovlp.leftShift;

				int32_t filteredPositions = 0;
				for (auto pos : curFilteredPos)
				{
					if (pos < ovlp.curBegin) continue;
					if (pos > ovlp.curEnd) break;
					++filteredPositions;
				}

				float normLen = std::max(ovlp.curRange(), 
										 ovlp.extRange()) - filteredPositions;
				float matchRate = (float)chainLength * 
								  _vertexIndex.getSampleRate() / normLen;
				matchRate = std::min(matchRate, 1.0f);
				//float repeatRate = (float)filteredPositions / ovlp.curRange();
				ovlp.seqDivergence = std::log(1 / matchRate) / kmerSize;
				ovlp.seqDivergence += _estimatorBias;

				float divThreshold = _maxDivergence;
				/*if (ovlp.curBegin < _maxJump || curLen - ovlp.curEnd < _maxJump ||
					ovlp.extBegin < _maxJump || extLen - ovlp.extEnd < _maxJump)
				{
					divThreshold += _badEndAdjustment;
				}*/
				if (ovlp.seqDivergence < divThreshold)
				{
					extOverlaps.push_back(ovlp);
				}

				//collecting overlap statistics
				size_t wnd = ovlp.curBegin / STAT_WND;
				assert(wnd < divStatWindows.size());
				if (ovlp.curRange() > divStatWindows[wnd].curRange())
				{
					divStatWindows[wnd] = ovlp;
				}

				//benchmarking divergence
				/*float alnDiff = kswAlign(fastaRec.sequence, ovlp.curBegin, ovlp.curRange(),
										 _seqContainer.getSeq(extId), ovlp.extBegin, 
										  ovlp.extRange(), 1, -2, 2, 1, false);
				fout << alnDiff << " " << ovlp.seqDivergence << std::endl;*/
				/*if (0.15 > alnDiff && ovlp.seqDivergence > 0.20)
				{
					kswAlign(fastaRec.sequence
								.substr(ovlp.curBegin, ovlp.curRange()),
							 _seqContainer.getSeq(extId)
								.substr(ovlp.extBegin, ovlp.extRange()),
							 1, -2, 2, 1, true);
					std::cout << alnDiff << " " << ovlp.seqDivergence << 
						" " << (float)filteredPositions / ovlp.curRange() << std::endl;
				}*/
			}
		}
		
		//selecting the best
		if (_onlyMaxExt)
		{
			const OverlapRange* maxOvlp = nullptr;
			for (const auto& ovlp : extOverlaps)
			{
				if (!maxOvlp || ovlp.score > maxOvlp->score)
				{
					maxOvlp = &ovlp;
				}
			}
			if (maxOvlp) detectedOverlaps.push_back(*maxOvlp);
		}
		else
		{
			std::vector<OverlapRange> primOverlaps;
			//sort by decreasing score
			std::sort(extOverlaps.begin(), extOverlaps.end(),
					  [](const OverlapRange& r1, const OverlapRange& r2)
					  {return r1.score > r2.score;});
			
			for (const auto& ovlp : extOverlaps)
			{
				bool isContained = false;
				for (const auto& prim : primOverlaps)
				{
					if (ovlp.containedBy(prim) && prim.score > ovlp.score)
					{
						isContained = true;
						break;
					}
				}
				if (!isContained)
				{
					primOverlaps.push_back(ovlp);
				}
			}
			for (const auto& ovlp : primOverlaps)
			{
				detectedOverlaps.push_back(ovlp);
			}
		}
	}

	timeDp += std::chrono::duration_cast<std::chrono::duration<float>>
				(std::chrono::system_clock::now() - timeStart).count();
	timeStart = std::chrono::system_clock::now();

	for (const auto& ovlp : divStatWindows)
	{
		if (ovlp.curRange() > 0)
		{
			divStats.add(ovlp.seqDivergence);
		}
	}
	return detectedOverlaps;
}

bool OverlapContainer::hasSelfOverlaps(FastaRecord::Id readId)
{
	this->lazySeqOverlaps(readId);
	if (!readId.strand()) readId = readId.rc();
	return _overlapIndex.find(readId).suggestChimeric;
}


std::vector<OverlapRange> 
	OverlapContainer::quickSeqOverlaps(FastaRecord::Id readId)
{
	bool suggestChimeric;
	const FastaRecord& record = _queryContainer.getRecord(readId);
	return _ovlpDetect.getSeqOverlaps(record, suggestChimeric, _divergenceStats);
}

const std::vector<OverlapRange>&
	OverlapContainer::lazySeqOverlaps(FastaRecord::Id readId)
{
	bool flipped = !readId.strand();
	if (flipped) readId = readId.rc();
	IndexVecWrapper wrapper;

	//upsert creates default value if it does not exist
	_overlapIndex.upsert(readId, 	
		[&wrapper](IndexVecWrapper& val)
			{wrapper = val;});
	if (wrapper.cached)
	{
		return !flipped ? *wrapper.fwdOverlaps : *wrapper.revOverlaps;
	}

	//otherwise, need to compute overlaps.
	//do it for forward strand to be distinct
	bool suggestChimeric;
	const FastaRecord& record = _queryContainer.getRecord(readId);
	std::vector<int> kmerPositions;
	auto overlaps = _ovlpDetect.getSeqOverlaps(record, suggestChimeric,
											   _divergenceStats);
	overlaps.shrink_to_fit();

	std::vector<OverlapRange> revOverlaps;
	revOverlaps.reserve(overlaps.size());
	for (const auto& ovlp : overlaps) revOverlaps.push_back(ovlp.complement());

	_overlapIndex.update_fn(readId,
		[&wrapper, &overlaps, &revOverlaps, &suggestChimeric, &flipped, this]
		(IndexVecWrapper& val)
		{
			if (!val.cached)
			{
				_indexSize += overlaps.size();
				*val.fwdOverlaps = std::move(overlaps);
				*val.revOverlaps = std::move(revOverlaps);
				val.suggestChimeric = suggestChimeric;
				val.cached = true;
			}
			wrapper = val;
		});

	return !flipped ? *wrapper.fwdOverlaps : *wrapper.revOverlaps;
}

void OverlapContainer::ensureTransitivity(bool onlyMaxExt)
{
	Logger::get().debug() << "Computing transitive closure for overlaps";
	
	std::vector<FastaRecord::Id> allSeqs;
	for (const auto& seqIt : _overlapIndex.lock_table()) 
	{
		allSeqs.push_back(seqIt.first);
		allSeqs.push_back(seqIt.first.rc());
	}

	int totalOverlaps = 0;
	for (const auto& seq : allSeqs)
	{
		const auto& curOvlps = this->unsafeSeqOverlaps(seq);
		totalOverlaps += curOvlps.size();
		std::vector<OverlapRange> ovlpsToAdd;
		for (auto& curOvlp : curOvlps)
		{
			auto& extOvlps = this->unsafeSeqOverlaps(curOvlp.extId);

			if (onlyMaxExt)
			{
				bool found = false;
				for (auto& extOvlp : extOvlps)
				{
					if (extOvlp.extId == curOvlp.curId)
					{
						if (curOvlp.score > extOvlp.score)
						{
							extOvlp = curOvlp.reverse();
						}
						found = true;
						break;
					}
				}
				if (!found)
				{
					ovlpsToAdd.push_back(curOvlp.reverse());
				}
			}
			else
			{
				ovlpsToAdd.push_back(curOvlp.reverse());
			}
		}
		for (const auto& ovlp : ovlpsToAdd)
		{
			this->unsafeSeqOverlaps(ovlp.curId).push_back(ovlp);
		}
	}
}


void OverlapContainer::findAllOverlaps()
{
	//Logger::get().info() << "Finding overlaps:";
	std::vector<FastaRecord::Id> allQueries;
	for (const auto& seq : _queryContainer.iterSeqs())
	{
		allQueries.push_back(seq.id);
	}

	std::mutex indexMutex;
	std::function<void(const FastaRecord::Id&)> indexUpdate = 
	[this, &indexMutex] (const FastaRecord::Id& seqId)
	{
		this->lazySeqOverlaps(seqId);	//automatically stores overlaps
	};
	processInParallel(allQueries, indexUpdate, 
					  Parameters::get().numThreads, true);
	this->ensureTransitivity(false);

	int numOverlaps = 0;
	for (const auto& seqOvlps : _overlapIndex.lock_table()) 
	{
		numOverlaps += seqOvlps.second.fwdOverlaps->size() * 2;
	}
	Logger::get().debug() << "Found " << numOverlaps << " overlaps";

	this->filterOverlaps();

	numOverlaps = 0;
	for (const auto& seqOvlps : _overlapIndex.lock_table()) 
	{
		numOverlaps += seqOvlps.second.fwdOverlaps->size() * 2;
	}
	Logger::get().debug() << "Left " << numOverlaps 
		<< " overlaps after filtering";
}

std::vector<OverlapRange>&
	OverlapContainer::unsafeSeqOverlaps(FastaRecord::Id seqId)
{
		FastaRecord::Id normId = seqId.strand() ? seqId : seqId.rc();
		_overlapIndex.insert(normId);	//ensure it's in the table
		IndexVecWrapper wrapper = _overlapIndex.find(normId);
		return seqId.strand() ? *wrapper.fwdOverlaps : 
								*wrapper.revOverlaps;
}

//TODO: potentially might become non-symmetric after filtering
void OverlapContainer::filterOverlaps()
{
	static const int MAX_ENDS_DIFF = Parameters::get().kmerSize;

	std::vector<FastaRecord::Id> seqIds;
	for (const auto& seq : _queryContainer.iterSeqs())
	{
		seqIds.push_back(seq.id);
	}

	std::function<void(const FastaRecord::Id& seqId)> filterParallel =
	[this] (const FastaRecord::Id& seqId)
	{
		auto& overlaps = this->unsafeSeqOverlaps(seqId);
		
		SetVec<OverlapRange*> overlapSets;
		for (auto& ovlp : overlaps) 
		{
			overlapSets.push_back(new SetNode<OverlapRange*>(&ovlp));
		}
		for (size_t i = 0; i < overlapSets.size(); ++i)
		{
			for (size_t j = 0; j < overlapSets.size(); ++j)
			{
				OverlapRange& ovlpOne = *overlapSets[i]->data;
				OverlapRange& ovlpTwo = *overlapSets[j]->data;

				if (ovlpOne.extId != ovlpTwo.extId) continue;
				int curDiff = ovlpOne.curRange() - ovlpOne.curIntersect(ovlpTwo);
				int extDiff = ovlpOne.extRange() - ovlpOne.extIntersect(ovlpTwo);

				if (curDiff < MAX_ENDS_DIFF && extDiff < MAX_ENDS_DIFF) 
				{
					unionSet(overlapSets[i], overlapSets[j]);
				}
			}
		}
		auto clusters = groupBySet(overlapSets);
		std::vector<OverlapRange> newOvlps;
		for (const auto& cluster : clusters)
		{
			OverlapRange* maxOvlp = nullptr;
			for (auto& ovlp : cluster.second)
			{
				if (!maxOvlp || ovlp->score > maxOvlp->score)
				{
					maxOvlp = ovlp;
				}
			}
			newOvlps.push_back(*maxOvlp);
		}
		overlaps = std::move(newOvlps);

		std::sort(overlaps.begin(), overlaps.end(), 
				  [](const OverlapRange& o1, const OverlapRange& o2)
				  {return o1.curBegin < o2.curBegin;});

	};
	processInParallel(seqIds, filterParallel, 
					  Parameters::get().numThreads, false);
}


void OverlapContainer::estimateOverlaperParameters(const std::vector<int>& kmerPositions, std::string& outFile)
{
	Logger::get().debug() << "Estimating k-mer identity bias";

	const int NEDEED_OVERLAPS = 1000;
	const int MAX_SEQS = 1000;

	std::vector<FastaRecord::Id> readsToCheck;
	std::remove(outFile.c_str());
	bool isSam = outFile.substr(outFile.find_last_of(".") + 1) == "sam";
	std::string samFile = isSam ? outFile : (outFile.substr(0, outFile.find_last_of(".")) + ".sam");
    std::ofstream fout;
    if (!isSam) fout.open(outFile);
    std::ofstream fsamout(samFile);
	std::vector<OverlapRange> readAlignments;
	std::vector<std::string> readSamAlignments;

	std::mutex indexMutex;

	std::function<void(const FastaRecord::Id&)> alignRead =
	[this, &indexMutex, &fout, &isSam, &fsamout, &readAlignments, &readSamAlignments]
	(const FastaRecord::Id& seqId)
	{
		auto seq_overlaps = this->quickSeqOverlaps(seqId);
		if (seq_overlaps.size() == 0) return;
        size_t maxKmers = 0;
        size_t bestOvlpId = 0;
		/////synchronized part
		indexMutex.lock();
        for (size_t j = 0; j < seq_overlaps.size(); ++j) {
            auto &ovlp = seq_overlaps[j];
            readAlignments.push_back(ovlp);
            if (ovlp.kmerMatches.size() > maxKmers) {
                maxKmers = ovlp.kmerMatches.size();
                bestOvlpId = j;
            }
        }
		indexMutex.unlock();
        if (maxKmers >= Parameters::get().minKmers){
            auto &ovlp = seq_overlaps[bestOvlpId];
            const FastaRecord fastaRec = _queryContainer.getRecord(seqId);
            if (_ovlpDetect._seqContainer.seqName(ovlp.extId).front()=='+') {
                size_t alignLen = std::min(ovlp.curLen - ovlp.curBegin-1, ovlp.extLen-ovlp.extBegin-1);
                //Logger::get().info() << alignLen << " " << ovlp.curBegin << " " << ovlp.extBegin << " "<<(ovlp.extEnd - ovlp.extBegin-1) << " " << (ovlp.curEnd - ovlp.curBegin-1);
                std::string s = kswAlign(_ovlpDetect._seqContainer.getSeq(ovlp.extId), _ovlpDetect._seqContainer.seqName(ovlp.extId), ovlp.extBegin,
                                        (ovlp.extEnd - ovlp.extBegin-1), fastaRec.sequence, fastaRec.description, ovlp.curBegin, (ovlp.curEnd - ovlp.curBegin-1),
                       /*match*/ 1, /*mm*/ -2, /*gap open*/ 2,
                        /*gap ext*/ 1, true);
		        indexMutex.lock();
                readSamAlignments.push_back(s);
		        indexMutex.unlock();
            }
        }
		/////
	};
    for (size_t i = 0; i < _queryContainer.iterSeqs().size(); ++i)
    {
        readsToCheck.push_back(_queryContainer.iterSeqs()[i].id);
    }
	processInParallel(readsToCheck, alignRead,
					  Parameters::get().numThreads, true);

    for (size_t j = 0; j < readAlignments.size(); ++j) {
        auto &ovlp = readAlignments[j];
        if (!isSam) {
            fout << "\tAln\t";
            ovlp.dump(fout, _queryContainer, _ovlpDetect._seqContainer);
            fout << "\n";
            for (auto &posPair : ovlp.kmerMatches) {
                fout << posPair.curPos << "\t" << posPair.extPos << "\t" << posPair.bonus << "\n";
            }
            fout << "\n";
        }
    }

    if (!Parameters::get().polishSam) {
		FastaRecord::Id refId = readAlignments[0].extId;
		std::string refName = _ovlpDetect._seqContainer.seqName(refId).substr(1);
		int32_t refLen = _ovlpDetect._seqContainer.seqLen(refId);
        fsamout << "@SQ\tSN:" << refName << "\tLN:" << refLen << "\n";
        fsamout << "@PG\tID:tandemMapper\tPN:tandemMapper\tVN:1.0\tCL:tandemmapper.py\n";
    }
    for (size_t j = 0; j < readSamAlignments.size(); ++j) {
        auto &s = readSamAlignments[j];
        fsamout << s;
    }

    fout.close();
    fsamout.close();
	Logger::get().info() << "Alignments saved to " << outFile;
    exit(0);

}


void OverlapContainer::setRelativeDivergenceThreshold(float relThreshold)
{
	_ovlpDetect._maxDivergence = _meanTrueOvlpDiv + relThreshold;
	Logger::get().debug() << "Max divergence threshold set to " 
		<< _ovlpDetect._maxDivergence;
}


void OverlapContainer::overlapDivergenceStats()
{
	std::vector<float> ovlpDivergence(_divergenceStats.divVec.begin(),
									  _divergenceStats.divVec.begin() + 
									  		_divergenceStats.vecSize);
	const int HIST_LENGTH = 100;
	const int HIST_HEIGHT = 20;
	const float HIST_MIN = 0;
	const float HIST_MAX = 0.5;
	const float mult = HIST_LENGTH / (HIST_MAX * 100);
	std::vector<int> histogram(HIST_LENGTH, 0);
	for (float d : ovlpDivergence)
	{
		if (HIST_MIN <= d && d < HIST_MAX) 
		{
			++histogram[int(d * mult * 100)];
		}
	}
	int histMax = 1;
	int threshold = _ovlpDetect._maxDivergence * mult * 100;
	for (int freq : histogram) histMax = std::max(histMax, freq);

	std::string histString = "\n";
	for (int height = HIST_HEIGHT - 1; height >= 0; --height)
	{
		histString += "    |";
		for (int i = 0; i < HIST_LENGTH; ++i)
		{
			if ((float)histogram[i] / histMax > (float)height / HIST_HEIGHT)
			{
				histString += '*';
			}
			else if (i != threshold)
			{
				histString += ' ';
			}
			else
			{
				histString += '|';
			}
		}
		histString += '\n';
	}
	histString += "    " + std::string(HIST_LENGTH,  '-') + "\n";
	std::string footer(HIST_LENGTH, ' ');
	for (int i = 0; i < 10; ++i)
	{
		size_t startPos = i * HIST_LENGTH / 10;
		auto s = std::to_string(i * 5) + "%";
		for (size_t j = 0; j < s.size(); ++j) footer[j + startPos] = s[j];
	}
	histString += "    " + footer + "\n";

	Logger::get().info() << "Median overlap divergence: " 
		<< quantile(ovlpDivergence, 50); 
	Logger::get().debug() << "Sequence divergence distribution: \n" << histString
		<< "\n    Q25 = " << std::setprecision(2)
		<< quantile(ovlpDivergence, 25) << ", Q50 = " 
		<< quantile(ovlpDivergence, 50)
		<< ", Q75 = " << quantile(ovlpDivergence, 75) << "\n"
		<< std::setprecision(6);
}


void OverlapContainer::buildIntervalTree()
{
	//Logger::get().debug() << "Building interval tree";
	std::vector<FastaRecord::Id> allSeqs;
	for (const auto& seqIt : _overlapIndex.lock_table()) 
	{
		allSeqs.push_back(seqIt.first);
		allSeqs.push_back(seqIt.first.rc());
	}

	for (const auto& seq : allSeqs)
	{
		std::vector<Interval<const OverlapRange*>> intervals;
		auto& overlaps = this->unsafeSeqOverlaps(seq);
		for (const auto& ovlp : overlaps)
		{
			intervals.emplace_back(ovlp.curBegin, ovlp.curEnd, &ovlp);
		}
		_ovlpTree[seq] = IntervalTree<const OverlapRange*>(intervals);
	}
}

std::vector<Interval<const OverlapRange*>> 
	OverlapContainer::getCoveringOverlaps(FastaRecord::Id seqId, 
								  int32_t start, int32_t end) const
{
	return _ovlpTree.at(seqId).findOverlapping(start, end);
}
