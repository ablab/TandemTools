//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include "repeat_graph.h"
#include "read_aligner.h"

//A simple class that assigns edges multiplicity based on the coverage
//and copmutes the mean coverage of all edges
class MultiplicityInferer
{
public:
	MultiplicityInferer(RepeatGraph& graph, ReadAligner& aligner,
						const SequenceContainer& asmSeqs, 
						const SequenceContainer& readSeqs):
		_graph(graph), _aligner(aligner), _asmSeqs(asmSeqs), 
		_readSeqs(readSeqs), _uniqueCovThreshold(0), _meanCoverage(0) {}

	void estimateCoverage();
	int  getMeanCoverage() const {return _meanCoverage;}
	void maskUnsupportedEdges();
	void removeUnsupportedEdges();
	void removeUnsupportedConnections();

	int  collapseHeterozygousLoops(bool removeAlternatives);
	int  collapseHeterozygousBulges(bool removeAlternatives);
	int  trimTips()
	{
		int totalShort = 0;
		int totalLong = 0;
		for (;;)
		{
			int iterShort = 0;
			int iterLong = 0;
			this->trimTipsIteration(iterShort, iterLong);
			totalShort += iterShort;
			totalLong += iterLong;
			if (iterShort + iterLong == 0) break;
		}

		Logger::get().debug() << "Clipped " << totalShort 
			<< " short and " << totalLong << " long tips";
		return totalShort + totalLong;
	}

	//coverage threshold for an edge to be considered "unique"
	int  getUniqueCovThreshold() const 	{return _uniqueCovThreshold;}

private:
	void trimTipsIteration(int& outShort, int& outLong);

	RepeatGraph& _graph;
	ReadAligner& _aligner;
	const SequenceContainer& _asmSeqs;
	const SequenceContainer& _readSeqs;
	int _uniqueCovThreshold; 
	int _meanCoverage;
};
