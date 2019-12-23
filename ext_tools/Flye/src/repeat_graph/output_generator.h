//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

//This module generates multiple types of output given the graph:
//it output edges sequences in FASTA format, the graph structure
//as dot or gfa. Also, output information about the unresolved
//repeats for the subsequent alasysis


#pragma once

#include "repeat_graph.h"
#include "graph_processing.h"

class OutputGenerator
{
public:
	OutputGenerator(RepeatGraph& graph, const ReadAligner& aligner,
				    const SequenceContainer& readSeqs):
		_graph(graph), _aligner(aligner), _readSeqs(readSeqs) {}

	void outputDot(const std::vector<UnbranchingPath>& paths, 
				   const std::string& filename);
	void outputGfa(const std::vector<UnbranchingPath>& paths, 
				   const std::string& filename);
	void outputFasta(const std::vector<UnbranchingPath>& paths, 
					 const std::string& filename);
	std::vector<FastaRecord> 
		generatePathSequences(const std::vector<UnbranchingPath>& paths) const;
private:

	RepeatGraph& _graph;
	const ReadAligner& _aligner;
	//const SequenceContainer& _asmSeqs;
	const SequenceContainer& _readSeqs;
};
