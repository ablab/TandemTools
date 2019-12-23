//(c) 2016-2017 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include "output_generator.h"
#include "../sequence/consensus_generator.h"
#include <iomanip>


//Generates FASTA from the given graph paths
std::vector<FastaRecord> OutputGenerator::
	generatePathSequences(const std::vector<UnbranchingPath>& paths) const
{
	std::vector<FastaRecord> contigSequences;

	for (auto& contig : paths)
	{
		//As each edge might correspond to multiple sequences,
		//we need to select them so as to minimize the
		//number of original contigs (that were used to build the graph)
		std::unordered_map<FastaRecord::Id, int> seqIdFreq;
		for (auto& edge : contig.path) 
		{
			std::unordered_set<FastaRecord::Id> edgeSeqIds;
			for (auto& seg: edge->seqSegments) 
			{
				edgeSeqIds.insert(seg.origSeqId);
			}
			for (auto& seqId : edgeSeqIds)
			{
				seqIdFreq[seqId] += 1;
			}
		}

		std::string nucSequence;
		//ContigPath contigPath;
		//contigPath.name = contig.name();
		//int32_t prevFlank = 0;
		//int32_t prevSubLength = 0;

		for (size_t i = 0; i < contig.path.size(); ++i) 
		{
			if (contig.path[i]->seqSegments.empty()) 
			{
				throw std::runtime_error("Edge without sequence");
			}

			//get the sequence with maximum frequency
			EdgeSequence* bestSegment = nullptr;
			for (auto& seg : contig.path[i]->seqSegments)
			{
				if (!bestSegment || 
					seqIdFreq[seg.origSeqId] > seqIdFreq[bestSegment->origSeqId])
				{
					bestSegment = &seg;
				}
			}
			if (bestSegment->seqLen == 0) continue;
			nucSequence += _graph.edgeSequences()
								.getSeq(bestSegment->edgeSeqId).str();
		}
		//contigParts.push_back(contigPath);
		contigSequences.push_back({DnaSequence(nucSequence), contig.name(), 
								  FastaRecord::ID_NONE});
	}

	return contigSequences;
}

void OutputGenerator::outputFasta(const std::vector<UnbranchingPath>& paths,
								  const std::string& filename)
{
	std::vector<UnbranchingPath> posStrandPaths;
	for (auto& path : paths)
	{
		if (path.id.strand()) posStrandPaths.push_back(path);
	}
	SequenceContainer::writeFasta(this->generatePathSequences(posStrandPaths), 
								  filename);
}

void OutputGenerator::outputGfa(const std::vector<UnbranchingPath>& paths,
							    const std::string& filename)
{
	auto sequences = this->generatePathSequences(paths);

	Logger::get().debug() << "Writing Gfa";
	FILE* fout = fopen(filename.c_str(), "w");
	if (!fout) throw std::runtime_error("Can't open " + filename);

	fprintf(fout, "H\tVN:Z:1.0\n");
	for (size_t i = 0; i < paths.size(); ++i)
	{
		if (!paths[i].id.strand()) continue;

		//size_t kmerCount = sequences[i].sequence.length() * paths[i].meanCoverage;
		fprintf(fout, "S\t%s\t%s\tdp:i:%d\n", paths[i].name().c_str(), 
				sequences[i].sequence.str().c_str(), (int)paths[i].meanCoverage);
	}

	for (auto& contigLeft : paths)
	{
		for (auto& contigRight : paths)
		{
			if (contigLeft.path.back()->nodeRight != 
				contigRight.path.front()->nodeLeft) continue;

			std::string leftSign = contigLeft.id.strand() ? "+" :"-";
			std::string leftName = contigLeft.nameUnsigned();

			std::string rightSign = contigRight.id.strand() ? "+" :"-";
			std::string rightName = contigRight.nameUnsigned();

			fprintf(fout, "L\t%s\t%s\t%s\t%s\t0M\n", leftName.c_str(), 
					leftSign.c_str(), rightName.c_str(), rightSign.c_str());
		}
	}
}

void OutputGenerator::outputDot(const std::vector<UnbranchingPath>& paths,
								const std::string& filename)
{
	Logger::get().debug() << "Writing Dot";

	std::ofstream fout(filename);
	if (!fout.is_open()) throw std::runtime_error("Can't open " + filename);

	fout << "digraph {\n";
	fout << "nodesep = 0.5;\n";
	fout << "node [shape = circle, label = \"\", height = 0.3];\n";
	
	///re-enumerating helper functions
	std::unordered_map<GraphNode*, int> nodeIds;
	int nextNodeId = 0;
	auto nodeToId = [&nodeIds, &nextNodeId](GraphNode* node)
	{
		if (!nodeIds.count(node))
		{
			nodeIds[node] = nextNodeId++;
		}
		return nodeIds[node];
	};

	for (auto& node : _graph.iterNodes())
	{
		if (node->isTelomere())
		{
			fout << "\"" << nodeToId(node) 
				<< "\" [style = \"filled\", fillcolor = \"grey\"];\n";
		}
	}

	//coloring repeat clusters
	const std::string COLORS[] = {"red", "darkgreen", "blue", "goldenrod", 
								  "cadetblue1", "darkorchid", "aquamarine1", 
								  "darkgoldenrod1", "deepskyblue1", 
								  "darkolivegreen3"};
	std::vector<GraphEdge*> dfsStack;
	std::unordered_set<GraphEdge*> visited;
	std::unordered_map<GraphEdge*, std::string> edgeColors;
	size_t colorId = 0;
	for (auto& edge: _graph.iterEdges())
	{
		if (!edge->isRepetitive() || visited.count(edge)) continue;
		dfsStack.push_back(edge);
		while(!dfsStack.empty())
		{
			auto curEdge = dfsStack.back(); 
			dfsStack.pop_back();
			if (visited.count(curEdge)) continue;
			edgeColors[curEdge] = COLORS[colorId];
			edgeColors[_graph.complementEdge(curEdge)] = COLORS[colorId];
			visited.insert(curEdge);
			visited.insert(_graph.complementEdge(curEdge));
			for (auto adjEdge: curEdge->adjacentEdges())
			{
				if (adjEdge->isRepetitive() && !visited.count(adjEdge))
				{
					dfsStack.push_back(adjEdge);
				}
			}
		}
		colorId = (colorId + 1) % (sizeof(COLORS) / sizeof(COLORS[0]));
	}

	for (auto& contig : paths)
	{
		std::stringstream lengthStr;
		if (contig.length < 5000)
		{
			lengthStr << std::fixed << std::setprecision(1) 
				<< (float)contig.length / 1000 << "k";
		}
		else
		{
			lengthStr << contig.length / 1000 << "k";
		}
		lengthStr << " " << contig.meanCoverage << "x";

		if (contig.repetitive)
		{
			std::string color = edgeColors[contig.path.front()];
			std::string direction = contig.path.front()->selfComplement ?
									", dir = both" : "";
			fout << "\"" << nodeToId(contig.path.front()->nodeLeft) 
				 << "\" -> \"" << nodeToId(contig.path.back()->nodeRight)
				 << "\" [label = \"id " << contig.id.signedId() << 
				 "\\l" << lengthStr.str() << "\", color = \"" 
				 << color << "\" " << ", penwidth = 3" << direction << "] ;\n";
		}
		else
		{
			fout << "\"" << nodeToId(contig.path.front()->nodeLeft)
				 << "\" -> \"" << nodeToId(contig.path.back()->nodeRight)
				 << "\" [label = \"id " << contig.id.signedId()
				 << "\\l" << lengthStr.str() << "\", color = \"black\"] ;\n";
		}
	}

	fout << "}\n";
}
