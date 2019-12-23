//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <deque>
#include <iomanip>
#include <cmath>

#include "../sequence/overlap.h"
#include "../sequence/vertex_index.h"
#include "../common/config.h"
#include "../common/disjoint_set.h"
#include "repeat_graph.h"
#include "graph_processing.h"


namespace
{
	struct Point2d
	{
		Point2d(FastaRecord::Id curId = FastaRecord::ID_NONE, 
				int32_t curPos = 0, 
			    FastaRecord::Id extId = FastaRecord::ID_NONE, 
				int32_t extPos = 0):
			curId(curId), curPos(curPos), extId(extId), extPos(extPos) {}
		
		FastaRecord::Id curId;
		int32_t curPos;
		FastaRecord::Id extId;
		int32_t extPos;
	};

	struct Point1d
	{
		Point1d(FastaRecord::Id seqId = FastaRecord::ID_NONE, int32_t pos = 0):
			seqId(seqId), pos(pos) {}
		
		FastaRecord::Id seqId;
		int32_t pos;
	};

}

bool GraphEdge::isTip() const
{
	return nodeLeft->inEdges.empty() || nodeRight->outEdges.empty();
}

std::unordered_set<GraphEdge*> GraphEdge::adjacentEdges()
{
	std::unordered_set<GraphEdge*> edges;
	for (auto& e: nodeLeft->inEdges) edges.insert(e);
	for (auto& e: nodeLeft->outEdges) edges.insert(e);
	for (auto& e: nodeRight->inEdges) edges.insert(e);
	for (auto& e: nodeRight->outEdges) edges.insert(e);
	edges.erase(this);
	return edges;
}

void RepeatGraph::build()
{
	//getting overlaps
	VertexIndex asmIndex(_asmSeqs, 
						 (int)Config::get("repeat_graph_kmer_sample"));
	asmIndex.countKmers(/*min freq*/ 1, /*genome size*/ 0);
	asmIndex.setRepeatCutoff(/*min freq*/ 1);
	asmIndex.buildIndex(/*min freq*/ 2);

	//float badEndAdj = (float)Config::get("repeat_graph_ovlp_end_adjust");
	OverlapDetector asmOverlapper(_asmSeqs, asmIndex, 
								  (int)Config::get("maximum_jump"), 
								  Parameters::get().minimumOverlap,
								  /*no overhang*/ 0, /*all overlaps*/ 0,
								  /*keep alignment*/ true, /*only max*/ false,
								  (float)Config::get("repeat_graph_ovlp_divergence"),
								  /* no bad end*/ 0, /*nucl alignment*/ true);

	OverlapContainer asmOverlaps(asmOverlapper, _asmSeqs);
	asmOverlaps.findAllOverlaps();
	asmOverlaps.buildIntervalTree();
	asmOverlaps.overlapDivergenceStats();

	this->getGluepoints(asmOverlaps);
	this->collapseTandems();
	this->initializeEdges(asmOverlaps);
	GraphProcessor proc(*this, _asmSeqs);
	proc.simplify();
	this->logEdges();
	this->updateEdgeSequences();
	//this->markChimericEdges();
}

void RepeatGraph::getGluepoints(OverlapContainer& asmOverlaps)
{
	//Convert local alignments to the set of gluing points
	//that are further used for repeat graph construction.
	//Note, that there are two steps of clustering here.
	//First, alignment endpoints are clustered wrt to their
	//projections within a single contig (X-coordinates). Then
	//each cluster is further split into subclusters based 
	//on the 2-nd projection of each clustered point (Y-coordinates).
	//Each subcluster will thus have a different Y-coordinate,
	//but they will share their X-coordinate and cluster ID
	//(this means they will be glued during repeat graph cosntruction)
	
	Logger::get().debug() << "Computing gluepoints";
	typedef SetNode<Point2d> SetPoint2d;
	std::unordered_map<FastaRecord::Id, SetVec<Point2d>> endpoints;

	//first, extract endpoints from all overlaps.
	//each point has X and Y coordinates (curSeq and extSeq)
	for (auto& seq : _asmSeqs.iterSeqs())
	{
		for (auto& ovlp : asmOverlaps.lazySeqOverlaps(seq.id))
		{
			endpoints[ovlp.curId]
				.push_back(new SetPoint2d(Point2d(ovlp.curId, ovlp.curBegin,
										  ovlp.extId, ovlp.extBegin)));
			endpoints[ovlp.curId]
				.push_back(new SetPoint2d(Point2d(ovlp.curId, ovlp.curEnd,
										  ovlp.extId, ovlp.extEnd)));
		}
	}

	//for each contig, cluster gluepoints that are close to each other
	//only cosider X coordinates for now
	for (auto& seqPoints : endpoints)
	{
		std::sort(seqPoints.second.begin(), seqPoints.second.end(),
				  [](const SetPoint2d* p1, const SetPoint2d* p2)
				  {return p1->data.curPos < p2->data.curPos;});

		for (size_t i = 0; i < seqPoints.second.size() - 1; ++i)
		{
			auto* p1 = seqPoints.second[i];
			auto* p2 = seqPoints.second[i + 1];
			if (abs(p1->data.curPos - p2->data.curPos) < _maxSeparation)
			{
				unionSet(p1, p2);
			}

		}
	}
	std::vector<SetPoint2d*> flatEndpoints;
	for (auto& seqPoints : endpoints)
	{
		flatEndpoints.insert(flatEndpoints.end(), seqPoints.second.begin(), 
							 seqPoints.second.end());
	}
	auto clusters = groupBySet(flatEndpoints);

	typedef SetNode<Point1d> SetPoint1d;
	std::unordered_map<FastaRecord::Id, SetVec<Point1d>> tempGluepoints;
	std::unordered_map<SetPoint1d*, SetPoint1d*> complements;

	//we will now split each cluster based on it's Y coordinates
	//and project these subgroups to the corresponding sequences
	for (auto& clustEndpoints : clusters)
	{
		//first, simply add projections for each point from the cluster
		FastaRecord::Id clustSeq = clustEndpoints.second.front().curId;
		if (!clustSeq.strand()) continue;	//only for forward strands

		std::vector<int32_t> positions;
		for (auto& ep : clustEndpoints.second) 
		{
			positions.push_back(ep.curPos);
		}
		int32_t clusterXpos = median(positions);

		std::vector<Point1d> clusterPoints;
		clusterPoints.emplace_back(clustSeq, clusterXpos);

		SetVec<Point2d> extCoords;
		for (auto& ep : clustEndpoints.second)
		{
			extCoords.push_back(new SetPoint2d(ep));
		}
		
		//Important part: extending set of gluing points
		//We need also add extra projections
		//for gluepoints that are inside overlaps
		//(handles situations with 'repeat hierarchy', when some
		//repeats are parts of the other bigger repeats)
		for (auto& interval : asmOverlaps.getCoveringOverlaps(clustSeq, 
											 clusterXpos - 1, clusterXpos + 1))
		{
			auto& ovlp = *interval.value;
			if (ovlp.curEnd - clusterXpos > _maxSeparation &&
				clusterXpos - ovlp.curBegin > _maxSeparation)
			{
				int32_t projectedPos = ovlp.project(clusterXpos);
				extCoords.push_back(new SetPoint2d(Point2d(clustSeq, clusterXpos,
											   	   ovlp.extId, 
											   	   projectedPos)));
			}
		}

		//Finally, cluster the projected points based on Y coordinates
		std::sort(extCoords.begin(), extCoords.end(),
				  [](const SetPoint2d* p1, const SetPoint2d* p2)
				  {return (p1->data.extId != p2->data.extId) ? 
				  		   p1->data.extId < p2->data.extId :
						   p1->data.extPos < p2->data.extPos;});
		for (size_t i = 0; i < extCoords.size() - 1; ++i)
		{
			auto* p1 = extCoords[i];
			auto* p2 = extCoords[i + 1];
			if (p1->data.extId == p2->data.extId &&
				abs(p1->data.extPos - p2->data.extPos) < _maxSeparation)
			{
				unionSet(p1, p2);
			}

		}
		auto extClusters = groupBySet(extCoords);
		//now, get coordinates for each cluster
		for (auto& extClust : extClusters)
		{
			std::vector<int32_t> positions;
			for (auto& ep : extClust.second) 
			{
				positions.push_back(ep.extPos);
			}
			int32_t clusterYpos = median(positions);

			FastaRecord::Id extSeq = extClust.second.front().extId;
			clusterPoints.emplace_back(extSeq, clusterYpos);
		}

		//We should now consider how newly generaetd clusters
		//are integrated with the existing ones, we might need to
		//merge some of them together
		std::vector<SetPoint1d*> toMerge;
		for (auto& clustPt : clusterPoints)
		{
			int32_t seqLen = _asmSeqs.seqLen(clustPt.seqId);
			Point1d complPt(clustPt.seqId.rc(), seqLen - clustPt.pos - 1);

			auto& seqGluepoints = tempGluepoints[clustPt.seqId];
			auto& complGluepoints = tempGluepoints[clustPt.seqId.rc()];

			//inserting into sorted vector
			auto cmp = [] (const SetPoint1d* gp, int32_t pos)
								{return gp->data.pos < pos;};
			size_t i = std::lower_bound(seqGluepoints.begin(), 
										seqGluepoints.end(),
										clustPt.pos, cmp) - seqGluepoints.begin();
			auto cmp2 = [] (int32_t pos, const SetPoint1d* gp)
								{return pos < gp->data.pos;};
			size_t ci = std::upper_bound(complGluepoints.begin(), 
										 complGluepoints.end(),
										 complPt.pos, cmp2) - complGluepoints.begin();
			if (!seqGluepoints.empty())
			{
				if (i > 0 && 
					clustPt.pos - seqGluepoints[i - 1]->data.pos < _maxSeparation)
				{
					toMerge.push_back(seqGluepoints[i - 1]);
				}
				if (i < seqGluepoints.size() && 
					seqGluepoints[i]->data.pos - clustPt.pos < _maxSeparation)
				{
					toMerge.push_back(seqGluepoints[i]);
				}
			}

			seqGluepoints.insert(seqGluepoints.begin() + i, 
								 new SetPoint1d(clustPt));
			auto fwdPtr = seqGluepoints[i];
			complGluepoints.insert(complGluepoints.begin() + ci,
								   new SetPoint1d(complPt));
			auto revPtr = complGluepoints[ci];

			complements[fwdPtr] = revPtr;
			complements[revPtr] = fwdPtr;
			toMerge.push_back(fwdPtr);
		}
		for (size_t i = 0; i < toMerge.size() - 1; ++i)
		{
			unionSet(toMerge[i], toMerge[i + 1]);
			unionSet(complements[toMerge[i]], 
					 complements[toMerge[i + 1]]);
		}
	}
	
	//Generating final gluepoints, we might need to additionally
	//split long clusters into parts (tandem repeats)
	size_t pointId = 0;
	std::unordered_map<SetPoint1d*, size_t> setToId;
	auto addConsensusPoint = [&setToId, this, &pointId]
		(const std::vector<SetPoint1d*>& group)
	{
		SetPoint1d* reprPoint = group.front();
		if (!setToId.count(findSet(reprPoint)))
		{
			setToId[findSet(reprPoint)] = pointId++;
		}
		int32_t clusterSize = group.back()->data.pos - 
							  group.front()->data.pos;

		//big cluster corresponding to a tandem repeat - 
		//split it into multiple short edges
		if (clusterSize > _maxSeparation)
		{
			_gluePoints[reprPoint->data.seqId]
				.emplace_back(setToId[findSet(reprPoint)],
							  reprPoint->data.seqId, group.front()->data.pos);

			int32_t repeats = std::floor(clusterSize / _maxSeparation);
			int32_t mode = clusterSize / repeats;
			for (int32_t i = 1; i < repeats; ++i)
			{
				int32_t pos = group.front()->data.pos + mode * i;
				_gluePoints[reprPoint->data.seqId]
					.emplace_back(setToId[findSet(reprPoint)],
								  reprPoint->data.seqId, pos);

			}

			_gluePoints[reprPoint->data.seqId]
				.emplace_back(setToId[findSet(reprPoint)],
							  reprPoint->data.seqId, group.back()->data.pos);
		}
		//"normal" endpoint - just take a consensus
		else
		{
			std::vector<int32_t> positions;
			for (auto& ep : group) 
			{
				positions.push_back(ep->data.pos);
			}
			int32_t clusterXpos = median(positions);

			_gluePoints[reprPoint->data.seqId]
				.emplace_back(setToId[findSet(reprPoint)],
							  reprPoint->data.seqId, clusterXpos);
		}

	};
	for (auto& seqGluepoints : tempGluepoints)
	{
		std::vector<SetPoint1d*> currentGroup;
		for (auto& gp : seqGluepoints.second)
		{
			if (currentGroup.empty() || 
				gp->data.pos - currentGroup.back()->data.pos < _maxSeparation)
			{
				currentGroup.push_back(gp);
			}
			else
			{
				addConsensusPoint(currentGroup);
				currentGroup.clear();
				currentGroup.push_back(gp);
			}
		}
		if (!currentGroup.empty())
		{
			addConsensusPoint(currentGroup);
		}
	}

	//ensure that coordinates ont forward and reverse contig copies are symmetric
	for (auto& seq : _asmSeqs.iterSeqs())
	{
		if (!seq.id.strand()) continue;

		auto& seqPoints = _gluePoints[seq.id];
		auto& complPoints = _gluePoints[seq.id.rc()];

		int32_t seqLen = _asmSeqs.seqLen(seq.id);
		for (size_t i = 0; i < seqPoints.size(); ++i)
		{
			complPoints[seqPoints.size() - i - 1].position = 
				seqLen - seqPoints[i].position - 1;
		}
	}

	//This time we are making sure again that
	//final glue points "valid" - have proper "projections"
	//against each overlap
	this->checkGluepointProjections(asmOverlaps);

	//add contig end points, if needed
	const int MAX_TIP = Parameters::get().minimumOverlap;
	for (auto& seq : _asmSeqs.iterSeqs())
	{
		if (!seq.id.strand()) continue;

		auto& seqPoints = _gluePoints[seq.id];
		auto& complPoints = _gluePoints[seq.id.rc()];

		if (seqPoints.empty() || seqPoints.front().position > MAX_TIP)
		{
			seqPoints.emplace(seqPoints.begin(), pointId++, 
							  seq.id, 0);
			complPoints.emplace_back(pointId++, seq.id.rc(),
							  		 _asmSeqs.seqLen(seq.id) - 1);
		}
		if (seqPoints.size() == 1 || 
			_asmSeqs.seqLen(seq.id) - seqPoints.back().position > MAX_TIP)
		{
			seqPoints.emplace_back(pointId++, seq.id, 
								   _asmSeqs.seqLen(seq.id) - 1);
			complPoints.emplace(complPoints.begin(), pointId++, 
							  	seq.id.rc(), 0);
		}
	}

	int numGluepoints = 0;
	for (auto& seqRec : _gluePoints) numGluepoints += seqRec.second.size();
	Logger::get().debug() << "Created " << numGluepoints << " gluepoints";
}

//This function ensures that all glue points have 
//proper projections against each overlap. Make
//as many iterations as needed
void RepeatGraph::checkGluepointProjections(const OverlapContainer& asmOverlaps)
{
	size_t MAX_ITER = 1000;
	for (size_t i = 0; i < MAX_ITER; ++i)
	{
		std::unordered_map<FastaRecord::Id, 
						   std::vector<GluePoint>> addedGluepoints;
		std::unordered_map<size_t, SetNode<size_t>*> mergedGluepoints;
		auto combinePts = [&mergedGluepoints](size_t idOne, size_t idTwo)
		{
			if (!mergedGluepoints.count(idOne))
			{
				mergedGluepoints[idOne] = new SetNode<size_t>(idOne);
			}
			if (!mergedGluepoints.count(idTwo))
			{
				mergedGluepoints[idTwo] = new SetNode<size_t>(idTwo);
			}
			unionSet(mergedGluepoints[idOne],
					 mergedGluepoints[idTwo]);
		};

		for (auto& gp : _gluePoints)
		{
			if (!gp.first.strand()) continue;

			for (size_t i = 0; i < gp.second.size(); ++i)
			{
				GluePoint pt = gp.second[i];
				GluePoint ptCompl = 
					_gluePoints[pt.seqId.rc()][gp.second.size() - i - 1];

				for (auto& interval : asmOverlaps
						.getCoveringOverlaps(gp.first, pt.position - 1,  
										   	 pt.position + 1))
				{
					auto& ovlp = *interval.value;
					auto& seqPoints = _gluePoints[ovlp.extId];

					int32_t projectedPos = ovlp.project(pt.position);
					bool isValid = false;

					auto cmp = [] (const GluePoint& gp, int32_t pos)
						{return gp.position < pos;};
					auto rBegin = std::lower_bound(seqPoints.begin(), seqPoints.end(),
										 		   projectedPos - _maxSeparation, cmp);
					auto rEnd = std::lower_bound(seqPoints.begin(), seqPoints.end(),
										 		 projectedPos + _maxSeparation, cmp);
					for (; rBegin != rEnd; ++rBegin)
					{
						size_t compPos = seqPoints.size() - (rBegin - seqPoints.begin()) - 1;
						GluePoint projCompl = _gluePoints[ovlp.extId.rc()][compPos];

						if (abs(rBegin->position - projectedPos) <=
							_maxSeparation) 
						{
							//make sure that projections have proper IDs
							if (pt.pointId != rBegin->pointId)
							{
								combinePts(pt.pointId, rBegin->pointId);
								combinePts(ptCompl.pointId, projCompl.pointId);
							}
							
							isValid = true;
						}
					}
					for (auto& otherPt : addedGluepoints[ovlp.extId])
					{
						if (abs(otherPt.position - projectedPos) <=
							_maxSeparation) isValid = true;
					}

					if (!isValid) 
					{
						int32_t seqLen = _asmSeqs.seqLen(ovlp.extId);
						GluePoint newPt(pt.pointId, ovlp.extId, projectedPos);
						GluePoint newPtCompl(ptCompl.pointId, ovlp.extId.rc(),
											 seqLen - newPt.position - 1);

						addedGluepoints[newPt.seqId].push_back(newPt);
						addedGluepoints[newPtCompl.seqId].push_back(newPtCompl);
					}
				}
			}
		}

		int totalAdded = 0;
		for (auto& ptVec : addedGluepoints)
		{
			totalAdded += ptVec.second.size();
			std::copy(ptVec.second.begin(), ptVec.second.end(),
					  std::back_inserter(_gluePoints[ptVec.first]));
			std::sort(_gluePoints[ptVec.first].begin(), 
					  _gluePoints[ptVec.first].end(),
					  [](const GluePoint& p1, const GluePoint& p2)
						{return p1.position < p2.position;});
		}
		Logger::get().debug() << "Added " << totalAdded 
			<< " gluepoint projections";

		for (auto& gp : _gluePoints)
		{
			for (auto& point : gp.second)
			{
				if (mergedGluepoints.count(point.pointId))
				{
					point.pointId = 
						findSet(mergedGluepoints[point.pointId])->data;
				}
			}
		}
		/*for (auto& point : mergedGluepoints)
		{
			Logger::get().debug() << "Merging " << point.first << " -> "
				<< findSet(point.second)->data;
		}*/
		for (auto& point : mergedGluepoints) delete point.second;

		if (!totalAdded) break;
	}
}

//Cleaning up some messy "tandem" repeats that are actually
//artifatcs of the alignment
void RepeatGraph::collapseTandems()
{
	std::unordered_map<size_t, std::unordered_set<size_t>> tandemLefts;
	std::unordered_map<size_t, std::unordered_set<size_t>> tandemRights;
	std::unordered_set<size_t> bigTandems;

	for (auto& seqPoints : _gluePoints)
	{
		size_t leftId = 0;
		size_t rightId = 0;
		while (rightId < seqPoints.second.size())
		{
			while (rightId < seqPoints.second.size() && 
				   seqPoints.second[leftId].pointId == 
				   		seqPoints.second[rightId].pointId)
			{
				++rightId;
			}

			if (rightId - leftId > 1)
			{
				size_t tandemId = seqPoints.second[leftId].pointId;
				if (seqPoints.second[rightId - 1].position - 
						seqPoints.second[leftId].position > 
						Parameters::get().minimumOverlap)
				{
					bigTandems.insert(tandemId);
					//bigTandems.insert(complPoints[tandemId]);
				}

				if (rightId < seqPoints.second.size())
				{
					tandemRights[tandemId]
						.insert(seqPoints.second[rightId].pointId);
				}
				else
				{
					tandemRights[tandemId].insert(-1);
				}
				if (leftId > 0)
				{
					tandemLefts[tandemId]
						.insert(seqPoints.second[leftId - 1].pointId);
				}
				else
				{
					tandemLefts[tandemId].insert(-1);
				}
			}
			leftId = rightId;
		}
	}

	int collapsedLeft = 0;
	int collapsedRight = 0;
	int collapsedBoth = 0;
	for (auto& seqPoints : _gluePoints)
	{
		std::vector<GluePoint> newPoints;
		size_t leftId = 0;
		size_t rightId = 0;
		while (rightId < seqPoints.second.size())
		{
			while (rightId < seqPoints.second.size() && 
				   seqPoints.second[leftId].pointId == 
				   		seqPoints.second[rightId].pointId)
			{
				++rightId;
			}

			size_t tandemId = seqPoints.second[leftId].pointId;
			//size_t complId = complPoints[tandemId];
			if (rightId - leftId == 1 || bigTandems.count(tandemId))
			{
				for (size_t i = leftId; i < rightId; ++i)
				{
					newPoints.push_back(seqPoints.second[i]);
				}
			}
			else	//see if we can collapse this tandem repeat
			{
				//making sure graph remains symmetric
				bool leftDetemined = tandemLefts[tandemId].size() == 1;
									 //tandemRights[complId].size() == 1;
				bool rightDetemined = tandemRights[tandemId].size() == 1;
									  //tandemLefts[complId].size() == 1;
				if (!leftDetemined && !rightDetemined)
				{
					for (size_t i = leftId; i < rightId; ++i)
					{
						newPoints.push_back(seqPoints.second[i]);
					}
				}
				else if(leftDetemined && !rightDetemined)
				{
					newPoints.push_back(seqPoints.second[rightId - 1]);
					++collapsedLeft;
				}
				else if (!leftDetemined && rightDetemined)
				{
					newPoints.push_back(seqPoints.second[leftId]);
					++collapsedRight;
				}
				else
				{
					int32_t newPos = (seqPoints.second[rightId - 1].position +
									  seqPoints.second[leftId].position) / 2;
					newPoints.push_back(seqPoints.second[leftId]);
					newPoints.back().position = newPos;
					++collapsedBoth;
				}
			}
			leftId = rightId;
		}
		seqPoints.second = newPoints;
	}

	Logger::get().debug() << "Tandems removed: " << collapsedLeft << " left, " 
		<< collapsedRight << " right, " << collapsedBoth << " both";
}

void RepeatGraph::initializeEdges(const OverlapContainer& asmOverlaps)
{
	Logger::get().debug() << "Initializing edges";

	typedef std::pair<GraphNode*, GraphNode*> NodePair;
	std::unordered_map<NodePair, std::vector<EdgeSequence>, 
					   pairhash> parallelSegments;
	std::unordered_map<NodePair, NodePair, pairhash> complEdges;

	std::unordered_map<size_t, GraphNode*> nodeIndex;
	auto idToNode = [&nodeIndex, this](size_t nodeId)
	{
		if (!nodeIndex.count(nodeId))
		{
			nodeIndex[nodeId] = this->addNode();
		}
		return nodeIndex[nodeId];
	};

	for (auto& seqEdgesPair : _gluePoints)
	{
		if (!seqEdgesPair.first.strand()) continue;
		if (seqEdgesPair.second.size() < 2) continue;

		FastaRecord::Id complId = seqEdgesPair.first.rc();

		if (seqEdgesPair.second.size() != _gluePoints[complId].size())
		{
			throw std::runtime_error("Graph is not symmetric");
		}

		for (size_t i = 0; i < seqEdgesPair.second.size() - 1; ++i)
		{
			GluePoint gpLeft = seqEdgesPair.second[i];
			GluePoint gpRight = seqEdgesPair.second[i + 1];

			size_t complPos = seqEdgesPair.second.size() - i - 2;
			GluePoint complLeft = _gluePoints[complId][complPos];
			GluePoint complRight = _gluePoints[complId][complPos + 1];

			GraphNode* leftNode = idToNode(gpLeft.pointId);
			GraphNode* rightNode = idToNode(gpRight.pointId);
			NodePair fwdPair = std::make_pair(leftNode, rightNode);

			GraphNode* complLeftNode = idToNode(complLeft.pointId);
			GraphNode* complRightNode = idToNode(complRight.pointId);
			NodePair revPair = std::make_pair(complLeftNode, complRightNode);

			int32_t seqLen = _asmSeqs.seqLen(gpLeft.seqId);
			parallelSegments[fwdPair].emplace_back(gpLeft.seqId, seqLen, 
												   gpLeft.position, 
							  					   gpRight.position);
			parallelSegments[revPair]
				.push_back(parallelSegments[fwdPair].back().complement());

			complEdges[fwdPair] = revPair;
			complEdges[revPair] = fwdPair;
		}
	}

	auto segIntersect = [] (const EdgeSequence& s, int32_t intBegin, 
							int32_t intEnd)
	{
		return std::max(std::min(intEnd, s.origSeqEnd) - 
						std::max(intBegin, s.origSeqStart), 0);
	};

	std::unordered_set<NodePair, pairhash> usedPairs;
	size_t singletonsFiltered = 0;
	for (auto& nodePairSeqs : parallelSegments)
	{
		if (usedPairs.count(nodePairSeqs.first)) continue;
		usedPairs.insert(complEdges[nodePairSeqs.first]);

		//creating set and building index
		SetVec<EdgeSequence*> segmentSets;
		typedef SetNode<EdgeSequence*> SetSegment;
		std::unordered_map<FastaRecord::Id, 
						   std::vector<SetSegment*>> segmentIndex;
		for (auto& seg : nodePairSeqs.second) 
		{
			segmentSets.push_back(new SetNode<EdgeSequence*>(&seg));
			segmentIndex[seg.origSeqId].push_back(segmentSets.back());
		}
		for (auto& seqSegments : segmentIndex)
		{
			std::sort(seqSegments.second.begin(), seqSegments.second.end(),
					  [](const SetSegment* s1, const SetSegment* s2)
					  {return s1->data->origSeqStart < s2->data->origSeqStart;});
		}


		//cluster segments based on their overlaps
		for (auto& setOne : segmentSets)
		{
			for (auto& interval : asmOverlaps
					.getCoveringOverlaps(setOne->data->origSeqId, 
										 setOne->data->origSeqStart,
										 setOne->data->origSeqEnd))
			{
				auto& ovlp = *interval.value;
				int32_t intersectOne = 
					segIntersect(*setOne->data, ovlp.curBegin, ovlp.curEnd);
				if (intersectOne <= 0) continue;

				auto& ss = segmentIndex[ovlp.extId];
				auto cmpBegin = [] (const SetSegment* s, int32_t pos)
								    {return s->data->origSeqStart < pos;};
				auto cmpEnd = [] (const SetSegment* s, int32_t pos)
								    {return s->data->origSeqEnd < pos;};
				auto startRange = std::lower_bound(ss.begin(), ss.end(),
												   ovlp.extBegin, cmpEnd);
				auto endRange = std::lower_bound(ss.begin(), ss.end(),
												 ovlp.extEnd, cmpBegin);
				if (endRange != ss.end()) ++endRange;
				for (;startRange != endRange; ++startRange)
				{
					auto* setTwo = *startRange;
					if (findSet(setOne) == findSet(setTwo)) continue;

					int32_t projStart = ovlp.project(setOne->data->origSeqStart);
					int32_t projEnd = ovlp.project(setOne->data->origSeqEnd);
					int32_t projIntersect =
						segIntersect(*setTwo->data, projStart, projEnd);

					if (projIntersect > setOne->data->seqLen / 2 && 
						projIntersect > setTwo->data->seqLen / 2)
					{
						unionSet(setOne, setTwo);
					}
				}
			}
		}
		auto edgeClusters = groupBySet(segmentSets);
		/*if (edgeClusters.size() > 0)
		{
			Logger::get().debug() << "Node with " << segmentSets.size() << " segments";
			Logger::get().debug() << "clusters: " << edgeClusters.size();
			for (auto& edgeClust : edgeClusters)
			{
				int sumLen = 0;
				for (auto s : edgeClust.second) sumLen += s->seqLen;
				Logger::get().debug() << "\tcl: " << edgeClust.second.size()
					<< " " << sumLen / edgeClust.second.size() << " "
					<< edgeClust.first;
				
				if (edgeClust.second.size() < 10)
				{
					for (auto s : edgeClust.second)
					{
						Logger::get().debug() << "\t\t" << _asmSeqs.seqName(s->origSeqId)
							<< " " << s->origSeqStart << " " << s->seqLen;

						//////////
						if (edgeClust.second.size() <= 3)
						{
							auto segOne = s;
							for (auto& interval : asmOverlaps
									.getCoveringOverlaps(segOne->origSeqId, segOne->origSeqStart,
														 segOne->origSeqEnd))
							{
								auto& ovlp = *interval.value;
								int32_t intersectOne = 
									segIntersect(*segOne, ovlp.curBegin, ovlp.curEnd);
								if (intersectOne < segOne->seqLen / 2) continue;

								auto& ss = segmentIndex[ovlp.extId];
								auto cmpBegin = [] (const SetSegment* s, int32_t pos)
													{return s->data->origSeqStart < pos;};
								auto cmpEnd = [] (const SetSegment* s, int32_t pos)
													{return s->data->origSeqEnd < pos;};
								auto startRange = std::lower_bound(ss.begin(), ss.end(),
																   ovlp.extBegin, cmpEnd);
								auto endRange = std::lower_bound(ss.begin(), ss.end(),
																 ovlp.extEnd, cmpBegin);
								if (endRange != ss.end()) ++endRange;
								for (;startRange != endRange; ++startRange)
								{
									auto* setTwo = *startRange;
									if (segOne->origSeqStart == setTwo->data->origSeqStart &&
										segOne->origSeqEnd == setTwo->data->origSeqEnd) continue;

									//projecting the interval endpoints
									//(overlap might be covering the actual segment)
									int32_t projStart = ovlp.project(segOne->origSeqStart);
									int32_t projEnd = ovlp.project(segOne->origSeqEnd);
									int32_t projIntersect =
										segIntersect(*setTwo->data, projStart, projEnd);
									
									if (projIntersect > 0)
									{
										Logger::get().debug() << "\t\t\t" 
											<< _asmSeqs.seqName(setTwo->data->origSeqId) << " "
											<< setTwo->data->origSeqStart << " " 
											<< setTwo->data->origSeqEnd << " "
											<< intersectOne << " " 
											<< projIntersect << " " << findSet(setTwo);
									}
								}
							}
						}
					}
					///////////
				}
				else
				{
					Logger::get().debug() << "\t\t\t...";
				}
			}
		}*/

		//add edge for each cluster
		std::vector<EdgeSequence> usedSegments;
		for (auto& edgeClust : edgeClusters)
		{
			//filtering segments that were not glued, but covered by overlaps
			if (edgeClusters.size() > 1 && edgeClust.second.size() == 1)
			{
				auto seg = edgeClust.second.front();
				bool covered = false;
				for (auto& interval : asmOverlaps
						.getCoveringOverlaps(seg->origSeqId, seg->origSeqStart,
										 	 seg->origSeqEnd))
				{
					auto& ovlp = *interval.value;
					int32_t intersect = 
						segIntersect(*seg, ovlp.curBegin, ovlp.curEnd);
					if (intersect == seg->seqLen) covered = true;
				}
				if (covered)
				{
					++singletonsFiltered;
					continue;
				}
			}

			//in case we have complement edges within the node pair
			auto& anySegment = *edgeClust.second.front();
			if (std::find(usedSegments.begin(), usedSegments.end(), anySegment) 
						  != usedSegments.end()) continue;

			GraphNode* leftNode = nodePairSeqs.first.first;
			GraphNode* rightNode = nodePairSeqs.first.second;
			GraphEdge newEdge(leftNode, rightNode, FastaRecord::Id(_nextEdgeId));
			for (auto& seg : edgeClust.second)
			{
				newEdge.seqSegments.push_back(*seg);
				usedSegments.push_back(seg->complement());
			}

			//check if it's self-complmenet
			bool selfComplement = std::find(usedSegments.begin(), 
						usedSegments.end(), anySegment) != usedSegments.end();
			newEdge.selfComplement = selfComplement;

			this->addEdge(std::move(newEdge));
			if (!selfComplement)
			{
				leftNode = complEdges[nodePairSeqs.first].first;
				rightNode = complEdges[nodePairSeqs.first].second;
				GraphEdge* complEdge = this->addEdge(GraphEdge(leftNode, rightNode, 
												FastaRecord::Id(_nextEdgeId + 1)));
				for (auto& seg : edgeClust.second)
				{
					complEdge->seqSegments.push_back(seg->complement());
				}
			}

			_nextEdgeId += 2;
		}
	}
	Logger::get().debug() << "Filtered " << singletonsFiltered 
		<< " singleton segments";
}


void RepeatGraph::logEdges()
{
	typedef std::pair<EdgeSequence*, GraphEdge*> SegEdgePair;
	std::unordered_map<FastaRecord::Id, 
					   std::vector<SegEdgePair>> sequenceEdges;
	for (auto& edge : this->iterEdges())
	{
		for (auto& segment : edge->seqSegments)
		{
			sequenceEdges[segment.origSeqId].push_back({&segment, edge});
		}
	}
	for (auto& seqEdgesPair : sequenceEdges)
	{
		std::sort(seqEdgesPair.second.begin(), seqEdgesPair.second.end(),
				  [](const SegEdgePair& s1, const SegEdgePair& s2)
				  	{return s1.first->origSeqStart < s2.first->origSeqStart;});
	}

	for (auto& seqEdgesPair : sequenceEdges)
	{
		if (!seqEdgesPair.first.strand()) continue;

		size_t begin = 0;
		while(begin < seqEdgesPair.second.size())
		{
			size_t end = begin + 1;
			while(end < seqEdgesPair.second.size() &&
				 	seqEdgesPair.second[begin].second->edgeId ==
				   	seqEdgesPair.second[end].second->edgeId) ++end;

			EdgeSequence* beginSeg = seqEdgesPair.second[begin].first;
			EdgeSequence* endSeg = seqEdgesPair.second[end - 1].first;
			GraphEdge* edge = seqEdgesPair.second[begin].second;

			std::string unique = edge->seqSegments.size() == 1 ? "*" : " ";
			std::string mult = end - begin == 1 ? "" : 
							   "(" + std::to_string(end - begin) + ")";
			Logger::get().debug() << unique << "\t" 
								  << edge->edgeId.signedId() << "\t" 
								  << _asmSeqs.seqName(beginSeg->origSeqId) << "\t"
								  << beginSeg->origSeqStart << "\t" 
								  << endSeg->origSeqEnd << "\t"
								  << endSeg->origSeqEnd - beginSeg->origSeqStart << "\t"
								  << mult;
			begin = end;
		}
	}
	Logger::get().debug() << "Total edges: " << _nextEdgeId / 2;
}

GraphPath RepeatGraph::complementPath(const GraphPath& path) const
{
	if (path.empty()) return {};

	GraphPath complEdges;
	for (auto itEdge = path.rbegin(); itEdge != path.rend(); ++itEdge)
	{
		complEdges.push_back(_idToEdge.at((*itEdge)->edgeId.rc()));
	}

	assert(!complEdges.empty());
	return complEdges;
}

GraphEdge* RepeatGraph::complementEdge(GraphEdge* edge) const
{
	return _idToEdge.at(edge->edgeId.rc());
}

GraphNode* RepeatGraph::complementNode(GraphNode* node) const
{
	if (!node->outEdges.empty())
	{
		return this->complementEdge(node->outEdges.front())->nodeRight;
	}
	else if(!node->inEdges.empty())
	{
		return this->complementEdge(node->inEdges.front())->nodeLeft;
	}
	return nullptr;
}

void RepeatGraph::storeGraph(const std::string& filename)
{
	size_t nextNodeId = 0;
	std::unordered_map<GraphNode*, size_t> nodeIds;
	for (auto& node : this->iterNodes())
	{
		nodeIds[node] = nextNodeId++;
	}

	std::ofstream fout(filename);
	if (!fout)
	{
		throw std::runtime_error("Can't open "  + filename);
	}

	for (auto& edge : this->iterEdges())
	{
		fout << "Edge\t" << edge->edgeId << "\t" 
			<< nodeIds[edge->nodeLeft] << "\t" << nodeIds[edge->nodeRight]
			<< "\t" << edge->repetitive << "\t" << edge->selfComplement 
			<< "\t" << edge->resolved << "\t" << edge->meanCoverage << "\n";

		for (auto& seg : edge->seqSegments)
		{
			fout << "\tSequence\t";
			seg.dump(fout, *_edgeSeqsContainer);
			fout << "\n";
		}
	}
}

void RepeatGraph::loadGraph(const std::string& filename)
{
	std::ifstream fin(filename);
	if (!fin)
	{
		throw std::runtime_error("Can't open "  + filename);
	}

	std::unordered_map<size_t, GraphNode*> idToNode;
	GraphEdge* currentEdge = 0;
	while(!fin.eof())
	{
		std::string buffer;
		fin >> buffer;
		if (buffer.empty()) continue;

		if (buffer == "Edge")
		{
			size_t leftNode = 0;
			size_t rightNode = 0;
			size_t edgeId = 0;
			fin >> edgeId >> leftNode >> rightNode;

			if (!idToNode.count(leftNode))
			{
				idToNode[leftNode] = this->addNode();
			}
			if (!idToNode.count(rightNode))
			{
				idToNode[rightNode] = this->addNode();
			}
			
			GraphEdge edge(idToNode[leftNode], idToNode[rightNode],
						   FastaRecord::Id(edgeId));
			fin >> edge.repetitive >> edge.selfComplement 
				>> edge.resolved >> edge.meanCoverage;
			currentEdge = this->addEdge(std::move(edge));
		}
		else if (buffer == "Sequence")
		{
			if (!currentEdge)std::runtime_error("Error parsing: " + filename);

			EdgeSequence seg;
			seg.parse(fin, *_edgeSeqsContainer);
			currentEdge->seqSegments.push_back(seg);
		}
		else throw std::runtime_error("Error parsing: " + filename);
	}
}


EdgeSequence RepeatGraph::addEdgeSequence(const DnaSequence& sequence, 
							 			  int32_t start, int32_t length,
							 			  const std::string& description)
{
	auto subSeq = sequence.substr(start, length);
	auto& newRec = _edgeSeqsContainer->addSequence(subSeq, description);
	return EdgeSequence(newRec.id, newRec.sequence.length());
}

void RepeatGraph::updateEdgeSequences()
{
	for (auto& edge : _graphEdges)
	{
		if (!edge->edgeId.strand()) continue;

		int num = 0;
		for (auto& edgeSeq : edge->seqSegments)
		{
			size_t len = edgeSeq.origSeqEnd - edgeSeq.origSeqStart;
			auto subSeq = _asmSeqs.getSeq(edgeSeq.origSeqId)
										.substr(edgeSeq.origSeqStart, len);
			std::string description = "edge_" + 
				std::to_string(edge->edgeId.signedId()) + 
				"_" + std::to_string(num++) + "_" +
				_asmSeqs.getRecord(edgeSeq.origSeqId).description + "_" +
				std::to_string(edgeSeq.origSeqStart) + "_" + 
				std::to_string(edgeSeq.origSeqEnd);
			auto& newRec = _edgeSeqsContainer->addSequence(subSeq, description);

			edgeSeq.edgeSeqId = newRec.id;
			edgeSeq.seqLen = newRec.sequence.length();
		}

		if (!edge->selfComplement)
		{
			GraphEdge* complEdge = this->complementEdge(edge);
			complEdge->seqSegments.clear();
			for (auto& seq : edge->seqSegments)
			{
				complEdge->seqSegments.push_back(seq.complement());
			}
		}
	}
	_edgeSeqsContainer->buildPositionIndex();
}

RepeatGraph::~RepeatGraph()
{
	std::unordered_set<GraphEdge*> toRemove;
	for (auto& edge : this->iterEdges()) toRemove.insert(edge);
	for (auto edge : toRemove) this->removeEdge(edge);
	for (auto node : _graphNodes) delete node;
	_graphNodes.clear();
}
