/* Copyright (C) 2015 Tenindra Abeywickrama
 *
 * This file is part of Road Network kNN Experimental Evaluation.
 *
 * Road Network kNN Experimental Evaluation is free software; you can
 * redistribute it and/or modify it under the terms of the GNU Affero 
 * General Public License as published by the Free Software Foundation; 
 * either version 3 of the License, or (at your option) any later version.
 *
 * Road Network kNN Experimental Evaluation is distributed in the hope 
 * that it will be useful, but WITHOUT ANY WARRANTY; without even the 
 * implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
 * PURPOSE. See the GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public 
 * License along with Road Network kNN Experimental Evaluation; see 
 * LICENSE.txt; if not, see <http://www.gnu.org/licenses/>.
 */

#include "TransformInputCommand.h"

#include "../processing/MortonList.h"
#include "../queue/BinaryMinHeap.h"
#include "../utility/StopWatch.h"
#include "../utility/utility.h"
#include "../utility/geometry.h"
#include "../common.h"

#include <unordered_set>
#include <cctype>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <limits>
#include <regex>
#include <unordered_map>
#include <assert.h>

void TransformInputCommand::execute(int argc, char* argv[])
{
    std::string method = "";
    std::string grFilePath = "";
    std::string coFilePath = "";
    std::string outputFilePrefix = "";
    std::string regionName = "";
    std::string subRegionName = "";
    int edgeWeightInflateFactor = 0, coordinateInflateFactor = 0;
    
   /*
     * Process Command Line Arguments
     */
    int opt;
    while ((opt = getopt (argc, argv, "m:e:n:o:w:c:r:s:")) != -1) {
        switch (opt) {
            case 'm':
                method = optarg;
                break;
            case 'e':
                grFilePath = optarg;
                break;
            case 'n':
                coFilePath = optarg;
                break;
            case 'o':
                outputFilePrefix = optarg;
                break;
            case 'w':
                edgeWeightInflateFactor = std::stoi(optarg);
                break;
            case 'c':
                coordinateInflateFactor = std::stoi(optarg);
                break;
            case 'r':
                regionName = optarg;
                break;
            case 's':
                subRegionName = optarg;
                break;
            default:
                std::cerr << "Unknown option(s) provided!\n\n";
                showCommandUsage(argv[0]);
                exit(1);
        }
    }  

    // Validate Command Line Arguments
    if (argc == 5 && method != "") {
        // This is 5 so that user can just enter method to find out what parameters are required for method
        this->showMethodUsage(method,argv[0]);
        exit(1);
    } 
    
    if (grFilePath == "" || coFilePath == "" || outputFilePrefix == "") {
        std::cerr << "Invalid argument(s)!\n\n";
        this->showCommandUsage(argv[0]);
        exit(1);
    }

    if (method == TRANSFORM_DIMACS) {
        if (argc < 11) {
            // Arguments: -e <edge file> -c <node file> -o <prefix for output files>
            std::cerr << "Too few arguments!\n\n";
            this->showMethodUsage(method,argv[0]);
            exit(1);
        }
        
        StopWatch sw;
        
        sw.start();
        this->transformDIMACSFiles(grFilePath,coFilePath,outputFilePrefix);
        sw.stop();
        
        double processingTimeMs = sw.getTimeMs();   
    
        std::cout << "\nDIMACS input files processed in " << processingTimeMs << "ms" << std::endl;
            
    } else if (method == TRANSFORM_TPQ) {
        if (argc < 13) {
            // Arguments: -e <edge file> -c <node file> -o <prefix for output files> -w <weight inflate factor> -c <coordinate inflate factor>
            std::cerr << "Too few arguments!\n\n";
            this->showMethodUsage(method,argv[0]);
            exit(1);
        }
        
        if (edgeWeightInflateFactor <= 0 || coordinateInflateFactor <= 0) {
            std::cerr << "Invalid argument(s)!\n\n";
            this->showCommandUsage(argv[0]);
            exit(1);
        }
        
        if (edgeWeightInflateFactor != coordinateInflateFactor) {
            std::cout << "Inflating weights and coordinates by different factor - this may cause many Euclidean corrections" << std::endl;
        }
        
        StopWatch sw;
        
        sw.start();
        this->transformTPQFiles(grFilePath,coFilePath,outputFilePrefix,edgeWeightInflateFactor,coordinateInflateFactor);
        sw.stop();
        
        double processingTimeMs = sw.getTimeMs();   
    
        std::cout << "\nTPQ input files processed in " << processingTimeMs << "ms" << std::endl;
            
    } else if (method == TRANSFORM_OSM_DIMACS) {
        if (argc < 15) {
            // Arguments: -e <edge file> -c <node file> -o <prefix for output files> -r <region name> -s <sub-region name>
            std::cerr << "Too few arguments!\n\n";
            this->showMethodUsage(method,argv[0]);
            exit(1);
        }
        
        if (regionName == "" || subRegionName == "") {
            std::cerr << "Invalid argument(s)!\n\n";
            this->showCommandUsage(argv[0]);
            exit(1);
        }
        
        StopWatch sw;
        
        sw.start();
        this->transformOSMDimacsFiles(grFilePath,coFilePath,outputFilePrefix,regionName,subRegionName);
        sw.stop();
        
        double processingTimeMs = sw.getTimeMs();   
    
        std::cout << "\nOSM DIMACS input files processed in " << processingTimeMs << "ms" << std::endl;
            
    } else {
        std::cerr << "Invalid query method!\n\n";
        this->showCommandUsage(argv[0]);
        exit(1);
    }
        
}

void TransformInputCommand::showCommandUsage(std::string programName)
{
    std::cerr << "Usage: " << programName << " -c " + constants::TRANSFORM_INPUT_CMD + " -m <query method>\n\n"
              << "Methods Options:\n"
              << utility::getFormattedUsageString(TRANSFORM_DIMACS,"Process DIMACS files") + "\n"
              << utility::getFormattedUsageString(TRANSFORM_TPQ,"Process TPQ files") + "\n"
              << utility::getFormattedUsageString(TRANSFORM_OSM_DIMACS,"Process DIMACS files created from OSM") + "\n";
}

void TransformInputCommand::showMethodUsage(std::string method, std::string programName)
{
    if (method == TRANSFORM_DIMACS ) {
        std::cerr << "Usage: " << programName << " -c " + constants::TRANSFORM_INPUT_CMD
                  << " -m " + TRANSFORM_DIMACS + " -e <graph edge file>"
                  << " -n <graph node file> \n-o <prefix for output files>\n";
    } else if (method == TRANSFORM_TPQ ) {
        std::cerr << "Usage: " << programName << " -c " + constants::TRANSFORM_INPUT_CMD
                  << " -m " + TRANSFORM_TPQ + " -e <graph edge file>"
                  << " -n <graph node file> \n-o <prefix for output files> -w <weight inflate factor>"
                  << " -c <coordinate inflate factor>\n";
    } else if (method == TRANSFORM_OSM_DIMACS ) {
        std::cerr << "Usage: " << programName << " -c " + constants::TRANSFORM_INPUT_CMD
                  << " -m " + TRANSFORM_OSM_DIMACS + " -e <graph edge file>"
                  << " -n <graph node file> \n-o <prefix for output files>-r <region name> -s <sub-region name>\n";
    } else {
        this->showCommandUsage(programName);
    }
}

void TransformInputCommand::addNode(NodeID node, Coordinate x, Coordinate y)
{
    if (node == nodes.size()) {
        this->nodes.push_back(node);
        this->coordinates.push_back(CoordinatePair(x,y));
    } else {
        std::cerr << "Nodes are not in order" << std::endl;
        std::exit(1);
    }
}

void TransformInputCommand::addEdge(NodeID source, NodeID target, EdgeWeight weight)
{
    if (weight > 0) {
        this->neighbours[source].push_back(NodeEdgeWeightPair(target,weight));
    } else {
        std::cerr << "Edge weight is zero" << std::endl;
        std::exit(1);
    }
}

bool TransformInputCommand::checkValidNode(NodeID node)
{
    if (node >= this->nodes.size()) {
        std::cerr << "Node (" << node << ") does not exist (too large)" << std::endl;
        exit(1);
    }
    return true;
}

bool TransformInputCommand::isValidEdge(NodeID source, NodeID target, int rawWeight, EdgeWeight weight)
{
    // We assume weight has been corrected for Euclidean distance
    
    bool isValidEdge = true;
    if (source == target) {
        // Remove self-loops
        ++this->numSelfLoops;
        isValidEdge = false;
    }
    
    if (isValidEdge && rawWeight < 0) {
        // Remove negative edges
        ++this->numNegativeEdgeWeights;
        isValidEdge = false;
    }
    
    // Remove parallel edges (i.e. multiple edges to same target)
    if (isValidEdge) {
        // Since weight it's positive if we reach here we cast to EdgeWeight
        for (std::size_t i = 0; i < this->neighbours[source].size(); ++i) {
            if (this->neighbours[source][i].first == target) {
                ++this->numParallelEdges;
                isValidEdge = false;
                break;
            }
        }
    }
    return isValidEdge;
}

EdgeWeight TransformInputCommand::correctEdgeWeightByEuclidDist(NodeID source, NodeID target, EdgeWeight weight)
{
    // For distance-based edge weights, edge weight cannot be 
    // smaller than euclidean distance between source and target
    if (this->edgeType == constants::DISTANCE_WEIGHTS) {
        // We ceil to make sure edge weight is always bigger than euclidean distance
        EdgeWeight minWeight = std::ceil(geometry::getEuclideanDist(this->coordinates[source].first,this->coordinates[source].second,
                                                                    this->coordinates[target].first,this->coordinates[target].second));
        if (minWeight == 0) {
            std::cerr << "There is an edge between two nodes that have same coordinates!" << std::endl;
            std::cerr << "Node " << source << " has the same coordinates as node " << target << std::endl;
            exit(1);
        }
        if (weight < minWeight) {
            this->weightAddedByEuclideanCorrections += (minWeight-weight);
            ++this->numEuclideanCorrections;
            if (weight*10 < minWeight) {
                //std::cout << "Corrected edge between node " << source << " and node " << target << " by over a factor of 10" << std::endl;
                this->numCorrectionsOverFactorTen++;
            }
            return minWeight;
        }
    }
    return weight;
}

void TransformInputCommand::fixShortestPathOverflow()
{
    StopWatch sw;
    sw.start();
    LongPathDistance estMaxSPDist = this->computePseudoDiameter(100);
    sw.stop();
    //std::cout << "Pseudo Diameter Calculated in " << sw.getTimeMs() << "ms" << std::endl;
    //std::cout << "Pseudo Diameter = " << estMaxSPDist << std::endl;
    
    // We adjust weights so that maximum possible shortest (i.e. one involving)
    // all edges fit in EdgeWeight (at it's smallest it is 32-bits)
    LongPathDistance maxPossibleSPDist = std::numeric_limits<EdgeWeight>::max();
    
    if (estMaxSPDist > maxPossibleSPDist) {
        this->wasDeflated = true;
        int deflateShift = 1;
        while ((estMaxSPDist >> deflateShift) > maxPossibleSPDist) {
            ++deflateShift;
        }
        this->deflationFactor = 1 << deflateShift;

        // We need to deflate X and Y coordinates first so Euclidean distance will match new weights
        for (std::size_t i = 0; i < this->nodes.size(); ++i) {
            this->coordinates[i].first = std::ceil(this->coordinates[i].first >> deflateShift);
            this->coordinates[i].second = std::ceil(this->coordinates[i].second>> deflateShift);
        }
        
        // Adjust all edge weights by deflation factor (but check Euclidean distances to be sure)
        unsigned long long newTotalEdgeWeight = 0;
        for (std::size_t i = 0; i < this->nodes.size(); ++i) {
            for (std::size_t j = 0; j < this->neighbours[i].size(); ++j) {
                this->neighbours[i][j].second = std::ceil(this->neighbours[i][j].second >> deflateShift);
                if (this->neighbours[i][j].second == 0) {
                    std::cerr << "Deflation has caused a zero edge weight!" << std::endl;
                    exit(1);
                }
                newTotalEdgeWeight += this->neighbours[i][j].second;
                EdgeWeight minWeight = std::ceil(geometry::getEuclideanDist(this->coordinates[i].first,this->coordinates[i].second,
                                                                            this->coordinates[this->neighbours[i][j].first].first,
                                                                            this->coordinates[this->neighbours[i][j].first].second));
                // Note: It's possible that if a Coordinate collision was caused by deflation the minWeight is
                // zero. That's fine, it will not affect the edge weight and we will remove non-unique Nodes later
                if (this->neighbours[i][j].second < minWeight) {
                    this->neighbours[i][j].second = minWeight;
                }
            }
        }
        
        this->oldTotalEdgeWeight = this->totalEdgeWeight;
        this->totalEdgeWeight = newTotalEdgeWeight;
    }
}

LongPathDistance TransformInputCommand::computePseudoDiameter(int numRounds)
{
    std::vector<bool> isNodeSettled(this->numNodes,false);
    BinaryMinHeap<LongPathDistance,NodeID> pqueue;
    
    LongPathDistance minDist, maxSPDist = 0, prevMaxSPDist = 0;
    NodeID minDistNodeID, adjNode, maxSPDistNode, maxSPDistNodeDegree = std::numeric_limits<NodeID>::max();
    
    // Initialize with priority queue with source node
    pqueue.insert(0,0);
    int i;
    maxSPDistNode = 0;
    for (i = 0; (i < numRounds || numRounds == -1) && pqueue.size() != 0; ++i) {
        maxSPDist = 0;
        while (pqueue.size() > 0) {
            // Extract and remove node with smallest distance from source
            // and mark it as "settled" so we do not inspect again
            minDist = pqueue.getMinKey();
            minDistNodeID = pqueue.extractMinElement();
            if (!isNodeSettled[minDistNodeID]) {
                isNodeSettled[minDistNodeID] = true;
                
                if (maxSPDist < minDist) {
                    maxSPDist = minDist;
                    maxSPDistNode = minDistNodeID;
                    maxSPDistNodeDegree = this->neighbours[minDistNodeID].size();
                } else if (maxSPDist == minDist && maxSPDistNodeDegree > this->neighbours[minDistNodeID].size()) {
                    // If the distances are equal we choose the node with the smaller degree
                    maxSPDistNode = minDistNodeID;
                    maxSPDistNodeDegree = this->neighbours[minDistNodeID].size();
                }

                for (std::size_t i = 0; i < this->neighbours[minDistNodeID].size(); ++i) {
                    adjNode = this->neighbours[minDistNodeID][i].first;
                    // Only update those we haven't already settled
                    if (!isNodeSettled[adjNode]) {
                        pqueue.insert(adjNode,minDist+this->neighbours[minDistNodeID][i].second);
                    }
                }
            }
        }    
        
        isNodeSettled.assign(this->numNodes,false);
        pqueue.clear();
        
        // Make the node with maximum shortest path the next node
        // unless we did not improve the max shortest path distance
        // (which case we stop irrespective of number of rounds)
        //std::cout << "New Pseudo Diameter Estimate = " << maxSPDist << std::endl;            
        if (prevMaxSPDist < maxSPDist) {
            pqueue.insert(maxSPDistNode,0);
            prevMaxSPDist = maxSPDist;
        }
    }
    //std::cout << "PseudoDiameter Calculation Rounds = " << i << std::endl;
    return prevMaxSPDist;
}

void TransformInputCommand::computeCoordinateRanges()
{
    Coordinate minX, minY, maxX, maxY;  
    
    minX = this->coordinates[0].first;
    maxX = this->coordinates[0].first;
    minY = this->coordinates[0].second;
    maxY = this->coordinates[0].second;    
    for (std::size_t i = 1; i < nodes.size(); ++i) {
        if (this->coordinates[i].first < minX) {
            minX = this->coordinates[i].first;
        }
        if (this->coordinates[i].first > maxX) {
            maxX = this->coordinates[i].first;
        }
        if (this->coordinates[i].second < minY) {
            minY = coordinates[i].second;
        }
        if (this->coordinates[i].second > maxY) {
            maxY = this->coordinates[i].second;
        }
    }

    // Take opportunity to note x and y translations required to make
    // origin (0,0) and hence use Morton numbers to check uniqueness later
    this->xTranslation = minX;
    this->yTranslation = minY;

}

void TransformInputCommand::computeCoordinateTranslations(std::vector<TransformInputCommand::Point>& coordinates, int& minX, int& minY)
{
    Coordinate testMinX = coordinates[0].x;
    Coordinate testMinY = coordinates[0].y;  
    for (std::size_t i = 1; i < nodes.size(); ++i) {
        if (coordinates[i].x < testMinX) {
            testMinX = coordinates[i].x;
        }
        if (coordinates[i].y < testMinY) {
            testMinY = coordinates[i].y;
        }
    }

    // Take opportunity to note x and y translations required to make
    // origin (0,0) and hence use Morton numbers to check uniqueness later
    minX = testMinX;
    minY = testMinY;
}

void TransformInputCommand::checkNonUniqueCoordinates()
{
    int nonUniqueCoordinates = 0;
    std::unordered_set<MortonNumber> zNumbers;    
    for (NodeID i = 0; i < this->nodes.size(); ++i) {
        MortonCode code(this->coordinates[i].first-this->xTranslation,this->coordinates[i].second-this->yTranslation,31);
        MortonNumber zNumber = code.getZNumber();
        if (zNumbers.find(zNumber) == zNumbers.end()) {
            zNumbers.insert(code.getZNumber());
        } else {
//             MortonCode testCode(zNumber,1);
//             Coordinate testX, testY;
//             int width;
//             testCode.decompose(31,testX,testY,width);
//             std::cout << "Node " << i << " has non-unique coordinates (" << testX+this->xTranslation << "," << testY+this->yTranslation << ")" << std::endl;
            this->nonUniqueNodes.insert(i);
            ++nonUniqueCoordinates;
        }
    }
    if (nonUniqueCoordinates > 0) {
        std::cout << "Found " << nonUniqueCoordinates << " non-unique coordinates" << std::endl;
    }
}

void TransformInputCommand::removeNonUniqueNodes()
{
    std::cerr << "This function doesn't work - removing nodes currently may disconnect entire subgraphs but we don't remove them" << std::endl;
    if (this->nonUniqueNodes.size() == 0) {
        // Nothing to do if there are no non-unique nodes
        return;
    }
    std::unordered_map<NodeID,NodeID> newToOldNodeID;
    std::unordered_map<NodeID,NodeID> oldToNewNodeID;
    std::vector<NodeID> newNodes;
    std::vector<std::vector<NodeEdgeWeightPair>> newNeighbours;
    std::vector<CoordinatePair> newCoordinates;
    NodeID newID = 0;
    for (NodeID i = 0; i < this->nodes.size(); ++i) {
        if (this->nonUniqueNodes.find(i) == this->nonUniqueNodes.end()) {
            // This node should not be removed so we give it a new renumbered ID
            newToOldNodeID[newID] = i;
            oldToNewNodeID[i] = newID;
            newNodes.push_back(newID);
            newCoordinates.push_back(std::move(this->coordinates[i]));
            ++newID;
            assert(newID == newNodes.size() && "Re-numbering has failed");
        } else {
            this->numNonUniqueNodesRemoved++;
            // This node's edges will be implicitly removed
            this->numNonUniqueNodeEdgesRemoved += this->neighbours[i].size();
        }
    }
    this->numNodes = newNodes.size();

    // We copy the old adjacency list for each node, excluding
    // neighbours that are non-unique (removed above)
    
    // Note: This might result in some nodes becoming disconnected, 
    // so we later call removeDisconnectedNodes() to remove them
    NodeID oldNodeID, neighbourNewID;
    this->numEdges = 0;
    for (NodeID i = 0; i < newNodes.size(); ++i) {
        assert(newToOldNodeID.find(i) != newToOldNodeID.end() && "New to old NodeID mapping doesn't exist");
        oldNodeID = newToOldNodeID[i];
        std::vector<NodeEdgeWeightPair> newAdjacencyList;
        for (std::size_t j = 0; j < this->neighbours[oldNodeID].size(); ++j) {
            if (this->nonUniqueNodes.find(this->neighbours[oldNodeID][j].first) == this->nonUniqueNodes.end()) {
                assert(oldToNewNodeID.find(this->neighbours[oldNodeID][j].first) != oldToNewNodeID.end() && "Old to new NodeID mapping doesn't exist");
                neighbourNewID = oldToNewNodeID[this->neighbours[oldNodeID][j].first];
                newAdjacencyList.push_back(NodeEdgeWeightPair(neighbourNewID,this->neighbours[oldNodeID][j].second));
            } else {
                this->numNonUniqueNodeEdgesRemoved++;
            }
        }
        this->numEdges += newAdjacencyList.size();
        newNeighbours.push_back(newAdjacencyList);
    }
    this->nodes = newNodes;
    this->neighbours = newNeighbours;
    this->coordinates = newCoordinates;
    assert(this->nodes.size() == this->neighbours.size());
    assert(this->neighbours.size() == this->coordinates.size());
}

void TransformInputCommand::removeIsolatedNodes()
{
    // Note: Removing disconnected nodes will not cause more disconnected
    // to be created (unlike removing non-unique) because by definition a
    // disconnected node is not present in any other nodes adjacency list
    std::unordered_map<NodeID,NodeID> newToOldNodeID;
    std::unordered_map<NodeID,NodeID> oldToNewNodeID;
    std::vector<NodeID> newNodes;
    std::vector<std::vector<NodeEdgeWeightPair>> newNeighbours;
    std::vector<CoordinatePair> newCoordinates;
    std::unordered_set<NodeID> disconnectedNodes;
    
    NodeID newID = 0;
    for (NodeID i = 0; i < this->nodes.size(); ++i) {
        if (this->neighbours[i].size() != 0) {
            // This is not a disconnected node
            newToOldNodeID[newID] = i;
            oldToNewNodeID[i] = newID;
            newNodes.push_back(newID);
            newCoordinates.push_back(std::move(this->coordinates[i]));
            ++newID;
            assert(newID == newNodes.size() && "Re-numbering has failed");
        } else {
            disconnectedNodes.insert(i);
            this->numDisconnectedNodesRemoved++;
            // Note: We do not need to count how many edges are removed
            // because we have already checked this node has no edges
        }
    }

    if (disconnectedNodes.size() == 0) {
        // Nothing to do, no disconnected nodes found
        return;
    }
    
    this->numNodes = newNodes.size();
    NodeID oldNodeID, neighbourNewID;
    for (NodeID i = 0; i < newNodes.size(); ++i) {
        assert(newToOldNodeID.find(i) != newToOldNodeID.end() && "New to old NodeID mapping doesn't exist");
        oldNodeID = newToOldNodeID[i];
        std::vector<NodeEdgeWeightPair> newAdjacencyList;
        for (std::size_t j = 0; j < this->neighbours[oldNodeID].size();  ++j) {
            if (disconnectedNodes.find(this->neighbours[oldNodeID][j].first) == disconnectedNodes.end()) {
                assert(oldToNewNodeID.find(this->neighbours[oldNodeID][j].first) != oldToNewNodeID.end() && "Old to new NodeID mapping doesn't exist");
                neighbourNewID = oldToNewNodeID[this->neighbours[oldNodeID][j].first];
                newAdjacencyList.push_back(NodeEdgeWeightPair(neighbourNewID,this->neighbours[oldNodeID][j].second));
            } else {
                assert(false && "Cannot reach here because disconnected node, by defintion, should not appear in any other node's adjacency list");
            }
        }
        newNeighbours.push_back(newAdjacencyList);
    }
    this->nodes = newNodes;
    this->neighbours = newNeighbours;
    this->coordinates = newCoordinates;
    assert(this->nodes.size() == this->neighbours.size());
    assert(this->neighbours.size() == this->coordinates.size());
}

void TransformInputCommand::removeIsolatedRegions()
{
    std::vector<std::unordered_set<NodeID>> nodeRegions;
    NodeID currentNode, neighbourNodeID;
    int unassignedNodes = this->numNodes;
    NodeID unassignedNode = 0;
    
    std::vector<bool> isNodeFound(this->numNodes,false);
    std::deque<NodeID> pqueue;
    std::unordered_set<NodeID> currentRegion;
    
    std::size_t biggestRegionSize = 0;
    
    while (unassignedNodes > 0) {
        // Clear previous iteration data
        this->numRegionsFound++;
        pqueue.clear();
        currentRegion.clear();
        
        pqueue.push_back(unassignedNode);
        while (pqueue.size() > 0) {
            currentNode = pqueue.front();
            pqueue.pop_front();
            if (!isNodeFound[currentNode]) {
                currentRegion.insert(currentNode);
                isNodeFound[currentNode] = true;
                unassignedNodes--;

                for (std::size_t i = 0; i < neighbours[currentNode].size(); ++i) {
                    neighbourNodeID = neighbours[currentNode][i].first;
                    // Only add those we haven't already seen
                    if (!isNodeFound[neighbourNodeID]) {
                        pqueue.push_back(neighbourNodeID);
                    }
                }
            }
        }
        
        if (currentRegion.size() > biggestRegionSize) {
            biggestRegionSize = currentRegion.size();
        }
        
        nodeRegions.push_back(std::move(currentRegion));
        
        if (unassignedNodes > 0) {
            // Find a new unassigned node (as root not of new BFS search)
            // (i.e. it is not present any of the region found so far)
            bool nodeFound = false;
            for (NodeID node = 0; node < this->numNodes; ++node) {
                if (!isNodeFound[node]) {
                    nodeFound = true;
                    unassignedNode = node;
                    break;
                }
            }
            assert(nodeFound && "There are unassigned nodes, but we could not find any one of them!");
        }
    
    }
    
    // Remove all isolated regions except the biggest region
    bool biggestRegionPreserved = false; // In case there are two regions with same biggest size
    for (std::size_t i = 0; i < nodeRegions.size(); ++i) {
        if (nodeRegions[i].size() != biggestRegionSize || biggestRegionPreserved) {
            if (nodeRegions[i].size() == biggestRegionSize) {
                biggestRegionPreserved = true;
            }
            for (NodeID node: nodeRegions[i]) {
                // Clear it's adjacency list so it will be deleted as isolated node
                // This is OK because we know all of the edge are only within this region
                this->numIsolatedRegionNodesRemoved++;
                this->numIsolatedRegionEdgesRemoved += this->neighbours[node].size();
                this->numEdges -= this->neighbours[node].size(); // this->numNodes will be decremented later
                this->neighbours[node].clear();
            }
        }
    }
    
}

void TransformInputCommand::outputStandardFormatFiles(std::string outputFilePrefix)
{
    std::string newCoFilePath = outputFilePrefix + "/" + this->networkName + "-" + this->edgeType + "." + constants::NODE_EXT;
    
    std::ofstream outputCoordFile(newCoFilePath, std::ios::out | std::ios::trunc);
    if (outputCoordFile.is_open()) {
        outputCoordFile << this->networkName << " " << this->edgeType << " " << this->numNodes << "\n";
        for (std::size_t i = 0; i < this->nodes.size(); ++i) {
            outputCoordFile << this->nodes[i] << " " << this->coordinates[i].first << " " << this->coordinates[i].second << "\n";
        }
    } else {
        std::cerr << "Cannot open coordinate output file " << newCoFilePath << std::endl;
    }
    
    std::string newGrFilePath = outputFilePrefix + "/" + networkName + "-" + edgeType + "." + constants::EDGE_EXT;
    
    std::ofstream outputGrFile(newGrFilePath, std::ios::out | std::ios::trunc);
    if (outputGrFile.is_open()) {
        outputGrFile << this->networkName << " " << this->edgeType << " " << this->numNodes << " " << this->numEdges << "\n";
        for (std::size_t i = 0; i < this->nodes.size(); ++i) {
            for (std::size_t j = 0; j < this->neighbours[i].size(); ++j) {
                outputGrFile << this->nodes[i] << " " << this->neighbours[i][j].first << " " << this->neighbours[i][j].second << "\n";
            }
        }
    } else {
        std::cerr << "Cannot open graph output file " << newGrFilePath << std::endl;
    }
}

void TransformInputCommand::printTransformationStatistics()
{
    std::cout << "\nGraph " << this->networkName << " processed with " << this->numNodes << " nodes and " << this->numEdges << " edges" << std::endl;
    std::cout << "Total edge weight = " << this->totalEdgeWeight << std::endl;
    std::cout << "Average edge weight = " << static_cast<double>(this->totalEdgeWeight)/this->numEdges << std::endl;

    std::cout << "\nStage 1: Validate Input File Graph" << std::endl;
    std::cout << "Number of self-loops removed = " << this->numSelfLoops << std::endl;
    std::cout << "Number of parallel edges removed = " << this->numParallelEdges << std::endl;        
    std::cout << "Number of negative weight edges removed = " << this->numNegativeEdgeWeights << std::endl;        
    std::cout << "Number of edges smaller than Euclidean distance corrected = " << this->numEuclideanCorrections << std::endl;        
    std::cout << "Number of edges smaller than Euclidean distance divided by 10 corrected = " << this->numCorrectionsOverFactorTen << std::endl;        
    if (this->numEuclideanCorrections != 0) {
        std::cout << "Average euclidean correction = " << static_cast<double>(this->weightAddedByEuclideanCorrections)/this->numEuclideanCorrections << std::endl;
    }

    std::cout << "\nStage 2: Correct Possible Shortest Path Distance Overflows" << std::endl;
    if (this->wasDeflated) {
        std::cout << "Sum of edge weights exceeded max possible shortest path distance value" << std::endl;
        std::cout << "Coordinate and EdgeWeight values were adjusted by deflation factor of " << this->deflationFactor << std::endl;
        std::cout << "Sum of edge weights was reduced from " << this->oldTotalEdgeWeight << " to " << this->totalEdgeWeight << std::endl;
    } else {
        std::cout << "No deflation performed, total of all edge weights below maximum possible shortest path distance value" << std::endl;
    }

    std::cout << "\nStage 3: Remove Non-Unique Nodes" << std::endl;
    std::cout << "Number of non-unique nodes removed = " << this->numNonUniqueNodesRemoved << std::endl;
    std::cout << "Number of non-unique node edges removed = " << this->numNonUniqueNodeEdgesRemoved << std::endl;

    std::cout << "\nStage 4: Disconnect Isolated Region Nodes" << std::endl;
    std::cout << "Number of isolated regions found  = " << this->numRegionsFound-1 << std::endl; 
    std::cout << "Number of isolated region nodes made disconnected (removed later as below) = " << this->numIsolatedRegionNodesRemoved << std::endl;
    std::cout << "Number of isolated region edges removed = " << this->numIsolatedRegionEdgesRemoved << std::endl;
    
    std::cout << "\nStage 5: Remove Disconnected Nodes" << std::endl;
    std::cout << "Number of disconnected nodes removed = " << this->numDisconnectedNodesRemoved << std::endl; 
    std::cout << "Number of disconnected node edges removed = " << this->numDisconnectedNodeEdgesRemoved << std::endl; 

}

void TransformInputCommand::transformDIMACSFiles(std::string grFilePath, std::string coFilePath, std::string outputFilePrefix)
{
    std::cout << "\nProcessing Coordinates File " << coFilePath << std::endl;
    
    std::string line = "", lineType = "", dummyField;
    NodeID sourceID, targetID;
    EdgeWeight weight;
    int x, y, rawWeight;
    std::string coordEdgeType;
    
    std::ifstream coordFile(coFilePath, std::ios::in);
    if (coordFile.is_open()) {
        // Note: Boost regular expression has to match whole input string
        // partial matches will not occur
        std::regex expression(".*\\-([dt])\\.([^.]+)\\.co$");
        std::smatch matches;
        if (std::regex_match(coFilePath,matches,expression)) {
            coordEdgeType = matches[1];
            this->networkName = matches[2];
        } else {
            std::cout << "Input coordinate file not named in expected DIMACS format!" << std::endl;
            exit(1);
        }
        
        if (coordEdgeType != constants::DISTANCE_WEIGHTS && coordEdgeType != constants::TIME_WEIGHTS) {
            std::cout << "Invalid edge weight type in coordinate file name" << std::endl;
            exit(1);
        }
        while (std::getline(coordFile, line))
        {
            std::stringstream ss(line);
            ss >> lineType;
            if (lineType == "p") {
                // Header line containing graph details
                ss >> dummyField >> dummyField >> dummyField >> this->numNodes;
            } else if (lineType == "v") {
                ss >> sourceID >> x >> y;
                
                // DIMACS nodes are numbered starting from 1 and they appear in order
                --sourceID;
                this->addNode(sourceID,x,y);
            }
        }
        
        if (this->numNodes == 0) {
            std::cerr << "File header information not found or does not contain total number of nodes" << std::endl;
            exit(1);
        }
        
        if (this->numNodes != this->nodes.size() || this->numNodes != this->coordinates.size()) {
            std::cerr << "Number of nodes read from file does not match number of nodes given in header" << std::endl;
            exit(1);
        }
        
    } else {
        std::cerr << "Cannot open coordinates file " << coFilePath << std::endl;
        exit(1);
    }    
    
    unsigned int edgesFileNumNodes, edgesAdded = 0;
    this->neighbours.resize(this->numNodes);
    
    std::ifstream grFile(grFilePath, std::ios::in);
    if (grFile.is_open()) {
        std::regex expression(".*\\-([dt])\\.([^.]+)\\.gr$");
        std::smatch matches;
        if (std::regex_match(grFilePath,matches,expression)) {
            this->edgeType = matches[1];
        } else {
            std::cout << "Input graph file not named in expected DIMACS format!" << std::endl;
            exit(1);
        }
        
        if (this->edgeType != constants::DISTANCE_WEIGHTS && this->edgeType != constants::TIME_WEIGHTS) {
            std::cout << "Invalid edge weight type in graph file name" << std::endl;
            exit(1);
        }

        if (this->edgeType == constants::DISTANCE_WEIGHTS) {
            std::cout << "Processing Distance-Weight Graph File " << grFilePath << std::endl;
        } else {
            std::cout << "Processing Time-Weight Graph File " << grFilePath << std::endl;
        }

        while (std::getline(grFile, line))
        {
            std::stringstream ss(line);
            ss >> lineType;
            if (lineType == "p") {
                // Header line containing graph details
                ss >> dummyField >> edgesFileNumNodes >> this->numEdges;
            } else if (lineType == "a") {
                ss >> sourceID >> targetID >> rawWeight;
                
                // DIMACS nodes are numbered starting from 1
                --sourceID;
                this->checkValidNode(sourceID);
                --targetID;
                this->checkValidNode(targetID);
                
                // Validate Edges
                weight = static_cast<EdgeWeight>(std::ceil(rawWeight));
                if (this->isValidEdge(sourceID,targetID,rawWeight,weight)) {
                    if (this->edgeType == constants::DISTANCE_WEIGHTS) {
                        // We only ensure weights are not smaller than Euclidean distances
                        // (correcting them if so) for distance edge weights
                        weight = this->correctEdgeWeightByEuclidDist(sourceID,targetID,weight);
                    }
                    this->addEdge(sourceID,targetID,weight);
                    this->totalEdgeWeight += weight;
                    ++edgesAdded;
                } else {
                    // We do not add invalid edges so we need to reduce number of edges
                    --numEdges;
                    // Note: This might result in some nodes becoming isolated (no neighbours), so we
                    // later call removeIsolatedNodes() to remove them
                }
            }
        }
        
        if (edgesFileNumNodes != this->numNodes) {
            std::cerr << "Number of nodes read from coordinates file does not match number of nodes given in edge file header" << std::endl;
            exit(1);
        }
        
        if (this->numEdges != edgesAdded) {
            std::cerr << "Number of edges read from file does not match number of edges given in header" << std::endl;
            exit(1);
        }
        
    } else {
        std::cerr << "Cannot open graph file " << grFilePath << std::endl;
        exit(1);
    }
    
    this->fixShortestPathOverflow();
    
    this->computeCoordinateRanges();
    
    this->checkNonUniqueCoordinates();
    
    //this->removeNonUniqueNodes();
    
    this->removeIsolatedRegions();
    
    this->removeIsolatedNodes();
    
    this->printTransformationStatistics();
    
    this->outputStandardFormatFiles(outputFilePrefix);
}

void TransformInputCommand::transformTPQFiles(std::string grFilePath, std::string coFilePath, std::string outputFilePrefix, int edgeWeightInflateFactor, int coordinateInflateFactor)
{
    std::cout << "\nProcessing Coordinates File " << coFilePath << std::endl;
    
    std::string line = "";
    NodeID sourceID, targetID, edgeID;
    EdgeWeight weight;
    Coordinate x, y;
    double rawX, rawY, rawWeight;
    
    std::ifstream coordFile(coFilePath, std::ios::in);
    if (coordFile.is_open()) {
        // Note: Boost regular expression has to match whole input string
        // partial matches will not occur
        std::regex expression(".*/([^/]+)\\.cnode$");
        std::smatch matches;
        if (std::regex_match(coFilePath,matches,expression)) {
            this->edgeType = constants::DISTANCE_WEIGHTS; // Edge are always distance with TPQ
            this->networkName = matches[1];
            // Make sure it is in uppercase (TPQ files are not usually)
            std::transform(this->networkName.begin(), this->networkName.end(), this->networkName.begin(), ::toupper);
        } else {
            std::cout << "Input coordinate file not named in expected TPQ format!" << std::endl;
            exit(1);
        }
        
        while (std::getline(coordFile, line))
        {
            std::stringstream ss(line);
            ss >> sourceID >> rawX >> rawY >> rawWeight;
            
            if (sourceID == this->nodes.size()) {
                x = static_cast<Coordinate>(rawX*coordinateInflateFactor);
                y = static_cast<Coordinate>(rawY*coordinateInflateFactor);
                this->addNode(numNodes,x,y);
                ++this->numNodes;
            } else {
                std::cerr << "Nodes are not in order in coordinates file" << std::endl;
            }
        }
        
    } else {
        std::cerr << "Cannot open coordinates file " << coFilePath << std::endl;
        exit(1);
    }    
    
    std::cout << "Processing Distance-Weight Graph File " << grFilePath << std::endl;

    this->neighbours.resize(this->numNodes);
    
    std::ifstream d(grFilePath, std::ios::in);
    if (d.is_open()) {
        
        while (std::getline(d, line))
        {
            std::stringstream ss(line);
            ss >> edgeID >> sourceID >> targetID >> rawWeight;

            this->checkValidNode(sourceID);
            this->checkValidNode(targetID);
            
            // Validate Edges
            weight = static_cast<EdgeWeight>(std::ceil(rawWeight*edgeWeightInflateFactor));
            if (this->isValidEdge(sourceID,targetID,rawWeight,weight)) {
                weight = this->correctEdgeWeightByEuclidDist(sourceID,targetID,weight);
                this->addEdge(sourceID,targetID,weight);
                this->totalEdgeWeight += weight;
                ++this->numEdges;
            }
            // TPQ data is undirected but only one edge is present in file so we insert reverse edge
            if (this->isValidEdge(targetID,sourceID,rawWeight,weight)) {
                this->addEdge(targetID,sourceID,weight);
                this->totalEdgeWeight += weight;
                ++this->numEdges;
            }
        }
        
    } else {
        std::cerr << "Cannot open graph file " << grFilePath << std::endl;
        exit(1);
    }

    this->fixShortestPathOverflow();
    
    this->computeCoordinateRanges();
    
    this->checkNonUniqueCoordinates();
    
    //this->removeNonUniqueNodes();
    
    this->removeIsolatedRegions();
    
    this->removeIsolatedNodes();
    
    this->printTransformationStatistics();    

    this->outputStandardFormatFiles(outputFilePrefix);
    
}

void TransformInputCommand::transformOSMDimacsFiles(std::string grFilePath, std::string coFilePath, std::string outputFilePrefix, std::string regionName, std::string subRegionName)
{
    std::cout << "\nProcessing Coordinates File " << coFilePath << std::endl;
    
    std::vector<Point> preNodeCoords;
    std::vector<std::vector<Edge>> preAdjLists;
    
    
    std::string line = "", lineType = "", dummyField;
    NodeID sourceID, targetID;
    EdgeWeight weight;
    int x, y, rawWeight;
    std::string coordEdgeType;
    
    // Preliminary Counts
    unsigned int preNodes = 0, preEdges = 0;
    
    std::ifstream coordFile(coFilePath, std::ios::in);
    if (coordFile.is_open()) {
        while (std::getline(coordFile, line))
        {
            std::stringstream ss(line);
            ss >> lineType;
            if (lineType == "p") {
                // Header line containing graph details
                ss >> dummyField >> dummyField >> dummyField >> preNodes;
                preNodeCoords.reserve(preNodes);
            } else if (lineType == "v") {
                ss >> sourceID >> x >> y;
                
                // DIMACS nodes are numbered starting from 1 and they appear in order
                --sourceID;
                if (sourceID != preNodeCoords.size()) {
                    std::cerr << "Node (" << sourceID+1 << ") node provided in order" << std::endl;
                    exit(1);
                }
                preNodeCoords.push_back(Point(x,y));
            }
        }
        
        if (preNodes == 0) {
            std::cerr << "File header information not found or does not contain total number of nodes" << std::endl;
            exit(1);
        }
        
        if (preNodes != preNodeCoords.size()) {
            std::cerr << "Number of nodes read from file does not match number of nodes given in header" << std::endl;
            exit(1);
        }
        
    } else {
        std::cerr << "Cannot open coordinates file " << coFilePath << std::endl;
        exit(1);
    }    
    
    std::cout << "\nProcessing Graph File " << grFilePath << std::endl;

    unsigned int edgesFileNumNodes, edgesAdded = 0, selfLoops = 0, parallelEdges = 0;
    preAdjLists.resize(preNodes);
    
    std::ifstream grFile(grFilePath, std::ios::in);
    if (grFile.is_open()) {
        while (std::getline(grFile, line))
        {
            std::stringstream ss(line);
            ss >> lineType;
            if (lineType == "p") {
                // Header line containing graph details
                ss >> dummyField >> edgesFileNumNodes >> preEdges;
            } else if (lineType == "a") {
                ss >> sourceID >> targetID >> rawWeight;
                
                // DIMACS nodes are numbered starting from 1
                --sourceID;
                --targetID;
                if (sourceID > preNodes || targetID > preNodes) {
                    std::cerr << "Invalid edge: Source = " << sourceID << ", Target = " << targetID << ", Weight = " << rawWeight << std::endl;
                    exit(1);
                }

                if (sourceID == targetID) {
                    // Remove self-loops
                    ++selfLoops;
                    --edgesFileNumNodes;
                    continue;
                }
    
                bool parallelEdgeFound = false;
                // Remove parallel edges (i.e. multiple edges to same target)
                for (std::size_t i = 0; i < preAdjLists[sourceID].size(); ++i) {
                    if (preAdjLists[sourceID][i].target == targetID) {
                        ++parallelEdges;
                        parallelEdgeFound = true;
                        break;
                    }
                }
                if (parallelEdgeFound) {
                    continue;
                }
                
                weight = static_cast<EdgeWeight>(std::ceil(rawWeight));
                preAdjLists[sourceID].push_back(Edge(targetID,weight));
                ++edgesAdded;
            }
        }
        if (edgesFileNumNodes != preNodes) {
            std::cerr << "Number of nodes read from coordinates file does not match number of nodes given in edge file header" << std::endl;
            exit(1);
        }
        if (edgesAdded != preEdges) {
            std::cerr << "Number of edges read from file does not match number of edges given in header" << std::endl;
            exit(1);
        }
    } else {
        std::cerr << "Cannot open graph file " << grFilePath << std::endl;
        exit(1);
    }
    
    std::cout << selfLoops << " self-loops removed in preliminary stage" << std::endl;
    std::cout << parallelEdges << " parallel edges removed in preliminary stage" << std::endl;

    std::cout << "\nCheck Undirected Edges" << std::endl;
    
    int directedEdges = 0;
    int directedEdgeWeights = 0;
    bool sourceFound;
    for (sourceID = 0; sourceID < preAdjLists.size(); ++sourceID) {
        for (std::size_t i = 0; i < preAdjLists[sourceID].size(); ++i) {
            targetID = preAdjLists[sourceID][i].target;
            sourceFound = false;
            for (std::size_t j = 0; j < preAdjLists[targetID].size(); ++j) {
                if (preAdjLists[targetID][j].target == sourceID) {
                    sourceFound = true;
                    if (preAdjLists[targetID][j].weight != preAdjLists[sourceID][i].weight) {
                        // Weights are different in each direction, choose the smaller weight
                        if (preAdjLists[targetID][j].weight > preAdjLists[sourceID][i].weight) {
                            preAdjLists[targetID][j].weight = preAdjLists[sourceID][i].weight;
                        } else {
                            preAdjLists[sourceID][i].weight = preAdjLists[targetID][j].weight;
                        }
                        directedEdgeWeights++;
                    }
                }
            }
            if (!sourceFound) {
                preAdjLists[targetID].push_back(Edge(sourceID,preAdjLists[sourceID][i].weight));
                directedEdges++;
            }
        }
    }
    std::cout << directedEdges << " directed edges made undirected" << std::endl;
    std::cout << directedEdgeWeights << " directed edges weight corrected" << std::endl;

    // Scan adjacency lists and if any node has the same coordinate as a neighbour 
    // merge them and updates edge to preserve connections

    int minX, minY, nonUniqueNeighbour = 0;
    this->computeCoordinateTranslations(preNodeCoords,minX,minY);
    MortonNumber sourceZNumber;
    
    NodeID sourceNeighbourID;
    bool merge;
    
    std::cout << "\nContracting Duplicate Coordinate Nodes"<< std::endl;
    
    std::unordered_map<MortonNumber,NodeID> mortonNumbers;

    for (sourceID = 0; sourceID < preAdjLists.size(); ++sourceID) {
        MortonCode code(preNodeCoords[sourceID].x-minX,preNodeCoords[sourceID].y-minY,31);
        sourceZNumber = code.getZNumber();
        if (mortonNumbers.find(sourceZNumber) == mortonNumbers.end()) {
            mortonNumbers[sourceZNumber] = sourceID;
        } else {
            nonUniqueNeighbour++;
            // If the coordinate is already possesed by some other node
            // then we merge together (deleting the current sourceID node)
            targetID = mortonNumbers[sourceZNumber];
            
            // Replaced the sourceID with the targetID wherever it occurs
            // Note: This may create self-loops and parallel edges (to be removed
            // when processed by above DIMACS function)
            for (std::size_t k = 0; k < preAdjLists[sourceID].size(); ++k) {
                preAdjLists[targetID].push_back(preAdjLists[sourceID][k]);
                sourceNeighbourID = preAdjLists[sourceID][k].target;
                for (std::size_t l = 0; l < preAdjLists[sourceNeighbourID].size(); ++l) {
                    if (preAdjLists[sourceNeighbourID][l].target == sourceID) {
                        preAdjLists[sourceNeighbourID][l].target = targetID;
                    }
                }
            }
            preAdjLists[sourceID].clear();
        }
    }
    std::cout << nonUniqueNeighbour << " co-located vertices removed in preliminary stage" << std::endl;
    
    NodeID newID = 0, newEdges = 0;
    std::vector<NodeID> newIDs(preAdjLists.size(),0);
    std::vector<bool> newIDSet(preAdjLists.size(),false);
    for (NodeID oldID = 0; oldID < preAdjLists.size(); ++oldID) {
        newEdges += preAdjLists[oldID].size();
        if (preAdjLists[oldID].size() != 0) {
            newIDs[oldID] = newID;
            newIDSet[oldID] = true;
            newID++;
        }
    }

    std::cout << "\nWriting DIMACS Formatted Coordinate File"<< std::endl;

    std::string newCoFilePath = outputFilePrefix + "/" + regionName + "-road-d." + subRegionName + ".co";
    std::ofstream outputCoordFile(newCoFilePath, std::ios::out | std::ios::trunc);
    if (outputCoordFile.is_open()) {
        outputCoordFile << "p aux sp co " << newID << "\n";
        for (NodeID oldID = 0; oldID < preAdjLists.size(); ++oldID) {
            if (preAdjLists[oldID].size() != 0) {
                outputCoordFile << "v " << newIDs[oldID]+1 << " " << preNodeCoords[oldID].x << " " << preNodeCoords[oldID].y << "\n";
            }
        }
    } else {
        std::cerr << "Cannot open coordinate output file " << newCoFilePath << std::endl;
    }
    
    std::cout << "\nWriting DIMACS Formatted Distance Graph File"<< std::endl;

    std::string newDgrFilePath = outputFilePrefix + "/" + regionName + "-road-d." + subRegionName + ".gr";
    std::ofstream outputDgrFile(newDgrFilePath, std::ios::out | std::ios::trunc);
    if (outputDgrFile.is_open()) {
        outputDgrFile << "p sp " << newID << " " << newEdges << "\n";
        for (NodeID oldID = 0; oldID < preAdjLists.size(); ++oldID) {
            for (std::size_t j = 0; j < preAdjLists[oldID].size(); ++j) {
                targetID = preAdjLists[oldID][j].target;
                if (!newIDSet[targetID]) {
                    std::cerr << "Node " << oldID << " has node " << targetID << " as a neighbour but it has no new node ID" << std::endl;
//                     std::exit(1);
                }
                // Ceil because it should be equal to larger than the euclidean distance
                weight = std::ceil(geometry::getEuclideanDist(preNodeCoords[oldID].x,preNodeCoords[oldID].y,
                                                              preNodeCoords[targetID].x,preNodeCoords[targetID].y));
                outputDgrFile << "a " << newIDs[oldID]+1 << " " << newIDs[targetID]+1 << " " << weight << "\n";
            }
        }
    } else {
        std::cerr << "Cannot open distance graph output file " << newDgrFilePath << std::endl;
    }    
    
    std::cout << "\nWriting DIMACS Formatted Time Graph File"<< std::endl;

    std::string newTgrFilePath = outputFilePrefix + "/" + regionName + "-road-t." + subRegionName + ".gr";
    std::ofstream outputTgrFile(newTgrFilePath, std::ios::out | std::ios::trunc);
    if (outputTgrFile.is_open()) {
        outputTgrFile << "p sp " << newID << " " << newEdges << "\n";
        for (NodeID oldID = 0; oldID < preAdjLists.size(); ++oldID) {
            for (std::size_t j = 0; j < preAdjLists[oldID].size(); ++j) {
                targetID = preAdjLists[oldID][j].target;
                // Ceil because it should be equal to larger than the euclidean distance
                outputTgrFile << "a " << newIDs[oldID]+1 << " " << newIDs[targetID]+1 << " " << preAdjLists[oldID][j].weight << "\n";
            }
        }
    } else {
        std::cerr << "Cannot open time graph output file " << newTgrFilePath << std::endl;
    }    
//     // Debug
//     for (sourceID = 0; sourceID < preAdjLists.size(); ++sourceID) {
//         for (std::size_t j = 0; j < preAdjLists[sourceID].size(); ++j) {
//             std::cout << sourceID << " " << preAdjLists[sourceID][j].target << " " << preAdjLists[sourceID][j].weight << "\n";
//         }
//     }

}
