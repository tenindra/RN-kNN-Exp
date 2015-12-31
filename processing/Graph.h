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

#ifndef _GRAPH_H
#define _GRAPH_H

#include "../common.h"

#include <vector>
#include <unordered_set>
#include <fstream>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/utility.hpp>

/*
 * Directed Graph Stored In-Memory
 */

class Graph {
    
    public:
        Graph() {}
        std::string getNetworkName();
        int getNumNodes();
        int getNumEdges();
        std::string getEdgeType();
        void setGraphData(std::string name, int numNodes, int numEdges, std::string edgeType);
        void insertNode(NodeID nodeID, Coordinate x, Coordinate y);
        void insertEdge (NodeID sourceID, NodeID targetID, EdgeWeight edgeWeight);
        void parseGraphFile(std::string grFilePath, std::string coFilePath);
        std::vector<NodeID> getNodesIDsVector();
        std::unordered_set<NodeID> getNodesIDsUset();
        Coordinate getMaxCoordinateRange();
        void getMinMaxCoordinates(Coordinate& minX, Coordinate& maxX, Coordinate& minY, Coordinate& maxY);
        NodeID getNeighbourAndEdgeWeightByPosition(NodeID u, EdgeID pos, EdgeWeight& edgeWeight);
        NodeID getNeighbourAndEdgeWeightByPosition(NodeID u, EdgeID pos, EdgeWeight& edgeWeight, NodeID& edgeIndex);
        void getCoordinates(NodeID u, Coordinate& x, Coordinate& y);
        void getTranslatedCoordinates(NodeID u, Coordinate& x, Coordinate& y, int xTranslation, int yTranslation);
        EuclideanDistance getEuclideanDistance(NodeID node1, NodeID node2);
        int getEdgeListStartIndex(NodeID node);
        int getEdgeListSize(NodeID node);
        
        // Graph Validation/Statistics Functions
        bool isConnected(NodeID node1, NodeID node2);
        bool isUndirectedGraph();
        bool isConnectedGraph();
        bool isContiguous();
        std::vector<int> countOutdegreeFrequencies();
        void printGraphStats();
        NodeID findRoadEnd(NodeID sourceNode, NodeID edgeIndex, int& outDegree, EdgeWeight& distanceToEnd);
        double getMaxGraphSpeedByEdge();
        double getMinGraphSpeedByEdge();
        
        // Convert Graph to Other Formats and Output
        void outputToTSVFile(std::string outputFilePath);
        void outputDIMACSFiles(std::string grFilePath, std::string coFilePath);
        void outputToDDSGFile(std::string outputFilePath);
        void outputZeroedCoordinatesToFile(std::string outputFilePath);
        Graph createSubGraph(std::string subgraphName, std::unordered_set<NodeID>& subgraphNodes);
        
        // INE Functions
        void parseObjectFile(std::string objectSetFilePath, std::string& setType, double& setDensity, int& setVariable, int& setSize);
        void parseObjectSet(std::vector<NodeID>& set);
        void resetAllObjects();
        bool isObject(NodeID nodeID);
        
        // Memory Usage Functions
        double computeIndexSize();
        double computeINEIndexSize();
        double computeMemoryUsage();
        
        // For each node in the graph (i.e. from 0 to numNodes), this vector
        // has the index of it's first edge in the edges graph. It also has one  The last
        // additional element with the size of the edges vector (so if the graph has
        // n nodes then this vector has n+1 elements).
        // Note: This is to stop us from falling off the edges vector when iterating.
        std::vector<unsigned int> firstEdgeIndex; 
        std::vector<NodeEdgeWeightPair> edges;
        std::vector<CoordinatePair> coordinates; // x is first, y is second (in pair)
        
    private:
        friend class boost::serialization::access;
        
        std::string networkName = "Unknown";
        unsigned int numNodes;
        unsigned int numEdges;
        std::string type;
        
        Coordinate minXCoord = 0, maxXCoord = 0, minYCoord = 0, maxYCoord = 0;
        NodeID minXNode = 0, maxXNode = 0, minYNode = 0, maxYNode = 0;
        double maxGraphSpeedByEdge = 0, minGraphSpeedByEdge = std::numeric_limits<double>::max();

        // Non-Serialized Members
        std::vector<bool> nodeObjects; // Used at run-time only
        
        // Boost Serialization
        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            ar & this->networkName;
            ar & this->numNodes;
            ar & this->numEdges;
            ar & this->type;
            ar & this->firstEdgeIndex;
            ar & this->edges;
            ar & this->coordinates;
            ar & this->minXCoord;
            ar & this->maxXCoord;
            ar & this->minYCoord;
            ar & this->maxYCoord;
            ar & this->minXNode;
            ar & this->maxXNode;
            ar & this->minYNode;
            ar & this->maxYNode;
            ar & this->maxGraphSpeedByEdge;
            ar & this->minGraphSpeedByEdge;
        }    
        
};

#endif // _GRAPH_H
