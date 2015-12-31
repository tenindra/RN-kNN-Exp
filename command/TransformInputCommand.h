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

#ifndef _TRANSFORMINPUTCOMMAND_H
#define _TRANSFORMINPUTCOMMAND_H

#include "Command.h"
#include "../common.h"

#include <vector>
#include <unordered_set>

class TransformInputCommand: public Command {

    public:
        void execute(int argc, char* argv[]);
        void showCommandUsage(std::string programName);
        void showMethodUsage(std::string method, std::string programName);
        
        // Input Types for Transform Command
        std::string const TRANSFORM_DIMACS = "dimacs";
        std::string const TRANSFORM_TPQ = "tpq";
        std::string const TRANSFORM_OSM_DIMACS = "osm_dimacs";
    
    private:
        struct Point {
            int x;
            int y;
            Point(int _x, int _y): x(_x), y(_y) {}
        };
        struct Edge {
            NodeID target;
            EdgeWeight weight;
            Edge(NodeID _target, EdgeWeight _weight): target(_target), weight(_weight) {}
        };
        
        void transformDIMACSFiles(std::string grFilePath, std::string coFilePath, std::string outputFilePrefix);
        void transformTPQFiles(std::string grFilePath, std::string coFilePath, std::string outputFilePrefix, int edgeWeightInflateFactor, int coordinateInflateFactor);
        void transformOSMDimacsFiles(std::string grFilePath, std::string coFilePath, std::string outputFilePrefix, std::string regionName, std::string subRegionName);
        void addNode(NodeID node, Coordinate x, Coordinate y);
        void addEdge(NodeID source, NodeID target, EdgeWeight weight);
        bool isValidEdge(NodeID source, NodeID target, int rawWeight, EdgeWeight weight);
        bool checkValidNode(NodeID node);
        EdgeWeight correctEdgeWeightByEuclidDist(NodeID source, NodeID target, EdgeWeight weight);
        void fixShortestPathOverflow();
        LongPathDistance computePseudoDiameter(int numRounds = 1);
        void computeCoordinateRanges();
        void computeCoordinateTranslations(std::vector<Point>& coordinates, int& minX, int& minY);
        void checkNonUniqueCoordinates();
        void outputStandardFormatFiles(std::string outputFilePrefix);
        void printTransformationStatistics();
        void removeNonUniqueNodes();
        void removeIsolatedNodes();
        void removeIsolatedRegions();
        
        // Graph Information
        unsigned int numNodes = 0, numEdges = 0;
        std::string networkName = "", edgeType = "";
        
        // Error Detected/Corrected Statistics
        int numSelfLoops = 0, numParallelEdges = 0, numNegativeEdgeWeights = 0, numEuclideanCorrections = 0, numCorrectionsOverFactorTen = 0;
        int numNonUniqueNodesRemoved = 0, numDisconnectedNodesRemoved = 0, numNonUniqueNodeEdgesRemoved = 0, numDisconnectedNodeEdgesRemoved = 0;
        int numRegionsFound = 0, numIsolatedRegionNodesRemoved = 0, numIsolatedRegionEdgesRemoved = 0;
        unsigned long long totalEdgeWeight = 0, weightAddedByEuclideanCorrections = 0;
        int xTranslation = 0, yTranslation = 0;
        bool wasDeflated = false;
        int deflationFactor = 0;
        unsigned long long oldTotalEdgeWeight = 0;

        // Graph
        std::vector<NodeID> nodes;
        std::vector<std::vector<NodeEdgeWeightPair>> neighbours;
        std::vector<CoordinatePair> coordinates;
        std::unordered_set<NodeID> nonUniqueNodes;
        
};

#endif // _TRANSFORMINPUTCOMMAND_H
