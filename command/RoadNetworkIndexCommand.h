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

#ifndef _ROADNETWORKINDEXCOMMAND_H
#define _ROADNETWORKINDEXCOMMAND_H

#include "Command.h"
#include "../processing/Graph.h"

class RoadNetworkIndexCommand: public Command {

    public:
        void execute(int argc, char* argv[]);
        void showCommandUsage(std::string programName);
        void showIndexTypeUsage(std::string index, std::string programName);
        
    private:
        void buildSILC(Graph& graph, int maxQuadtreeLeafSize, std::string idxOutputFile, std::string statsOutputFile);
        void buildALT(Graph& graph, int numLandmarks, std::string landmarkMethod, std::string idxOutputFile, std::string statsOutputFile);
        void buildGtree(Graph& graph, int fanout, std::size_t maxLeafSize, std::string idxOutputFile, std::string statsOutputFile);
        void buildRouteOverlay(Graph& graph, int fanout, int levels, std::string idxOutputFile, std::string statsOutputFile);
        void buildPHL(std::string txtGrFile, std::string idxOutputFile, std::string statsOutputFile, std::string networkName, int numNodes, int numEdges);
        void buildCH(std::string graphFilePath, std::string coordinateFilePath, std::string bgrIntFilePath, 
                     std::string bcoIntFilePath, std::string chFilePath, std::string noFilePath, 
                     std::string statsOutputFile, std::string networkName, int numNodes, int numEdges);
        void buildIntermediaryFiles(std::string method, std::string bgrFile, std::string outputEdgeFile, 
                                    std::string outputNodeFIle, std::string& networkName, int& numNodes, 
                                    int& numEdges);
        void buildTNR(int gridSize, std::string bgrIntFilePath, std::string bcoIntFilePath, std::string chFilePath, std::string tnrFilePath, 
                      std::string statsOutputFile, std::string networkName, int numNodes, int numEdges);
};

#endif // _ROADNETWORKINDEXCOMMAND_H