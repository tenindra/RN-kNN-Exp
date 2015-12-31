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

#ifndef _EXPERIMENTSCOMMAND_H
#define _EXPERIMENTSCOMMAND_H

#include "Command.h"
#include "../processing/Graph.h"

#include <unordered_map>

class ExperimentsCommand: public Command {

    public:
        void execute(int argc, char* argv[]);
        void showCommandUsage(std::string programName);
        void showPhaseUsage(std::string phase, std::string programName);
        
    private:
        void buildIndexes(std::string bgrFileName, std::string parameters, std::string filePathPrefix, std::string statsOutputFile);
        void buildExternalIndexes(std::string bgrFileName, std::string parameters, std::string filePathPrefix, std::string statsOutputFile);
        std::unordered_map<std::string,std::string> getParameters(std::string parameters);
        void buildGtree(Graph& graph, int fanout, std::size_t maxLeafSize, std::string idxOutputFile, std::string statsOutputFile, std::vector<std::string> specialFields);
        void buildRouteOverlay(Graph& graph, int fanout, int levels, std::string idxOutputFile, std::string statsOutputFile, std::vector<std::string> specialFields);
        void buildSILC(Graph& graph, int maxQuadtreeLeafSize, std::string idxOutputFile, std::string statsOutputFile, std::vector<std::string> specialFields);
        void buildPHL(std::string networkName, int numNodes, int numEdges, std::string dataOutputFile, std::string idxOutputFile, std::string statsOutputFile, std::vector<std::string> specialFields);
        void buildCH(std::string networkName, int numNodes, int numEdges, std::string graphFilePath, std::string coordinateFilePath, std::string bgrIntFilePath, std::string bcoIntFilePath,
                     std::string chFilePath, std::string noFilePath, std::string statsOutputFile, std::vector<std::string> specialFields);
        void buildTNR(std::string networkName, int numNodes, int numEdges, std::string bgrIntFilePath, std::string bcoIntFilePath, std::string chFilePath, 
                      std::string tnrFilePath, int gridSize, std::string statsOutputFile, std::vector<std::string> specialFields);
        void buildALT(Graph& graph, int numLandmarks, std::string selectionMethod, std::string idxOutputFile, std::string statsOutputFile, std::vector<std::string> specialFields);
        void buildObjectIndexes(std::string bgrFileName, std::string parameters, std::size_t numSets, std::string objDensities, 
                                std::string objTypes, std::string objVariable, std::string filePathPrefix, std::string statsOutputFile, std::string queryNodeFile = "", unsigned int numQueryNodes = 0);
        void generateObjectSets(Graph& graph, std::size_t numSets, std::vector<double> objDensities, 
                                std::vector<std::string> objTypes, std::vector<int> objVariable, std::string filePathPrefix, std::string statsOutputFile, 
                                std::vector<std::string>& parameterKeys, std::vector<std::string>& parameterValues, std::string queryNodeFile = "", unsigned int numQueryNodes = 0);
        void buildOccurenceLists(std::string gtreeIdxFile, std::size_t numSets, std::vector<double> objDensities, 
                                 std::vector<std::string> objTypes, std::vector<int> objVariable, std::string filePathPrefix, std::string statsOutputFile, 
                                 std::vector<std::string>& parameterKeys, std::vector<std::string>& parameterValues);
        void buildAssociationDirectories(std::string routeOverlayIdxFile, std::size_t numSets, std::vector<double> objDensities, 
                                       std::vector<std::string> objTypes, std::vector<int> objVariable, std::string filePathPrefix, std::string statsOutputFile, 
                                       std::vector<std::string>& parameterKeys, std::vector<std::string>& parameterValues);
        void buildQuadtrees(Graph& graph, std::vector<int>& maxLeafSizes, std::size_t numSets, std::vector<double> objDensities, 
                            std::vector<std::string> objTypes, std::vector<int> objVariable, std::string filePathPrefix, std::string statsOutputFile);
        void buildRtrees(Graph& graph, std::vector<int>& branchFactors, std::size_t numSets, std::vector<double> objDensities, 
                         std::vector<std::string> objTypes, std::vector<int> objVariable, std::string filePathPrefix, std::string statsOutputFile);
        void runQueries(std::string bgrFileName, std::string queryNodeFile, std::string kValues, std::string parameters, std::size_t numSets, 
                        std::string objDensities, std::string objTypes, std::string objVariable, std::string filePathPrefix, std::string statsOutputFile);
        void runINEQueries(Graph& graph, std::vector<NodeID>& queryNodes, std::vector<int>& kValues, std::size_t numSets, 
                           std::vector<double> objDensities, std::vector<std::string> objTypes, std::vector<int> objVariable, std::string filePathPrefix, 
                           std::string statsOutputFile, std::vector<std::string> specialFields = {});
        void runGtreeQueries(Graph& graph, std::string gtreeIdxFile, std::vector<NodeID>& queryNodes, std::vector<int>& kValues, 
                             std::size_t numSets, std::vector<double> objDensities, std::vector<std::string> objTypes, std::vector<int> objVariable, std::string filePathPrefix, 
                             std::string statsOutputFile, bool verifyKNN, std::vector<std::string>& parameterKeys, std::vector<std::string>& parameterValues, 
                             std::vector<std::string> specialFields = {});
        void runROADQueries(Graph& graph, std::string routeOverlayIdxFile, std::vector<NodeID>& queryNodes, std::vector<int>& kValues, 
                            std::size_t numSets, std::vector<double> objDensities, std::vector<std::string> objTypes, std::vector<int> objVariable, std::string filePathPrefix, 
                            std::string statsOutputFile, bool verifyKNN, std::vector<std::string>& parameterKeys, std::vector<std::string>& parameterValues,
                            std::vector<std::string> specialFields = {});
        void runSILCBasedQueries(Graph& graph, std::vector<std::string>& methods, std::string silcIdxFile, std::string suppIdxFile, 
                                 std::vector<int>& maxLeafSizes, std::vector<NodeID>& queryNodes, std::vector<int>& kValues, std::size_t numSets, 
                                 std::vector<double> objDensities, std::vector<std::string> objTypes, std::vector<int> objVariable, std::string filePathPrefix, 
                                 std::string statsOutputFile, bool verifyKNN, std::vector<std::string> specialFields = {});
        void runIERQueries(Graph& graph, std::string index, std::string idxFilePath, std::vector<int>& branchFactors, 
                           std::vector<NodeID>& queryNodes, std::vector<int>& kValues, std::size_t numSets, 
                           std::vector<double> objDensities, std::vector<std::string> objTypes, std::vector<int> objVariable, std::string filePathPrefix, 
                           std::string statsOutputFile, bool verifyKNN, std::vector<std::string> specialFields = {});
        void runSingleMethodQueries(std::string method, std::string bgrFileName, std::string queryNodeFile, std::string kValues, 
                                    std::string parameters, std::size_t numSets, std::string objDensities, std::string objTypes, std::string objVariable, 
                                    std::string filePathPrefix, std::string statsOutputFile);
        void runINEQueriesByDynamicGraph(Graph& graph, std::string dynBgrFileName, std::vector<NodeID>& queryNodes, std::vector<int>& kValues, std::size_t numSets, 
                                         std::vector<double> objDensities, std::vector<std::string> objTypes, std::vector<int> objVariable, std::string filePathPrefix, 
                                         std::string statsOutputFile, bool verifyKNN, std::vector<std::string> specialFields = {});
        void runRealWorldPOIQueries(std::string bgrFileName, std::string queryNodeFile, std::string kValues, std::string parameters, 
                                    std::string filePathPrefix, std::string statsOutputFile, std::string rwPOISetListFile);
        void buildRealWorldObjIndexes(std::string bgrFileName, std::string parameters, std::string filePathPrefix,
                                      std::string statsOutputFile, std::string rwPOISetListFile);
};

#endif // _EXPERIMENTSCOMMAND_H