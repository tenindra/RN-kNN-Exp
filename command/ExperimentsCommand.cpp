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

#include "ExperimentsCommand.h"

#include "../processing/DynamicGraph.h"
#include "../processing/Gtree.h"
#include "../processing/ROAD.h"
#include "../processing/MortonList.h"
#include "../processing/INE.h"
#include "../processing/IER.h"
#include "../processing/SetGenerator.h"
#include "../processing/ALT.h"
#include "../processing/ShortestPathWrapper.h"
#include "../tuple/IndexTuple.h"
#include "../tuple/ObjectIndexTuple.h"
#include "../tuple/KnnQueryTuple.h"
#include "../common.h"
#include "../utility/StopWatch.h"
#include "../utility/Statistics.h"
#include "../utility/utility.h"
#include "../utility/serialization.h"

#include <cstdio>

void ExperimentsCommand::execute(int argc, char* argv[])
{
    std::string experiment = "";
    std::string bgrFilePath = "";
    std::string filePathPrefix = "";
    std::string statsOutputFile = "";
    std::string rwPOISetListFile = "";
    std::string parameters = "";
    std::string objDensities = "";
    std::string objTypes = "";
    std::string objVariable = "";
    unsigned int numSets = 0;
    std::string queryNodeFile = "";
    std::string kValues = "";
    std::string method = "";
    unsigned int numPoints = 0;
    
    /*
     * Process Command Line Arguments
     */
    int opt;
    while ((opt = getopt (argc, argv, "e:g:p:f:s:n:d:t:q:k:m:v:l:r:")) != -1) {
        switch (opt) {
            case 'e':
                experiment = optarg;
                break;
            case 'g':
                bgrFilePath = optarg;
                break;
            case 'p':
                parameters = optarg;
                break;
            case 'f':
                filePathPrefix = optarg;
                break;
            case 's':
                statsOutputFile = optarg;
                break;
            case 'n':
                numSets = std::stoul(optarg);
                break;
            case 'd':
                objDensities = optarg;
                break;
            case 't':
                objTypes = optarg;
                break;
            case 'q':
                queryNodeFile = optarg;
                break;
            case 'k':
                kValues = optarg;
                break;
            case 'm':
                method = optarg;
                break;
            case 'v':
                objVariable = optarg;
                break;
            case 'l':
                numPoints = std::stoul(optarg);
                break;
            case 'r':
                rwPOISetListFile = optarg;
                break;
            default:
                std::cerr << "Unknown option(s) provided!\n\n";
                showCommandUsage(argv[0]);
                exit(1);
        }
    }

    // Validate Command Line Arguments
    if (argc == 5) {
        // This is 5 so that user can just enter method to find out what parameters are required for method
        this->showPhaseUsage(experiment,argv[0]);
        exit(1);
    } 
    
    if (experiment == "") {
        std::cerr << "Invalid argument(s)!\n\n";
        this->showCommandUsage(argv[0]);
        exit(1);
    }
    
    if (experiment == constants::EXP_BUILD_INDEXES) {
        if (argc < 13) {
            // Arguments: -g <binary graph file> -p <parameters> -f <index output file path prefix> -s <stats output file>
            std::cerr << "Too few arguments!\n\n";
            this->showPhaseUsage(experiment,argv[0]);
            exit(1);
        }
    
        if (bgrFilePath == "" || parameters == "" || filePathPrefix == "" || statsOutputFile == "") {
            std::cerr << "Invalid argument(s)!\n\n";
            this->showPhaseUsage(experiment,argv[0]);
            exit(1);
        }
        
        this->buildIndexes(bgrFilePath,parameters,filePathPrefix,statsOutputFile);
        this->buildExternalIndexes(bgrFilePath,parameters,filePathPrefix,statsOutputFile);
        
    } else if (experiment == constants::EXP_BUILD_OBJ_INDEXES) {
        if (argc < 21) {
            // Arguments: -g <binary graph file>  -p <parameters> -n <num sets> -d <list of object densities> 
            // -t <list of object types> -v <list of some object variable> -f <index output file path prefix> 
            // -s <stats output file>  [-q <query node file if generating query nodes related to object set> -l <num points for query node set]]
            std::cerr << "Too few arguments!\n\n";
            this->showPhaseUsage(experiment,argv[0]);
            exit(1);
        }
    
        if (bgrFilePath == "" || parameters == "" || filePathPrefix == "" || statsOutputFile == "" 
            || numSets == 0 || objDensities == "" || objTypes == "" || objVariable == "") {
            std::cerr << "Invalid argument(s)!\n\n";
            this->showPhaseUsage(experiment,argv[0]);
            exit(1);
        }
        
        this->buildObjectIndexes(bgrFilePath,parameters,numSets,objDensities,objTypes,objVariable,filePathPrefix,statsOutputFile,queryNodeFile,numPoints);
        
    } else if (experiment == constants::EXP_RUN_KNN) {
        if (argc < 25) {
            // Arguments: -g <binary graph file> -q <query node file> -k <k values> -p <parameters> -n <num sets> -d <list of object densities> -t <list of object types> 
            // -v <list of some object variable> -f <index output file path prefix> -s <stats output file>
            std::cerr << "Too few arguments!\n\n";
            this->showPhaseUsage(experiment,argv[0]);
            exit(1);
        }
    
        if (bgrFilePath == "" || parameters == "" || filePathPrefix == "" || statsOutputFile == "" 
            || numSets == 0 || objDensities == "" || objTypes == "" || queryNodeFile == "" || kValues == "" || objVariable == "") {
            std::cerr << "Invalid argument(s)!\n\n";
            this->showPhaseUsage(experiment,argv[0]);
            exit(1);
        }
        
        this->runQueries(bgrFilePath,queryNodeFile,kValues,parameters,numSets,objDensities,objTypes,objVariable,filePathPrefix,statsOutputFile);

    } else if (experiment == constants::EXP_RUN_KNN_OPTIMIZATIONS) {
        if (argc < 27) {
            // Arguments: -m <method> -g <binary graph file> -q <query node file> 
            // -k <k values> -p <parameters> -n <num sets> -d <list of object densities> -t <list of object types> -v <list of some object variable>
            // -f <index output file path prefix> -s <stats output file>
            std::cerr << "Too few arguments!\n\n";
            this->showPhaseUsage(experiment,argv[0]);
            exit(1);
        }
    
        if (bgrFilePath == "" || parameters == "" || filePathPrefix == "" || statsOutputFile == "" 
            || numSets == 0 || objDensities == "" || objTypes == "" || queryNodeFile == "" 
            || kValues == "" || method == "") {
            std::cerr << "Invalid argument(s)!\n\n";
            this->showPhaseUsage(experiment,argv[0]);
            exit(1);
        }
        
        this->runSingleMethodQueries(method,bgrFilePath,queryNodeFile,kValues,parameters,numSets,objDensities,objTypes,objVariable,filePathPrefix,statsOutputFile);

    } else if (experiment == constants::EXP_RUN_KNN_RW_POI) {
        if (argc < 19) {
            // Arguments: -g <binary graph file> -q <query node file> -k <k values> -p <parameters> 
            // -f <index output file path prefix> -s <stats output file> -r <file listing real-world POI set files name>
            std::cerr << "Too few arguments!\n\n";
            this->showPhaseUsage(experiment,argv[0]);
            exit(1);
        }
    
        if (bgrFilePath == "" || parameters == "" || filePathPrefix == "" || statsOutputFile == "" 
            || rwPOISetListFile == "" || queryNodeFile == "" || kValues == "") {
            std::cerr << "Invalid argument(s)!\n\n";
            this->showPhaseUsage(experiment,argv[0]);
            exit(1);
        }
        
        this->runRealWorldPOIQueries(bgrFilePath,queryNodeFile,kValues,parameters,filePathPrefix,statsOutputFile,rwPOISetListFile);

    } else if (experiment == constants::EXP_BUILD_RW_POI_OBJ_INDEXES) {
        if (argc < 15) {
            // Arguments: -g <binary graph file>  -p <parameters> -f <index output file path prefix> 
            // -s <stats output file> -r <file listing real-world POI set files name>
            std::cerr << "Too few arguments!\n\n";
            this->showPhaseUsage(experiment,argv[0]);
            exit(1);
        }
    
        if (bgrFilePath == "" || parameters == "" || filePathPrefix == "" || statsOutputFile == "" 
            || rwPOISetListFile == "") {
            std::cerr << "Invalid argument(s)!\n\n";
            this->showPhaseUsage(experiment,argv[0]);
            exit(1);
        }
        
        this->buildRealWorldObjIndexes(bgrFilePath,parameters,filePathPrefix,statsOutputFile,rwPOISetListFile);
        
    } else {
        std::cerr << "Invalid experimental step!\n\n";
        this->showCommandUsage(argv[0]);
        exit(1);
    }
    
}

void ExperimentsCommand::showCommandUsage(std::string programName)
{
    std::cerr << "Usage: " << programName << " -c " + constants::EXPERIMENTS_CMD + " -e <experimental step>\n\n"
              << "Steps:\n"
              << utility::getFormattedUsageString(constants::EXP_BUILD_INDEXES,"1. Build all indexes") + "\n"
              << utility::getFormattedUsageString(constants::EXP_BUILD_OBJ_INDEXES,"2. Generate object sets and build object indexes") + "\n"
              << utility::getFormattedUsageString(constants::EXP_RUN_KNN,"3. Run all standard kNN query experiments") + "\n"
              << utility::getFormattedUsageString(constants::EXP_RUN_KNN_OPTIMIZATIONS,"4. Run individual method kNN query experiments") + "\n"
              << utility::getFormattedUsageString(constants::EXP_RUN_KNN_RW_POI,"5. Build object indexes for real-world POIs") + "\n"
              << utility::getFormattedUsageString(constants::EXP_RUN_KNN_RW_POI,"6. Run kNN querying on real-world POIs experiments") + "\n";
}

void ExperimentsCommand::showPhaseUsage(std::string method, std::string programName)
{
    if (method == constants::EXP_BUILD_INDEXES ) {
        std::cerr << "Usage: " << programName << " -c " + constants::EXPERIMENTS_CMD
                  << " -e " + constants::EXP_BUILD_INDEXES + " -g <binary graph file>\n"
                  << "-p <parameter key-value pairs> -f <index output file path prefix> -s <stats output file>\n";
    } else if (method == constants::EXP_BUILD_OBJ_INDEXES ) {
        std::cerr << "Usage: " << programName << " -c " + constants::EXPERIMENTS_CMD
                  << " -e " + constants::EXP_BUILD_OBJ_INDEXES + " -g <binary graph file>\n"
                  << "-p <parameter key-value pairs> -n <num sets> -d <list of object set densities or partitions>\n"
                  << "-t <list of object types>  -v <list of some object variable> -f <index output file path prefix> -s <stats output file>\n";
    } else if (method == constants::EXP_RUN_KNN ) {
        std::cerr << "Usage: " << programName << " -c " + constants::EXPERIMENTS_CMD
                  << " -e " + constants::EXP_RUN_KNN + " -g <binary graph file>\n"
                  << "-q <query node file> -k <k values> -p <parameters> -n <num sets> -d <list of object set densities or partitions>\n"
                  << "-t <list of object types>  -v <list of some object variable> -f <index output file path prefix> -s <stats output file>\n";
    } else if (method == constants::EXP_RUN_KNN_OPTIMIZATIONS ) {
        std::cerr << "Usage: " << programName << " -c " + constants::EXPERIMENTS_CMD
                  << " -e " + constants::EXP_RUN_KNN_OPTIMIZATIONS + " -m <method>\n"
                  << "-g <binary graph file> -q <query node file> -k <k values> -p <parameters>\n"
                  << "-n <num sets> -d <list of object densities> -t <list of object types>\n"
                  << " -v <list of some object variable> -f <index output file path prefix> -s <stats output file>\n";
    } else if (method == constants::EXP_RUN_KNN_RW_POI ) {
        std::cerr << "Usage: " << programName << " -c " + constants::EXPERIMENTS_CMD
                  << " -e " + constants::EXP_RUN_KNN_RW_POI + " -g <binary graph file>\n"
                  << "-q <query node file> -k <k values> -p <parameters> -r <rw POI set list file>\n"
                  << "-f <index output file path prefix> -s <stats output file>\n";
    } else if (method == constants::EXP_BUILD_RW_POI_OBJ_INDEXES ) {
        std::cerr << "Usage: " << programName << " -c " + constants::EXPERIMENTS_CMD
                  << " -e " + constants::EXP_RUN_KNN_RW_POI + " -g <binary graph file>\n"
                  << "-p <parameters> -r <rw POI set list file>\n"
                  << "-f <index output file path prefix> -s <stats output file>\n";
    } else {
        std::cerr << "Invalid experiment phase!" << std::endl;
        this->showCommandUsage(programName);
    }
}

std::unordered_map<std::string, std::string> ExperimentsCommand::getParameters(std::string parameters)
{
    std::unordered_map<std::string,std::string> parameterMap;
    std::vector<std::string> pairs = utility::splitByDelim(parameters,',');
    
    for (std::size_t i = 0; i < pairs.size(); ++i) {
        std::vector<std::string> pair = utility::splitByDelim(pairs[i],'=');
        if (pair.size() == 2) {
//             std::cout << "Key = " << pair[0] << std::endl;
//             std::cout << "Value = " << pair[1] << std::endl;
            parameterMap[pair[0]] = pair[1];
        } else {
            std::cerr << "Invalid key-value pair in parameter string" << std::endl;
            exit(1);
        }
    }
    return parameterMap;
}


void ExperimentsCommand::buildIndexes(std::string bgrFileName, std::string parameters, std::string filePathPrefix, std::string statsOutputFile)
{
    std::unordered_map<std::string,std::string> parameterMap = this->getParameters(parameters);
    
    bool buildGtreeIdx = parameterMap["gtree"] == "1";
    bool buildRoadIdx = parameterMap["road"] == "1";
    bool buildSILCIdx = parameterMap["silc"] == "1";
    bool buildALTIdx = parameterMap["alt"] == "1";
    
    // Find all additional fields we need to add to index stats tuples
    std::vector<std::string> specialFields;
    int field = 0;
    while (true) {
        std::string key = "special_field_" + std::to_string(field);
        if (parameterMap.find(key) != parameterMap.end()) {
            specialFields.push_back(parameterMap[key]);
            ++field;
        } else {
            break;
        }
    }    
    
    std::string idxOutputFile = "";
    std::vector<std::string> parameterKeys, parameterValues;

    Graph graph = serialization::getIndexFromBinaryFile<Graph>(bgrFileName);
    
    // Put all of these into functions (except input validation) so indexes get released after serialization
    
    if (buildGtreeIdx) {
        int fanout = std::stoi(parameterMap["gtree_fanout"]);
        std::size_t maxLeafSize = std::stoi(parameterMap["gtree_maxleafsize"]);
        if (fanout < 2 || maxLeafSize < 32) {
            std::cerr << "Invalid Gtree parameters!\n";
            exit(1);    
        }
        parameterKeys.clear();
        parameterValues.clear();
        parameterKeys.push_back("fanout");
        parameterValues.push_back(parameterMap["gtree_fanout"]);
        parameterKeys.push_back("maxleafsize");
        parameterValues.push_back(parameterMap["gtree_maxleafsize"]);
        idxOutputFile = filePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_GTREE_CMD,parameterKeys,parameterValues);
        this->buildGtree(graph,fanout,maxLeafSize,idxOutputFile,statsOutputFile,specialFields);
    }
    
    if (buildRoadIdx) {
        int fanout = std::stoi(parameterMap["road_fanout"]);
        std::size_t levels = std::stoi(parameterMap["road_levels"]);
        if (fanout < 2 || levels < 2) {
            std::cerr << "Invalid Route Overlay parameters!\n";
            exit(1);    
        }
        parameterKeys.clear();
        parameterValues.clear();
        parameterKeys.push_back("fanout");
        parameterValues.push_back(parameterMap["road_fanout"]);
        parameterKeys.push_back("levels");
        parameterValues.push_back(parameterMap["road_levels"]);
        idxOutputFile = filePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_ROUTEOVERLAY_CMD,parameterKeys,parameterValues);
        this->buildRouteOverlay(graph,fanout,levels,idxOutputFile,statsOutputFile,specialFields);
    }
    
    if (buildSILCIdx) {
        if (graph.getEdgeType() != constants::TIME_WEIGHTS) {
            int maxQuadtreeLeafSize = std::stoi(parameterMap["silc_maxquadtreeleafsize"]);
            if (maxQuadtreeLeafSize < 1) {
                std::cerr << "Invalid SILC parameters!\n";
                exit(1);    
            }
            // Note: madQuadtreeLeafSize doesn't affect structure of any Morton list, only the intermediary
            // quadtree used to build Morton lists that we discard anyway so we don't include in file name
            parameterKeys.clear();
            parameterValues.clear();
            idxOutputFile = filePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_SILC_CMD,parameterKeys,parameterValues);
            this->buildSILC(graph,maxQuadtreeLeafSize,idxOutputFile,statsOutputFile,specialFields);
        } else {
            std::cout << "SILC does not support kNN queries on travel time road network graphs!" << std::endl;
        }
    }

    if (buildALTIdx) {
        int numLandmarks = std::stoi(parameterMap["alt_numlandmarks"]);
        std::string selectionMethod = parameterMap["alt_landmarktype"];
        if (numLandmarks < 1 || selectionMethod == "") {
            std::cerr << "Invalid ALT parameters!\n";
            exit(1);    
        }
        parameterKeys.clear();
        parameterValues.clear();
        parameterKeys.push_back("landmarks");
        parameterValues.push_back(parameterMap["alt_numlandmarks"]);
        parameterKeys.push_back("type");
        parameterValues.push_back(parameterMap["alt_landmarktype"]);
        idxOutputFile = filePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_ALT_CMD,parameterKeys,parameterValues);
        this->buildALT(graph,numLandmarks,selectionMethod,idxOutputFile,statsOutputFile,specialFields);
    }
    
    std::cout << "Index building complete!" << std::endl;
}

void ExperimentsCommand::buildExternalIndexes(std::string bgrFileName, std::string parameters, std::string filePathPrefix, std::string statsOutputFile)
{
    std::unordered_map<std::string,std::string> parameterMap = this->getParameters(parameters);
    
    bool buildPHLIdx = parameterMap["phl"] == "1";
    bool buildCHIdx = parameterMap["ch"] == "1";
    bool buildTNRIdx = parameterMap["tnr"] == "1";
    
    // Find all additional fields we need to add to index stats tuples
    std::vector<std::string> specialFields;
    int field = 0;
    while (true) {
        std::string key = "special_field_" + std::to_string(field);
        if (parameterMap.find(key) != parameterMap.end()) {
            specialFields.push_back(parameterMap[key]);
            ++field;
        } else {
            break;
        }
    }    
    
    // Put all of these into functions (except input validation) so indexes get released after serialization
    std::string idxOutputFile = "";
    std::vector<std::string> parameterKeys, parameterValues;

    // Some additional file paths to create intermediary and additional files for some indexes
    std::string graphDataOutputFile = "", coordDataOutputFile = "", bgrIntFile = "", bcoIntFile = "", chFile = "";
    
    std::string networkName = "";
    int numNodes = 0, numEdges = 0;
    if (buildCHIdx || buildPHLIdx || buildTNRIdx) {
        // Use scope to destroy graph data structure to have extra memory for methods that don't need it
        Graph tempGraph = serialization::getIndexFromBinaryFile<Graph>(bgrFileName);
        networkName = tempGraph.getNetworkName();
        numNodes = tempGraph.getNumNodes();
        numEdges = tempGraph.getNumEdges();
        if (buildCHIdx) {
            // Create appropriate text graph files for input
            parameterKeys.clear();
            parameterValues.clear();
            graphDataOutputFile = filePathPrefix + "/data/" + utility::constructIndexFileName(networkName,"ddsg",parameterKeys,parameterValues);
            coordDataOutputFile = filePathPrefix + "/data/" + utility::constructIndexFileName(networkName,"zero_coords",parameterKeys,parameterValues);
            tempGraph.outputToDDSGFile(graphDataOutputFile);
            tempGraph.outputZeroedCoordinatesToFile(coordDataOutputFile);
        }
        if (buildPHLIdx) {
            parameterKeys.clear();
            parameterValues.clear();
            graphDataOutputFile = filePathPrefix + "/data/" + utility::constructIndexFileName(networkName,"tsv",parameterKeys,parameterValues);
            tempGraph.outputToTSVFile(graphDataOutputFile);
        }
    }
    
    if (buildCHIdx) {
        parameterKeys.clear();
        parameterValues.clear();
        graphDataOutputFile = filePathPrefix + "/data/" + utility::constructIndexFileName(networkName,"ddsg",parameterKeys,parameterValues);
        coordDataOutputFile = filePathPrefix + "/data/" + utility::constructIndexFileName(networkName,"zero_coords",parameterKeys,parameterValues);
        bgrIntFile = filePathPrefix + "/indexes/" + utility::constructIndexFileName(networkName,"ext_bin",parameterKeys,parameterValues);
        bcoIntFile = filePathPrefix + "/indexes/" + utility::constructIndexFileName(networkName,"ext_co_bin",parameterKeys,parameterValues);
        chFile = filePathPrefix + "/indexes/" + utility::constructIndexFileName(networkName,"hcn",parameterKeys,parameterValues);
        idxOutputFile = filePathPrefix + "/indexes/" + utility::constructIndexFileName(networkName,constants::IDX_CH_CMD,parameterKeys,parameterValues);
        this->buildCH(networkName,numNodes,numEdges,graphDataOutputFile,coordDataOutputFile,bgrIntFile,bcoIntFile,idxOutputFile,chFile,statsOutputFile,specialFields);
        
        // Delete intermediary text files as they are not needed
        std::remove(graphDataOutputFile.c_str());
        std::remove(coordDataOutputFile.c_str());
    }
    
    if (buildTNRIdx) {
        int gridSize = std::stoi(parameterMap["tnr_gridsize"]);
        if (gridSize < 16) {
            std::cerr << "Invalid TNR grid size!\n";
            exit(1);
        }
        parameterKeys.clear();
        parameterValues.clear();
        bgrIntFile = filePathPrefix + "/indexes/" + utility::constructIndexFileName(networkName,"ext_bin",parameterKeys,parameterValues);
        bcoIntFile = filePathPrefix + "/indexes/" + utility::constructIndexFileName(networkName,"ext_co_bin",parameterKeys,parameterValues);
        chFile = filePathPrefix + "/indexes/" + utility::constructIndexFileName(networkName,constants::IDX_CH_CMD,parameterKeys,parameterValues);
        parameterKeys.push_back("gridsize");
        parameterValues.push_back(parameterMap["tnr_gridsize"]);
        idxOutputFile = filePathPrefix + "/indexes/" + utility::constructIndexFileName(networkName,constants::IDX_TNR_CMD,parameterKeys,parameterValues);
        this->buildTNR(networkName,numNodes,numEdges,bgrIntFile,bcoIntFile,chFile,idxOutputFile,gridSize,statsOutputFile,specialFields);
    }

    if (buildPHLIdx) {
        parameterKeys.clear();
        parameterValues.clear();
        graphDataOutputFile = filePathPrefix + "/data/" + utility::constructIndexFileName(networkName,"tsv",parameterKeys,parameterValues);
        idxOutputFile = filePathPrefix + "/indexes/" + utility::constructIndexFileName(networkName,constants::IDX_PHL_CMD,parameterKeys,parameterValues);
        this->buildPHL(networkName,numNodes,numEdges,graphDataOutputFile,idxOutputFile,statsOutputFile,specialFields);

        // Delete intermediary text file as it is not needed
        std::remove(graphDataOutputFile.c_str());
    }

    std::cout << "External index building complete!" << std::endl;
}

void ExperimentsCommand::runQueries(std::string bgrFileName, std::string queryNodeFile, std::string kValues, 
                                    std::string parameters, std::size_t numSets, std::string objDensities, 
                                    std::string objTypes, std::string objVariable, std::string filePathPrefix, std::string statsOutputFile)
{
    Graph graph = serialization::getIndexFromBinaryFile<Graph>(bgrFileName);
    std::unordered_map<std::string,std::string> parameterMap = this->getParameters(parameters);

    std::vector<NodeID> queryNodes = utility::getPointSetFromFile(queryNodeFile);
    if (queryNodes.size() == 0) {
        std::cerr << "No query points were provided!\n";
        exit(1);    
    }
    if (numSets == 0) {
        std::cerr << "No object sets were provided!\n";
        exit(1);    
    }

    bool verifykNN = parameterMap["verify"] == "1";
    bool queryINE = parameterMap["ine"] == "1";
    bool queryGtree = parameterMap["gtree"] == "1";
    bool queryRoad = parameterMap["road"] == "1";
    bool queryIER = parameterMap["ier"] == "1";
    bool querySILC = parameterMap["silc"] == "1";
    bool queryDistBrws = parameterMap["dist_brws"] == "1";
    bool queryIERPHL = parameterMap["ier_phl"] == "1";
    
    std::vector<std::string> strObjDensitiesVec = utility::splitByDelim(objDensities,',');
    std::vector<double> objDensitiesVec;
    for(std::size_t i = 0; i < strObjDensitiesVec.size(); ++i) {
        double density = std::stod(strObjDensitiesVec[i]);
        if (density > 0) {
            objDensitiesVec.push_back(density);
        } else {
            std::cerr << "Invalid density in list provided!\n";
            exit(1);    
        }
    }
    std::vector<std::string> objTypesVec = utility::splitByDelim(objTypes,',');
    std::vector<std::string> strKValuesVec = utility::splitByDelim(kValues,',');
    std::vector<int> kValuesVec;
    for(std::size_t i = 0; i < strKValuesVec.size(); ++i) {
        int k = std::stoi(strKValuesVec[i]);
        if (k > 0) {
            kValuesVec.push_back(k);
        } else {
            std::cerr << "Invalid k value in list provided!\n";
            exit(1);    
        }
    }
    std::vector<std::string> strObjVariables = utility::splitByDelim(objVariable,',');
    std::vector<int> objVariableVec;
    for(std::size_t i = 0; i < strObjVariables.size(); ++i) {
        int variable = std::stoi(strObjVariables[i]);
        if (variable > 0) {
            objVariableVec.push_back(variable);
        } else {
            std::cerr << "Invalid variable in list provided (must be greater than zero)!\n";
            exit(1);    
        }
    }
    if (objDensitiesVec.size() == 0 || objTypesVec.size() == 0 || objVariableVec.size() == 0 || kValuesVec.size() == 0) {
        std::cerr << "Not enough densities or types provided!\n";
        exit(1);    
    }
    
    if (queryINE) {
        this->runINEQueries(graph,queryNodes,kValuesVec,numSets,objDensitiesVec,objTypesVec,objVariableVec,filePathPrefix,statsOutputFile);
    }
    
    std::string gtreeIdxFile, roadIdxFile, silcIdxFile, juncIdxFile, phlIdxFile;
    std::vector<std::string> parameterKeys, parameterValues;
    
    if (queryRoad) {
        int fanout = std::stoi(parameterMap["road_fanout"]);
        std::size_t levels = std::stoi(parameterMap["road_levels"]);
        if (fanout < 2 || levels < 2) {
            std::cerr << "Invalid Route Overlay parameters!\n";
            exit(1);    
        }
        parameterKeys.clear();
        parameterValues.clear();
        parameterKeys.push_back("fanout");
        parameterValues.push_back(parameterMap["road_fanout"]);
        parameterKeys.push_back("levels");
        parameterValues.push_back(parameterMap["road_levels"]);
        roadIdxFile = filePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_ROUTEOVERLAY_CMD,parameterKeys,parameterValues);        
        this->runROADQueries(graph,roadIdxFile,queryNodes,kValuesVec,numSets,objDensitiesVec,objTypesVec,objVariableVec,filePathPrefix,statsOutputFile,verifykNN,parameterKeys,parameterValues);
    }
    
    if (queryGtree) {
        int fanout = std::stoi(parameterMap["gtree_fanout"]);
        std::size_t maxLeafSize = std::stoi(parameterMap["gtree_maxleafsize"]);
        if (fanout < 2 || maxLeafSize < 32) {
            std::cerr << "Invalid Gtree parameters!\n";
            exit(1);    
        }
        parameterKeys.clear();
        parameterValues.clear();
        parameterKeys.push_back("fanout");
        parameterValues.push_back(parameterMap["gtree_fanout"]);
        parameterKeys.push_back("maxleafsize");
        parameterValues.push_back(parameterMap["gtree_maxleafsize"]);
        gtreeIdxFile = filePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_GTREE_CMD,parameterKeys,parameterValues);
        this->runGtreeQueries(graph,gtreeIdxFile,queryNodes,kValuesVec,numSets,objDensitiesVec,objTypesVec,objVariableVec,filePathPrefix,statsOutputFile,verifykNN,parameterKeys,parameterValues);
    }
    
    if (queryIER) {
        std::vector<int> branchFactors = utility::getIngetersFromStringList(parameterMap["rtree_branchfactor"],':',0);
        if (branchFactors.size() == 0) {
            std::cerr << "No R-tree branch factors provided!\n";
            exit(1);
        }
        parameterKeys.clear();
        parameterValues.clear();
        parameterKeys.push_back("fanout");
        parameterValues.push_back(parameterMap["gtree_fanout"]);
        parameterKeys.push_back("maxleafsize");
        parameterValues.push_back(parameterMap["gtree_maxleafsize"]);
        gtreeIdxFile = filePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_GTREE_CMD,parameterKeys,parameterValues);
        this->runIERQueries(graph,constants::GTREE_SPDIST_QUERY,gtreeIdxFile,branchFactors,queryNodes,kValuesVec,numSets,objDensitiesVec,objTypesVec,objVariableVec,filePathPrefix,statsOutputFile,verifykNN);
    }
    
    if (queryIERPHL) {
        std::vector<int> branchFactors = utility::getIngetersFromStringList(parameterMap["rtree_branchfactor"],':',0);
        if (branchFactors.size() == 0) {
            std::cerr << "No R-tree branch factors provided!\n";
            exit(1);
        }
        parameterKeys.clear();
        parameterValues.clear();
        phlIdxFile = filePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_PHL_CMD,parameterKeys,parameterValues);
        if (utility::fileExists(phlIdxFile)) {
            // It's possible PHL file does not exist for larger road networks
            // so we need to check in this case (but junc index file will exist)
            this->runIERQueries(graph,constants::PHL_SPDIST_QUERY,phlIdxFile,branchFactors,queryNodes,kValuesVec,numSets,objDensitiesVec,objTypesVec,objVariableVec,filePathPrefix,statsOutputFile,verifykNN);
        }
    }
    
    if (queryDistBrws) {
        
        std::vector<int> branchFactors = utility::getIngetersFromStringList(parameterMap["rtree_branchfactor"],':',0);
        if (branchFactors.size() != 0) {
            std::vector<std::string> methods;
            methods.push_back(constants::DB_RTREE_KNN_QUERY);
            parameterKeys.clear();
            parameterValues.clear();
            silcIdxFile = filePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_SILC_CMD,parameterKeys,parameterValues);
            juncIdxFile = filePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_JUNC_CMD,parameterKeys,parameterValues);
            if (utility::fileExists(silcIdxFile) && utility::fileExists(juncIdxFile)) {
                // It's possible SILC file does not exist for larger road networks
                // so we need to check in this case (but junc index file will exist)
                this->runSILCBasedQueries(graph,methods,silcIdxFile,juncIdxFile,branchFactors,queryNodes,kValuesVec,numSets,objDensitiesVec,objTypesVec,objVariableVec,filePathPrefix,statsOutputFile,verifykNN);        
            }
        }
//         // Original DistBrws 
//         std::vector<int> maxLeafSizes = utility::getIngetersFromStringList(parameterMap["db_quadtree_maxleafsize"],':',0);
//         if (maxLeafSizes.size() != 0) {
//             std::vector<std::string> methods;
//             methods.push_back(constants::OPT_DISTBRWS_KNN_QUERY);
//             parameterKeys.clear();
//             parameterValues.clear();
//             silcIdxFile = filePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_SILC_CMD,parameterKeys,parameterValues);
//             juncIdxFile = filePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_JUNC_CMD,parameterKeys,parameterValues);
//             this->runSILCBasedQueries(graph,methods,silcIdxFile,juncIdxFile,maxLeafSizes,queryNodes,kValuesVec,numSets,objDensitiesVec,objTypesVec,objVariableVec,filePathPrefix,statsOutputFile,verifykNN);        
//         } else {
//             // If they are not provided then this we are using a dataset without a SILC index
//             //std::cout << "No object hierarchy quadtree max leaf sizes provided!\n";
//         }
    }

    std::cout << "Query testing complete!" << std::endl;
}

void ExperimentsCommand::buildGtree(Graph& graph, int fanout, std::size_t maxLeafSize, std::string idxOutputFile, std::string statsOutputFile, std::vector<std::string> specialFields)
{
    StopWatch sw;
    
    sw.start();
    Gtree gtree(graph.getNetworkName(),graph.getNumNodes(),graph.getNumEdges(),fanout,maxLeafSize);
    gtree.buildGtree(graph);
    sw.stop();

    double processingTimeMs = sw.getTimeMs();
    // Note: G-tree uses Graph to search for object within leaf (only edges and edge weights same parts as INE)
    double memoryUsage = gtree.computeIndexSize() + graph.computeINEIndexSize();
    
    IndexTuple stats(gtree.getNetworkName(),gtree.getNumNodes(),gtree.getNumEdges(),
                     constants::IDX_GTREE_CMD,processingTimeMs,memoryUsage);
    stats.addSupplementaryFields("max_leaf_size",std::to_string(maxLeafSize));
    stats.addSupplementaryFields("fanout",std::to_string(fanout));
    stats.addSupplementaryFields("levels",std::to_string(gtree.getNumLevels()));
    int numTreeNodes = gtree.getTreeSize(), numBorders = gtree.getNumBorders(), b2bRelationships = gtree.getBorderToBorderRelationships(), avgPathCost = gtree.getAvgPathCost();
    stats.addSupplementaryFields("num_tree_nodes",std::to_string(numTreeNodes));
    stats.addSupplementaryFields("num_borders",std::to_string(numBorders));
    stats.addSupplementaryFields("num_b2b_rel",std::to_string(b2bRelationships));
    stats.addSupplementaryFields("avg_path_cost",std::to_string(avgPathCost));
    stats.setAdditionalFields(specialFields);

//     std::cout << stats.getMultilineTupleString();

    serialization::outputIndexToBinaryFile<Gtree>(gtree,idxOutputFile);

    this->outputCommandStats(statsOutputFile,stats.getTupleString());
    
    std::cout << "Gtree index successfully created!" << std::endl;
}

void ExperimentsCommand::buildRouteOverlay(Graph& graph, int fanout, int levels, std::string idxOutputFile, std::string statsOutputFile, std::vector<std::string> specialFields)
{
    StopWatch sw;
    sw.start();
    ROAD road(graph.getNetworkName(),graph.getNumNodes(),graph.getNumEdges(),fanout,levels);
    road.buildRouteOverlay(graph);
    sw.stop();
    
    double processingTimeMs = sw.getTimeMs();   
    double memoryUsage = road.computeIndexSize();
    
    IndexTuple stats(road.getNetworkName(),road.getNumNodes(),road.getNumEdges(),
                     constants::IDX_ROUTEOVERLAY_CMD,processingTimeMs,memoryUsage);
    stats.addSupplementaryFields("fanout",std::to_string(fanout));
    stats.addSupplementaryFields("levels",std::to_string(levels));
    stats.addSupplementaryFields("real_levels",std::to_string(road.getRealNumLevels()));
    int numTreeNodes = road.getRnetTreeSize(), numBorders = road.getNumBorders(), b2bRelationships = road.getBorderToBorderRelationships();
    stats.addSupplementaryFields("num_tree_nodes",std::to_string(numTreeNodes));
    stats.addSupplementaryFields("num_borders",std::to_string(numBorders));
    stats.addSupplementaryFields("num_b2b_rel",std::to_string(b2bRelationships));
    stats.setAdditionalFields(specialFields);

//     std::cout << stats.getMultilineTupleString();

    serialization::outputIndexToBinaryFile<ROAD>(road,idxOutputFile);

    this->outputCommandStats(statsOutputFile,stats.getTupleString());
    
    std::cout << "Route Overlay index successfully created!" << std::endl;
}

void ExperimentsCommand::buildSILC(Graph& graph, int maxQuadtreeLeafSize, std::string idxOutputFile, std::string statsOutputFile, std::vector<std::string> specialFields)
{
    StopWatch sw;
    
    sw.start();
    Quadtree quadtree(maxQuadtreeLeafSize);
    quadtree.buildFromGraph(graph);
    sw.stop();
    double quadtreeConstructionTime = sw.getTimeMs();    

    sw.reset();
    sw.start();
    SILCPathOracle silc;
    silc.buildPathOracle(quadtree,graph);
    sw.stop();

    double processingTimeMs = quadtreeConstructionTime + sw.getTimeMs();
    // Note: SILC uses Graph to get adj node and weight associated with colour (only edges and edge weights same parts as INE)
    double memoryUsage = silc.computeIndexSize() + graph.computeINEIndexSize();

    IndexTuple stats(silc.getNetworkName(),silc.getNumNodes(),silc.getNumEdges(),constants::IDX_SILC_CMD,processingTimeMs,memoryUsage);
    stats.setAdditionalFields(specialFields);

//     std::cout << stats.getMultilineTupleString();

    serialization::outputIndexToBinaryFile<SILCPathOracle>(silc,idxOutputFile);
    
    this->outputCommandStats(statsOutputFile,stats.getTupleString());

    std::cout << "SILC index successfully created!" << std::endl;
}

void ExperimentsCommand::buildPHL(std::string networkName, int numNodes, int numEdges, std::string dataOutputFile, std::string idxOutputFile, std::string statsOutputFile, std::vector<std::string> specialFields)
{
    PrunedHighwayLabeling phl;
    phl.ConstructLabel(dataOutputFile.c_str());
    
    double processingTimeMs = phl.getConstructionTime()*1000; // convert to ms
    double memoryUsage = phl.computeIndexSize();

    IndexTuple stats(networkName,numNodes,numEdges,constants::IDX_PHL_CMD,processingTimeMs,memoryUsage);
    stats.setAdditionalFields(specialFields);
    this->outputCommandStats(statsOutputFile,stats.getTupleString());

    phl.StoreLabel(idxOutputFile.c_str());
    
    std::cout << "PHL index successfully created!" << std::endl;
}

void ExperimentsCommand::buildCH(std::string networkName, int numNodes, int numEdges, std::string graphFilePath, std::string coordinateFilePath, std::string bgrIntFilePath, 
                                 std::string bcoIntFilePath, std::string chFilePath, std::string noFilePath, std::string statsOutputFile, 
                                 std::vector<std::string> specialFields)
{
    ShortestPathWrapper spw;
    spw.buildIntermediaryGraph(graphFilePath,coordinateFilePath,bgrIntFilePath,bcoIntFilePath);
    std::cout << "Intermediary graph indexes successfully created!" << std::endl;
    
    spw.buildCH(bgrIntFilePath,chFilePath,noFilePath);
    
    double processingTimeMs = spw.getConstructionTimeMs();
    double memoryUsage = spw.getIndexSizeMB();

    IndexTuple stats(networkName,numNodes,numEdges,constants::IDX_CH_CMD,processingTimeMs,memoryUsage);
    stats.setAdditionalFields(specialFields);
    this->outputCommandStats(statsOutputFile,stats.getTupleString());

    std::cout << "CH index successfully created!" << std::endl;
}

void ExperimentsCommand::buildTNR(std::string networkName, int numNodes, int numEdges, std::string bgrIntFilePath, std::string bcoIntFilePath, std::string chFilePath, 
                                  std::string tnrFilePath, int gridSize, std::string statsOutputFile, std::vector<std::string> specialFields)
{
    ShortestPathWrapper spw;
    spw.buildTNR(bgrIntFilePath,bcoIntFilePath,chFilePath,tnrFilePath,gridSize);

    double processingTimeMs = spw.getConstructionTimeMs();
    double memoryUsage = spw.getIndexSizeMB();
    
    IndexTuple stats(networkName,numNodes,numEdges,constants::IDX_TNR_CMD,processingTimeMs,memoryUsage);
    stats.setAdditionalFields(specialFields);
    stats.addSupplementaryFields("gridsize",std::to_string(gridSize));
    this->outputCommandStats(statsOutputFile,stats.getTupleString());

    std::cout << "TNR index successfully created!" << std::endl;
}

void ExperimentsCommand::buildALT(Graph& graph, int numLandmarks, std::string selectionMethod, std::string idxOutputFile, 
                                  std::string statsOutputFile, std::vector<std::string> specialFields)
{
    StopWatch sw;
    bool validMethod;
    
    sw.start();
    ALT alt(graph.getNetworkName(),graph.getNumNodes(),graph.getNumEdges());
    LANDMARK_TYPE landmarkType = alt.getLandmarkType(selectionMethod,validMethod);
    if (!validMethod) {
        std::cerr << "Invalid landmark selection method!" << std::endl;
        std::exit(1);
    }
    alt.buildALT(graph,landmarkType,numLandmarks);
    sw.stop();
    
    double processingTimeMs = sw.getTimeMs();    
    double memoryUsage = alt.computeIndexSize();

    IndexTuple stats(alt.getNetworkName(),alt.getNumNodes(),alt.getNumEdges(),constants::IDX_ALT_CMD,processingTimeMs,memoryUsage);
    stats.setAdditionalFields(specialFields);
    stats.addSupplementaryFields("landmarks",std::to_string(numLandmarks));
    stats.addSupplementaryFields("type",selectionMethod);

//     std::cout << stats.getMultilineTupleString();

    serialization::outputIndexToBinaryFile<ALT>(alt,idxOutputFile);
    
    this->outputCommandStats(statsOutputFile,stats.getTupleString());

    std::cout << "ALT index successfully created!" << std::endl;
}

void ExperimentsCommand::buildObjectIndexes(std::string bgrFileName, std::string parameters, std::size_t numSets, std::string objDensities, 
                                            std::string objTypes, std::string objVariable, std::string filePathPrefix, std::string statsOutputFile, std::string queryNodeFile, unsigned int numQueryNodes)
{
    Graph graph = serialization::getIndexFromBinaryFile<Graph>(bgrFileName);
    std::unordered_map<std::string,std::string> parameterMap = this->getParameters(parameters);

    bool generateObjSets = parameterMap["obj_sets"] == "1";
    bool buildOccList = parameterMap["occ_list"] == "1";
    bool buildAssocDir = parameterMap["assoc_dir"] == "1";
    bool buildQuadtree = parameterMap["quadtree"] == "1";
    bool buildRtree = parameterMap["rtree"] == "1";

    std::vector<std::string> strDensities = utility::splitByDelim(objDensities,',');
    std::vector<double> densities;
    for(std::size_t i = 0; i < strDensities.size(); ++i) {
        double density = std::stod(strDensities[i]);
        if (density > 0) {
            densities.push_back(density);
        } else {
            std::cerr << "Invalid density in list provided!\n";
            exit(1);
        }
    }
    std::vector<std::string> strVariables = utility::splitByDelim(objVariable,',');
    std::vector<int> variables;
    for(std::size_t i = 0; i < strVariables.size(); ++i) {
        int variable = std::stoi(strVariables[i]);
        if (variable > 0) {
            variables.push_back(variable);
        } else {
            std::cerr << "Invalid variable in list provided (must be greater than zero)!\n";
            exit(1);    
        }
    }
    std::vector<std::string> types = utility::splitByDelim(objTypes,',');
    if (densities.size() == 0 || types.size() == 0 || variables.size() == 0) {
        std::cerr << "Not enough densities or types provided!\n";
        exit(1);
    }

    std::string gtreeIdxFile, roadIdxFile;
    std::vector<std::string> parameterKeys, parameterValues;
    
    if (generateObjSets) {
        // Generate and Output Object Sets
        this->generateObjectSets(graph,numSets,densities,types,variables,filePathPrefix,statsOutputFile,parameterKeys,parameterValues,queryNodeFile,numQueryNodes);
    }
    
    if (buildOccList) {
        int fanout = std::stoi(parameterMap["gtree_fanout"]);
        std::size_t maxLeafSize = std::stoi(parameterMap["gtree_maxleafsize"]);
        if (fanout < 2 || maxLeafSize < 32) {
            std::cerr << "Invalid Gtree parameters!\n";
            exit(1);
        }
        parameterKeys.clear();
        parameterValues.clear();
        parameterKeys.push_back("fanout");
        parameterValues.push_back(parameterMap["gtree_fanout"]);
        parameterKeys.push_back("maxleafsize");
        parameterValues.push_back(parameterMap["gtree_maxleafsize"]);
        gtreeIdxFile = filePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_GTREE_CMD,parameterKeys,parameterValues);
        this->buildOccurenceLists(gtreeIdxFile,numSets,densities,types,variables,filePathPrefix,statsOutputFile,parameterKeys,parameterValues);
    }
    
    if (buildAssocDir) {
        int fanout = std::stoi(parameterMap["road_fanout"]);
        std::size_t levels = std::stoi(parameterMap["road_levels"]);
        if (fanout < 2 || levels < 2) {
            std::cerr << "Invalid Route Overlay parameters!\n";
            exit(1);
        }
        parameterKeys.clear();
        parameterValues.clear();
        parameterKeys.push_back("fanout");
        parameterValues.push_back(parameterMap["road_fanout"]);
        parameterKeys.push_back("levels");
        parameterValues.push_back(parameterMap["road_levels"]);
        roadIdxFile = filePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_ROUTEOVERLAY_CMD,parameterKeys,parameterValues);        
        this->buildAssociationDirectories(roadIdxFile,numSets,densities,types,variables,filePathPrefix,statsOutputFile,parameterKeys,parameterValues);
    }
    
    if (buildRtree) {
        std::vector<int> branchFactors = utility::getIngetersFromStringList(parameterMap["rtree_branchfactor"],':',0);
        if (branchFactors.size() == 0) {
            std::cerr << "No R-tree branch factors provided!\n";
            exit(1);
        }
        this->buildRtrees(graph,branchFactors,numSets,densities,types,variables,filePathPrefix,statsOutputFile);        
    }

    if (buildQuadtree) {
        // Build Quadtrees for required by both SILC and DistBrws (but only build unique quadtrees)
        std::unordered_set<int> uniqueMaxLeafSizes; 
        std::vector<int> silcMaxLeafSizes = utility::getIngetersFromStringList(parameterMap["silc_quadtree_maxleafsize"],':',0);
        std::vector<int> dbMaxLeafSizes = utility::getIngetersFromStringList(parameterMap["db_quadtree_maxleafsize"],':',0);
        uniqueMaxLeafSizes.insert(silcMaxLeafSizes.begin(),silcMaxLeafSizes.end());
        uniqueMaxLeafSizes.insert(dbMaxLeafSizes.begin(),dbMaxLeafSizes.end());
        std::vector<int> maxLeafSizes(uniqueMaxLeafSizes.begin(),uniqueMaxLeafSizes.end());
        if (maxLeafSizes.size() == 0) {
            std::cerr << "No object hierarchy quadtree max leaf sizes provided!\n";
            exit(1);
        }
        this->buildQuadtrees(graph,maxLeafSizes,numSets,densities,types,variables,filePathPrefix,statsOutputFile);
    }
}

void ExperimentsCommand::generateObjectSets(Graph& graph, std::size_t numSets, std::vector<double> objDensities, 
                                            std::vector<std::string> objTypes, std::vector<int> objVariable, std::string filePathPrefix, std::string statsOutputFile, 
                                            std::vector<std::string>& parameterKeys, std::vector<std::string>& parameterValues, std::string queryNodeFile, unsigned int numQueryNodes)
{
    SetGenerator sg;
    int numNodes = graph.getNumNodes(), setSize;
    std::vector<NodeID> sampleSet;
    std::string objSetOutputFile = "";
    StopWatch sw;
    double totalINETime, totalINEMemory;
    int totalSets = numSets;
    
    // Generate and Output Object Sets
    for (std::size_t i = 0; i < objTypes.size(); ++i) {
        if (objTypes[i] == constants::RAND_OBJ_SET) {
            for (std::size_t j = 0; j < objDensities.size(); ++j) {
                for (std::size_t l = 0; l < objVariable.size(); ++l) {
                    setSize = std::ceil(numNodes*objDensities[j]);
                    totalINETime = 0;
                    totalINEMemory = 0;
                    for (std::size_t k = 0; k < numSets; ++k) {
                        if (objDensities[j] != 1) {
                            sampleSet = sg.generateRandomSampleSet(numNodes,setSize);
                        } else {
                            sampleSet = graph.getNodesIDsVector();
                        }

                        if (objDensities[j] != 1 || k == 0) {
                            // If density is one we only need to create one object set
                            objSetOutputFile = filePathPrefix + "/obj_indexes/" + utility::constructObjsectSetFileName(graph.getNetworkName(),objTypes[i],objDensities[j],objVariable[l],k);
                            utility::writeSampleSet(objSetOutputFile,graph.getNetworkName(),objTypes[i],objDensities[j],setSize,objVariable[l],sampleSet);
                        }
                        // We only include the size of the input to the other object indexes and the time taken
                        // to load load this information into the graph as the baselines for INE
                        graph.resetAllObjects();
                        sw.reset();
                        sw.start();
                        graph.parseObjectSet(sampleSet);
                        sw.stop();
                        totalINETime += sw.getTimeUs();
                        totalINEMemory += static_cast<double>(sizeof(NodeID)*sampleSet.size())/(1024*1024);
                        
                        if (objDensities[j] == 1) {
                            totalSets = 1;
                            // If density is 1 then no don't iterate more than once
                            break;
                        } else {
                            totalSets = numSets;
                        }
                    }
                    
                    // Average stats over all object sets and then write to stats file
                    double processingTimeMs = totalINETime/totalSets;
                    double memoryUsage = totalINEMemory/totalSets*1024; // in KB
                    
                    ObjectIndexTuple stats(graph.getNetworkName(),graph.getNumNodes(),graph.getNumEdges(),
                                            constants::INE_KNN_QUERY,processingTimeMs,memoryUsage,objTypes[i],objDensities[j],objVariable[l],setSize);

                    this->outputCommandStats(statsOutputFile,stats.getTupleString());
                }
            }
        } else if (objTypes[i] == constants::PARTITION_OBJ_SET) {
            // If all densities are greater than 1 then we will use densities
            // as the number of partitions, otherwise we calculate the number  
            // of partitions using the density and cluster size
            int numPartitions, clusterSize;
            bool fractionalDensityFound = false;
            for (std::size_t j = 0; j < objDensities.size(); ++j) {
                if (objDensities[j] < 1) {
                    fractionalDensityFound = true;
                    break;
                }
            }
            for (std::size_t j = 0; j < objDensities.size(); ++j) {
                for (std::size_t l = 0; l < objVariable.size(); ++l) {
                    clusterSize = objVariable[l];
                    if (fractionalDensityFound) {
                        setSize = objDensities[j]*graph.getNumNodes();
                        numPartitions = objDensities[j]*graph.getNumNodes()/clusterSize;
                    } else {
                        // Otherwise objDensities[j] represents the number cluster to create
                        setSize = objDensities[j]*clusterSize;
                        numPartitions = objDensities[j];
                    }
                    if (numPartitions < 1 || numPartitions > graph.getNumNodes()) {
                        std::cerr << "Number of partitions (" << numPartitions << ") must be greater than or equal to 1 and less than or equal to the total number of nodes" << std::endl;
                        std::exit(1);
                    }
                    totalINETime = 0;
                    totalINEMemory = 0;
                    for (std::size_t k = 0; k < numSets; ++k) {
                        sampleSet = sg.generateRandomClusteredSampleSet(graph,0.05,setSize,clusterSize);
                        setSize = sampleSet.size();

                        objSetOutputFile = filePathPrefix + "/obj_indexes/" + utility::constructObjsectSetFileName(graph.getNetworkName(),objTypes[i],objDensities[j],objVariable[l],k);
                        utility::writeSampleSet(objSetOutputFile,graph.getNetworkName(),objTypes[i],objDensities[j],setSize,objVariable[l],sampleSet);

                        // We only include the size of the input to the other object indexes as the baseline 
                        // as there's no realy baseline time to measure like in the input graph case
                        totalINETime += 0;
                        totalINEMemory += static_cast<double>(sizeof(NodeID)*sampleSet.size())/(1024*1024);
                    }
                    
                    // Average stats over 100 object sets and then write to stats file
                    double processingTimeMs = totalINETime/totalSets;
                    double memoryUsage = totalINEMemory/totalSets;
                    
                    ObjectIndexTuple stats(graph.getNetworkName(),graph.getNumNodes(),graph.getNumEdges(),
                                           objTypes[i],processingTimeMs,memoryUsage,objTypes[i],objDensities[j],objVariable[l],setSize);

                    this->outputCommandStats(statsOutputFile,stats.getTupleString());
                }
            }
        } else if (objTypes[i] == constants::MINND_OBJ_SET || objTypes[i] == constants::MINMAXND_OBJ_SET) {
            assert(queryNodeFile != "" && numQueryNodes != 0);
            int numBuckets = objVariable.size()+1; // Note: Input doesn't include "0th bucket" which contains query nodes
            SetGenerator sg;
            std::vector<std::vector<NodeID>> buckets = sg.getNodeBuckets(graph,numBuckets,"powersof2");

            // Query Nodes are selected from the 0th bucket (closest to centre)
            std::vector<NodeID> queryNodes;
            if (buckets[0].size() >= numQueryNodes) {
                std::vector<NodeID> bucketZeroNodes = sg.generateMinMaxNetworkDistSampleSet(buckets,0,0,numQueryNodes);
                queryNodes.insert(queryNodes.begin(),bucketZeroNodes.begin(),bucketZeroNodes.end());
            } else if (buckets[0].size() >= numQueryNodes/2) {
                queryNodes.insert(queryNodes.begin(),buckets[0].begin(),buckets[0].end());
                while (queryNodes.size() < numQueryNodes) {
                    for (auto it = buckets[0].begin(); it != buckets[0].end() && queryNodes.size() < numQueryNodes; ++it) {
                        queryNodes.push_back(*it);
                    }
                }
            } else {
                std::cerr << "Unabled to create query set of size " << numQueryNodes << std::endl;
                std::cout << "Nodes Available = " << buckets[0].size() << std::endl;
                std::cout << "Please use fewer buckets" <<  std::endl;
                std::exit(1);
            }
            utility::writeSampleSet(queryNodeFile,graph.getNetworkName(),objTypes[i],0,numQueryNodes,0,queryNodes);

            for (std::size_t j = 0; j < objDensities.size(); ++j) {
                setSize = std::ceil(numNodes*objDensities[j]);
                for (std::size_t l = 0; l < objVariable.size(); ++l) {
                    totalINETime = 0;
                    totalINEMemory = 0;
                    for (std::size_t k = 0; k < numSets; ++k) {
                        if (objTypes[i] == constants::MINND_OBJ_SET) {
                            sampleSet = sg.generateMinNetworkDistSampleSet(buckets,objVariable[l],setSize);
                        } else if (objTypes[i] == constants::MINMAXND_OBJ_SET) {
                             sampleSet = sg.generateMinMaxNetworkDistSampleSet(buckets,objVariable[l],objVariable[l],setSize); // Get objects only from current set
                        }

                        objSetOutputFile = filePathPrefix + "/obj_indexes/" + utility::constructObjsectSetFileName(graph.getNetworkName(),objTypes[i],objDensities[j],objVariable[l],k);
                        utility::writeSampleSet(objSetOutputFile,graph.getNetworkName(),objTypes[i],objDensities[j],setSize,objVariable[l],sampleSet);

                        // We only include the size of the input to the other object indexes as the baseline 
                        // as there's no realy baseline time to measure like in the input graph case
                        totalINETime += 0;
                        totalINEMemory += static_cast<double>(sizeof(NodeID)*sampleSet.size())/(1024*1024);
                    }
                    
                    // Average stats over 100 object sets and then write to stats file
                    double processingTimeMs = totalINETime/totalSets;
                    double memoryUsage = totalINEMemory/totalSets;
                    
                    ObjectIndexTuple stats(graph.getNetworkName(),graph.getNumNodes(),graph.getNumEdges(),
                                           objTypes[i],processingTimeMs,memoryUsage,objTypes[i],objDensities[j],objVariable[l],setSize);

                    this->outputCommandStats(statsOutputFile,stats.getTupleString());
                }
            }
        } else {
            std::cerr << "Invalid object set type provided!\n";
            exit(1);
        }
    }
    std::cout << "Object sets successfully generated for " << graph.getNetworkName() << std::endl;
}


void ExperimentsCommand::buildOccurenceLists(std::string gtreeIdxFile, std::size_t numSets, std::vector<double> objDensities, 
                                             std::vector<std::string> objTypes, std::vector<int> objVariable, std::string filePathPrefix, std::string statsOutputFile, 
                                             std::vector<std::string>& parameterKeys, std::vector<std::string>& parameterValues)
{
    Gtree gtree = serialization::getIndexFromBinaryFile<Gtree>(gtreeIdxFile);
    
    int setSize, setVariable;
    std::string setType;
    double setDensity;
    std::string objSetFile, objIdxOutputFile;
    std::vector<NodeID> objectNodes;
    
    StopWatch sw;
    double totalTime = 0, totalMemory = 0;
    int totalSets = numSets;
    for (std::size_t i = 0; i < objTypes.size(); ++i) {
        for (std::size_t j = 0; j < objDensities.size(); ++j) {
            for (std::size_t l = 0; l < objVariable.size(); ++l) {
                totalTime = 0;
                totalMemory = 0;
                
                for (std::size_t k = 0; k < numSets; ++k) {
                    objSetFile = filePathPrefix + "/obj_indexes/" + utility::constructObjsectSetFileName(gtree.getNetworkName(),objTypes[i],objDensities[j],objVariable[l],k);
                    std::vector<NodeID> objectNodes = utility::getPointSetFromFile(objSetFile,setType,setDensity,setSize,setVariable);
                    //assert(objTypes[i] == setType && objDensities[j] == setDensity);

                    sw.start();
                    OccurenceList occList(setType,setDensity,setVariable,setSize);
                    for (auto objIt = objectNodes.begin(); objIt != objectNodes.end(); ++objIt) {
                        // Find the Gtree leaf index for this object and add it to occurence list
                        // for that leaf (create list if it doesn't exist)
                        int leafIdx = gtree.getLeafIndex(*objIt);
                        occList.addLeafOccurence(leafIdx,*objIt);
                        
                        // Propagate this to occurence lists of parents of leaf node
                        int childIdx = leafIdx;
                        int parentIdx = gtree.getParentIndex(leafIdx);
                        while (parentIdx != -1) {
                            occList.addParentOccurence(parentIdx, childIdx);
                            
                            // Go up a level (until we reach root)
                            childIdx = parentIdx;
                            parentIdx = gtree.getParentIndex(childIdx);
                        }
                    }
                    sw.stop();
                    totalTime += sw.getTimeUs();
                    totalMemory += occList.computeIndexSize()*1024; // in KB

                    if (objDensities[j] != 1 || k == 0 || objTypes[i] != constants::RAND_OBJ_SET) {
                        // If density is one we only need to create one object index                    
                        objIdxOutputFile = filePathPrefix + "/obj_indexes/" + utility::constructObjectIndexFileName(gtree.getNetworkName(),constants::OBJ_IDX_GTREE,objTypes[i],objDensities[j],objVariable[l],k,parameterKeys,parameterValues);
                        serialization::outputIndexToBinaryFile<OccurenceList>(occList,objIdxOutputFile);
                    }
                    
                    if (objDensities[j] == 1 && objTypes[i] == constants::RAND_OBJ_SET) {
                        totalSets = 1;
                        // If density is 1 then no don't iterate more than once
                        break;
                    } else {
                        totalSets = numSets;
                    }
                    
                }
                
                // Average stats over 100 object sets and then write to stats file
                double processingTimeMs = totalTime/totalSets;
                double memoryUsage = totalMemory/totalSets;
                
                ObjectIndexTuple stats(gtree.getNetworkName(),gtree.getNumNodes(),gtree.getNumEdges(),
                                    constants::OBJ_IDX_GTREE,processingTimeMs,memoryUsage,objTypes[i],objDensities[j],objVariable[l],setSize);

//                 std::cout << stats.getMultilineTupleString();

                this->outputCommandStats(statsOutputFile,stats.getTupleString());
            }
        }
    }
    std::cout << "Occurrence lists successfully created for " << gtree.getNetworkName() << std::endl;
}

void ExperimentsCommand::buildAssociationDirectories(std::string routeOverlayIdxFile, std::size_t numSets, std::vector<double> objDensities, 
                                                   std::vector<std::string> objTypes, std::vector<int> objVariable, std::string filePathPrefix, std::string statsOutputFile, 
                                                   std::vector<std::string>& parameterKeys, std::vector<std::string>& parameterValues)
{
    ROAD road = serialization::getIndexFromBinaryFile<ROAD>(routeOverlayIdxFile);
    
    int setSize, setVariable;
    std::string setType;
    double setDensity;
    std::string objSetFile, objIdxOutputFile;
    std::vector<NodeID> objectNodes;
    
    StopWatch sw;
    double totalTime, totalMemory;
    int totalSets = numSets;
    for (std::size_t i = 0; i < objTypes.size(); ++i) {
        for (std::size_t j = 0; j < objDensities.size(); ++j) {
            for (std::size_t l = 0; l < objVariable.size(); ++l) {
                totalTime = 0;
                totalMemory = 0;
                
                for (std::size_t k = 0; k < numSets; ++k) {
                    objSetFile = filePathPrefix + "/obj_indexes/" + utility::constructObjsectSetFileName(road.getNetworkName(),objTypes[i],objDensities[j],objVariable[l],k);
                    std::vector<NodeID> objectNodes = utility::getPointSetFromFile(objSetFile,setType,setDensity,setSize,setVariable);
                    //assert(objTypes[i] == setType && objDensities[j] == setDensity);

                    sw.start();
                    AssociationDirectory assocDir(setType,setDensity,setVariable,setSize,road.getRnetTreeSize());
                    for (auto objIt = objectNodes.begin(); objIt != objectNodes.end(); ++objIt) {
                        // Find the ROAD leaf Rnet for this object and add entry into association
                        // directory entry for that leaf then propagant up Rnet hierarchy
                        assocDir.addObject(*objIt);
                        
                        for (int rnetIdx: road.routeOverlay[*objIt].leafIdxs) {
                            assocDir.addRnetAssociation(rnetIdx);
                        
                            // Propagate this to occurence lists of parents of leaf node
                            int parentRnetIdx = road.getParentRnet(rnetIdx);
                            while (parentRnetIdx != -1) {
                                assocDir.addRnetAssociation(parentRnetIdx);
                                
                                // Go up a level (until we reach root)
                                parentRnetIdx = road.getParentRnet(parentRnetIdx);
                            }
                        }
                    }
                    sw.stop();
                    totalTime += sw.getTimeUs();
                    totalMemory += assocDir.computeIndexSize()*1024; // in KB
                    
                    if (objDensities[j] != 1 || k == 0 || objTypes[i] != constants::RAND_OBJ_SET) {
                        // If density is one we only need to create one object index                    
                        objIdxOutputFile = filePathPrefix + "/obj_indexes/" + utility::constructObjectIndexFileName(road.getNetworkName(),constants::OBJ_IDX_ROAD,objTypes[i],objDensities[j],objVariable[l],k,parameterKeys,parameterValues);
                        serialization::outputIndexToBinaryFile<AssociationDirectory>(assocDir,objIdxOutputFile);
                    }
                    if (objDensities[j] == 1 && objTypes[i] == constants::RAND_OBJ_SET) {
                        totalSets = 1;
                        // If density is 1 then no don't iterate more than once
                        break;
                    } else {
                        totalSets = numSets;
                    }
                    
                }
                
                // Average stats over 100 object sets and then write to stats file
                double processingTimeMs = totalTime/totalSets;
                double memoryUsage = totalMemory/totalSets;
                
                ObjectIndexTuple stats(road.getNetworkName(),road.getNumNodes(),road.getNumEdges(),
                                    constants::OBJ_IDX_ROAD,processingTimeMs,memoryUsage,objTypes[i],objDensities[j],objVariable[l],setSize);

//                 std::cout << stats.getMultilineTupleString();

                this->outputCommandStats(statsOutputFile,stats.getTupleString());
            }
        }
    }
    std::cout << "Association directories successfully created for " << road.getNetworkName() << std::endl;
}

void ExperimentsCommand::buildQuadtrees(Graph& graph, std::vector<int>& maxLeafSizes, std::size_t numSets, std::vector<double> objDensities, 
                                        std::vector<std::string> objTypes, std::vector<int> objVariable, std::string filePathPrefix, std::string statsOutputFile)
{
    int setSize, setVariable;
    std::string setType;
    double setDensity;
    std::string objSetFile, objIdxOutputFile;
    std::vector<NodeID> objectNodes;
    std::vector<std::string> parameterKeys, parameterValues;
    
    StopWatch sw;
    double totalTime, totalMemory;
    int totalSets = numSets;
    for (std::size_t i = 0; i < objTypes.size(); ++i) {
        for (std::size_t j = 0; j < objDensities.size(); ++j) {
            for (std::size_t m = 0; m < objVariable.size(); ++m) {
                for (std::size_t l = 0; l < maxLeafSizes.size(); ++l) {
                    parameterKeys.clear();
                    parameterValues.clear();
                    parameterKeys.push_back("maxleafsize");
                    parameterValues.push_back(std::to_string(maxLeafSizes[l]));
                    totalTime = 0;
                    totalMemory = 0;
                    
                    for (std::size_t k = 0; k < numSets; ++k) {
                        objSetFile = filePathPrefix + "/obj_indexes/" + utility::constructObjsectSetFileName(graph.getNetworkName(),objTypes[i],objDensities[j],objVariable[m],k);
                        std::vector<NodeID> objectNodes = utility::getPointSetFromFile(objSetFile,setType,setDensity,setSize,setVariable);
                        //assert(objTypes[i] == setType && objDensities[j] == setDensity);

                        sw.start();
                        SimpleQuadtree qt(maxLeafSizes[l]);
                        qt.buildFromPointSet(graph,objectNodes,setType,setDensity,setVariable,setSize);
                        sw.stop();
                        totalTime += sw.getTimeUs();
                        totalMemory += qt.computeIndexSize()*1024; // in KB
                        
                        if (objDensities[j] != 1 || k == 0 || objTypes[i] != constants::RAND_OBJ_SET) {
                            // If density is one we only need to create one object index                    
                            objIdxOutputFile = filePathPrefix + "/obj_indexes/" + utility::constructObjectIndexFileName(graph.getNetworkName(),constants::OBJ_IDX_QUADTREE,objTypes[i],objDensities[j],objVariable[m],k,parameterKeys,parameterValues);
                            serialization::outputIndexToBinaryFile<SimpleQuadtree>(qt,objIdxOutputFile);
                        }
                        
                    
                        if (objDensities[j] == 1 && objTypes[i] == constants::RAND_OBJ_SET) {
                            totalSets = 1;
                            // If density is 1 then no don't iterate more than once
                            break;
                        } else {
                            totalSets = numSets;
                        }
                    }
                
                    // Average stats over 100 object sets and then write to stats file
                    double processingTimeMs = totalTime/totalSets;
                    double memoryUsage = totalMemory/totalSets;
                    
                    ObjectIndexTuple stats(graph.getNetworkName(),graph.getNumNodes(),graph.getNumEdges(),
                                        constants::OBJ_IDX_QUADTREE,processingTimeMs,memoryUsage,objTypes[i],objDensities[j],objVariable[m],setSize);
                    stats.addSupplementaryFields("obj_qt_max_leaf_size",std::to_string(maxLeafSizes[l]));

//                     std::cout << stats.getMultilineTupleString();

                    this->outputCommandStats(statsOutputFile,stats.getTupleString());
                }
            }
        }
    }
    std::cout << "Quadtrees successfully created for " << graph.getNetworkName() << std::endl;
}

void ExperimentsCommand::buildRtrees(Graph& graph, std::vector<int>& branchFactors, std::size_t numSets, std::vector<double> objDensities, 
                                     std::vector<std::string> objTypes, std::vector<int> objVariable, std::string filePathPrefix, std::string statsOutputFile)
{
    int setSize, setVariable;
    std::string setType;
    double setDensity;
    std::string objSetFile, objIdxOutputFile;
    std::vector<NodeID> objectNodes;
    std::vector<std::string> parameterKeys, parameterValues;

    StopWatch sw;
    double totalTime, totalMemory;
    int totalSets = numSets;
    for (std::size_t i = 0; i < objTypes.size(); ++i) {
        for (std::size_t j = 0; j < objDensities.size(); ++j) {
            for (std::size_t m = 0; m < objVariable.size(); ++m) {
                for (std::size_t l = 0; l < branchFactors.size(); ++l) {
                    parameterKeys.clear();
                    parameterValues.clear();
                    parameterKeys.push_back("branchfactor");
                    parameterValues.push_back(std::to_string(branchFactors[l]));
                    totalTime = 0;
                    totalMemory = 0;

                    for (std::size_t k = 0; k < numSets; ++k) {
                        objSetFile = filePathPrefix + "/obj_indexes/" + utility::constructObjsectSetFileName(graph.getNetworkName(),objTypes[i],objDensities[j],objVariable[m],k);
                        std::vector<NodeID> objectNodes = utility::getPointSetFromFile(objSetFile,setType,setDensity,setSize,setVariable);
                        //assert(objTypes[i] == setType && objDensities[j] == setDensity);

                        sw.start();
                        std::vector<CoordinatePair> objectCoords;
                        for (std::size_t i = 0; i < objectNodes.size(); ++i) {
                            CoordinatePair objectCoordPair;
                            graph.getCoordinates(objectNodes[i],objectCoordPair.first,objectCoordPair.second);
                            objectCoords.push_back(objectCoordPair);
                        }
                        StaticRtree rtree(branchFactors[l],setType,setDensity,setVariable,setSize);
                        rtree.bulkLoad(objectNodes,objectCoords);
                        sw.stop();
                        totalTime += sw.getTimeUs();
                        totalMemory += rtree.computeIndexSize()*1024; // in KB

                        if (objDensities[j] != 1 || k == 0 || objTypes[i] != constants::RAND_OBJ_SET) {
                            // If density is one we only need to create one object index                    
                            objIdxOutputFile = filePathPrefix + "/obj_indexes/" + utility::constructObjectIndexFileName(graph.getNetworkName(),constants::OBJ_IDX_RTREE,objTypes[i],objDensities[j],objVariable[m],k,parameterKeys,parameterValues);
                            serialization::outputIndexToBinaryFile<StaticRtree>(rtree,objIdxOutputFile);
                        }

                        if (objDensities[j] == 1 && objTypes[i] == constants::RAND_OBJ_SET) {
                            totalSets = 1;
                            // If density is 1 then no don't iterate more than once
                            break;
                        } else {
                            totalSets = numSets;
                        }
                    }

                    // Average stats over 100 object sets and then write to stats file
                    double processingTimeMs = totalTime/totalSets;
                    double memoryUsage = totalMemory/totalSets;
                    
                    ObjectIndexTuple stats(graph.getNetworkName(),graph.getNumNodes(),graph.getNumEdges(),
                                        constants::OBJ_IDX_RTREE,processingTimeMs,memoryUsage,objTypes[i],objDensities[j],objVariable[m],setSize);
                    stats.addSupplementaryFields("rtree_branch_factor",std::to_string(branchFactors[l]));

//                     std::cout << stats.getMultilineTupleString();

                    this->outputCommandStats(statsOutputFile,stats.getTupleString());
                }
            }
        }
    }
    std::cout << "Rtrees successfully created for " << graph.getNetworkName() << std::endl;
}

void ExperimentsCommand::runINEQueries(Graph& graph, std::vector<NodeID>& queryNodes, std::vector<int>& kValues, std::size_t numSets, 
                                       std::vector<double> objDensities, std::vector<std::string> objTypes, std::vector<int> objVariable, std::string filePathPrefix, 
                                       std::string statsOutputFile, std::vector<std::string> specialFields)
{
    std::vector<NodeID> kNNs;
    std::vector<EdgeWeight> kNNDistances;
    std::vector<NodeID> localQueryNodes;

    std::string setType;
    double setDensity;
    int setSize, setVariable;

    StopWatch sw;
    double totalQueryTime;
    int totalQueries = numSets*queryNodes.size();
    localQueryNodes = queryNodes;
    int totalObjects;

    INE ine;
    Statistics knnStats;
    
    for (std::size_t i = 0; i < objTypes.size(); ++i) {
        if (objTypes[i] == constants::MINND_OBJ_SET || objTypes[i] == constants::MINMAXND_OBJ_SET) {
            // Reduce number of queries for INE over long distances or it will take days
            numSets = std::ceil(static_cast<double>(numSets)/10);
            SetGenerator sg;
            localQueryNodes = sg.generateSampleSubet(queryNodes, std::ceil(static_cast<double>(queryNodes.size())/5));
            totalQueries = numSets*localQueryNodes.size();
            if (numSets == 0) {
                std::cerr << "Must use at least 10 object sets for min. obj. distance experiments" << std::endl;
                std::exit(1);
            }
        }
        for (std::size_t j = 0; j < objDensities.size(); ++j) {
            for (std::size_t m = 0; m < objVariable.size(); ++m) {
                for (std::size_t k = 0; k < kValues.size(); ++k) {
                    totalQueryTime = 0;
                    totalObjects = 0;
#if defined(COLLECT_STATISTICS)
                    knnStats.clear();
#endif
                    for (std::size_t l = 0; l < numSets; ++l) {
                        std::string objSetFile = filePathPrefix + "/obj_indexes/" + utility::constructObjsectSetFileName(graph.getNetworkName(),objTypes[i],objDensities[j],objVariable[m],l);
                        std::vector<NodeID> objectNodes = utility::getPointSetFromFile(objSetFile,setType,setDensity,setSize,setVariable);
                        totalObjects += setSize;
            
                        graph.resetAllObjects();
                        graph.parseObjectSet(objectNodes);
                        
                        for (auto queryNodeIt = localQueryNodes.begin(); queryNodeIt != localQueryNodes.end(); ++queryNodeIt) {
                            kNNs.clear();
                            kNNDistances.clear();
                            kNNs.reserve(kValues[k]);
                            kNNDistances.reserve(kValues[k]);
                            sw.reset();
                            sw.start();
                            ine.getKNNs(graph,kValues[k],*queryNodeIt,kNNs,kNNDistances);
                            sw.stop();
                            totalQueryTime += sw.getTimeUs();
#if defined(COLLECT_STATISTICS)
                            knnStats.mergeStatistics(ine.stats);
#endif
                        }        
                    }

                    double queryTimeMs = totalQueryTime/totalQueries;
#if defined(COLLECT_STATISTICS)
                    knnStats.normalizeStatistics(totalQueries);
#endif
                    
                    // Collect stats and return to output to file
                    kNNs.clear(); // Clear so that we don't pass last executed queries results to stats tuple
                    kNNDistances.clear();
                    KnnQueryTuple stats(graph.getNetworkName(),graph.getNumNodes(),graph.getNumEdges(),totalQueries,constants::INE_KNN_QUERY,
                                        kValues[k],queryTimeMs,objTypes[i],objDensities[j],objVariable[m],static_cast<int>(totalObjects/numSets),kNNs,kNNDistances);
                    stats.setAdditionalFields(specialFields);
#if defined(COLLECT_STATISTICS)
                    knnStats.populateTupleFields(stats,0);
#endif
                    this->outputCommandStats(statsOutputFile,stats.getTupleString());
                }
            }
        }
    }
    std::cout << "INE kNN queries successfully executed for " << graph.getNetworkName() << std::endl;
}

void ExperimentsCommand::runGtreeQueries(Graph& graph, std::string gtreeIdxFile, std::vector<NodeID>& queryNodes, 
                                         std::vector<int>& kValues, std::size_t numSets, std::vector<double> objDensities, 
                                         std::vector<std::string> objTypes, std::vector<int> objVariable, std::string filePathPrefix, std::string statsOutputFile, 
                                         bool verifyKNN,  std::vector<std::string>& parameterKeys, std::vector<std::string>& parameterValues,
                                         std::vector<std::string> specialFields)
{
    std::vector<NodeID> kNNs, ineKNNs;
    std::vector<EdgeWeight> kNNDistances, ineKNNDistances;

    std::string objSetType;
    double objSetDensity;
    int objSetSize, objSetVariable;
    
    StopWatch sw;
    double totalQueryTime;
    int totalQueries = numSets*queryNodes.size();
    int totalObjects;
    
    INE ine;
    Statistics knnStats;
    std::string message;
    
    Gtree gtree = serialization::getIndexFromBinaryFile<Gtree>(gtreeIdxFile);
#if defined(GTREE_STL_HASH_TABLE_DIST_MATRIX) || defined(GTREE_GOOGLE_DENSEHASH_DIST_MATRIX)
    // We create hash-table versions of distance matrix for different testing
    sw.reset();
    sw.start();
    gtree.populateUnorderedMapDistanceMatrix();
    sw.stop();
//     std::cout << "Hash-Table Distance Matrix Population Time: " << sw.getTimeMs() << "ms" << std::endl;
//     std::cout << "Hash-Table Distance Matrix Memory Usage: " << gtree.computeDistanceMatrixMemoryUsage() << "MB" << std::endl;    
#endif
    
    for (std::size_t i = 0; i < objTypes.size(); ++i) {
        for (std::size_t j = 0; j < objDensities.size(); ++j) {
            for (std::size_t m = 0; m < objVariable.size(); ++m) {
                for (std::size_t k = 0; k < kValues.size(); ++k) {
                    totalQueryTime = 0;
                    totalObjects = 0;
#if defined(COLLECT_STATISTICS)
                    knnStats.clear();
#endif
                    for (std::size_t l = 0; l < numSets; ++l) {
                        std::string objIdxFilePath = filePathPrefix + "/obj_indexes/" + utility::constructObjectIndexFileName(gtree.getNetworkName(),constants::OBJ_IDX_GTREE,objTypes[i],objDensities[j],objVariable[m],l,parameterKeys,parameterValues);
                        OccurenceList occList = serialization::getIndexFromBinaryFile<OccurenceList>(objIdxFilePath);
                        totalObjects += occList.getObjSetSize();
            
                        if (verifyKNN) {
                            graph.resetAllObjects();
                            std::string objSetFile = filePathPrefix + "/obj_indexes/" + utility::constructObjsectSetFileName(graph.getNetworkName(),objTypes[i],objDensities[j],objVariable[m],l);
                            graph.parseObjectFile(objSetFile,objSetType,objSetDensity,objSetVariable,objSetSize);
                        }
                        
                        for (auto queryNodeIt = queryNodes.begin(); queryNodeIt != queryNodes.end(); ++queryNodeIt) {
                            kNNs.clear();
                            kNNDistances.clear();
                            kNNs.reserve(kValues[k]);
                            kNNDistances.reserve(kValues[k]);
                            sw.reset();
                            sw.start();
                            gtree.getKNNs(occList,kValues[k],*queryNodeIt,kNNs,kNNDistances,graph);
                            sw.stop();
                            totalQueryTime += sw.getTimeUs();
#if defined(COLLECT_STATISTICS)
                            knnStats.mergeStatistics(gtree.stats);
#endif
                            
                            if (verifyKNN) {
                                ineKNNs.clear();
                                ineKNNDistances.clear();
                                ine.getKNNs(graph,kValues[k],*queryNodeIt,ineKNNs,ineKNNDistances);
                                if (!utility::verifyKNN(ineKNNs,ineKNNDistances,kNNs,kNNDistances,false,kValues[k],message,true)) {
                                    std::cout << "Verfication failed for Gtree on object index " << objIdxFilePath << " for query node " << *queryNodeIt << " with k = " << kValues[k] << std::endl;
                                    std::cout << "Message: " << message << std::endl;
                                }
                            }
                        }        
                    }

                    double queryTimeMs = totalQueryTime/totalQueries;
#if defined(COLLECT_STATISTICS)
                    knnStats.normalizeStatistics(totalQueries);
#endif

                    // Collect stats and return to output to file
                    kNNs.clear(); // Clear so that we don't pass last executed queries results to stats tuple
                    kNNDistances.clear();
                    KnnQueryTuple stats(gtree.getNetworkName(),gtree.getNumNodes(),gtree.getNumEdges(),totalQueries,constants::GTREE_KNN_QUERY,
                                        kValues[k],queryTimeMs,objTypes[i],objDensities[j],objVariable[m],static_cast<int>(totalObjects/numSets),kNNs,kNNDistances);
                    stats.setAdditionalFields(specialFields);
                    stats.addSupplementaryFields("max_leaf_size",std::to_string(gtree.getMaxLeafSize()));
                    stats.addSupplementaryFields("fanout",std::to_string(gtree.getFanout()));
#if defined(COLLECT_STATISTICS)
                    knnStats.populateTupleFields(stats,0);
#endif
                    this->outputCommandStats(statsOutputFile,stats.getTupleString());
                }
            }
        }
    }
    std::cout << "Gtree kNN queries successfully executed for " << gtree.getNetworkName() << std::endl;
    if (verifyKNN) {
        std::cout << "Gtree kNN verification completed for " << graph.getNetworkName() << std::endl;
    }
}

void ExperimentsCommand::runROADQueries(Graph& graph, std::string routeOverlayIdxFile, std::vector<NodeID>& queryNodes, 
                                        std::vector<int>& kValues, std::size_t numSets, std::vector<double> objDensities, 
                                        std::vector<std::string> objTypes, std::vector<int> objVariable, std::string filePathPrefix, std::string statsOutputFile, 
                                        bool verifyKNN,  std::vector<std::string>& parameterKeys, std::vector<std::string>& parameterValues,
                                        std::vector<std::string> specialFields)
{
    std::vector<NodeID> kNNs, ineKNNs;
    std::vector<EdgeWeight> kNNDistances, ineKNNDistances;

    std::string objSetType;
    double objSetDensity;
    int objSetSize, objSetVariable;
    
    StopWatch sw;
    double totalQueryTime;
    int totalQueries = numSets*queryNodes.size();
    int totalObjects;
    
    INE ine;
    Statistics knnStats;
    std::string message;

    ROAD road = serialization::getIndexFromBinaryFile<ROAD>(routeOverlayIdxFile);
    
    for (std::size_t i = 0; i < objTypes.size(); ++i) {
        for (std::size_t j = 0; j < objDensities.size(); ++j) {
            for (std::size_t m = 0; m < objVariable.size(); ++m) {
                for (std::size_t k = 0; k < kValues.size(); ++k) {
                    totalQueryTime = 0;
                    totalObjects = 0;
#if defined(COLLECT_STATISTICS)
                    knnStats.clear();
#endif
                    for (std::size_t l = 0; l < numSets; ++l) {
                        std::string objIdxFilePath = filePathPrefix + "/obj_indexes/" + utility::constructObjectIndexFileName(road.getNetworkName(),constants::OBJ_IDX_ROAD,objTypes[i],objDensities[j],objVariable[m],l,parameterKeys,parameterValues);
                        AssociationDirectory assocDir = serialization::getIndexFromBinaryFile<AssociationDirectory>(objIdxFilePath);
                        totalObjects += assocDir.getObjSetSize();
                        
                        if (verifyKNN) {
                            graph.resetAllObjects();
                            std::string objSetFile = filePathPrefix + "/obj_indexes/" + utility::constructObjsectSetFileName(graph.getNetworkName(),objTypes[i],objDensities[j],objVariable[m],l);
                            graph.parseObjectFile(objSetFile,objSetType,objSetDensity,objSetVariable,objSetSize);
                        }
                        
                        for (auto queryNodeIt = queryNodes.begin(); queryNodeIt != queryNodes.end(); ++queryNodeIt) {
                            kNNs.clear();
                            kNNDistances.clear();
                            kNNs.reserve(kValues[k]);
                            kNNDistances.reserve(kValues[k]);
                            sw.reset();
                            sw.start();
                            road.getKNNs(assocDir,kValues[k],*queryNodeIt,kNNs,kNNDistances);
                            sw.stop();
                            totalQueryTime += sw.getTimeUs();
#if defined(COLLECT_STATISTICS)
                            knnStats.mergeStatistics(road.stats);
#endif

                            if (verifyKNN) {
                                ineKNNs.clear();
                                ineKNNDistances.clear();
                                ine.getKNNs(graph,kValues[k],*queryNodeIt,ineKNNs,ineKNNDistances);
                                if (!utility::verifyKNN(ineKNNs,ineKNNDistances,kNNs,kNNDistances,false,kValues[k],message,true)) {
                                    std::cout << "Verfication failed for ROAD on object index " << objIdxFilePath << " for query node " << *queryNodeIt << " with k = " << kValues[k] << std::endl;
                                    std::cout << "Message: " << message << std::endl;
                                }
                            }
                        }        
                    }

                    double queryTimeMs = totalQueryTime/totalQueries;  
#if defined(COLLECT_STATISTICS)
                    knnStats.normalizeStatistics(totalQueries);
#endif

                    // Collect stats and return to output to file
                    kNNs.clear(); // Clear so that we don't pass last executed queries results to stats tuple
                    kNNDistances.clear();
                    KnnQueryTuple stats(road.getNetworkName(),road.getNumNodes(),road.getNumEdges(),totalQueries,constants::ROAD_KNN_QUERY,
                                        kValues[k],queryTimeMs,objTypes[i],objDensities[j],objVariable[m],static_cast<int>(totalObjects/numSets),kNNs,kNNDistances);
                    stats.setAdditionalFields(specialFields);
                    stats.addSupplementaryFields("levels",std::to_string(road.getNumLevels()));
                    stats.addSupplementaryFields("fanout",std::to_string(road.getFanout()));
#if defined(COLLECT_STATISTICS)
                    knnStats.populateTupleFields(stats,0);
#endif
                    this->outputCommandStats(statsOutputFile,stats.getTupleString());
                }
            }
        }
    }
    std::cout << "ROAD kNN queries successfully executed for " << road.getNetworkName() << std::endl;
    if (verifyKNN) {
        std::cout << "ROAD kNN verification completed for " << road.getNetworkName() << std::endl;
    }
}

void ExperimentsCommand::runSILCBasedQueries(Graph& graph, std::vector<std::string>& methods, std::string silcIdxFile, std::string suppIdxFile, 
                                 std::vector<int>& maxLeafSizes, std::vector<NodeID>& queryNodes, std::vector<int>& kValues, std::size_t numSets, 
                                 std::vector<double> objDensities, std::vector<std::string> objTypes, std::vector<int> objVariable, std::string filePathPrefix, 
                                 std::string statsOutputFile, bool verifyKNN, std::vector<std::string> specialFields)
{
    std::vector<NodeID> kNNs, ineKNNs;
    std::vector<EdgeWeight> kNNDistances, ineKNNDistances;
    
    std::string objSetType;
    double objSetDensity;
    int objSetSize, objSetVariable;
    
    std::vector<std::string> parameterKeys, parameterValues;
    
    StopWatch sw;
    double totalQueryTime;
    int totalQueries = numSets*queryNodes.size();
    int totalObjects;
    
    INE ine;
    Statistics knnStats;
    std::string message, desc;    
    
//     SILCPathOracle silc = serialization::getIndexFromBinaryFile<SILCPathOracle>(silcIdxFile);  
//     Junction junc = serialization::getIndexFromBinaryFile<Junction>(suppIdxFile);
    SILCPathOracle silc;
    serialization::populateIndexFromBinaryFile<SILCPathOracle>(silcIdxFile,silc);
    Junction junc;
    serialization::populateIndexFromBinaryFile<Junction>(suppIdxFile,junc);
    
    for (std::size_t x = 0; x < methods.size(); ++x) {
        if (methods[x] == constants::DB_RTREE_KNN_QUERY) {
            desc = "DistBrws using Rtree";
        } else if (methods[x] == constants::SILC_KNN_QUERY) {
            desc = "SILC";
        } else if (methods[x] == constants::DISTBRWS_KNN_QUERY) {
            desc = "DistBrws";
        } else if (methods[x] == constants::OPT_SILC_KNN_QUERY) {
            desc = "Optimised SILC";
        } else if (methods[x] == constants::OPT_DISTBRWS_KNN_QUERY) {
            desc = "Optimised DistBrws";
        } else {
            std::cerr << "Invalid SILC kNN method provided!" << std::endl;
            std::exit(1);
        }
        
        for (std::size_t i = 0; i < objTypes.size(); ++i) {
            for (std::size_t m = 0; m < maxLeafSizes.size(); ++m) {
                for (std::size_t j = 0; j < objDensities.size(); ++j) {
                for (std::size_t n = 0; n < objVariable.size(); ++n) {
                        for (std::size_t k = 0; k < kValues.size(); ++k) {
                            totalQueryTime = 0;
                            totalObjects = 0;
#if defined(COLLECT_STATISTICS)
                            knnStats.clear();
#endif
                            for (std::size_t l = 0; l < numSets; ++l) {
                                std::string objIdxFilePath, objSetFilePath;
                                if (methods[x] == constants::DB_RTREE_KNN_QUERY) {
                                    parameterKeys.clear();
                                    parameterValues.clear();
                                    parameterKeys.push_back("branchfactor");
                                    parameterValues.push_back(std::to_string(maxLeafSizes[m])); // maxLeafSizes in this case holds the branch factor(s)
                                    objIdxFilePath = filePathPrefix + "/obj_indexes/" + utility::constructObjectIndexFileName(graph.getNetworkName(),constants::OBJ_IDX_RTREE,objTypes[i],objDensities[j],objVariable[n],l,parameterKeys,parameterValues);
                                    StaticRtree rtree = serialization::getIndexFromBinaryFile<StaticRtree>(objIdxFilePath);
                                    totalObjects += rtree.getObjSetSize();
                                    
                                    if (verifyKNN) {
                                        graph.resetAllObjects();
                                        std::string objSetFile = filePathPrefix + "/obj_indexes/" + utility::constructObjsectSetFileName(graph.getNetworkName(),objTypes[i],objDensities[j],objVariable[n],l);
                                        graph.parseObjectFile(objSetFile,objSetType,objSetDensity,objSetVariable,objSetSize);
                                    }
                                    
                                    for (auto queryNodeIt = queryNodes.begin(); queryNodeIt != queryNodes.end(); ++queryNodeIt) {
                                        kNNs.clear();
                                        kNNDistances.clear();
                                        kNNs.reserve(kValues[k]);
                                        kNNDistances.reserve(kValues[k]);
                                        sw.reset();
                                        sw.start();
                                        silc.getKNNsByDistanceBrowsingViaRtree(rtree,kValues[k],*queryNodeIt,kNNs,kNNDistances,graph,junc);
                                        sw.stop();
                                        totalQueryTime += sw.getTimeUs();
#if defined(COLLECT_STATISTICS)
                                        knnStats.mergeStatistics(silc.stats);
#endif
                                        if (verifyKNN) {
                                            ineKNNs.clear();
                                            ineKNNDistances.clear();
                                            ine.getKNNs(graph,kValues[k],*queryNodeIt,ineKNNs,ineKNNDistances);
                                            // We need to fully refine all kNN distance in order to verify
                                            silc.findRealKNNDistances(graph,*queryNodeIt,kNNs,kNNDistances);
                                            if (!utility::verifyKNN(ineKNNs,ineKNNDistances,kNNs,kNNDistances,false,kValues[k],message,true)) {
                                                std::cout << "Verfication failed for " << desc << " on object index " << objIdxFilePath << " for query node " << *queryNodeIt << " with k = " << kValues[k] << std::endl;
                                                std::cout << "Message: " << message << std::endl;
                                            }
                                        }
                                    }     

                                } else {
                                    parameterKeys.clear();
                                    parameterValues.clear();
                                    parameterKeys.push_back("maxleafsize");
                                    parameterValues.push_back(std::to_string(maxLeafSizes[m]));
                                    objIdxFilePath = filePathPrefix + "/obj_indexes/" + utility::constructObjectIndexFileName(graph.getNetworkName(),constants::OBJ_IDX_QUADTREE,objTypes[i],objDensities[j],objVariable[n],l,parameterKeys,parameterValues);
                                    SimpleQuadtree objHierarchy = serialization::getIndexFromBinaryFile<SimpleQuadtree>(objIdxFilePath);
                                    totalObjects += objHierarchy.getObjSetSize();
                                    
                                    if (verifyKNN) {
                                        graph.resetAllObjects();
                                        std::string objSetFile = filePathPrefix + "/obj_indexes/" + utility::constructObjsectSetFileName(graph.getNetworkName(),objTypes[i],objDensities[j],objVariable[n],l);
                                        graph.parseObjectFile(objSetFile,objSetType,objSetDensity,objSetVariable,objSetSize);
                                    }

                                    for (auto queryNodeIt = queryNodes.begin(); queryNodeIt != queryNodes.end(); ++queryNodeIt) {
                                        kNNs.clear();
                                        kNNDistances.clear();
                                        kNNs.reserve(kValues[k]);
                                        kNNDistances.reserve(kValues[k]);
                                        if (methods[x] == constants::SILC_KNN_QUERY) {
                                            sw.reset();
                                            sw.start();
                                            silc.getKNNs(objHierarchy,kValues[k],*queryNodeIt,kNNs,kNNDistances,graph);
                                            sw.stop();
                                        } else if (methods[x] == constants::DISTBRWS_KNN_QUERY) {
                                            sw.reset();
                                            sw.start();
                                            silc.getKNNsByDistanceBrowsing(objHierarchy,kValues[k],*queryNodeIt,kNNs,kNNDistances,graph);
                                            sw.stop();
                                        } else if (methods[x] == constants::OPT_SILC_KNN_QUERY) {
                                            sw.reset();
                                            sw.start();
                                            silc.getKNNsByOptimisedSILC(objHierarchy,kValues[k],*queryNodeIt,kNNs,kNNDistances,graph,junc);
                                            sw.stop();
                                        } else if (methods[x] == constants::OPT_DISTBRWS_KNN_QUERY) {
                                            sw.reset();
                                            sw.start();
                                            silc.getKNNsByOptimisedDistanceBrowsing(objHierarchy,kValues[k],*queryNodeIt,kNNs,kNNDistances,graph,junc);
                                            sw.stop();
                                        }
                                        totalQueryTime += sw.getTimeUs();
#if defined(COLLECT_STATISTICS)
                                        knnStats.mergeStatistics(silc.stats);
#endif
                                        if (verifyKNN) {
                                            ineKNNs.clear();
                                            ineKNNDistances.clear();
                                            ine.getKNNs(graph,kValues[k],*queryNodeIt,ineKNNs,ineKNNDistances);
                                            // We need to fully refine all kNN distance in order to verify
                                            silc.findRealKNNDistances(graph,*queryNodeIt,kNNs,kNNDistances);
                                            if (!utility::verifyKNN(ineKNNs,ineKNNDistances,kNNs,kNNDistances,false,kValues[k],message,true)) {
                                                std::cout << "Verfication failed for " << desc << " on object index " << objIdxFilePath << " for query node " << *queryNodeIt << " with k = " << kValues[k] << std::endl;
                                                std::cout << "Message: " << message << std::endl;
                                            }
                                        }
                                    }                                        
                                    
                                } 
                            }

                            double queryTimeMs = totalQueryTime/totalQueries;
#if defined(COLLECT_STATISTICS)
                            knnStats.normalizeStatistics(totalQueries);
#endif
                            
                            // Collect stats and return to output to file
                            kNNs.clear(); // Clear so that we don't pass last executed queries results to stats tuple
                            kNNDistances.clear();
                            KnnQueryTuple stats(silc.getNetworkName(),silc.getNumNodes(),silc.getNumEdges(),totalQueries,methods[x],
                                                kValues[k],queryTimeMs,objTypes[i],objDensities[j],objVariable[n],static_cast<int>(totalObjects/numSets),kNNs,kNNDistances);
                            stats.setAdditionalFields(specialFields);
                            if (methods[x] == constants::DB_RTREE_KNN_QUERY) {
                                stats.addSupplementaryFields("obj_rt_branch_factor",std::to_string(maxLeafSizes[m]));
                            } else {
                                stats.addSupplementaryFields("obj_qt_max_leaf_size",std::to_string(maxLeafSizes[m]));
                            }
#if defined(COLLECT_STATISTICS)
                            knnStats.populateTupleFields(stats,0);
#endif
                            this->outputCommandStats(statsOutputFile,stats.getTupleString());
                        }
                    }
                }
            }
        }

        std::cout << desc << " kNN queries successfully executed for " << silc.getNetworkName() << std::endl;
        if (verifyKNN) {
            std::cout << desc << " kNN verification completed for " << silc.getNetworkName() << std::endl;
        }
    }
}

void ExperimentsCommand::runIERQueries(Graph& graph, std::string index, std::string idxFilePath, std::vector<int>& branchFactors, 
                                       std::vector<NodeID>& queryNodes, std::vector<int>& kValues, std::size_t numSets, 
                                       std::vector<double> objDensities, std::vector<std::string> objTypes, std::vector<int> objVariable, std::string filePathPrefix,
                                       std::string statsOutputFile, bool verifyKNN, std::vector<std::string> specialFields)
{
    std::vector<NodeID> kNNs, ineKNNs;
    std::vector<EdgeWeight> kNNDistances, ineKNNDistances;
    std::vector<NodeID> localQueryNodes;

    std::string objSetType;
    double objSetDensity;
    int objSetSize, objSetVariable;
    
    std::vector<std::string> parameterKeys, parameterValues;
    
    StopWatch sw;
    double totalQueryTime;
    int totalQueries = numSets*queryNodes.size();
    int totalObjects;
    
    INE ine;
    Statistics knnStats;
    std::string message;
    
    IER ier;
    
    Gtree gtree;
    SILCPathOracle silc;
    PrunedHighwayLabeling phl;
    std::string desc = "", method = "";
    localQueryNodes = queryNodes;
    if (index == constants::GTREE_SPDIST_QUERY) {
        gtree = serialization::getIndexFromBinaryFile<Gtree>(idxFilePath);
        method = constants::IER_GTREE_KNN_QUERY;
        desc = " using G-tree";
    } else if (index == constants::SILC_SPDIST_QUERY) {
        silc = serialization::getIndexFromBinaryFile<SILCPathOracle>(idxFilePath);
        method = constants::IER_SILC_KNN_QUERY;
        desc = " using SILC";
    } else if (index == constants::PHL_SPDIST_QUERY) {
        phl.LoadLabel(idxFilePath.c_str());
        method = constants::IER_PHL_KNN_QUERY;
        desc = " using PHL";
    } else {
        method = constants::IER_DIJKSTRA_KNN_QUERY;
        desc = " using Dijkstra";
        // Reduce number of queries for IER-Dijk or it will take days
        numSets = std::ceil(static_cast<double>(numSets)/10);
        SetGenerator sg;
        localQueryNodes = sg.generateSampleSubet(queryNodes, std::ceil(static_cast<double>(queryNodes.size())/5));
        totalQueries = numSets*localQueryNodes.size();
    }
    
    for (std::size_t i = 0; i < objTypes.size(); ++i) {
        for (std::size_t m = 0; m < branchFactors.size(); ++m) {
            parameterKeys.clear();
            parameterValues.clear();
            parameterKeys.push_back("branchfactor");
            parameterValues.push_back(std::to_string(branchFactors[m]));
            for (std::size_t j = 0; j < objDensities.size(); ++j) {
                for (std::size_t n = 0; n < objVariable.size(); ++n) {
                    for (std::size_t k = 0; k < kValues.size(); ++k) {
                        totalQueryTime = 0;
                        totalObjects = 0;
#if defined(COLLECT_STATISTICS)
                        knnStats.clear();
#endif
                        for (std::size_t l = 0; l < numSets; ++l) {
                            std::string objIdxFilePath = filePathPrefix + "/obj_indexes/" + utility::constructObjectIndexFileName(graph.getNetworkName(),constants::OBJ_IDX_RTREE,objTypes[i],objDensities[j],objVariable[n],l,parameterKeys,parameterValues);
                            StaticRtree rtree = serialization::getIndexFromBinaryFile<StaticRtree>(objIdxFilePath);
                            totalObjects += rtree.getObjSetSize();
                            
                            if (verifyKNN) {
                                graph.resetAllObjects();
                                std::string objSetFile = filePathPrefix + "/obj_indexes/" + utility::constructObjsectSetFileName(graph.getNetworkName(),objTypes[i],objDensities[j],objVariable[n],l);
                                graph.parseObjectFile(objSetFile,objSetType,objSetDensity,objSetVariable,objSetSize);
                            }
                            
                            for (auto queryNodeIt = localQueryNodes.begin(); queryNodeIt != localQueryNodes.end(); ++queryNodeIt) {
                                kNNs.clear();
                                kNNDistances.clear();
                                kNNs.reserve(kValues[k]);
                                kNNDistances.reserve(kValues[k]);
                                if (index == constants::GTREE_SPDIST_QUERY) {
                                    if (graph.getEdgeType() == constants::DISTANCE_WEIGHTS) {
                                        sw.reset();
                                        sw.start();
                                        ier.getKNNsByGtree(rtree,kValues[k],*queryNodeIt,kNNs,kNNDistances,graph,gtree);
                                        sw.stop();
                                    } else {
                                        //assert(graph.getEdgeType() == constants::TIME_WEIGHTS);
                                        sw.reset();
                                        sw.start();
                                        ier.getKNNsByGtreeTravelTimes(rtree,kValues[k],*queryNodeIt,kNNs,kNNDistances,graph,gtree);
                                        sw.stop();
                                    }
                                } else if (index == constants::SILC_SPDIST_QUERY) {
                                    sw.reset();
                                    sw.start();
                                    ier.getKNNsBySILC(rtree,kValues[k],*queryNodeIt,kNNs,kNNDistances,graph,silc);
                                    sw.stop();
                                } else if (index == constants::PHL_SPDIST_QUERY) {
                                    if (graph.getEdgeType() == constants::DISTANCE_WEIGHTS) {
                                        sw.reset();
                                        sw.start();
                                        ier.getKNNsByPHL(rtree,kValues[k],*queryNodeIt,kNNs,kNNDistances,graph,phl);
                                        sw.stop();
                                    } else {
                                        //assert(graph.getEdgeType() == constants::TIME_WEIGHTS);
                                        sw.reset();
                                        sw.start();
                                        ier.getKNNsByPHLTravelTimes(rtree,kValues[k],*queryNodeIt,kNNs,kNNDistances,graph,phl);
                                        sw.stop();
                                    }
                                } else {
                                    if (graph.getEdgeType() == constants::DISTANCE_WEIGHTS) {
                                        sw.reset();
                                        sw.start();
                                        ier.getKNNsByDijkstra(rtree,kValues[k],*queryNodeIt,kNNs,kNNDistances,graph);
                                        sw.stop();
                                    } else {
                                        //assert(graph.getEdgeType() == constants::TIME_WEIGHTS);
                                        sw.reset();
                                        sw.start();
                                        ier.getKNNsByDijkstraTravelTimes(rtree,kValues[k],*queryNodeIt,kNNs,kNNDistances,graph);
                                        sw.stop();
                                    }
                                }
                                totalQueryTime += sw.getTimeUs();
#if defined(COLLECT_STATISTICS)
                                knnStats.mergeStatistics(ier.stats);
#endif
                                if (verifyKNN) {
                                    ineKNNs.clear();
                                    ineKNNDistances.clear();
                                    ine.getKNNs(graph,kValues[k],*queryNodeIt,ineKNNs,ineKNNDistances);
                                    // Note: We are marking Distance Browsing kNN as unsorted until equal max priority problem in algorithm is resolved
                                    if (!utility::verifyKNN(ineKNNs,ineKNNDistances,kNNs,kNNDistances,false,kValues[k],message,true)) {
                                        std::cout << "Verfication failed for IER" << desc << " on object index " << objIdxFilePath << " for query node " << *queryNodeIt << " with k = " << kValues[k] << std::endl;
                                        std::cout << "Message: " << message << std::endl;
                                    }
                                }
                            }
                        }

                        double ierQueryTimeMs = totalQueryTime/totalQueries;
#if defined(COLLECT_STATISTICS)
                        knnStats.normalizeStatistics(totalQueries);
#endif

                        // Collect stats and return to output to file
                        kNNs.clear(); // Clear so that we don't pass last executed queries results to stats tuple
                        kNNDistances.clear();
                        KnnQueryTuple stats(graph.getNetworkName(),graph.getNumNodes(),graph.getNumEdges(),totalQueries,method,
                                            kValues[k],ierQueryTimeMs,objTypes[i],objDensities[j],objVariable[n],static_cast<int>(totalObjects/numSets),kNNs,kNNDistances);
                        stats.setAdditionalFields(specialFields);
                        stats.addSupplementaryFields("rtree_branch_factor",std::to_string(branchFactors[m]));
#if defined(COLLECT_STATISTICS)
                        knnStats.populateTupleFields(stats,0);
#endif
                        this->outputCommandStats(statsOutputFile,stats.getTupleString());
                    }
                }
            }
        }
    }
    std::cout << "IER" << desc << " kNN queries successfully executed for " << graph.getNetworkName() << std::endl;
    if (verifyKNN) {
        std::cout << "IER" << desc << " kNN verification completed for " << graph.getNetworkName() << std::endl;
    }
}

void ExperimentsCommand::runSingleMethodQueries(std::string method, std::string bgrFileName, 
                                                std::string queryNodeFile, std::string kValues, std::string parameters, 
                                                std::size_t numSets, std::string objDensities, std::string objTypes, std::string objVariable, 
                                                std::string filePathPrefix, std::string statsOutputFile)
{
    Graph graph = serialization::getIndexFromBinaryFile<Graph>(bgrFileName);
    std::unordered_map<std::string,std::string> parameterMap = this->getParameters(parameters);

    std::vector<NodeID> queryNodes = utility::getPointSetFromFile(queryNodeFile);

    bool verifykNN = parameterMap["verify"] == "1";
    
    // Find all additional fields we need to add to kNN stats tuples
    std::vector<std::string> specialFields;
    int field = 0;
    while (true) {
        std::string key = "special_field_" + std::to_string(field);
        if (parameterMap.find(key) != parameterMap.end()) {
            specialFields.push_back(parameterMap[key]);
            ++field;
        } else {
            break;
        }
    }
    
    std::vector<std::string> strObjDensitiesVec = utility::splitByDelim(objDensities,',');
    std::vector<double> objDensitiesVec;
    for(std::size_t i = 0; i < strObjDensitiesVec.size(); ++i) {
        double density = std::stod(strObjDensitiesVec[i]);
        if (density > 0) {
            objDensitiesVec.push_back(density);
        } else {
            std::cerr << "Invalid density in list provided!\n";
            exit(1);    
        }
    }
    std::vector<std::string> objTypesVec = utility::splitByDelim(objTypes,',');
    std::vector<std::string> strKValuesVec = utility::splitByDelim(kValues,',');
    std::vector<int> kValuesVec;
    for(std::size_t i = 0; i < strKValuesVec.size(); ++i) {
        int k = std::stoi(strKValuesVec[i]);
        if (k > 0) {
            kValuesVec.push_back(k);
        } else {
            std::cerr << "Invalid k value in list provided!\n";
            exit(1);    
        }
    }
    std::vector<std::string> strObjVariables = utility::splitByDelim(objVariable,',');
    std::vector<int> objVariableVec;
    for(std::size_t i = 0; i < strObjVariables.size(); ++i) {
        int variable = std::stoi(strObjVariables[i]);
        if (variable > 0) {
            objVariableVec.push_back(variable);
        } else {
            std::cerr << "Invalid variable in list provided (must be greater than zero)!\n";
            exit(1);    
        }
    }
    if (objDensitiesVec.size() == 0 || objTypesVec.size() == 0 || objVariableVec.size() == 0 || kValuesVec.size() == 0) {
        std::cerr << "Not enough densities or types provided!\n";
        exit(1);    
    }
    
    std::string gtreeIdxFile, roadIdxFile, silcIdxFile, juncIdxFile, phlIdxFile, altIdxFile;
    std::vector<std::string> parameterKeys, parameterValues;
    
    if (method == constants::INE_KNN_QUERY) {
        this->runINEQueries(graph,queryNodes,kValuesVec,numSets,objDensitiesVec,objTypesVec,objVariableVec,filePathPrefix,statsOutputFile,specialFields);
    } else if (method == "bad_ine") {
        // Note: This is to compare INE when non-ideal data structures are chosen
        std::string dynBgrFileName = filePathPrefix + "/indexes/" + graph.getNetworkName() + "_dynamic.bin";
        this->runINEQueriesByDynamicGraph(graph,dynBgrFileName,queryNodes,kValuesVec,numSets,objDensitiesVec,objTypesVec,objVariableVec,filePathPrefix,statsOutputFile,verifykNN,specialFields);
    } else if (method == constants::GTREE_KNN_QUERY) {
        int fanout = std::stoi(parameterMap["gtree_fanout"]);
        std::size_t maxLeafSize = std::stoi(parameterMap["gtree_maxleafsize"]);
        if (fanout < 2 || maxLeafSize < 32) {
            std::cerr << "Invalid Gtree parameters!\n";
            exit(1);    
        }
        parameterKeys.clear();
        parameterValues.clear();
        parameterKeys.push_back("fanout");
        parameterValues.push_back(parameterMap["gtree_fanout"]);
        parameterKeys.push_back("maxleafsize");
        parameterValues.push_back(parameterMap["gtree_maxleafsize"]);
        gtreeIdxFile = filePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_GTREE_CMD,parameterKeys,parameterValues);
        this->runGtreeQueries(graph,gtreeIdxFile,queryNodes,kValuesVec,numSets,objDensitiesVec,objTypesVec,objVariableVec,filePathPrefix,statsOutputFile,verifykNN,parameterKeys,parameterValues,specialFields);
    } else if (method == constants::ROAD_KNN_QUERY) {
        int fanout = std::stoi(parameterMap["road_fanout"]);
        std::size_t levels = std::stoi(parameterMap["road_levels"]);
        if (fanout < 2 || levels < 2) {
            std::cerr << "Invalid Route Overlay parameters!\n";
            exit(1);
        }
        parameterKeys.clear();
        parameterValues.clear();
        parameterKeys.push_back("fanout");
        parameterValues.push_back(parameterMap["road_fanout"]);
        parameterKeys.push_back("levels");
        parameterValues.push_back(parameterMap["road_levels"]);
        roadIdxFile = filePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_ROUTEOVERLAY_CMD,parameterKeys,parameterValues);        
        this->runROADQueries(graph,roadIdxFile,queryNodes,kValuesVec,numSets,objDensitiesVec,objTypesVec,objVariableVec,filePathPrefix,statsOutputFile,verifykNN,parameterKeys,parameterValues,specialFields);
    } else if (method == constants::SILC_KNN_QUERY) {
        std::vector<int> maxLeafSizes = utility::getIngetersFromStringList(parameterMap["silc_quadtree_maxleafsize"],':',0);
        if (maxLeafSizes.size() == 0) {
            std::cerr << "No object hierarchy quadtree max leaf sizes provided " << method << "!\n";
            exit(1);
        }
        parameterKeys.clear();
        parameterValues.clear();
        silcIdxFile = filePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_SILC_CMD,parameterKeys,parameterValues);
        juncIdxFile = filePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_JUNC_CMD,parameterKeys,parameterValues);
        std::vector<std::string> methods;
        methods.push_back(constants::SILC_KNN_QUERY);
        this->runSILCBasedQueries(graph,methods,silcIdxFile,juncIdxFile,maxLeafSizes,queryNodes,kValuesVec,numSets,objDensitiesVec,objTypesVec,objVariableVec,filePathPrefix,statsOutputFile,verifykNN,specialFields);
    } else if (method == constants::DISTBRWS_KNN_QUERY) {
        std::vector<int> maxLeafSizes = utility::getIngetersFromStringList(parameterMap["db_quadtree_maxleafsize"],':',0);
        if (maxLeafSizes.size() == 0) {
            std::cerr << "No object hierarchy quadtree max leaf sizes provided " << method << "!\n";
            exit(1);
        }
        parameterKeys.clear();
        parameterValues.clear();
        silcIdxFile = filePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_SILC_CMD,parameterKeys,parameterValues);
        juncIdxFile = filePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_JUNC_CMD,parameterKeys,parameterValues);
        std::vector<std::string> methods;
        methods.push_back(constants::DISTBRWS_KNN_QUERY);
        this->runSILCBasedQueries(graph,methods,silcIdxFile,juncIdxFile,maxLeafSizes,queryNodes,kValuesVec,numSets,objDensitiesVec,objTypesVec,objVariableVec,filePathPrefix,statsOutputFile,verifykNN,specialFields);
    } else if (method == constants::OPT_SILC_KNN_QUERY) {
        std::vector<int> maxLeafSizes = utility::getIngetersFromStringList(parameterMap["silc_quadtree_maxleafsize"],':',0);
        if (maxLeafSizes.size() == 0) {
            std::cerr << "No object hierarchy quadtree max leaf sizes provided " << method << "!\n";
            exit(1);
        }
        parameterKeys.clear();
        parameterValues.clear();
        silcIdxFile = filePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_SILC_CMD,parameterKeys,parameterValues);
        juncIdxFile = filePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_JUNC_CMD,parameterKeys,parameterValues);
        std::vector<std::string> methods;
        methods.push_back(constants::OPT_SILC_KNN_QUERY);
        this->runSILCBasedQueries(graph,methods,silcIdxFile,juncIdxFile,maxLeafSizes,queryNodes,kValuesVec,numSets,objDensitiesVec,objTypesVec,objVariableVec,filePathPrefix,statsOutputFile,verifykNN,specialFields);
    } else if (method == constants::OPT_DISTBRWS_KNN_QUERY) {
        std::vector<int> maxLeafSizes = utility::getIngetersFromStringList(parameterMap["db_quadtree_maxleafsize"],':',0);
        if (maxLeafSizes.size() == 0) {
            std::cerr << "No object hierarchy quadtree max leaf sizes provided for " << method << "!\n";
            exit(1);
        }
        parameterKeys.clear();
        parameterValues.clear();
        silcIdxFile = filePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_SILC_CMD,parameterKeys,parameterValues);
        juncIdxFile = filePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_JUNC_CMD,parameterKeys,parameterValues);
        std::vector<std::string> methods;
        methods.push_back(constants::OPT_DISTBRWS_KNN_QUERY);
        this->runSILCBasedQueries(graph,methods,silcIdxFile,juncIdxFile,maxLeafSizes,queryNodes,kValuesVec,numSets,objDensitiesVec,objTypesVec,objVariableVec,filePathPrefix,statsOutputFile,verifykNN,specialFields);
    } else if (method == constants::DB_RTREE_KNN_QUERY) {
        std::vector<int> branchFactors = utility::getIngetersFromStringList(parameterMap["rtree_branchfactor"],':',0);
        if (branchFactors.size() == 0) {
            std::cerr << "No R-tree branch factors provided " << method << "!\n";
            exit(1);
        }
        parameterKeys.clear();
        parameterValues.clear();
        silcIdxFile = filePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_SILC_CMD,parameterKeys,parameterValues);
        juncIdxFile = filePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_JUNC_CMD,parameterKeys,parameterValues);
        std::vector<std::string> methods;
        methods.push_back(constants::DB_RTREE_KNN_QUERY);
        this->runSILCBasedQueries(graph,methods,silcIdxFile,juncIdxFile,branchFactors,queryNodes,kValuesVec,numSets,objDensitiesVec,objTypesVec,objVariableVec,filePathPrefix,statsOutputFile,verifykNN,specialFields);
    } else if (method == constants::IER_DIJKSTRA_KNN_QUERY) {
        std::vector<int> branchFactors = utility::getIngetersFromStringList(parameterMap["rtree_branchfactor"],':',0);
        if (branchFactors.size() == 0) {
            std::cerr << "No R-tree branch factors provided " << method << "!\n";
            exit(1);
        }
        parameterKeys.clear();
        parameterValues.clear();
        this->runIERQueries(graph,constants::DIJKSTRA_SPDIST_QUERY,gtreeIdxFile,branchFactors,queryNodes,kValuesVec,numSets,objDensitiesVec,objTypesVec,objVariableVec,filePathPrefix,statsOutputFile,verifykNN,specialFields);
    } else if (method == constants::IER_GTREE_KNN_QUERY) {
        std::vector<int> branchFactors = utility::getIngetersFromStringList(parameterMap["rtree_branchfactor"],':',0);
        if (branchFactors.size() == 0) {
            std::cerr << "No R-tree branch factors provided " << method << "!\n";
            exit(1);
        }
        parameterKeys.clear();
        parameterValues.clear();
        parameterKeys.push_back("fanout");
        parameterValues.push_back(parameterMap["gtree_fanout"]);
        parameterKeys.push_back("maxleafsize");
        parameterValues.push_back(parameterMap["gtree_maxleafsize"]);
        gtreeIdxFile = filePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_GTREE_CMD,parameterKeys,parameterValues);
        this->runIERQueries(graph,constants::GTREE_SPDIST_QUERY,gtreeIdxFile,branchFactors,queryNodes,kValuesVec,numSets,objDensitiesVec,objTypesVec,objVariableVec,filePathPrefix,statsOutputFile,verifykNN,specialFields);
    } else if (method == constants::IER_SILC_KNN_QUERY) {
        std::vector<int> branchFactors = utility::getIngetersFromStringList(parameterMap["rtree_branchfactor"],':',0);
        if (branchFactors.size() == 0) {
            std::cerr << "No R-tree branch factors provided " << method << "!\n";
            exit(1);
        }
        parameterKeys.clear();
        parameterValues.clear();
        silcIdxFile = filePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_SILC_CMD,parameterKeys,parameterValues);
        if (utility::fileExists(silcIdxFile)) {
            this->runIERQueries(graph,constants::SILC_SPDIST_QUERY,silcIdxFile,branchFactors,queryNodes,kValuesVec,numSets,objDensitiesVec,objTypesVec,objVariableVec,filePathPrefix,statsOutputFile,verifykNN,specialFields);
        } else {
            std::cerr << "SILC index not built for " << graph.getNetworkName() << ", unable to execute IER-SILC\n";
        }
    } else if (method == constants::IER_PHL_KNN_QUERY) {
        std::vector<int> branchFactors = utility::getIngetersFromStringList(parameterMap["rtree_branchfactor"],':',0);
        if (branchFactors.size() == 0) {
            std::cerr << "No R-tree branch factors provided " << method << "!\n";
            exit(1);
        }
        parameterKeys.clear();
        parameterValues.clear();
        phlIdxFile = filePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_PHL_CMD,parameterKeys,parameterValues);
        if (utility::fileExists(phlIdxFile)) {
            this->runIERQueries(graph,constants::PHL_SPDIST_QUERY,phlIdxFile,branchFactors,queryNodes,kValuesVec,numSets,objDensitiesVec,objTypesVec,objVariableVec,filePathPrefix,statsOutputFile,verifykNN,specialFields);
        } else {
            std::cerr << "PHL index not built for " << graph.getNetworkName() << ", unable to execute IER-PHL\n";
        }
    } else if (method == constants::IER_TNR_KNN_QUERY) {
        std::vector<int> branchFactors = utility::getIngetersFromStringList(parameterMap["rtree_branchfactor"],':',0);
        if (branchFactors.size() == 0) {
            std::cerr << "No R-tree branch factors provided " << method << "!\n";
            exit(1);
        }
        int gridSize = std::stoi(parameterMap["tnr_gridsize"]);
        if (gridSize < 16) {
            std::cerr << "Invalid TNR grid size!\n";
            exit(1);
        }
        parameterKeys.clear();
        parameterValues.clear();
        std::string bgrIntFile = filePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),"ext_bin",parameterKeys,parameterValues);
        std::string bcoIntFile = filePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),"ext_co_bin",parameterKeys,parameterValues);
        std::string chIdxFile = filePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_CH_CMD,parameterKeys,parameterValues);
        parameterKeys.push_back("gridsize");
        parameterValues.push_back(parameterMap["tnr_gridsize"]);
        std::string tnrIdxFile = filePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_TNR_CMD,parameterKeys,parameterValues);        
        ShortestPathWrapper spw;
        spw.runIERQueries(graph,method,bgrIntFile,bcoIntFile,chIdxFile,tnrIdxFile,branchFactors,queryNodes,kValuesVec,numSets,objDensitiesVec,objTypesVec,objVariableVec,filePathPrefix,statsOutputFile,verifykNN,specialFields);
    } else if (method == constants::IER_CH_KNN_QUERY) {
        std::vector<int> branchFactors = utility::getIngetersFromStringList(parameterMap["rtree_branchfactor"],':',0);
        if (branchFactors.size() == 0) {
            std::cerr << "No R-tree branch factors provided " << method << "!\n";
            exit(1);
        }
        parameterKeys.clear();
        parameterValues.clear();
        std::string bgrIntFile = filePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),"ext_bin",parameterKeys,parameterValues);
        std::string bcoIntFile = filePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),"ext_co_bin",parameterKeys,parameterValues);
        std::string chIdxFile = filePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_CH_CMD,parameterKeys,parameterValues);
        std::string tnrIdxFile = "";     
        ShortestPathWrapper spw;
        spw.runIERQueries(graph,method,bgrIntFile,bcoIntFile,chIdxFile,tnrIdxFile,branchFactors,queryNodes,kValuesVec,numSets,objDensitiesVec,objTypesVec,objVariableVec,filePathPrefix,statsOutputFile,verifykNN,specialFields);
    } else {
        std::cerr << "Could not recognise method - check kNN query command" << std::endl;
        exit(1);
    }

}

void ExperimentsCommand::runINEQueriesByDynamicGraph(Graph& graph, std::string dynBgrFileName, std::vector<NodeID>& queryNodes, std::vector<int>& kValues, std::size_t numSets, 
                                                     std::vector<double> objDensities, std::vector<std::string> objTypes, std::vector<int> objVariable, std::string filePathPrefix, 
                                                     std::string statsOutputFile, bool verifyKNN, std::vector<std::string> specialFields)
{
    std::vector<NodeID> kNNs, ineKNNs;
    std::vector<EdgeWeight> kNNDistances, ineKNNDistances;

    std::string objSetType;
    double objSetDensity;
    int objSetSize, setVariable;
    
    StopWatch sw;
    double totalQueryTime;
    int totalQueries = numSets*queryNodes.size();
    int totalObjects;
    
    INE ine;
    std::string message;
    Statistics knnStats;
    
    DynamicGraph dynGraph = serialization::getIndexFromBinaryFile<DynamicGraph>(dynBgrFileName);
    
    for (std::size_t i = 0; i < objTypes.size(); ++i) {
        for (std::size_t j = 0; j < objDensities.size(); ++j) {
            for (std::size_t m = 0; m < objVariable.size(); ++m) {
                for (std::size_t k = 0; k < kValues.size(); ++k) {
                    totalQueryTime = 0;
                    totalObjects = 0;
#if defined(COLLECT_STATISTICS)
                    knnStats.clear();
#endif
                    for (std::size_t l = 0; l < numSets; ++l) {
                        std::string objSetFile = filePathPrefix + "/obj_indexes/" + utility::constructObjsectSetFileName(graph.getNetworkName(),objTypes[i],objDensities[j],objVariable[m],l);
                        std::vector<NodeID> objectNodes = utility::getPointSetFromFile(objSetFile,objSetType,objSetDensity,objSetSize,setVariable);
                        totalObjects += objSetSize;
            
                        dynGraph.resetAllObjects();
                        dynGraph.parseObjectSet(objectNodes);
            
                        if (verifyKNN) {
                            graph.resetAllObjects();
                            graph.parseObjectSet(objectNodes);
                        }
                        
                        for (auto queryNodeIt = queryNodes.begin(); queryNodeIt != queryNodes.end(); ++queryNodeIt) {
                            kNNs.clear();
                            kNNDistances.clear();
                            kNNs.reserve(kValues[k]);
                            kNNDistances.reserve(kValues[k]);
                            sw.reset();
                            sw.start();
                            ine.getKNNsByDynamicGraph(dynGraph,kValues[k],*queryNodeIt,kNNs,kNNDistances);
                            sw.stop();
                            totalQueryTime += sw.getTimeUs();
#if defined(COLLECT_STATISTICS)
                            knnStats.mergeStatistics(ine.stats);
#endif
                            if (verifyKNN) {
                                ineKNNs.clear();
                                ineKNNDistances.clear();
                                ine.getKNNs(graph,kValues[k],*queryNodeIt,ineKNNs,ineKNNDistances);
                                if (!utility::verifyKNN(ineKNNs,ineKNNDistances,kNNs,kNNDistances,false,kValues[k],message,true)) {
                                    std::cout << "Verfication failed for INE on DynamicGraph on object set file " << objSetFile << " for query node " << *queryNodeIt << " with k = " << kValues[k] << std::endl;
                                    std::cout << "Message: " << message << std::endl;
                                }
                            }
                        }        
                    }

                    double queryTimeMs = totalQueryTime/totalQueries;
#if defined(COLLECT_STATISTICS)
                    knnStats.normalizeStatistics(totalQueries);
#endif

                    // Collect stats and return to output to file
                    kNNs.clear(); // Clear so that we don't pass last executed queries results to stats tuple
                    kNNDistances.clear();
                    KnnQueryTuple stats(dynGraph.getNetworkName(),dynGraph.getNumNodes(),dynGraph.getNumEdges(),0,constants::INE_KNN_QUERY,
                                        kValues[k],queryTimeMs,objTypes[i],objDensities[j],objVariable[m],static_cast<int>(totalObjects/numSets),kNNs,kNNDistances);
                    stats.setAdditionalFields(specialFields);
#if defined(COLLECT_STATISTICS)
                    knnStats.populateTupleFields(stats,0);
#endif
                    this->outputCommandStats(statsOutputFile,stats.getTupleString());
                }
            }
        }
    }
    std::cout << "INE kNN on DynamicGraph graph queries successfully executed for " << dynGraph.getNetworkName() << std::endl;
    if (verifyKNN) {
        std::cout << "INE kNN on DynamicGraph graph verification completed for " << dynGraph.getNetworkName() << std::endl;
    }
}

void ExperimentsCommand::runRealWorldPOIQueries(std::string bgrFileName, std::string queryNodeFile, std::string kValues, 
                                                std::string parameters, std::string filePathPrefix, std::string statsOutputFile, 
                                                std::string rwPOISetListFile)
{
    Graph graph = serialization::getIndexFromBinaryFile<Graph>(bgrFileName);
    std::unordered_map<std::string,std::string> parameterMap = this->getParameters(parameters);

    std::vector<NodeID> queryNodes = utility::getPointSetFromFile(queryNodeFile);
    if (queryNodes.size() == 0) {
        std::cerr << "No query points were provided!\n";
        exit(1);
    }

    std::vector<std::string> strKValuesVec = utility::splitByDelim(kValues,',');
    std::vector<int> kValuesVec;
    for(std::size_t i = 0; i < strKValuesVec.size(); ++i) {
        int k = std::stoi(strKValuesVec[i]);
        if (k > 0) {
            kValuesVec.push_back(k);
        } else {
            std::cerr << "Invalid k value in list provided!\n";
            exit(1);
        }
    }
    
    // Read POI set list file to extract file name for real-world datasets
    std::string fileName, setName;
    std::vector<std::string> fileNames, setNames;
    std::ifstream rwSetListFS(rwPOISetListFile, std::ios::in);
    if (rwSetListFS.is_open()) {
        while (rwSetListFS >> fileName >> setName)  {
            fileNames.push_back(fileName);
            setNames.push_back(setName);
        }
    } else {
        std::cerr << "Cannot open real world set list file " << rwPOISetListFile << std::endl;
        exit(1);
    }
    
    bool queryINE = parameterMap["ine"] == "1";
    bool queryGtree = parameterMap["gtree"] == "1";
    bool queryRoad = parameterMap["road"] == "1";
    bool queryIER = parameterMap["ier"] == "1";
    bool querySILC = parameterMap["silc"] == "1";
    bool queryDistBrws = parameterMap["dist_brws"] == "1";
    bool queryIERPHL = parameterMap["ier_phl"] == "1";
    
    std::vector<NodeID> kNNs;
    std::vector<EdgeWeight> kNNDistances;

    StopWatch sw;
    double totalQueryTime = 0;
    int totalQueries = queryNodes.size();
    std::string setType;
    double setDensity;
    int setSize, setVariable;
    
    if (queryINE) {
        INE ine;
        for (std::size_t i = 0; i < kValuesVec.size(); ++i) {
            kNNs.reserve(kValuesVec[i]);
            kNNDistances.reserve(kValuesVec[i]);
            for (std::size_t j = 0; j < fileNames.size(); ++j) {
                std::string objSetFile = filePathPrefix + "/real_world_pois/" + graph.getNetworkName() + "/" +fileNames[j] + ".txt";
                std::vector<NodeID> objectNodes = utility::getPointSetFromFile(objSetFile,setType,setDensity,setSize,setVariable);
                
                graph.resetAllObjects();
                graph.parseObjectSet(objectNodes);
                
                sw.reset();
                sw.start();
                for (auto queryNodeIt = queryNodes.begin(); queryNodeIt != queryNodes.end(); ++queryNodeIt) {
                    kNNs.clear();
                    kNNDistances.clear();
                    ine.getKNNs(graph,kValuesVec[i],*queryNodeIt,kNNs,kNNDistances);
                }
                sw.stop();
                totalQueryTime = sw.getTimeUs();
                double queryTimeMs = totalQueryTime/totalQueries;
                
                // Collect stats and return to output to file
                kNNs.clear(); // Clear so that we don't pass last executed queries results to stats tuple
                kNNDistances.clear();
                KnnQueryTuple stats(graph.getNetworkName(),graph.getNumNodes(),graph.getNumEdges(),totalQueries,constants::INE_KNN_QUERY,
                                    kValuesVec[i],queryTimeMs,setNames[j],setDensity,setVariable,setSize,kNNs,kNNDistances);
                this->outputCommandStats(statsOutputFile,stats.getTupleString());
            }
        }
    }
    
    std::string gtreeIdxFile, roadIdxFile, silcIdxFile, juncIdxFile, phlIdxFile;
    std::vector<std::string> parameterKeys, parameterValues;
    
    if (queryRoad) {
        int fanout = std::stoi(parameterMap["road_fanout"]);
        std::size_t levels = std::stoi(parameterMap["road_levels"]);
        if (fanout < 2 || levels < 2) {
            std::cerr << "Invalid Route Overlay parameters!\n";
            exit(1);
        }
        parameterKeys.clear();
        parameterValues.clear();
        parameterKeys.push_back("fanout");
        parameterValues.push_back(parameterMap["road_fanout"]);
        parameterKeys.push_back("levels");
        parameterValues.push_back(parameterMap["road_levels"]);
        roadIdxFile = filePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_ROUTEOVERLAY_CMD,parameterKeys,parameterValues);
        ROAD road = serialization::getIndexFromBinaryFile<ROAD>(roadIdxFile);
        for (std::size_t i = 0; i < kValuesVec.size(); ++i) {
            kNNs.reserve(kValuesVec[i]);
            kNNDistances.reserve(kValuesVec[i]);
            for (std::size_t j = 0; j < fileNames.size(); ++j) {
                std::string objIdxFilePath = filePathPrefix + "/real_world_pois/" + graph.getNetworkName() + "/" +fileNames[j] + "." + constants::OBJ_IDX_ROAD;
                AssociationDirectory assocDir = serialization::getIndexFromBinaryFile<AssociationDirectory>(objIdxFilePath);
                
                sw.reset();
                sw.start();
                for (auto queryNodeIt = queryNodes.begin(); queryNodeIt != queryNodes.end(); ++queryNodeIt) {
                    kNNs.clear();
                    kNNDistances.clear();
                    road.getKNNs(assocDir,kValuesVec[i],*queryNodeIt,kNNs,kNNDistances);
                }
                sw.stop();
                totalQueryTime = sw.getTimeUs();
                double queryTimeMs = totalQueryTime/totalQueries;
                
                // Collect stats and return to output to file
                kNNs.clear(); // Clear so that we don't pass last executed queries results to stats tuple
                kNNDistances.clear();
                KnnQueryTuple stats(graph.getNetworkName(),graph.getNumNodes(),graph.getNumEdges(),totalQueries,constants::ROAD_KNN_QUERY,
                                    kValuesVec[i],queryTimeMs,setNames[j],assocDir.getObjSetDensity(),assocDir.getObjSetVariable(),assocDir.getObjSetSize(),kNNs,kNNDistances);
                this->outputCommandStats(statsOutputFile,stats.getTupleString());
            }
        }
    }

    if (queryGtree) {
        int fanout = std::stoi(parameterMap["gtree_fanout"]);
        std::size_t maxLeafSize = std::stoi(parameterMap["gtree_maxleafsize"]);
        if (fanout < 2 || maxLeafSize < 32) {
            std::cerr << "Invalid Gtree parameters!\n";
            exit(1);
        }
        parameterKeys.clear();
        parameterValues.clear();
        parameterKeys.push_back("fanout");
        parameterValues.push_back(parameterMap["gtree_fanout"]);
        parameterKeys.push_back("maxleafsize");
        parameterValues.push_back(parameterMap["gtree_maxleafsize"]);
        gtreeIdxFile = filePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_GTREE_CMD,parameterKeys,parameterValues);
        Gtree gtree = serialization::getIndexFromBinaryFile<Gtree>(gtreeIdxFile);
        for (std::size_t i = 0; i < kValuesVec.size(); ++i) {
            kNNs.reserve(kValuesVec[i]);
            kNNDistances.reserve(kValuesVec[i]);
            for (std::size_t j = 0; j < fileNames.size(); ++j) {
                std::string objIdxFilePath = filePathPrefix + "/real_world_pois/" + graph.getNetworkName() + "/" +fileNames[j] + "." + constants::OBJ_IDX_GTREE;
                OccurenceList occList = serialization::getIndexFromBinaryFile<OccurenceList>(objIdxFilePath);
                
                sw.reset();
                sw.start();
                for (auto queryNodeIt = queryNodes.begin(); queryNodeIt != queryNodes.end(); ++queryNodeIt) {
                    kNNs.clear();
                    kNNDistances.clear();
                    gtree.getKNNs(occList,kValuesVec[i],*queryNodeIt,kNNs,kNNDistances,graph);
                }
                sw.stop();
                totalQueryTime = sw.getTimeUs();
                double queryTimeMs = totalQueryTime/totalQueries;
                
                // Collect stats and return to output to file
                kNNs.clear(); // Clear so that we don't pass last executed queries results to stats tuple
                kNNDistances.clear();
                KnnQueryTuple stats(graph.getNetworkName(),graph.getNumNodes(),graph.getNumEdges(),totalQueries,constants::GTREE_KNN_QUERY,
                                    kValuesVec[i],queryTimeMs,setNames[j],occList.getObjSetDensity(),occList.getObjSetVariable(),occList.getObjSetSize(),kNNs,kNNDistances);
                this->outputCommandStats(statsOutputFile,stats.getTupleString());
            }
        }
    }
    
    std::string otherIdxFile;

    if (queryIER) {
        std::vector<int> branchFactors = utility::getIngetersFromStringList(parameterMap["rtree_branchfactor"],':',0);
        if (branchFactors.size() == 0) {
            std::cerr << "No R-tree branch factors provided!\n";
            exit(1);
        }
        parameterKeys.clear();
        parameterValues.clear();
        parameterKeys.push_back("fanout");
        parameterValues.push_back(parameterMap["gtree_fanout"]);
        parameterKeys.push_back("maxleafsize");
        parameterValues.push_back(parameterMap["gtree_maxleafsize"]);
        gtreeIdxFile = filePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_GTREE_CMD,parameterKeys,parameterValues);
        IER ier;
        Gtree gtree = serialization::getIndexFromBinaryFile<Gtree>(gtreeIdxFile);
        for (std::size_t i = 0; i < kValuesVec.size(); ++i) {
            kNNs.reserve(kValuesVec[i]);
            kNNDistances.reserve(kValuesVec[i]);
            for (std::size_t j = 0; j < fileNames.size(); ++j) {
                std::string objIdxFilePath = filePathPrefix + "/real_world_pois/" + graph.getNetworkName() + "/" +fileNames[j] + "." + constants::OBJ_IDX_RTREE;
                StaticRtree rtree = serialization::getIndexFromBinaryFile<StaticRtree>(objIdxFilePath);
                
                if (graph.getEdgeType() == constants::DISTANCE_WEIGHTS) {
                    sw.reset();
                    sw.start();
                    for (auto queryNodeIt = queryNodes.begin(); queryNodeIt != queryNodes.end(); ++queryNodeIt) {
                        kNNs.clear();
                        kNNDistances.clear();
                        ier.getKNNsByGtree(rtree,kValuesVec[i],*queryNodeIt,kNNs,kNNDistances,graph,gtree);
                    }
                    sw.stop();
                    totalQueryTime = sw.getTimeUs();
                } else {
                    sw.reset();
                    sw.start();
                    for (auto queryNodeIt = queryNodes.begin(); queryNodeIt != queryNodes.end(); ++queryNodeIt) {
                        kNNs.clear();
                        kNNDistances.clear();
                        ier.getKNNsByGtreeTravelTimes(rtree,kValuesVec[i],*queryNodeIt,kNNs,kNNDistances,graph,gtree);
                    }
                    sw.stop();
                    totalQueryTime = sw.getTimeUs();
                }
                
                double queryTimeMs = totalQueryTime/totalQueries;
                
                // Collect stats and return to output to file
                kNNs.clear(); // Clear so that we don't pass last executed queries results to stats tuple
                kNNDistances.clear();
                KnnQueryTuple stats(graph.getNetworkName(),graph.getNumNodes(),graph.getNumEdges(),totalQueries,constants::IER_GTREE_KNN_QUERY,
                                    kValuesVec[i],queryTimeMs,setNames[j],rtree.getObjSetDensity(),rtree.getObjSetVariable(),rtree.getObjSetSize(),kNNs,kNNDistances);
                stats.addSupplementaryFields("rtree_branch_factor",std::to_string(rtree.getBranchFactor()));
                this->outputCommandStats(statsOutputFile,stats.getTupleString());
            }
        }
    }

    if (queryIERPHL) {
        std::vector<int> branchFactors = utility::getIngetersFromStringList(parameterMap["rtree_branchfactor"],':',0);
        if (branchFactors.size() == 0) {
            std::cerr << "No R-tree branch factors provided!\n";
            exit(1);
        }
        parameterKeys.clear();
        parameterValues.clear();
        phlIdxFile = filePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_PHL_CMD,parameterKeys,parameterValues);
        if (utility::fileExists(phlIdxFile)) {
            IER ier;
            PrunedHighwayLabeling phl;
            phl.LoadLabel(phlIdxFile.c_str());
            for (std::size_t i = 0; i < kValuesVec.size(); ++i) {
                kNNs.reserve(kValuesVec[i]);
                kNNDistances.reserve(kValuesVec[i]);
                for (std::size_t j = 0; j < fileNames.size(); ++j) {
                    std::string objIdxFilePath = filePathPrefix + "/real_world_pois/" + graph.getNetworkName() + "/" +fileNames[j] + "." + constants::OBJ_IDX_RTREE;
                    StaticRtree rtree = serialization::getIndexFromBinaryFile<StaticRtree>(objIdxFilePath);
                    
                    if (graph.getEdgeType() == constants::DISTANCE_WEIGHTS) {
                        sw.reset();
                        sw.start();
                        for (auto queryNodeIt = queryNodes.begin(); queryNodeIt != queryNodes.end(); ++queryNodeIt) {
                            kNNs.clear();
                            kNNDistances.clear();
                            ier.getKNNsByPHL(rtree,kValuesVec[i],*queryNodeIt,kNNs,kNNDistances,graph,phl);
                        }
                        sw.stop();
                        totalQueryTime = sw.getTimeUs();
                    } else {
                        sw.reset();
                        sw.start();
                        for (auto queryNodeIt = queryNodes.begin(); queryNodeIt != queryNodes.end(); ++queryNodeIt) {
                            kNNs.clear();
                            kNNDistances.clear();
                            ier.getKNNsByPHLTravelTimes(rtree,kValuesVec[i],*queryNodeIt,kNNs,kNNDistances,graph,phl);
                        }
                        sw.stop();
                        totalQueryTime = sw.getTimeUs();
                    }
                    
                    double queryTimeMs = totalQueryTime/totalQueries;
                    
                    // Collect stats and return to output to file
                    kNNs.clear(); // Clear so that we don't pass last executed queries results to stats tuple
                    kNNDistances.clear();
                    KnnQueryTuple stats(graph.getNetworkName(),graph.getNumNodes(),graph.getNumEdges(),totalQueries,constants::IER_PHL_KNN_QUERY,
                                        kValuesVec[i],queryTimeMs,setNames[j],rtree.getObjSetDensity(),rtree.getObjSetVariable(),rtree.getObjSetSize(),kNNs,kNNDistances);
                    stats.addSupplementaryFields("rtree_branch_factor",std::to_string(rtree.getBranchFactor()));
                    this->outputCommandStats(statsOutputFile,stats.getTupleString());
                }
            }
        }
    }
    
    if (queryDistBrws) {
        std::vector<int> branchFactors = utility::getIngetersFromStringList(parameterMap["rtree_branchfactor"],':',0);
        if (branchFactors.size() == 0) {
            std::cerr << "No R-tree branch factors provided!\n";
            exit(1);
        }
        parameterKeys.clear();
        parameterValues.clear();
        silcIdxFile = filePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_SILC_CMD,parameterKeys,parameterValues);
        juncIdxFile = filePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_JUNC_CMD,parameterKeys,parameterValues);
        if (utility::fileExists(silcIdxFile) && utility::fileExists(juncIdxFile) && graph.getEdgeType() == constants::DISTANCE_WEIGHTS) {
            // It's possible SILC file does not exist for larger road networks
            // so we need to check in this case (but junc index file will exist)
            // Note: SILC-based kNN method do not support travel-time graphs (unlike SILC shortest path)
            SILCPathOracle silc = serialization::getIndexFromBinaryFile<SILCPathOracle>(silcIdxFile);
            Junction junc = serialization::getIndexFromBinaryFile<Junction>(juncIdxFile);
            for (std::size_t i = 0; i < kValuesVec.size(); ++i) {
                kNNs.reserve(kValuesVec[i]);
                kNNDistances.reserve(kValuesVec[i]);
                for (std::size_t j = 0; j < fileNames.size(); ++j) {
                    // Note: We don't bother setting branch factor in file name for real-world POIs
                    std::string objIdxFilePath = filePathPrefix + "/real_world_pois/" + graph.getNetworkName() + "/" +fileNames[j] + "." + constants::OBJ_IDX_RTREE;
                    StaticRtree rtree = serialization::getIndexFromBinaryFile<StaticRtree>(objIdxFilePath);
                    
                    sw.reset();
                    sw.start();
                    for (auto queryNodeIt = queryNodes.begin(); queryNodeIt != queryNodes.end(); ++queryNodeIt) {
                        kNNs.clear();
                        kNNDistances.clear();
                        silc.getKNNsByDistanceBrowsingViaRtree(rtree,kValuesVec[i],*queryNodeIt,kNNs,kNNDistances,graph,junc);
                    }
                    sw.stop();
                    totalQueryTime = sw.getTimeUs();
                    
                    double queryTimeMs = totalQueryTime/totalQueries;
                    
                    // Collect stats and return to output to file
                    kNNs.clear(); // Clear so that we don't pass last executed queries results to stats tuple
                    kNNDistances.clear();
                    KnnQueryTuple stats(graph.getNetworkName(),graph.getNumNodes(),graph.getNumEdges(),totalQueries,constants::DB_RTREE_KNN_QUERY,
                                        kValuesVec[i],queryTimeMs,setNames[j],rtree.getObjSetDensity(),rtree.getObjSetVariable(),rtree.getObjSetSize(),kNNs,kNNDistances);
                    stats.addSupplementaryFields("obj_rt_branch_factor",std::to_string(rtree.getBranchFactor()));
                    this->outputCommandStats(statsOutputFile,stats.getTupleString());
                }
            }
        }
    }
    
    std::cout << "Real-world POI set query testing complete!" << std::endl;

}

void ExperimentsCommand::buildRealWorldObjIndexes(std::string bgrFileName, std::string parameters, 
                                                  std::string filePathPrefix, std::string statsOutputFile, 
                                                  std::string rwPOISetListFile)
{
    Graph graph = serialization::getIndexFromBinaryFile<Graph>(bgrFileName);
    std::unordered_map<std::string,std::string> parameterMap = this->getParameters(parameters);

    bool buildOccList = parameterMap["occ_list"] == "1";
    bool buildAssocDir = parameterMap["assoc_dir"] == "1";
    bool buildQuadtree = parameterMap["quadtree"] == "1";
    bool buildRtree = parameterMap["rtree"] == "1";

    // Read POI set list file to extract file name for real-world datasets
    std::string fileName, setName;
    std::vector<std::string> fileNames, setNames;
    std::ifstream rwSetListFS(rwPOISetListFile, std::ios::in);
    if (rwSetListFS.is_open()) {
        while (rwSetListFS >> fileName >> setName)  {
            fileNames.push_back(fileName);
            setNames.push_back(setName);
        }
    } else {
        std::cerr << "Cannot open real world set list file " << rwPOISetListFile << std::endl;
        exit(1);
    }
    
    std::string gtreeIdxFile, roadIdxFile;
    std::vector<std::string> parameterKeys, parameterValues;
    StopWatch sw;
    double totalTime = 0, totalMemory = 0;
    int setSize, setVariable, totalSetSize;
    std::string setType;
    double setDensity;
    
    if (buildOccList) {
        int fanout = std::stoi(parameterMap["gtree_fanout"]);
        std::size_t maxLeafSize = std::stoi(parameterMap["gtree_maxleafsize"]);
        if (fanout < 2 || maxLeafSize < 32) {
            std::cerr << "Invalid Gtree parameters!\n";
            exit(1);
        }
        parameterKeys.clear();
        parameterValues.clear();
        parameterKeys.push_back("fanout");
        parameterValues.push_back(parameterMap["gtree_fanout"]);
        parameterKeys.push_back("maxleafsize");
        parameterValues.push_back(parameterMap["gtree_maxleafsize"]);
        gtreeIdxFile = filePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_GTREE_CMD,parameterKeys,parameterValues);
        Gtree gtree = serialization::getIndexFromBinaryFile<Gtree>(gtreeIdxFile);
        totalTime = 0;
        totalMemory = 0;
        totalSetSize = 0;
        for (std::size_t j = 0; j < fileNames.size(); ++j) {
            std::string objSetFile = filePathPrefix + "/real_world_pois/" + graph.getNetworkName() + "/" +fileNames[j] + ".txt";
            std::vector<NodeID> objectNodes = utility::getPointSetFromFile(objSetFile,setType,setDensity,setSize,setVariable);
            
            sw.start();
            OccurenceList occList(setType,setDensity,setVariable,setSize);
            for (auto objIt = objectNodes.begin(); objIt != objectNodes.end(); ++objIt) {
                // Find the Gtree leaf index for this object and add it to occurence list
                // for that leaf (create list if it doesn't exist)
                int leafIdx = gtree.getLeafIndex(*objIt);
                occList.addLeafOccurence(leafIdx,*objIt);
                
                // Propagate this to occurence lists of parents of leaf node
                int childIdx = leafIdx;
                int parentIdx = gtree.getParentIndex(leafIdx);
                while (parentIdx != -1) {
                    occList.addParentOccurence(parentIdx, childIdx);
                    
                    // Go up a level (until we reach root)
                    childIdx = parentIdx;
                    parentIdx = gtree.getParentIndex(childIdx);
                }
            }
            sw.stop();
            totalTime += sw.getTimeUs();
            totalMemory += occList.computeIndexSize(); // in MB
            totalSetSize += occList.getObjSetSize();
            
            // Serialize object index to output file
            std::string objIdxOutputFile = filePathPrefix + "/real_world_pois/" + graph.getNetworkName() + "/" + fileNames[j] + "." + constants::OBJ_IDX_GTREE;
            serialization::outputIndexToBinaryFile<OccurenceList>(occList,objIdxOutputFile);
            
            // Collect stats and return to output to file
            //ObjectIndexTuple stats(gtree.getNetworkName(),gtree.getNumNodes(),gtree.getNumEdges(),
            //                    constants::OBJ_IDX_GTREE,totalTime,totalMemory,setNames[j],setDensity,setVariable,setSize);

            //this->outputCommandStats(statsOutputFile,stats.getTupleString());
        }
        // Collect stats and return to output to file
        ObjectIndexTuple stats(gtree.getNetworkName(),gtree.getNumNodes(),gtree.getNumEdges(),
                            constants::OBJ_IDX_GTREE,totalTime,totalMemory,"all_rw_sets",setDensity,setVariable,totalSetSize);

        this->outputCommandStats(statsOutputFile,stats.getTupleString());
    }
    
    if (buildAssocDir) {
        int fanout = std::stoi(parameterMap["road_fanout"]);
        std::size_t levels = std::stoi(parameterMap["road_levels"]);
        if (fanout < 2 || levels < 2) {
            std::cerr << "Invalid Route Overlay parameters!\n";
            exit(1);
        }
        parameterKeys.clear();
        parameterValues.clear();
        parameterKeys.push_back("fanout");
        parameterValues.push_back(parameterMap["road_fanout"]);
        parameterKeys.push_back("levels");
        parameterValues.push_back(parameterMap["road_levels"]);
        roadIdxFile = filePathPrefix + "/indexes/" + utility::constructIndexFileName(graph.getNetworkName(),constants::IDX_ROUTEOVERLAY_CMD,parameterKeys,parameterValues);
        ROAD road = serialization::getIndexFromBinaryFile<ROAD>(roadIdxFile);
        totalTime = 0;
        totalMemory = 0;
        totalSetSize = 0;
        for (std::size_t j = 0; j < fileNames.size(); ++j) {
            std::string objSetFile = filePathPrefix + "/real_world_pois/" + graph.getNetworkName() + "/" +fileNames[j] + ".txt";
            std::vector<NodeID> objectNodes = utility::getPointSetFromFile(objSetFile,setType,setDensity,setSize,setVariable);
            
            sw.start();
            AssociationDirectory assocDir(setType,setDensity,setVariable,setSize,road.getRnetTreeSize());
            for (auto objIt = objectNodes.begin(); objIt != objectNodes.end(); ++objIt) {
                // Find the ROAD leaf Rnet for this object and add entry into association
                // directory entry for that leaf then propagant up Rnet hierarchy
                assocDir.addObject(*objIt);
                
                for (int rnetIdx: road.routeOverlay[*objIt].leafIdxs) {
                    assocDir.addRnetAssociation(rnetIdx);
                
                    // Propagate this to occurence lists of parents of leaf node
                    int parentRnetIdx = road.getParentRnet(rnetIdx);
                    while (parentRnetIdx != -1) {
                        assocDir.addRnetAssociation(parentRnetIdx);
                        
                        // Go up a level (until we reach root)
                        parentRnetIdx = road.getParentRnet(parentRnetIdx);
                    }
                }
            }
            sw.stop();
            totalTime += sw.getTimeUs();
            totalMemory += assocDir.computeIndexSize(); // in MB
            totalSetSize += assocDir.getObjSetSize();
            
            // Serialize object index to output file
            std::string objIdxOutputFile = filePathPrefix + "/real_world_pois/" + graph.getNetworkName() + "/" + fileNames[j] + "." + constants::OBJ_IDX_ROAD;
            serialization::outputIndexToBinaryFile<AssociationDirectory>(assocDir,objIdxOutputFile);
            
            // Collect stats and return to output to file
            //ObjectIndexTuple stats(road.getNetworkName(),road.getNumNodes(),road.getNumEdges(),
            //                    constants::OBJ_IDX_ROAD,totalTime,totalMemory,setNames[j],setDensity,setVariable,setSize);

            //this->outputCommandStats(statsOutputFile,stats.getTupleString());
        }
        // Collect stats and return to output to file
        ObjectIndexTuple stats(road.getNetworkName(),road.getNumNodes(),road.getNumEdges(),
                            constants::OBJ_IDX_ROAD,totalTime,totalMemory,"all_rw_sets",setDensity,setVariable,totalSetSize);

        this->outputCommandStats(statsOutputFile,stats.getTupleString());
    }
    
    if (buildRtree) {
        std::vector<int> branchFactors = utility::getIngetersFromStringList(parameterMap["rtree_branchfactor"],':',0);
        if (branchFactors.size() == 0) {
            std::cerr << "No R-tree branch factors provided!\n";
            exit(1);
        }
        for (std::size_t l = 0; l < branchFactors.size(); ++l) {
            parameterKeys.clear();
            parameterValues.clear();
            parameterKeys.push_back("branchfactor");
            parameterValues.push_back(std::to_string(branchFactors[l]));
            totalTime = 0;
            totalMemory = 0;
            totalSetSize = 0;
            for (std::size_t j = 0; j < fileNames.size(); ++j) {
                // Note: We don't bother setting branch factor in file name for real-world POIs
                std::string objSetFile = filePathPrefix + "/real_world_pois/" + graph.getNetworkName() + "/" +fileNames[j] + ".txt";
                std::vector<NodeID> objectNodes = utility::getPointSetFromFile(objSetFile,setType,setDensity,setSize,setVariable);
                
                sw.start();
                std::vector<CoordinatePair> objectCoords;
                for (std::size_t i = 0; i < objectNodes.size(); ++i) {
                    CoordinatePair objectCoordPair;
                    graph.getCoordinates(objectNodes[i],objectCoordPair.first,objectCoordPair.second);
                    objectCoords.push_back(objectCoordPair);
                }
                StaticRtree rtree(branchFactors[l],setType,setDensity,setVariable,setSize);
                rtree.bulkLoad(objectNodes,objectCoords);
                sw.stop();
                totalTime += sw.getTimeUs();
                totalMemory += rtree.computeIndexSize(); // in MB
                totalSetSize += rtree.getObjSetSize();

                // Serialize object index to output file
                std::string objIdxOutputFile = filePathPrefix + "/real_world_pois/" + graph.getNetworkName() + "/" + fileNames[j] + "." + constants::OBJ_IDX_RTREE;
                serialization::outputIndexToBinaryFile<StaticRtree>(rtree,objIdxOutputFile);
                
                // Collect stats and return to output to file
                //ObjectIndexTuple stats(graph.getNetworkName(),graph.getNumNodes(),graph.getNumEdges(),
                //                    constants::OBJ_IDX_RTREE,totalTime,totalMemory,setNames[j],setDensity,setVariable,setSize);
                //stats.addSupplementaryFields("rtree_branch_factor",std::to_string(branchFactors[l]));

                //this->outputCommandStats(statsOutputFile,stats.getTupleString());
            }
            // Collect stats and return to output to file
            ObjectIndexTuple stats(graph.getNetworkName(),graph.getNumNodes(),graph.getNumEdges(),
                                constants::OBJ_IDX_RTREE,totalTime,totalMemory,"all_rw_sets",setDensity,setVariable,totalSetSize);
            stats.addSupplementaryFields("rtree_branch_factor",std::to_string(branchFactors[l]));

            this->outputCommandStats(statsOutputFile,stats.getTupleString());
        }
    }

}
