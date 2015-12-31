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

#include "RoadNetworkIndexCommand.h"

#include "../processing/MortonList.h"
#include "../processing/Quadtree.h"
#include "../processing/ALT.h"
#include "../processing/Gtree.h"
#include "../processing/ROAD.h"
#include "../processing/pruned_highway_labeling.h"
#include "../processing/ShortestPathWrapper.h"
#include "../tuple/IndexTuple.h"
#include "../common.h"
#include "../utility/utility.h"
#include "../utility/StopWatch.h"
#include "../utility/serialization.h"

#include <iostream>

void RoadNetworkIndexCommand::execute(int argc, char* argv[])
{
    std::string indexType = "";
    std::string bgrFilePath = "";
    std::string indexOutputFilePath = "";
    std::string statsOutputFilePath = "";
    std::string dataFilePrefix = "";
    std::string suppIdxOuputFilePrefix = "";
    int maxQuadTreeLeafItems = 0; // SILC
    unsigned int numLandmarks = 0; // ALT
    std::string landmarkMethod = ""; // ALT
    unsigned int fanout = 0; // Gtree, ROAD
    unsigned int maxLeafSize = 0; // Gtree
    unsigned int numLevels = 0; // ROAD
    unsigned int gridSize = 0; // TNR
    
    /*
     * Process Command Line Arguments
     */
    int opt;
    while ((opt = getopt (argc, argv, "i:b:o:s:q:n:l:f:t:h:d:a:g:")) != -1) {
        switch (opt) {
            case 'i':
                indexType = optarg;
                break;
            case 'b':
                bgrFilePath = optarg;
                break;
            case 'o':
                indexOutputFilePath = optarg;
                break;
            case 's':
                statsOutputFilePath = optarg;
                break;
            case 'q':
                maxQuadTreeLeafItems = std::stoi(optarg);
                break;
            case 'n':
                numLandmarks = std::stoul(optarg);
                break;
            case 'l':
                landmarkMethod = optarg;
                break;
            case 'f':
                fanout = std::stoul(optarg);
                break;
            case 't':
                maxLeafSize = std::stoul(optarg);
                break;
            case 'h':
                numLevels = std::stoul(optarg);
                break;
            case 'd':
                dataFilePrefix = optarg;
                break;
            case 'a':
                suppIdxOuputFilePrefix = optarg;
                break;
            case 'g':
                gridSize = std::stoul(optarg);
                break;
            default:
                std::cerr << "Unknown option(s) provided!\n\n";
                showCommandUsage(argv[0]);
                exit(1);
        }
    }   
    
    if (argc == 5) {
        // This is 5 so that user can just enter method to find out what parameters are required for method
        this->showIndexTypeUsage(indexType,argv[0]);
        exit(1);
    } else if (argc < 11) {
        // Arguments: -i <index type> -b <binary graph file> -o <index output file> -s <stats file>
        std::cerr << "Too few main arguments!\n\n";
        this->showIndexTypeUsage(indexType,argv[0]);
        exit(1);
    }
    
    if (indexType == "" || bgrFilePath == "" || indexOutputFilePath == ""
        || statsOutputFilePath == "") {
        std::cerr << "Invalid main argument(s)!\n\n";
        this->showCommandUsage(argv[0]);
        exit(1);
    }
    
    std::string networkName = "";
    int numNodes = 0, numEdges = 0;
    
    if (indexType == constants::IDX_SILC_CMD) {
        if (argc < 13) {
            // Arguments: -i <index type> -b <binary graph file> -o <index output file> -s <stats file> -q <max leaf items in quadtree>
            std::cerr << "Too few arguments!\n\n";
            this->showIndexTypeUsage(indexType,argv[0]);
            exit(1);
        }

        if (maxQuadTreeLeafItems == 0) {
            std::cerr << "Invalid maximum quadtree leaf size!\n";
            this->showIndexTypeUsage(indexType,argv[0]);
            exit(1);
        }

        Graph graph = serialization::getIndexFromBinaryFile<Graph>(bgrFilePath);
        this->buildSILC(graph,maxQuadTreeLeafItems,indexOutputFilePath,statsOutputFilePath);
    } else if (indexType == constants::IDX_ALT_CMD) {
        if (argc < 15) {
            // Arguments: -i <index type> -b <binary graph file> -o <index output file> -s <stats file> -n <num landmarks> -l <landmark method>
            std::cerr << "Too few arguments!\n\n";
            this->showIndexTypeUsage(indexType,argv[0]);
            exit(1);
        }

        if (numLandmarks == 0 || landmarkMethod == "") {
            std::cerr << "Invalid number of landmarks or selection method!\n";
            this->showIndexTypeUsage(indexType,argv[0]);
            exit(1);
        }

        Graph graph = serialization::getIndexFromBinaryFile<Graph>(bgrFilePath);
        this->buildALT(graph,numLandmarks,landmarkMethod,indexOutputFilePath,statsOutputFilePath);
    } else if (indexType == constants::IDX_GTREE_CMD) {
        if (argc < 15) {
            // Arguments: -i <index type> -b <binary graph file> -o <index output file> -s <stats file> -f <fanout> -t <max leaf size>
            std::cerr << "Too few arguments!\n\n";
            this->showIndexTypeUsage(indexType,argv[0]);
            exit(1);
        }

        if (fanout < 2 || maxLeafSize == 0) {
            std::cerr << "Invalid fanout or max leaf size (tau)!\n";
            this->showIndexTypeUsage(indexType,argv[0]);
            exit(1);
        }

        Graph graph = serialization::getIndexFromBinaryFile<Graph>(bgrFilePath);
        this->buildGtree(graph,fanout,maxLeafSize,indexOutputFilePath,statsOutputFilePath);
    } else if (indexType == constants::IDX_ROUTEOVERLAY_CMD) {
        if (argc < 15) {
            // Arguments: -i <index type> -b <binary graph file> -o <index output file> -s <stats file> -f <fanout> -h <number of levels>
            std::cerr << "Too few arguments!\n\n";
            this->showIndexTypeUsage(indexType,argv[0]);
            exit(1);
        }

        if (fanout < 2 || numLevels == 0) {
            std::cerr << "Invalid fanout or number of levels (height)!\n";
            this->showIndexTypeUsage(indexType,argv[0]);
            exit(1);
        }

        Graph graph = serialization::getIndexFromBinaryFile<Graph>(bgrFilePath);
        this->buildRouteOverlay(graph,fanout,numLevels,indexOutputFilePath,statsOutputFilePath);
    } else if (indexType == constants::IDX_PHL_CMD) {
        if (argc < 13) {
            // Arguments: -i <index type> -b <binary graph file> -o <index output file> -e <intermediary data output file> -s <stats file>
            std::cerr << "Too few arguments!\n\n";
            this->showIndexTypeUsage(indexType,argv[0]);
            exit(1);
        }

        if (dataFilePrefix == "") {
            std::cerr << "Invalid intermediary data file prefix!\n";
            this->showIndexTypeUsage(indexType,argv[0]);
            exit(1);
        }

        std::string dataFilePath = dataFilePrefix + ".tsv";
        this->buildIntermediaryFiles(indexType,bgrFilePath,dataFilePath,"",networkName,numNodes,numEdges);
        this->buildPHL(dataFilePath,indexOutputFilePath,statsOutputFilePath,networkName,numNodes,numEdges);
        // Delete intermediary text file as it is not needed
        std::remove(dataFilePath.c_str());
    } else if (indexType == constants::IDX_CH_CMD) {
        if (argc < 15) {
            // Arguments: -i <index type> -b <binary graph file> -o <ch index file> -d <data output file prefix> 
            // -a <supplementary index file output prefix> -s <stats file>
            std::cerr << "Too few arguments!\n\n";
            this->showIndexTypeUsage(indexType,argv[0]);
            exit(1);
        }

        if (dataFilePrefix == "" || suppIdxOuputFilePrefix == "") {
            std::cerr << "Invalid intermediary data or supplementary index file prefix!\n";
            this->showIndexTypeUsage(indexType,argv[0]);
            exit(1);
        }

        std::string txtEdgeFilePath = dataFilePrefix + ".ddsg";
        std::string txtNodeFilePath = dataFilePrefix + ".zero_co";
        this->buildIntermediaryFiles(indexType,bgrFilePath,txtEdgeFilePath,txtNodeFilePath,networkName,numNodes,numEdges);
        
        std::string bgrIntFile = suppIdxOuputFilePrefix  + ".ext_bin";
        std::string bcoIntFile = suppIdxOuputFilePrefix  + ".ext_co_bin";
        std::string hcnFile = suppIdxOuputFilePrefix  + ".hcn";
        this->buildCH(txtEdgeFilePath,txtNodeFilePath,bgrIntFile,bcoIntFile,indexOutputFilePath,hcnFile,statsOutputFilePath,networkName,numNodes,numEdges);
        
        // Delete intermediary text file as it is not needed
        std::remove(txtEdgeFilePath.c_str());
        std::remove(txtNodeFilePath.c_str());
    } else if (indexType == constants::IDX_TNR_CMD) {
        if (argc < 17) {
            // Arguments: -i <index type> -b <binary graph file> -o <ch index file> -d <data output file prefix> 
            // -a <supplementary index file output prefix> -s <stats file> -g <grid size>
            std::cerr << "Too few arguments!\n\n";
            this->showIndexTypeUsage(indexType,argv[0]);
            exit(1);
        }

        if (dataFilePrefix == "" || suppIdxOuputFilePrefix == "" || gridSize < 16) {
            std::cerr << "Invalid intermediary data or supplementary index file prefix or grid size!\n";
            this->showIndexTypeUsage(indexType,argv[0]);
            exit(1);
        }

        std::string bgrIntFile = suppIdxOuputFilePrefix  + ".ext_bin";
        std::string bcoIntFile = suppIdxOuputFilePrefix  + ".ext_co_bin";
        std::string chFile = suppIdxOuputFilePrefix  + "." + constants::IDX_CH_CMD;
        this->buildIntermediaryFiles(indexType,bgrFilePath,"","",networkName,numNodes,numEdges);
        this->buildTNR(gridSize,bgrIntFile,bcoIntFile,chFile,indexOutputFilePath,statsOutputFilePath,networkName,numNodes,numEdges);
    } else {
        std::cerr << "Invalid set type!\n\n";
        this->showCommandUsage(argv[0]);
        exit(1);
    }
    
}

void RoadNetworkIndexCommand::showCommandUsage(std::string programName)
{
    std::cerr << "Usage: " << programName << " -c " + constants::NETWORK_IDX_CMD + " -i <index type>\n\n"
              << "Index Type Options:\n"
              << utility::getFormattedUsageString(constants::IDX_SILC_CMD,"Create SILC index") << "\n"
              << utility::getFormattedUsageString(constants::IDX_ALT_CMD,"Create ALT index") << "\n"
              << utility::getFormattedUsageString(constants::IDX_GTREE_CMD,"Create Gtree index") << "\n"
              << utility::getFormattedUsageString(constants::IDX_ROUTEOVERLAY_CMD,"Create Route Overlay index") << "\n"
              << utility::getFormattedUsageString(constants::IDX_PHL_CMD,"Create PHL index") << "\n"
              << utility::getFormattedUsageString(constants::IDX_CH_CMD,"Create Contraction Hierarchies index") << "\n"
              << utility::getFormattedUsageString(constants::IDX_TNR_CMD,"Create Transit Node Routing index") << "\n";
}

void RoadNetworkIndexCommand::showIndexTypeUsage(std::string indexType, std::string programName)
{
    if (indexType == constants::IDX_SILC_CMD) {
        std::cerr << "Usage: " << programName << " -c " + constants::NETWORK_IDX_CMD
                  << " -i " + constants::IDX_SILC_CMD + " -b <binary graph file>\n"
                  << "-o <index output file> -s <stats file> -q <max leaf items in quadtree>\n";
    } else if (indexType == constants::IDX_ALT_CMD) {
        std::cerr << "Usage: " << programName << " -c " + constants::NETWORK_IDX_CMD
                  << " -i " + constants::IDX_ALT_CMD + " -b <binary graph file>\n"
                  << "-o <index output file> -s <stats file> -n <num landmarks>\n"
                  << "-l <landmark method>\n";
    } else if (indexType == constants::IDX_GTREE_CMD) {
        std::cerr << "Usage: " << programName << " -c " + constants::NETWORK_IDX_CMD
                  << " -i " + constants::IDX_GTREE_CMD + " -b <binary graph file>\n"
                  << "-o <index output file> -s <stats file> -f <fanout>\n"
                  << "-t <max leaf size (tau)>\n";
    } else if (indexType == constants::IDX_ROUTEOVERLAY_CMD) {
        std::cerr << "Usage: " << programName << " -c " + constants::NETWORK_IDX_CMD
                  << " -i " + constants::IDX_ROUTEOVERLAY_CMD + " -b <binary graph file>\n"
                  << "-o <index output file> -s <stats file> -f <fanout>\n"
                  << "-h <num levels (height)>\n";
    } else if (indexType == constants::IDX_PHL_CMD) {
        std::cerr << "Usage: " << programName << " -c " + constants::NETWORK_IDX_CMD
                  << " -i " + constants::IDX_PHL_CMD + " -b <binary graph file>\n"
                  << "-o <index output file> -s <stats file> -d <intermediary data output prefix>\n";
    } else if (indexType == constants::IDX_CH_CMD) {
        std::cerr << "Usage: " << programName << " -c " + constants::NETWORK_IDX_CMD
                  << " -i " + constants::IDX_CH_CMD + " -b <binary graph file>\n"
                  << "-o <index output file> -s <stats file> -d <intermediary data output prefix>\n"
                  << "-a <additional index output prefix>\n";
    } else if (indexType == constants::IDX_TNR_CMD) {
        std::cerr << "Usage: " << programName << " -c " + constants::NETWORK_IDX_CMD
                  << " -i " + constants::IDX_TNR_CMD + " -b <binary graph file>\n"
                  << "-o <index output file> -s <stats file> -d <intermediary data output prefix>\n"
                  << "-a <additional index output prefix> -g <grid size>\n";
    } else {
        std::cerr << "Invalid index type!" << std::endl;
        this->showCommandUsage(programName);
    }
}

void RoadNetworkIndexCommand::buildSILC(Graph& graph, int maxQuadtreeLeafSize, std::string idxOutputFile, std::string statsOutputFile)
{
    /*
     * Build SILC Index
     */
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
    
    /*
     * Serialize to Binary File
     */    
    //serialization::outputIndexToBinaryFile<Quadtree>(quadtree,quadtreeOutputFilePath); // Not needed again
    serialization::outputIndexToBinaryFile<SILCPathOracle>(silc,idxOutputFile);
    
    /*
     * Collect Stats and Output
     */
    // Note: SILC uses Graph to get adj node and weight associated with colour (only edges and edge weights same parts as INE)
    double indexSizeMB = silc.computeIndexSize() + graph.computeINEIndexSize();;

    IndexTuple stats(silc.getNetworkName(),silc.getNumNodes(),silc.getNumEdges(),constants::IDX_SILC_CMD,processingTimeMs,indexSizeMB);

    std::cout << stats.getMultilineTupleString();

    this->outputCommandStats(statsOutputFile,stats.getTupleString());

    std::cout << "\nNote: Quadtree Construction Time = " << quadtreeConstructionTime << "ms\n" << std::endl;
    
    std::cout << "Binary SILC index successfully created!" << std::endl;
}

void RoadNetworkIndexCommand::buildALT(Graph& graph, int numLandmarks, std::string landmarkMethod, 
                                       std::string idxOutputFile, std::string statsOutputFile)
{
    bool validMethod;
    ALT alt(graph.getNetworkName(),graph.getNumNodes(),graph.getNumEdges());
    LANDMARK_TYPE landmarkType = alt.getLandmarkType(landmarkMethod,validMethod);
    if (!validMethod) {
        std::cerr << "Invalid landmark selection method!" << std::endl;
        std::exit(1);
    }
    
    /*
     * Build ALT Index
     */
    StopWatch sw;
    sw.start();
    alt.buildALT(graph,landmarkType,numLandmarks);
    sw.stop();
    double processingTimeMs = sw.getTimeMs();    

    /*
     * Serialize to Binary File
     */    
    //serialization::outputIndexToBinaryFile<Quadtree>(quadtree,quadtreeOutputFilePath); // Not needed again
    serialization::outputIndexToBinaryFile<ALT>(alt,idxOutputFile);
    
    /*
     * Collect Stats and Output
     */
    double indexSizeMB = alt.computeIndexSize();

    IndexTuple stats(alt.getNetworkName(),alt.getNumNodes(),alt.getNumEdges(),constants::IDX_ALT_CMD,processingTimeMs,indexSizeMB);

    std::cout << stats.getMultilineTupleString();

    this->outputCommandStats(statsOutputFile,stats.getTupleString());

    std::cout << "Binary ALT index successfully created!" << std::endl;
}

void RoadNetworkIndexCommand::buildGtree(Graph& graph, int fanout, std::size_t maxLeafSize, std::string idxOutputFile, std::string statsOutputFile)
{
    /*
     * Building Gtree Index
     */
    StopWatch sw;
    
    sw.start();
    Gtree gtree(graph.getNetworkName(),graph.getNumNodes(),graph.getNumEdges(),fanout,maxLeafSize);
    gtree.buildGtree(graph);
    sw.stop();
    
    /*
     * Serialize to Binary File
     */    
    serialization::outputIndexToBinaryFile<Gtree>(gtree,idxOutputFile);

    /*
     * Collect Stats and Output
     */
    double processingTimeMs = sw.getTimeMs();
    // Note: G-tree uses Graph to search for object within leaf (only edges and edge weights same parts as INE)
    double memoryUsage = gtree.computeIndexSize() + graph.computeINEIndexSize();
    
    IndexTuple stats(gtree.getNetworkName(),gtree.getNumNodes(),gtree.getNumEdges(),
                     constants::IDX_GTREE_CMD,processingTimeMs,memoryUsage);
    stats.addSupplementaryFields("max_leaf_size",std::to_string(maxLeafSize));
    stats.addSupplementaryFields("fanout",std::to_string(fanout));
    stats.addSupplementaryFields("levels",std::to_string(gtree.getNumLevels()));

    std::cout << stats.getMultilineTupleString();

    this->outputCommandStats(statsOutputFile,stats.getTupleString());
    
    std::cout << "Gtree index successfully created!" << std::endl;

}

void RoadNetworkIndexCommand::buildRouteOverlay(Graph& graph, int fanout, int levels, std::string idxOutputFile, std::string statsOutputFile)
{
    /*
     * Building ROAD Index
     */
    StopWatch sw;
    
    sw.start();
    ROAD road(graph.getNetworkName(),graph.getNumNodes(),graph.getNumEdges(),fanout,levels);
    road.buildRouteOverlay(graph);
    sw.stop();
    
    /*
     * Serialize to Binary File
     */    
    serialization::outputIndexToBinaryFile<ROAD>(road,idxOutputFile);

    /*
     * Collect Stats and Output
     */
    double processingTimeMs = sw.getTimeMs();   
    double memoryUsage = road.computeIndexSize();
    
    IndexTuple stats(road.getNetworkName(),road.getNumNodes(),road.getNumEdges(),
                     constants::IDX_ROUTEOVERLAY_CMD,processingTimeMs,memoryUsage);
    stats.addSupplementaryFields("fanout",std::to_string(fanout));
    stats.addSupplementaryFields("levels",std::to_string(levels));

    std::cout << stats.getMultilineTupleString();

    this->outputCommandStats(statsOutputFile,stats.getTupleString());
    
    std::cout << "Route Overlay index successfully created!" << std::endl;

}

void RoadNetworkIndexCommand::buildPHL(std::string txtGrFile, std::string idxOutputFile, std::string statsOutputFile, std::string networkName, int numNodes, int numEdges)
{
    PrunedHighwayLabeling phl;
    phl.ConstructLabel(txtGrFile.c_str());
    
    double processingTimeMs = phl.getConstructionTime()*1000; // convert to ms
    double memoryUsage = phl.computeIndexSize();

    IndexTuple stats(networkName,numNodes,numEdges,constants::IDX_PHL_CMD,processingTimeMs,memoryUsage);

    std::cout << stats.getMultilineTupleString();

    phl.StoreLabel(idxOutputFile.c_str());
    
    this->outputCommandStats(statsOutputFile,stats.getTupleString());
    
    std::cout << "PHL index successfully created!" << std::endl;
}

void RoadNetworkIndexCommand::buildCH(std::string graphFilePath, std::string coordinateFilePath, std::string bgrIntFilePath,
                                      std::string bcoIntFilePath, std::string chFilePath, std::string noFilePath, 
                                      std::string statsOutputFile, std::string networkName, int numNodes, int numEdges)
{
    ShortestPathWrapper spw;

    spw.buildIntermediaryGraph(graphFilePath,coordinateFilePath,bgrIntFilePath,bcoIntFilePath);

    std::cout << "Intermediary graph indexes successfully created!" << std::endl;
    
    spw.buildCH(bgrIntFilePath,chFilePath,noFilePath);
    
    double processingTimeMs = spw.getConstructionTimeMs();
    double memoryUsage = spw.getIndexSizeMB();

    IndexTuple stats(networkName,numNodes,numEdges,constants::IDX_CH_CMD,processingTimeMs,memoryUsage);

    std::cout << stats.getMultilineTupleString();

    this->outputCommandStats(statsOutputFile,stats.getTupleString());

    std::cout << "CH index successfully created!" << std::endl;
}

void RoadNetworkIndexCommand::buildTNR(int gridSize, std::string bgrIntFilePath, std::string bcoIntFilePath,
                                       std::string chFilePath, std::string tnrFilePath, std::string statsOutputFile, 
                                       std::string networkName, int numNodes, int numEdges)
{
    ShortestPathWrapper spw;
    spw.buildTNR(bgrIntFilePath,bcoIntFilePath,chFilePath,tnrFilePath,gridSize);

    double processingTimeMs = spw.getConstructionTimeMs();
    double memoryUsage = spw.getIndexSizeMB();
    
    IndexTuple stats(networkName,numNodes,numEdges,constants::IDX_TNR_CMD,processingTimeMs,memoryUsage);

    std::cout << stats.getMultilineTupleString();
    
    this->outputCommandStats(statsOutputFile,stats.getTupleString());

    std::cout << "TNR index successfully created!" << std::endl;}


void RoadNetworkIndexCommand::buildIntermediaryFiles(std::string method, std::string bgrFile, std::string outputEdgeFile, std::string outputNodeFIle, 
                                                     std::string& networkName, int& numNodes, int& numEdges)
{
    Graph graph = serialization::getIndexFromBinaryFile<Graph>(bgrFile);
    networkName = graph.getNetworkName();
    numNodes = graph.getNumNodes();
    numEdges = graph.getNumEdges();
    if (method == constants::IDX_CH_CMD) {
        graph.outputToDDSGFile(outputEdgeFile);
        graph.outputZeroedCoordinatesToFile(outputNodeFIle);
    } else if (method == constants::IDX_PHL_CMD) {
        graph.outputToTSVFile(outputEdgeFile);
    }
}
