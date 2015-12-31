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

#include "KnnQueryCommand.h"

#include "../utility/StopWatch.h"
#include "../utility/utility.h"
#include "../utility/serialization.h"
#include "../processing/INE.h"
#include "../processing/Gtree.h"
#include "../processing/ROAD.h"
#include "../processing/MortonList.h"
#include "../processing/SimpleQuadtree.h"
#include "../processing/IER.h"
#include "../processing/StaticRtree.h"

void KnnQueryCommand::execute(int argc, char* argv[])
{
    std::string method = "";
    std::string queryPointsFilePath = "";
    std::string resultsFilePath = "";
    std::string statsFilePath = "";
    unsigned int k = 0;
    
    // Specific Variables
    std::string objSetFilePath; // INE, Optimise SILC
    std::string bgrFilePath = ""; // INE, Gtree, SILC, Distance Browsing
    std::string idxFilePath = ""; // Full Index, Gtree, ROAD, SILC, Distance Browsing
    std::string addIdxFilePath = ""; // Optimised SILC, Optimised Distance Browsing
    std::string objIdxFilePath = ""; // Gtree, ROAD, SILC, Distance Browsing
    
    /*
     * Process Command Line Arguments
     */
    int opt;
    while ((opt = getopt (argc, argv, "m:k:q:r:s:i:g:b:j:a:")) != -1) {
        switch (opt) {
            case 'm':
                method = optarg;
                break;
            case 'k':
                k = std::stoi(optarg);
                break;
            case 'q':
                queryPointsFilePath = optarg;
                break;
            case 'r':
                resultsFilePath = optarg;
                break;
            case 's':
                statsFilePath = optarg;
                break;
            case 'i':
                objSetFilePath = optarg;
                break;
            case 'g':
                bgrFilePath = optarg;
                break;
            case 'b':
                idxFilePath = optarg;
                break;
            case 'j':
                objIdxFilePath = optarg;
                break;
            case 'a':
                addIdxFilePath = optarg;
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
        this->showMethodUsage(method,argv[0]);
        exit(1);
    } 
    
    if (queryPointsFilePath == "" || k <= 0 || resultsFilePath == "" 
        || statsFilePath == "") {
        std::cerr << "Invalid argument(s)!\n\n";
        this->showCommandUsage(argv[0]);
        exit(1);
    }

    /*
     * Parse Query Point File
     */
    std::vector<NodeID> queryNodesIDs = utility::getPointSetFromFile(queryPointsFilePath);
    std::vector<KnnQueryTuple> resultTuples;
    
    std::cout << "Queries will be executed for " << queryNodesIDs.size() << " query point(s)" << std::endl;
    
    if (method == constants::INE_KNN_QUERY) {
        if (argc < 17) {
            // Arguments: -k <number of nearest neighbours -i <object set file 1> [-i <object set file 2> ... ]
            // -q <query point set file> -g <binary graph file> -r <results file> -s <stats file>
            std::cerr << "Too few arguments!\n\n";
            this->showMethodUsage(method,argv[0]);
            exit(1);
        }

        if (objSetFilePath == "" || bgrFilePath == "") {
            std::cerr << "Invalid argument(s)!\n\n";
            this->showMethodUsage(method,argv[0]);
            exit(1);
        }

        resultTuples = this->executeINEQuery(k,queryNodesIDs,objSetFilePath,bgrFilePath);

    } else if (method == constants::GTREE_KNN_QUERY) {
        if (argc < 19) {
            //Arguments: -k <number of nearest neighbours -b <index binary file> -j <object index file>
            //-g <binary graph file> -q <query point set file> -r <results file> -s <stats file>
            std::cerr << "Too few arguments!\n\n";
            this->showMethodUsage(method,argv[0]);
            exit(1);
        }
        
        if (idxFilePath == "" || objIdxFilePath == "" || bgrFilePath == "") {
            std::cerr << "Invalid argument(s)!\n\n";
            this->showMethodUsage(method,argv[0]);
            exit(1);
        }

        resultTuples = this->executeGtreeQuery(k,queryNodesIDs,idxFilePath,objIdxFilePath,bgrFilePath);

    } else if (method == constants::ROAD_KNN_QUERY) {
        if (argc < 17) {
            //Arguments: -k <number of nearest neighbours -b <index binary file> -j <object index file>
            //-q <query point set file> -r <results file> -s <stats file>
            std::cerr << "Too few arguments!\n\n";
            this->showMethodUsage(method,argv[0]);
            exit(1);
        }
        
        if (idxFilePath == "" || objIdxFilePath == "") {
            std::cerr << "Invalid argument(s)!\n\n";
            this->showMethodUsage(method,argv[0]);
            exit(1);
        }

        resultTuples = this->executeROADQuery(k,queryNodesIDs,idxFilePath,objIdxFilePath);

    } else if (method == constants::SILC_KNN_QUERY) {
        if (argc < 19) {
            //Arguments: -k <number of nearest neighbours -b <index binary file> -j <object index file>
            //-g <binary graph file> -q <query point set file> -r <results file> -s <stats file>
            std::cerr << "Too few arguments!\n\n";
            this->showMethodUsage(method,argv[0]);
            exit(1);
        }
        
        if (idxFilePath == "" || objIdxFilePath == "" || bgrFilePath == "") {
            std::cerr << "Invalid argument(s)!\n\n";
            this->showMethodUsage(method,argv[0]);
            exit(1);
        }

        resultTuples = this->executeSILCQuery(method,k,queryNodesIDs,idxFilePath,objIdxFilePath,bgrFilePath);

    } else if (method == constants::OPT_SILC_KNN_QUERY) {
        if (argc < 21) {
            //Arguments: -k <number of nearest neighbours -b <index binary file> -j <object index file>
            //-g <binary graph file> -q <query point set file> -r <results file> -s <stats file>
            //-a <junc index file> -i <object set file path>
            std::cerr << "Too few arguments!\n\n";
            this->showMethodUsage(method,argv[0]);
            exit(1);
        }
        
        if (idxFilePath == "" || objIdxFilePath == "" || bgrFilePath == "" || addIdxFilePath == "" || objSetFilePath == "") {
            std::cerr << "Invalid argument(s)!\n\n";
            this->showMethodUsage(method,argv[0]);
            exit(1);
        }

        resultTuples = this->executeSILCQuery(method,k,queryNodesIDs,idxFilePath,objIdxFilePath,bgrFilePath,addIdxFilePath,objSetFilePath);

    } else if (method == constants::DISTBRWS_KNN_QUERY) {
        if (argc < 19) {
            //Arguments: -k <number of nearest neighbours -b <index binary file> -j <object index file>
            //-g <binary graph file> -q <query point set file> -r <results file> -s <stats file>
            std::cerr << "Too few arguments!\n\n";
            this->showMethodUsage(method,argv[0]);
            exit(1);
        }
        
        if (idxFilePath == "" || objIdxFilePath == "" || bgrFilePath == "") {
            std::cerr << "Invalid argument(s)!\n\n";
            this->showMethodUsage(method,argv[0]);
            exit(1);
        }

        resultTuples = this->executeSILCQuery(method,k,queryNodesIDs,idxFilePath,objIdxFilePath,bgrFilePath);

    } else if (method == constants::OPT_DISTBRWS_KNN_QUERY) {
        if (argc < 21) {
            //Arguments: -k <number of nearest neighbours -b <index binary file> -j <object index file>
            //-g <binary graph file> -q <query point set file> -r <results file> -s <stats file>
            //-a <junc index file> -i <object set file path>
            std::cerr << "Too few arguments!\n\n";
            this->showMethodUsage(method,argv[0]);
            exit(1);
        }
        
        if (idxFilePath == "" || objIdxFilePath == "" || bgrFilePath == "" || addIdxFilePath == "" || objSetFilePath == "") {
            std::cerr << "Invalid argument(s)!\n\n";
            this->showMethodUsage(method,argv[0]);
            exit(1);
        }

        resultTuples = this->executeSILCQuery(method,k,queryNodesIDs,idxFilePath,objIdxFilePath,bgrFilePath,addIdxFilePath,objSetFilePath);

    } else if (method == constants::IER_DIJKSTRA_KNN_QUERY) {
        if (argc < 17) {
            //Arguments: -k <number of nearest neighbours -j <object index file>
            //-g <binary graph file> -q <query point set file> -r <results file> -s <stats file>
            std::cerr << "Too few arguments!\n\n";
            this->showMethodUsage(method,argv[0]);
            exit(1);
        }
        
        if (objIdxFilePath == "" || bgrFilePath == "") {
            std::cerr << "Invalid argument(s)!\n\n";
            this->showMethodUsage(method,argv[0]);
            exit(1);
        }

        resultTuples = this->executeIERQuery(method,k,queryNodesIDs,objIdxFilePath,bgrFilePath,idxFilePath);
    } else if (method == constants::IER_GTREE_KNN_QUERY) {
        if (argc < 19) {
            //Arguments: -k <number of nearest neighbours -j <object index file>
            //-g <binary graph file> -q <query point set file> -r <results file>
            //-s <stats file> -b <index binary file>
            std::cerr << "Too few arguments!\n\n";
            this->showMethodUsage(method,argv[0]);
            exit(1);
        }
        
        if (objIdxFilePath == "" || bgrFilePath == "" || idxFilePath == "") {
            std::cerr << "Invalid argument(s)!\n\n";
            this->showMethodUsage(method,argv[0]);
            exit(1);
        }

        resultTuples = this->executeIERQuery(method,k,queryNodesIDs,objIdxFilePath,bgrFilePath,idxFilePath);
    } else if (method == constants::DB_RTREE_KNN_QUERY) {
        if (argc < 19) {
            //Arguments: -k <number of nearest neighbours -j <object index file>
            //-g <binary graph file> -q <query point set file> -r <results file>
            //-s <stats file> -b <index binary file>
            std::cerr << "Too few arguments!\n\n";
            this->showMethodUsage(method,argv[0]);
            exit(1);
        }
        
        if (objIdxFilePath == "" || bgrFilePath == "" || idxFilePath == "" || addIdxFilePath == "") {
            std::cerr << "Invalid argument(s)!\n\n";
            this->showMethodUsage(method,argv[0]);
            exit(1);
        }

        resultTuples = this->executeIERQuery(method,k,queryNodesIDs,objIdxFilePath,bgrFilePath,idxFilePath,addIdxFilePath);
    } else if (method == constants::IER_PHL_KNN_QUERY) {
        if (argc < 19) {
            //Arguments: -k <number of nearest neighbours -j <object index file>
            //-g <binary graph file> -q <query point set file> -r <results file>
            //-s <stats file> -b <index binary file>
            std::cerr << "Too few arguments!\n\n";
            this->showMethodUsage(method,argv[0]);
            exit(1);
        }
        
        if (objIdxFilePath == "" || bgrFilePath == "" || idxFilePath == "") {
            std::cerr << "Invalid argument(s)!\n\n";
            this->showMethodUsage(method,argv[0]);
            exit(1);
        }

        resultTuples = this->executeIERQuery(method,k,queryNodesIDs,objIdxFilePath,bgrFilePath,idxFilePath);
    } else {
        std::cerr << "Invalid query method!\n\n";
        this->showCommandUsage(argv[0]);
        exit(1);
    }
    
    // Write kNN results/stats to output/files
    for (auto it = resultTuples.begin(); it != resultTuples.end(); ++it) {
        std::cout << it->getMultilineTupleString();
        this->outputCommandStats(statsFilePath,it->getTupleString());
        this->outputCommandStats(resultsFilePath,it->getResultsString());        
    }
    
}

void KnnQueryCommand::showCommandUsage(std::string programName)
{
    std::cerr << "Usage: " << programName << " -c " + constants::QUERY_KNN_CMD + " -m <query method>\n\n"
              << "Methods Options:\n"
              << utility::getFormattedUsageString(constants::INE_KNN_QUERY,"Incremental Network Expansion") + "\n"
              << utility::getFormattedUsageString(constants::GTREE_KNN_QUERY,"kNN using Gtree") + "\n"
              << utility::getFormattedUsageString(constants::ROAD_KNN_QUERY,"kNN using ROAD") + "\n"
              << utility::getFormattedUsageString(constants::SILC_KNN_QUERY,"kNN using SILC incremental search") + "\n"
              << utility::getFormattedUsageString(constants::OPT_SILC_KNN_QUERY,"kNN using SILC with junction optimisation") + "\n"
              << utility::getFormattedUsageString(constants::DISTBRWS_KNN_QUERY,"kNN using SILC Distace Browsing") + "\n"
              << utility::getFormattedUsageString(constants::OPT_DISTBRWS_KNN_QUERY,"kNN using SILC Distace Browsing with junction optimisation") + "\n"
              << utility::getFormattedUsageString(constants::DB_RTREE_KNN_QUERY,"kNN using SILC Distance Browsing by Euclidean NNs") + "\n"
              << utility::getFormattedUsageString(constants::IER_DIJKSTRA_KNN_QUERY,"kNN using IER using Dijkstra") + "\n"
              << utility::getFormattedUsageString(constants::IER_GTREE_KNN_QUERY,"kNN using IER using Gtree") + "\n"
              << utility::getFormattedUsageString(constants::IER_PHL_KNN_QUERY,"kNN using IER using PHL") + "\n";
}

void KnnQueryCommand::showMethodUsage(std::string method, std::string programName)
{
    if (method == constants::INE_KNN_QUERY ) {
        std::cerr << "Usage: " << programName << " -c " + constants::QUERY_KNN_CMD
                  << "-m " + constants::INE_KNN_QUERY + " -k <num nearest neighbours>\n"
                  << "-i <object set file> -q <query point file> -g <binary graph file>\n"
                  << "-r <knn results file> -s <query stats file>\n";
    } else if (method == constants::GTREE_KNN_QUERY ) {
        std::cerr << "Usage: " << programName << " -c " + constants::QUERY_KNN_CMD 
                  << "-m " + constants::GTREE_KNN_QUERY + " -k <num nearest neighbours>\n"
                  << "-b <gtree index binary file> -j <object index file> -g <binary graph file>\n"
                  << "-q <query point file> -r <knn results file> -s <query stats file>\n";
    } else if (method == constants::GTREE_KNN_QUERY ) {
        std::cerr << "Usage: " << programName << " -c " + constants::QUERY_KNN_CMD 
                  << "-m " + constants::ROAD_KNN_QUERY + " -k <num nearest neighbours>\n"
                  << "-b <road index binary file> -j <object index file> -q <query point file>\n"
                  << "-r <knn results file> -s <query stats file>\n";
    } else if (method == constants::SILC_KNN_QUERY ) {
        std::cerr << "Usage: " << programName << " -c " + constants::QUERY_KNN_CMD 
                  << "-m " + constants::SILC_KNN_QUERY + " -k <num nearest neighbours>\n"
                  << "-b <SILC index binary file> -j <object quadtree file> -g <binary graph file>\n"
                  << "-q <query point file> -r <knn results file> -s <query stats file>\n";
    } else if (method == constants::OPT_SILC_KNN_QUERY ) {
        std::cerr << "Usage: " << programName << " -c " + constants::QUERY_KNN_CMD 
                  << "-m " + constants::OPT_SILC_KNN_QUERY + " -k <num nearest neighbours>\n"
                  << "-b <SILC index binary file> -j <object quadtree file> -g <binary graph file>\n"
                  << "-q <query point file> -r <knn results file> -s <query stats file>\n"
                  << "-a <junc index file> -i <object set file path>\n";
    } else if (method == constants::DISTBRWS_KNN_QUERY ) {
        std::cerr << "Usage: " << programName << " -c " + constants::QUERY_KNN_CMD 
                  << "-m " + constants::DISTBRWS_KNN_QUERY + " -k <num nearest neighbours>\n"
                  << "-b <SILC index binary file> -j <object quadtree file> -g <binary graph file>\n"
                  << "-q <query point file> -r <knn results file> -s <query stats file>\n";
    } else if (method == constants::OPT_DISTBRWS_KNN_QUERY ) {
        std::cerr << "Usage: " << programName << " -c " + constants::QUERY_KNN_CMD 
                  << "-m " + constants::OPT_DISTBRWS_KNN_QUERY + " -k <num nearest neighbours>\n"
                  << "-b <SILC index binary file> -j <object quadtree file> -g <binary graph file>\n"
                  << "-q <query point file> -r <knn results file> -s <query stats file>\n"
                  << "-a <junc index file> -i <object set file path>\n";
    } else if (method == constants::IER_DIJKSTRA_KNN_QUERY ) {
        std::cerr << "Usage: " << programName << " -c " + constants::QUERY_KNN_CMD 
                  << "-m " + constants::IER_DIJKSTRA_KNN_QUERY + " -k <num nearest neighbours>\n"
                  << "-j <object index file> -q <query point file> -g <binary graph file>\n"
                  << "-r <knn results file> -s <query stats file>\n";
    } else if (method == constants::IER_GTREE_KNN_QUERY ) {
        std::cerr << "Usage: " << programName << " -c " + constants::QUERY_KNN_CMD 
                  << "-m " + method + " -k <num nearest neighbours>\n"
                  << "-j <object index file> -q <query point file> -g <binary graph file>\n"
                  << "-b <index binary file> -r <knn results file> -s <query stats file>\n";
    } else {
        std::cerr << "Invalid query method!" << std::endl;
        this->showCommandUsage(programName);
    }
}

std::vector<KnnQueryTuple> KnnQueryCommand::executeINEQuery(int k, std::vector<NodeID>& queryNodesIDs, 
                                                            std::string objSetFilePath, std::string bgrFilePath)
{
    // Store results in vector of tuples
    std::vector<KnnQueryTuple> resultTuples;
    
    /*
     * Deserialize Graph from Binary File
     */
    Graph graph = serialization::getIndexFromBinaryFile<Graph>(bgrFilePath);
//     int lastindex = bgrFilePath.find_last_of("."); 
//     std::string dynBgrFilePath = bgrFilePath.substr(0, lastindex) + "_dynamic.bin"; 
//     DynamicGraph graph = serialization::getIndexFromBinaryFile<DynamicGraph>(dynBgrFilePath);

    int nodes = graph.getNumNodes();
    int edges = graph.getNumEdges();
    std::string networkName = graph.getNetworkName();
    INE ine;
    
    /*
     * Execute knn query for each query point for each object set and each query point
     */
    StopWatch sw;
    double queryTimeMs, objSetDensity;
    int objSetSize, objSetVariable;
    std:: string objSetType;
    std::vector<NodeID> kNNs;
    std::vector<EdgeWeight> kNNDistances;
    kNNs.reserve(k);
    kNNDistances.reserve(k);
    
    graph.resetAllObjects();
    graph.parseObjectFile(objSetFilePath,objSetType,objSetDensity,objSetVariable,objSetSize);

    for (auto query_node_it = queryNodesIDs.begin(); query_node_it != queryNodesIDs.end(); ++query_node_it) {
        kNNs.clear();
        kNNDistances.clear();
        sw.start();
        ine.getKNNs(graph,k,*query_node_it,kNNs,kNNDistances);
        sw.stop();
        queryTimeMs = sw.getTimeMs();
        sw.reset();
        
        // Collect stats and return to output to file
        KnnQueryTuple stats(networkName,nodes,edges,*query_node_it,constants::INE_KNN_QUERY,
                            k,queryTimeMs,objSetType,objSetDensity,objSetVariable,objSetSize,kNNs,kNNDistances);
        resultTuples.push_back(stats);
    }        
    
    return resultTuples;
    
}

std::vector<KnnQueryTuple> KnnQueryCommand::executeGtreeQuery(int k, std::vector<NodeID>& queryNodesIDs, 
                                                              std::string idxFilePath, std::string objIdxFilePath,
                                                              std::string bgrFilePath)
{
    // Store results in vector of tuples
    std::vector<KnnQueryTuple> resultTuples;

    /*
     * Deserialize Gtree Index and Occurence List from Binary File
     */
    Gtree gtree = serialization::getIndexFromBinaryFile<Gtree>(idxFilePath);

    Graph graph = serialization::getIndexFromBinaryFile<Graph>(bgrFilePath);
    OccurenceList occList = serialization::getIndexFromBinaryFile<OccurenceList>(objIdxFilePath);
    int nodes = gtree.getNumNodes();
    int edges = gtree.getNumEdges();
    std::string networkName = gtree.getNetworkName();
    double objSetDensity = occList.getObjSetDensity();
    int objSetSize = occList.getObjSetSize(), objSetVariable = occList.getObjSetVariable();
    std::string objSetType = occList.getObjSetType();
    
    /*
     * Execute knn query for each query point for each object set and each query point
     */
    StopWatch sw;
    double queryTimeMs;
    std::vector<NodeID> kNNs;
    std::vector<EdgeWeight> kNNDistances;
    kNNs.reserve(k);
    kNNDistances.reserve(k);
        
    for (auto it = queryNodesIDs.begin(); it != queryNodesIDs.end(); ++it) {
        kNNs.clear();
        kNNDistances.clear();
        sw.start();
        gtree.getKNNs(occList,k,*it,kNNs,kNNDistances,graph); // This will move return value to kNNs
        sw.stop();
        queryTimeMs = sw.getTimeMs();
        sw.reset();
        // Collect stats and return to output to file
        KnnQueryTuple stats(networkName,nodes,edges,*it,constants::GTREE_KNN_QUERY,
                            k,queryTimeMs,objSetType,objSetDensity,objSetVariable,objSetSize,kNNs,kNNDistances);
        stats.addSupplementaryFields("max_leaf_size",std::to_string(gtree.getMaxLeafSize()));
        stats.addSupplementaryFields("fanout",std::to_string(gtree.getFanout()));
        resultTuples.push_back(stats);
    }        
    
    return resultTuples;
}

std::vector<KnnQueryTuple> KnnQueryCommand::executeROADQuery(int k, std::vector<NodeID>& queryNodesIDs, 
                                                             std::string idxFilePath, std::string objIdxFilePath)
{
    // Store results in vector of tuples
    std::vector<KnnQueryTuple> resultTuples;

    /*
     * Deserialize Gtree Index and Occurence List from Binary File
     */
    ROAD road = serialization::getIndexFromBinaryFile<ROAD>(idxFilePath);
    AssociationDirectory assocDir = serialization::getIndexFromBinaryFile<AssociationDirectory>(objIdxFilePath);
    int nodes = road.getNumNodes();
    int edges = road.getNumEdges();
    std::string networkName = road.getNetworkName();
    double objSetDensity = assocDir.getObjSetDensity();
    int objSetSize = assocDir.getObjSetSize(), objSetVariable = assocDir.getObjSetVariable();
    std::string objSetType = assocDir.getObjSetType();
    
    /*
     * Execute knn query for each query point for each object set and each query point
     */
    StopWatch sw;
    double queryTimeMs;
    std::vector<NodeID> kNNs;
    std::vector<EdgeWeight> kNNDistances;
    kNNs.reserve(k);
    kNNDistances.reserve(k);
        
    for (auto it = queryNodesIDs.begin(); it != queryNodesIDs.end(); ++it) {
        kNNs.clear();
        kNNDistances.clear();
        sw.start();
        road.getKNNs(assocDir,k,*it,kNNs,kNNDistances); // This will move return value to kNNs
        sw.stop();
        queryTimeMs = sw.getTimeMs();
        sw.reset();
        // Collect stats and return to output to file
        KnnQueryTuple stats(networkName,nodes,edges,*it,constants::ROAD_KNN_QUERY,
                            k,queryTimeMs,objSetType,objSetDensity,objSetVariable,objSetSize,kNNs,kNNDistances);
        stats.addSupplementaryFields("levels",std::to_string(road.getNumLevels()));
        stats.addSupplementaryFields("fanout",std::to_string(road.getFanout()));
        resultTuples.push_back(stats);
    }        
    
    return resultTuples;
}

std::vector<KnnQueryTuple> KnnQueryCommand::executeSILCQuery(std::string method, int k, std::vector<NodeID>& queryNodesIDs, 
                                                             std::string idxFilePath, std::string objIdxFilePath,
                                                             std::string bgrFilePath, std::string addIdxFile, 
                                                             std::string objSetFilePath)
{
    // Store results in vector of tuples
    std::vector<KnnQueryTuple> resultTuples;
    
    /*
     * Deserialize SILC Index and Occurence List from Binary File
     */
    Graph graph = serialization::getIndexFromBinaryFile<Graph>(bgrFilePath);
    SILCPathOracle silc = serialization::getIndexFromBinaryFile<SILCPathOracle>(idxFilePath);
    Junction junc;
    if (method == constants::OPT_SILC_KNN_QUERY || method == constants::OPT_DISTBRWS_KNN_QUERY) {
        junc = serialization::getIndexFromBinaryFile<Junction>(addIdxFile);
    }
    int nodes = silc.getNumNodes();
    int edges = silc.getNumEdges();
    std::string networkName = silc.getNetworkName();
    double objSetDensity;
    int objSetSize, objSetVariable;
    std::string objSetType;
    SimpleQuadtree objHierarchy = serialization::getIndexFromBinaryFile<SimpleQuadtree>(objIdxFilePath);
    objSetDensity = objHierarchy.getObjSetDensity();
    objSetSize = objHierarchy.getObjSetSize();
    objSetType = objHierarchy.getObjSetType();
    objSetVariable = objHierarchy.getObjSetVariable();
    
    /*
     * Execute knn query for each query point for each object set and each query point
     */
    StopWatch sw;
    double queryTimeMs;
    std::vector<NodeID> kNNs;
    std::vector<EdgeWeight> kNNDistances;
    kNNs.reserve(k);
    kNNDistances.reserve(k);
        
    for (auto it = queryNodesIDs.begin(); it != queryNodesIDs.end(); ++it) {
        kNNs.clear();
        kNNDistances.clear();
        if (method == constants::SILC_KNN_QUERY) {
            sw.start();
            silc.getKNNs(objHierarchy,k,*it,kNNs,kNNDistances,graph); // This will move return value to kNNs
            sw.stop();
        } else if (method == constants::DISTBRWS_KNN_QUERY) {
            sw.start();
            silc.getKNNsByDistanceBrowsing(objHierarchy,k,*it,kNNs,kNNDistances,graph);
            sw.stop();
        } else if (method == constants::OPT_SILC_KNN_QUERY) {
            sw.start();
            silc.getKNNsByOptimisedSILC(objHierarchy,k,*it,kNNs,kNNDistances,graph,junc);
            sw.stop();
        } else if (method == constants::OPT_DISTBRWS_KNN_QUERY) {
            sw.start();
            silc.getKNNsByOptimisedDistanceBrowsing(objHierarchy,k,*it,kNNs,kNNDistances,graph,junc);
            sw.stop();
        } else {
            std::cerr << "Invalid version of SILC kNN query!" << std::endl;
            exit(1);
        }
        queryTimeMs = sw.getTimeMs();
        sw.reset();
        // Collect stats and return to output to file
        KnnQueryTuple stats(networkName,nodes,edges,*it,method,k,
                            queryTimeMs,objSetType,objSetDensity,objSetVariable,
                            objSetSize,kNNs,kNNDistances);
        stats.addSupplementaryFields("obj_qt_max_leaf_size",std::to_string(objHierarchy.getMaxLeafItems()));
        resultTuples.push_back(stats);
    }        
    
    return resultTuples;
}

std::vector<KnnQueryTuple> KnnQueryCommand::executeIERQuery(std::string method, int k, std::vector<NodeID>& queryNodesIDs, 
                                                            std::string objIdxFilePath, std::string bgrFilePath, 
                                                            std::string idxFilePath, std::string addIdxFile)
{
    // Store results in vector of tuples
    std::vector<KnnQueryTuple> resultTuples;

    /*
     * Deserialize Grap, Rtree and (potentially) other indexes from binary files
     */
    Graph graph = serialization::getIndexFromBinaryFile<Graph>(bgrFilePath);
    StaticRtree rtree = serialization::getIndexFromBinaryFile<StaticRtree>(objIdxFilePath);
    Gtree gtree;
    SILCPathOracle silc;
    Junction junc;
    PrunedHighwayLabeling phl;
    if (method == constants::IER_GTREE_KNN_QUERY) {
        gtree = serialization::getIndexFromBinaryFile<Gtree>(idxFilePath);
    } else if (method == constants::DB_RTREE_KNN_QUERY) {
        silc = serialization::getIndexFromBinaryFile<SILCPathOracle>(idxFilePath);
        junc = serialization::getIndexFromBinaryFile<Junction>(addIdxFile);
    } else if (method == constants::IER_PHL_KNN_QUERY) {
        phl.LoadLabel(idxFilePath.c_str());
    }

    int nodes = graph.getNumNodes();
    int edges = graph.getNumEdges();
    std::string networkName = graph.getNetworkName();
    IER ier;

    /*
     * Execute knn query for each query point for each object set and each query point
     */
    StopWatch sw;
    double queryTimeMs, objSetDensity = rtree.getObjSetDensity();
    int objSetSize = rtree.getObjSetSize(), objSetVariable = rtree.getObjSetVariable();
    std::string objSetType = rtree.getObjSetType();
    std::vector<NodeID> kNNs;
    std::vector<EdgeWeight> kNNDistances;
    kNNs.reserve(k);
    kNNDistances.reserve(k);

    for (auto query_node_it = queryNodesIDs.begin(); query_node_it != queryNodesIDs.end(); ++query_node_it) {
        kNNs.clear();
        kNNDistances.clear();
        sw.start();
        if (method == constants::IER_DIJKSTRA_KNN_QUERY) {
            ier.getKNNsByDijkstra(rtree,k,*query_node_it,kNNs,kNNDistances,graph);
        } else if (method == constants::IER_GTREE_KNN_QUERY) {
            ier.getKNNsByGtree(rtree,k,*query_node_it,kNNs,kNNDistances,graph,gtree);
        } else if (method == constants::DB_RTREE_KNN_QUERY) {
            silc.getKNNsByDistanceBrowsingViaRtree(rtree,k,*query_node_it,kNNs,kNNDistances,graph,junc);
        } else if (method == constants::IER_PHL_KNN_QUERY) {
            ier.getKNNsByPHL(rtree,k,*query_node_it,kNNs,kNNDistances,graph,phl);
        } else {
            std::cerr << "Invalid version of IER kNN query!" << std::endl;
            exit(1);
        }
        sw.stop();
        queryTimeMs = sw.getTimeMs();
        sw.reset();

        // Collect stats and return to output to file
        KnnQueryTuple stats(networkName,nodes,edges,*query_node_it,method,
                            k,queryTimeMs,objSetType,objSetDensity,objSetVariable,objSetSize,kNNs,kNNDistances);
        resultTuples.push_back(stats);
    }

    return resultTuples;

}
