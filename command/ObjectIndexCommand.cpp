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

#include "ObjectIndexCommand.h"

#include "../processing/Gtree.h"
#include "../processing/ROAD.h"
#include "../processing/SimpleQuadtree.h"
#include "../processing/StaticRtree.h"
#include "../utility/StopWatch.h"
#include "../tuple/ObjectIndexTuple.h"
#include "../utility/utility.h"
#include "../utility/serialization.h"

#include <iostream>

void ObjectIndexCommand::execute(int argc, char* argv[])
{
    std::string indexType = "";
    std::string objSetFilePath = "";
    std::string binaryIndexFilePath = "";
    std::string outputFilePath = "";
    std::string statsFilePath = "";
    int maxQuadTreeLeafSize = 0;
    int branchFactor = 2;
    
    /*
     * Process Command Line Arguments
     */
    int opt;
    while ((opt = getopt (argc, argv, "t:i:b:o:s:l:f:")) != -1) {
        switch (opt) {
            case 't':
                indexType = optarg;
                break;
            case 'i':
                objSetFilePath = optarg;
                break;
            case 'b':
                binaryIndexFilePath = optarg;
                break;
            case 'o':
                outputFilePath = optarg;
                break;
            case 's':
                statsFilePath = optarg;
                break;
            case 'l':
                maxQuadTreeLeafSize = std::stoi(optarg);
                break;
            case 'f':
                branchFactor = std::stoi(optarg);
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
        // Arguments: -t type -i <object set file> -b <binary index file> -o <index output file> -s <stats file>
        std::cerr << "Too few arguments!\n\n";
        this->showIndexTypeUsage(indexType,argv[0]);
        return;
    }
    
    if (indexType == "" || objSetFilePath == "" || binaryIndexFilePath == ""
        || outputFilePath == "" || statsFilePath == "") {
        std::cerr << "Invalid argument(s)!\n\n";
        this->showCommandUsage(argv[0]);
        exit(1);
    }
    
    if (indexType == constants::OBJ_IDX_GTREE) {
        if (argc < 11) {
            // Arguments: -t type -i <object set file> -b <binary index file> -o <index output file> -s <stats file>
            std::cerr << "Too few arguments!\n\n";
            this->showIndexTypeUsage(indexType,argv[0]);
            exit(1);
        }
        this->createGtreeOccurenceList(objSetFilePath,binaryIndexFilePath,outputFilePath,statsFilePath);
    } else if (indexType == constants::OBJ_IDX_ROAD) {
        if (argc < 11) {
            // Arguments: -t type -i <object set file> -b <binary index file> -o <index output file> -s <stats file>
            std::cerr << "Too few arguments!\n\n";
            this->showIndexTypeUsage(indexType,argv[0]);
            exit(1);
        }
        this->createROADAssociationDirectory(objSetFilePath,binaryIndexFilePath,outputFilePath,statsFilePath);
    } else if (indexType == constants::OBJ_IDX_QUADTREE) {
        if (argc < 13) {
            // Arguments: -t type -i <object set file> -b <binary index file> -l <max quadtree leaf size> -o <index output file> -s <stats file>
            std::cerr << "Too few arguments!\n\n";
            this->showIndexTypeUsage(indexType,argv[0]);
            exit(1);
        }

        if (maxQuadTreeLeafSize <= 0) {
            std::cerr << "Invalid maximum quadtree leaf size!\n";
            this->showIndexTypeUsage(indexType,argv[0]);
            exit(1);
        }
        this->createQuadtreeObjectHierarchy(objSetFilePath,binaryIndexFilePath,outputFilePath,statsFilePath,maxQuadTreeLeafSize);
    } else if (indexType == constants::OBJ_IDX_RTREE) {
        if (argc < 13) {
            // Arguments: -t type -i <object set file> -b <binary index file> -f <branch factor> -o <index output file> -s <stats file>
            std::cerr << "Too few arguments!\n\n";
            this->showIndexTypeUsage(indexType,argv[0]);
            exit(1);
        }

        if (branchFactor < 2) {
            std::cerr << "Invalid branch factor!\n";
            this->showIndexTypeUsage(indexType,argv[0]);
            exit(1);
        }
        this->createRtree(objSetFilePath,binaryIndexFilePath,outputFilePath,statsFilePath,branchFactor);
    } else {
        std::cerr << "Invalid set type!\n\n";
        this->showCommandUsage(argv[0]);
        exit(1);
    }
    
}

void ObjectIndexCommand::showCommandUsage(std::string programName)
{
    std::cerr << "Usage: " << programName << " -c " + constants::OBJ_IDX_CMD + " -t <index type>\n\n"
              << "Index Type Options:\n"
              << utility::getFormattedUsageString(constants::OBJ_IDX_GTREE,"Create object index for Gtree queries") << "\n"
              << utility::getFormattedUsageString(constants::OBJ_IDX_ROAD,"Create association directory for ROAD queries") << "\n"
              << utility::getFormattedUsageString(constants::OBJ_IDX_QUADTREE,"Create quadtree on object set for SILC queries") << "\n"
              << utility::getFormattedUsageString(constants::OBJ_IDX_RTREE,"Create Rtree on object set for IER queries") << "\n";
}

void ObjectIndexCommand::showIndexTypeUsage(std::string indexType, std::string programName)
{
    if (indexType == constants::OBJ_IDX_GTREE) {
        std::cerr << "Usage: " << programName << " -c " + constants::OBJ_IDX_CMD
                  << " -t " + constants::OBJ_IDX_GTREE + " -i <object set file>\n"
                  << "-b <binary gtree index file> -o <object index output file> -s <stats file>\n";
    } else if (indexType == constants::OBJ_IDX_ROAD) {
        std::cerr << "Usage: " << programName << " -c " + constants::OBJ_IDX_CMD
                  << " -t " + constants::OBJ_IDX_ROAD + " -i <object set file>\n"
                  << "-b <binary road index file> -o <object index output file> -s <stats file>\n";
    } else if (indexType == constants::OBJ_IDX_QUADTREE) {
        std::cerr << "Usage: " << programName << " -c " + constants::OBJ_IDX_CMD
                  << " -t " + constants::OBJ_IDX_QUADTREE + " -i <object set file>\n"
                  << "-b <binary graph file> -l <max quadtree leaf size> -o <object index output file>\n"
                  << " -s <stats file>\n";
    } else if (indexType == constants::OBJ_IDX_RTREE) {
        std::cerr << "Usage: " << programName << " -c " + constants::OBJ_IDX_CMD
                  << " -t " + constants::OBJ_IDX_RTREE + " -i <object set file>\n"
                  << "-b <binary graph file> -f <branch factor> -o <object index output file>\n"
                  << " -s <stats file>\n";
    } else {
        std::cerr << "Invalid index type!" << std::endl;
        this->showCommandUsage(programName);
    }
}

void ObjectIndexCommand::createGtreeOccurenceList(std::string objSetFilePath, std::string binaryIndexFilePath, 
                                                  std::string outputFilePath, std::string statsFilePath)
{
    /*
     * Create Gtree Occurence List
     */
    
    StopWatch sw;
    std::string setType;
    double setDensity;
    int setSize, setVariable;

    // Load index and object set file
    Gtree gtree = serialization::getIndexFromBinaryFile<Gtree>(binaryIndexFilePath);
    std::vector<NodeID> objectNodes = utility::getPointSetFromFile(objSetFilePath,setType,setDensity,setSize,setVariable);
        
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
    
    /*
     * Serialize to Binary File
     */
    serialization::outputIndexToBinaryFile<OccurenceList>(occList,outputFilePath);

    /*
     * Collect Stats and Output
     */
    double processingTimeMs = sw.getTimeMs();   
    double memoryUsage = occList.computeIndexSize();
    
    ObjectIndexTuple stats(gtree.getNetworkName(),gtree.getNumNodes(),gtree.getNumEdges(),
                    constants::OBJ_IDX_GTREE,processingTimeMs,memoryUsage,setType,setDensity,setVariable,setSize);

    std::cout << stats.getMultilineTupleString();

    this->outputCommandStats(statsFilePath,stats.getTupleString());
    
    std::cout << "Gtree Occurence List successfully created!" << std::endl;
}

void ObjectIndexCommand::createROADAssociationDirectory(std::string objSetFilePath, std::string binaryIndexFilePath, 
                                                        std::string outputFilePath, std::string statsFilePath)
{
    /*
     * Create ROAD Associate Directory
     */

    StopWatch sw;
    std::string setType;
    double setDensity;
    int setSize, setVariable;

    // Load index and object set file
    ROAD road = serialization::getIndexFromBinaryFile<ROAD>(binaryIndexFilePath);
    std::vector<NodeID> objectNodes = utility::getPointSetFromFile(objSetFilePath,setType,setDensity,setSize,setVariable);
    
    sw.start();
    AssociationDirectory assocDir(setType,setDensity,setVariable,setSize,road.getRnetTreeSize());
    for (auto objIt = objectNodes.begin(); objIt != objectNodes.end(); ++objIt) {
        // Find the ROAD leaf Rnet for this object and add entry into association
        // directory entry for that leaf then propagant up Rnet hierarchy
        assocDir.addObject(*objIt);
//         std::vector<int> leafRnets = road.getLeafRnets(*objIt);
        
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
    
    /*
     * Serialize to Binary File
     */    
    serialization::outputIndexToBinaryFile<AssociationDirectory>(assocDir,outputFilePath);

    /*
     * Collect Stats and Output
     */
    double processingTimeMs = sw.getTimeMs();   
    double memoryUsage = assocDir.computeIndexSize();
    
    ObjectIndexTuple stats(road.getNetworkName(),road.getNumNodes(),road.getNumEdges(),
                    constants::OBJ_IDX_ROAD,processingTimeMs,memoryUsage,setType,setDensity,setVariable,setSize);

    std::cout << stats.getMultilineTupleString();

    this->outputCommandStats(statsFilePath,stats.getTupleString());
    
    std::cout << "ROAD Association Directory successfully created!" << std::endl;
}

void ObjectIndexCommand::createQuadtreeObjectHierarchy(std::string objSetFilePath, std::string bgrFilePath, 
                                                       std::string outputFilePath, std::string statsFilePath, 
                                                       int maxQuadTreeLeafSize) {
    /*
     * Create Quadtree Object Hierarchy
     */

    StopWatch sw;
    std::string setType;
    double setDensity;
    int setSize, setVariable;
        
    // Load index and object set file
    Graph graph = serialization::getIndexFromBinaryFile<Graph>(bgrFilePath);
    std::vector<NodeID> objectNodes = utility::getPointSetFromFile(objSetFilePath,setType,setDensity,setSize,setVariable);
    
    sw.start();
    SimpleQuadtree qt(maxQuadTreeLeafSize);
    qt.buildFromPointSet(graph,objectNodes,setType,setDensity,setVariable,setSize);
    sw.stop();
    
    /*
     * Serialize to Binary File
     */    
    serialization::outputIndexToBinaryFile<SimpleQuadtree>(qt,outputFilePath);

    /*
     * Collect Stats and Output
     */
    double processingTimeMs = sw.getTimeMs();   
    double memoryUsage = qt.computeIndexSize();
    
    ObjectIndexTuple stats(graph.getNetworkName(),graph.getNumNodes(),graph.getNumEdges(),
                    constants::OBJ_IDX_QUADTREE,processingTimeMs,memoryUsage,
                    setType,setDensity,setVariable,setSize);
    stats.addSupplementaryFields("max_leaf_size",std::to_string(maxQuadTreeLeafSize));

    std::cout << stats.getMultilineTupleString();

    this->outputCommandStats(statsFilePath,stats.getTupleString());
    
    std::cout << "Quadtree Object Hierarchy successfully created!" << std::endl;
}

void ObjectIndexCommand::createRtree(std::string objSetFilePath, std::string bgrFilePath, 
                                     std::string outputFilePath, std::string statsFilePath, 
                                     int branchFactor)
{
    /*
     * Create Rtree
     */

    StopWatch sw;
    std::string setType;
    double setDensity;
    int setSize, setVariable;
        
    // Load index and object set file
    Graph graph = serialization::getIndexFromBinaryFile<Graph>(bgrFilePath);
    std::vector<NodeID> objectNodes = utility::getPointSetFromFile(objSetFilePath,setType,setDensity,setSize,setVariable);
    
    sw.start();
    // Note: We include time to build rtree input vector because it 
    // requires retrieval of coordinates from graph (Quadtree does this inside itself)
    std::vector<CoordinatePair> objectCoords;
    for (std::size_t i = 0; i < objectNodes.size(); ++i) {
        CoordinatePair objectCoordPair;
        graph.getCoordinates(objectNodes[i],objectCoordPair.first,objectCoordPair.second);
        objectCoords.push_back(objectCoordPair);
    }    
    StaticRtree rtree(branchFactor,setType,setDensity,setVariable,setSize);
    rtree.bulkLoad(objectNodes,objectCoords);
    sw.stop();
    
    /*
     * Serialize to Binary File
     */    
    serialization::outputIndexToBinaryFile<StaticRtree>(rtree,outputFilePath);

    /*
     * Collect Stats and Output
     */
    double processingTimeMs = sw.getTimeMs();   
    double memoryUsage = rtree.computeIndexSize();
    
    ObjectIndexTuple stats(graph.getNetworkName(),graph.getNumNodes(),graph.getNumEdges(),
                    constants::OBJ_IDX_RTREE,processingTimeMs,memoryUsage,
                    setType,setDensity,setVariable,setSize);
    stats.addSupplementaryFields("branch_factor",std::to_string(branchFactor));

    std::cout << stats.getMultilineTupleString();

    this->outputCommandStats(statsFilePath,stats.getTupleString());
    
    std::cout << "Rtree successfully created!" << std::endl;
}
