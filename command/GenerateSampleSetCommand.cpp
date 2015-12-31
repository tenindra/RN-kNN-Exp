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

#include "GenerateSampleSetCommand.h"

#include "../processing/Graph.h"
#include "../processing/SetGenerator.h"
#include "../processing/StaticRtree.h"
#include "../utility/StopWatch.h"
#include "../utility/utility.h"
#include "../utility/serialization.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <boost/algorithm/string/iter_find.hpp>
#include <boost/algorithm/string/finder.hpp>
#include <boost/algorithm/string.hpp>
#include <set>

void GenerateSampleSetCommand::execute(int argc, char* argv[])
{
    std::string setType = "";
    double setDensity = 0;
    int setVariable = 0;
    int setSize = 0;
    std::string bgrFilePath = "";
    std::string outputFilePath = "";
    std::string status = "ignore";
    std::string poiFilePath = "";
    int outdegree = 0;
    int numPartitions = 0;
    double chanceHD = 0;
    double chanceLD = 0;
    double highDensity = 0;
    double lowDensity = 0;

    /*
     * Process Command Line Arguments
     */
    int opt;
    while ((opt = getopt (argc, argv, "t:i:o:d:n:s:e:p:v:f:c:h:l:z:")) != -1) {
        switch (opt) {
            case 't':
                setType = optarg;
                break;
            case 'i':
                bgrFilePath = optarg;
                break;
            case 'o':
                outputFilePath = optarg;
                break;
            case 'd':
                setDensity = std::stod(optarg);
                break;
            case 'v':
                setVariable = std::stoi(optarg);
                break;
            case 'n':
                setSize = std::stoi(optarg);
                break;
            case 's':
                status = optarg;
                break;
            case 'e':
                outdegree = std::stoi(optarg);
                break;
            case 'p':
                numPartitions = std::stoi(optarg);
                break;
            case 'f':
                poiFilePath = optarg;
                break;
            case 'c':
                chanceHD = std::stod(optarg);
                break;
            case 'z':
                chanceLD = std::stod(optarg);
                break;
            case 'h':
                highDensity = std::stod(optarg);
                break;
            case 'l':
                lowDensity = std::stod(optarg);
                break;
            default:
                std::cerr << "Unknown option(s) provided!\n\n";
                showCommandUsage(argv[0]);
                exit(1);
        }
    }

    // Validate Command Line Arguments
    if (argc == 5 && setType != "") {
        // This is 5 so that user can just enter method to find out what parameters are required for method
        this->showSetTypeUsage(setType,argv[0]);
        exit(1);
    }    
    
    if (setType == "" || bgrFilePath == "" || outputFilePath == "") {
        std::cerr << "Invalid argument(s)!\n\n";
        this->showCommandUsage(argv[0]);
        exit(1);
    }

    /*
     * Serialize Graph from Binary File
     */
    Graph graph = serialization::getIndexFromBinaryFile<Graph>(bgrFilePath);

    /*
     * Sample Graph to Generate Object Set
     */          
    int numNodes = graph.getNumNodes();
    SetGenerator sg;
    StopWatch sw;

    if (setType == constants::RAND_OBJ_SET) {
        /*
         * Create Randomly Selected Object Set
         */
        
        // Process Set Type Command-Line Arguments
        if (argc < 11) {
            // Arguments: -c sample_set -t random -i <binary graph file> -o <output file> 
            // (at least one of -d <density> or -n <size>) [-s <ignore|include|exclude> [-e <outdegree of nodes to include or exclude>]] [-v <set variable - currently not used>]
            std::cerr << "Too few arguments!\n\n";
            this->showSetTypeUsage(setType,argv[0]);
            exit(1);
        } else if ((setDensity <= 0 || setDensity >= 1) && setSize <= 0) {
            std::cerr << "Invalid argument(s)!\n\n";
            this->showSetTypeUsage(setType,argv[0]);
            exit(1);
        }
        
        // Determine how many objects must be generated if density was provided
        if (setSize == 0) {
            setSize = std::ceil(numNodes*setDensity);
        } else {
            setDensity = setSize/numNodes;
        }

        // Generat set and write to file
        sw.start();
        std::vector<NodeID> sampleSet;
        std::string message = "";
        if (status == "include") {
            if (outdegree <= 0) {
                std::cerr << "Invalid out degree, must be greater than zero\n\n";
                this->showSetTypeUsage(setType,argv[0]);
                exit(1);            
            }
            sampleSet = sg.generateRandomSampleSet(graph,setSize,1,outdegree);
            message = " including only nodes with degree " + std::to_string(outdegree);
        } else if (status == "exclude") {
            if (outdegree <= 0) {
                std::cerr << "Invalid out degree, must be greater than zero\n\n";
                this->showSetTypeUsage(setType,argv[0]);
                exit(1);            
            }
            sampleSet = sg.generateRandomSampleSet(graph,setSize,0,outdegree);
            message = " excluding all nodes with degree " + std::to_string(outdegree);
        } else if (status == "ignore") {
            sampleSet = sg.generateRandomSampleSet(graph,setSize,-1);
            message = " including nodes with any degree";
        } else {
            std::cerr << "-s option argument does not have a valid value, it must be either \"include\", \"exclude\" or \"ignore\" if provided at all\n\n";
            this->showSetTypeUsage(setType,argv[0]);
            exit(1);            
        }
        sw.stop();
        utility::writeSampleSet(outputFilePath,graph.getNetworkName(),setType,setDensity,setSize,setVariable,sampleSet);
        std::cout << "Random sample set of size " << sampleSet.size() << message << " successfully created in " << sw.getTimeMs() << "ms!" << std::endl;
        
    } else if (setType == constants::PARTITION_OBJ_SET) {
        /*
         * Create Object Set Sampling Each Object from Separate Graph Partition
         */
        int clusterSize = setVariable;
        
        // Process Set Type Command-Line Arguments
        if (argc < 13) {
            // Arguments: -c sample_set -t partition -i <binary graph file> -o <output file> 
            // -p <number of partitions> -h <high density> -l <low density> -c <proportion of high density partitions> 
            // -z <proportion of low density partitions from remaining non-high density partitions>
            // [-n <object set size limit> or -d <object set density limit>]
            std::cerr << "Too few arguments!\n\n";
            this->showSetTypeUsage(setType,argv[0]);
            exit(1);
        } else if (numPartitions == 0 || highDensity <= 0 || lowDensity <= 0 || chanceHD <= 0 || chanceLD <= 0) {
            std::cerr << "Invalid argument(s)!\n\n";
            this->showSetTypeUsage(setType,argv[0]);
            exit(1);
        }
        
        int numObjects = 0;
        if (setSize != 0) {
            numObjects = setSize;
        } else if (setDensity != 0) {
            numObjects = std::ceil(graph.getNumNodes()*setDensity);
        }
        
        // Generat set and write to file
        sw.start();
        std::vector<std::unordered_set<NodeID>> partitionsNodeSets;
        std::vector<std::vector<NodeID>> partitionsNodes;
        sg.generatePartitions(graph,numPartitions,partitionsNodeSets,partitionsNodes);
        std::vector<NodeID> sampleSet = sg.generatePartitionSampleSet(graph,chanceHD,chanceLD,highDensity,lowDensity,partitionsNodeSets,partitionsNodes,numObjects);
        sw.stop();
        utility::writeSampleSet(outputFilePath,graph.getNetworkName(),setType,setDensity,setSize,clusterSize,sampleSet);
        std::cout << "Sample set of size " << sampleSet.size() << " successfully created using " << numPartitions 
                  << " partitions and a cluster size of " << clusterSize << " in " << sw.getTimeMs() << "ms!" << std::endl;
        
    } else if (setType == constants::RAND_PAIRS_SET) {
        /*
         * Create Randomly Selected Pair Set
         */

        // Process Set Type Command-Line Arguments
        if (argc < 11) {
            // Arguments: -c sample_set -t random_pairs -i <binary graph file> -o 
            // <output file> -n <num pairs>
            std::cerr << "Too few arguments!\n\n";
            this->showSetTypeUsage(setType,argv[0]);
            exit(1);
        } else if (setSize <= 0) {
            std::cerr << "Invalid argument(s)!\n\n";
            this->showSetTypeUsage(setType,argv[0]);
            exit(1);
        }
        
        // Generat pair set and write to file
        sw.start();
        std::vector<NodePair> pairSet = sg.generateRandomPairSet(graph,setSize);
        sw.stop();
        utility::writePairSet(outputFilePath,graph.getNetworkName(),setType,setSize,pairSet);
        std::cout << "Sample set with " << pairSet.size() << " pairs of nodes successfully created in " << sw.getTimeMs() << "ms!" << std::endl;
        
    } else if (setType == constants::POI_OBJECT_SET) {
        /*
         * Create Object Sets from Formatted POI File
         */

        // Process Set Type Command-Line Arguments
        if (argc < 11) {
            // Arguments: -c sample_set -i <binary graph file> -f <poi file> -o <output file path> 
            std::cerr << "Too few arguments!\n\n";
            this->showSetTypeUsage(setType,argv[0]);
            exit(1);
        } else if (poiFilePath == "") {
            std::cerr << "Invalid argument, no POI file not provided!\n\n";
            this->showSetTypeUsage(setType,argv[0]);
            exit(1);
        }
        std::unordered_map<std::string,std::vector<Point>> objectSets;
        
        std::string line, type, id, name, countryCode, alternateNames, xCoord, yCoord, tags;
        std::ifstream poiFile(poiFilePath, std::ios::in);
        
        if (poiFile.is_open()) {
            
            // gisgraphy POI file
            while(std::getline(poiFile, line))
            {
                std::stringstream linestream(line);
                std::getline(linestream, type, '\t');
                std::getline(linestream, id, '\t');
                std::getline(linestream, name, '\t');
                std::getline(linestream, countryCode, '\t');
                std::getline(linestream, alternateNames, '\t');
                std::getline(linestream, xCoord, '\t');
                std::getline(linestream, yCoord, '\t');
                std::getline(linestream, tags, '\t');
//                 std::cout << "Type: " << type << "\n";
//                 std::cout << "ID: " << id << "\n";
//                 std::cout << "Name: " << name << "\n";
//                 std::cout << "Country Code: " << countryCode << "\n";
//                 std::cout << "Alternative Names: " << alternateNames << "\n";
//                 std::cout << "Point: (" << xCoord << "," << yCoord << ")\n";
//                 std::cout << "Tags: " << tags << "\n";
                
                std::vector<std::string> setNames;
                boost::iter_split(setNames, tags, boost::first_finder("__", boost::is_iequal()));
                
                for (std::string setName: setNames) {
                    std::replace(setName.begin(), setName.end(), '_', ' ');
                    boost::trim(setName);
                    std::replace(setName.begin(), setName.end(), ' ', '_');
                    boost::to_lower(setName);
                    if (setName != "") {
                        if (objectSets.find(setName) == objectSets.end()) {
                            objectSets.emplace(setName,std::vector<Point>());
                        }
                        objectSets[setName].push_back(Point(std::stod(xCoord)*1000000,std::stod(yCoord)*1000000));
                    }
                }
            }

//             // OsmPoiPbf output file
//             std::unordered_map<std::string,std::string> poiTypeNames {
//                 {"1","Alcohol"}, {"2","Eatery"}, {"3","School"}, {"4","College"}, 
//                 {"5","Medical_Centre"}, {"6","Dentist"}, {"7","Pharmacy"}, {"8","Hospital"}, 
//                 {"9","Cinema"}, {"10","Theatre"}, {"11","Courthouse"}, {"12","Fire_Station"}, 
//                 {"13","Police"}, {"14","Post Office"}, {"15","Library"}, {"16","Fuel"}, 
//                 {"17","Bus_Stop"}, {"18","Bank"}, {"19","Place_of_Worship"}, {"20","Parking"}, 
//                 {"21","University"}, {"23","Convenience_Store"}, {"25","Supermarket"}, 
//                 {"26","Caravan_Site"}, {"27","Accommodation"}, {"28","Park"},
//                 {"29","Sports_Ground"}, {"30","Swimming_Pool"}, {"32","Aerodrome"}
//             };
// 
//             sw.start();
//             while(std::getline(poiFile, line))
//             {
//                 std::stringstream linestream(line);
//                 std::getline(linestream, type, '|');
//                 std::getline(linestream, id, '|');
//                 std::getline(linestream, yCoord, '|');
//                 std::getline(linestream, xCoord, '|');
//                 std::getline(linestream, name, '|');
// //                 std::cout << "Type: " << type << "\n";
// //                 std::cout << "ID: " << id << "\n";
// //                 std::cout << "Point: (" << xCoord << "," << yCoord << ")\n";
// //                 std::cout << "Name: " << name << "\n";
//                 
//                 if (poiTypeNames.find(type) != poiTypeNames.end()) {
//                     std::string setName = poiTypeNames[type];
//                     if (objectSets.find(setName) == objectSets.end()) {
//                         objectSets.emplace(setName,std::vector<Point>());
//                     }
//                     objectSets[setName].push_back(Point(std::stod(xCoord)*1000000,std::stod(yCoord)*1000000));
//                 }                
//             }
            
            std::set<std::pair<int,std::string>> objectSetsBySize;
            for (auto it = objectSets.begin(); it != objectSets.end(); ++it) {
                objectSetsBySize.insert(std::make_pair(it->second.size(),it->first));
            }
            
            std::vector<NodeID> graphNodes = graph.getNodesIDsVector();
            StaticRtree rtree(8,"",1,1,graph.getNumNodes());
            rtree.bulkLoad(graphNodes,graph.coordinates);
            std::unordered_set<NodeID> inSet;
            std::vector<NodeID> sampleSet;
            std::vector<NodeID> kNNs;
            std::vector<EuclideanDistanceSqrd> kNNDistSqrd;
            NodeID locationID;
            EdgeWeight dist;
    
            for (auto it = objectSetsBySize.rbegin(); it != objectSetsBySize.rend(); ++it) {
                if (it->first > 2000) {
                    sampleSet.clear();
                    inSet.clear();
                    for (auto innerIt = objectSets[it->second].begin(); innerIt != objectSets[it->second].end(); ++innerIt) {            
                        kNNs.clear();
                        kNNDistSqrd.clear();
                        rtree.getKNNs(1,innerIt->x,innerIt->y,kNNs,kNNDistSqrd);
                        locationID = kNNs[0];
                        dist = std::ceil(std::sqrt(kNNDistSqrd[0]));
                        if (dist < 4000 && inSet.find(locationID) == inSet.end()) {
                            sampleSet.push_back(locationID);
                            inSet.insert(locationID);
                        }
                    }
                    std::cout << it->second << ": Original Size = " << it->first << ", Projected Size = " << sampleSet.size() << "\n";
                    setType = "real";
                    setSize = sampleSet.size();
                    setDensity = setSize/graphNodes.size();
                    utility::writeSampleSet(outputFilePath+"/"+it->second+".txt",graph.getNetworkName(),setType,setDensity,setSize,1,sampleSet);
                }
            }
            sw.stop();
            
            std::cout << "Sample set with " << sampleSet.size() << " nodes successfully created in " << sw.getTimeMs() << "ms!" << std::endl;
        
        } else {
            std::cerr << "Cannot open POI file " << poiFilePath << std::endl;
            exit(1);
        }    
        poiFile.close();
        
    } else if (setType == constants::CLUSTER_OBJ_SET) {
        /*
         * Create Object Set Sampling Randomly and Growing Clusters
         */
        int clusterSize = setVariable;
        
        // Process Set Type Command-Line Arguments
        if (argc < 13) {
            // Arguments: -c sample_set -t partition -i <binary graph file> -o <output file> 
            // (at least one of -d <density> or -p <number of clusters> or -n <size of object set>) -v <size of each cluster> -h <cluster probability>
            std::cerr << "Too few arguments!\n\n";
            this->showSetTypeUsage(setType,argv[0]);
            exit(1);
        } else if ((numPartitions == 0 && setDensity == 0 && setSize == 0) || clusterSize <= 0 || chanceHD <= 0) {
            std::cerr << "Invalid argument(s)!\n\n";
            this->showSetTypeUsage(setType,argv[0]);
            exit(1);
        }
        
        int numObjects = 0;
        if (setDensity != 0) {
            numObjects = std::ceil(graph.getNumNodes()*setDensity);
        } else if (setSize != 0) {
            numObjects = setSize;
        } else {
            numObjects = numPartitions*clusterSize;
        }
        
        // Generat set and write to file
        sw.start();
        std::vector<NodeID> sampleSet = sg.generateRandomClusteredSampleSet(graph,chanceHD,numObjects,clusterSize);
        sw.stop();
        utility::writeSampleSet(outputFilePath,graph.getNetworkName(),setType,setDensity,setSize,clusterSize,sampleSet);
        std::cout << "Sample set of size " << sampleSet.size() << " successfully created using " << numPartitions 
                  << " partitions and a cluster size of " << clusterSize << " in " << sw.getTimeMs() << "ms!" << std::endl;
        
    } else {
        std::cerr << "Invalid set type!\n\n";
        this->showCommandUsage(argv[0]);
        exit(1);
    }

}

void GenerateSampleSetCommand::showCommandUsage(std::string programName)
{
    std::cerr << "Usage: " << programName << " -c " + constants::SAMPLE_SET_CMD + " -t <set type>\n\n"
              << "Set Type Options:\n"
              << utility::getFormattedUsageString(constants::RAND_OBJ_SET,"Randomly sample graph nodes") << "\n"
              << utility::getFormattedUsageString(constants::PARTITION_OBJ_SET,"Sample from graph partitions") << "\n"
              << utility::getFormattedUsageString(constants::RAND_PAIRS_SET,"Randomly sample pairs of graph nodes") << "\n"
              << utility::getFormattedUsageString(constants::POI_OBJECT_SET,"Create object sets from formatted POI file") << "\n"
              << utility::getFormattedUsageString(constants::CLUSTER_OBJ_SET,"Randomly sample cluster centres and grow them") << "\n";
}

void GenerateSampleSetCommand::showSetTypeUsage(std::string setType, std::string programName)
{
    if (setType == constants::RAND_OBJ_SET) {
        std::cerr << "Usage: " << programName << " -c " + constants::SAMPLE_SET_CMD + 
                     " -t " +constants::RAND_OBJ_SET + " -i <binary graph file>"
                  << " -o <set output file> (one option and argument from group 1)\n "
                  << "[-s <ignore|include|exclude> [-e <outdegree of nodes to include or exclude>]]\n\n"
                  << "Group 1:\n"
                  << utility::getFormattedUsageString("-d <fraction>","Density of sampled set") << "\n"
                  << utility::getFormattedUsageString("-n <size>","Size of sampled set") << "\n";
    } else if (setType == constants::PARTITION_OBJ_SET) {
        std::cerr << "Usage: " << programName << " -c " + constants::SAMPLE_SET_CMD + 
                     " -t " +constants::PARTITION_OBJ_SET + " -i <binary graph file>"
                  << " -o <set output file> (one option and argument from group 1) -v <size of each cluster>\n"
                  << "-h <high density> -l <low density> -c <fraction of high-density partitions> \n"
                  << "-z <fraction lower-density partitions (out of those not high density)>\n"
                  << "[-n <object size limit> or -d <object set density limit>]\n\n"
                  << "Group 1:\n"
                  << utility::getFormattedUsageString("-p <parts>","Number of clusters") << "\n";
    } else if (setType == constants::RAND_PAIRS_SET) {
        std::cerr << "Usage: " << programName << " -c " + constants::SAMPLE_SET_CMD + 
                     " -t " + constants::RAND_PAIRS_SET + " -i <binary graph file>"
                  << " -o <set output file> -n <number of pairs>\n";
    } else if (setType == constants::POI_OBJECT_SET) {
        std::cerr << "Usage: " << programName << " -c " + constants::SAMPLE_SET_CMD + 
                     " -t " + constants::POI_OBJECT_SET + " -i <binary graph file>"
                  << " -o <output file path> -f <input POI file>\n";
    } else if (setType == constants::CLUSTER_OBJ_SET) {
        std::cerr << "Usage: " << programName << " -c " + constants::SAMPLE_SET_CMD + 
                     " -t " +constants::CLUSTER_OBJ_SET + " -i <binary graph file>"
                  << " -o <set output file> (one option and argument from group 1) -v <size of each cluster> -c <cluster probability>\n\n"
                  << "Group 1:\n"
                  << utility::getFormattedUsageString("-d <fraction>","Density of sampled set") << "\n"
                  << utility::getFormattedUsageString("-p <parts>","Number of clusters") << "\n"
                  << utility::getFormattedUsageString("-n <size>","Size of sampled set") << "\n";
    } else {
        this->showCommandUsage(programName);
    }
}
