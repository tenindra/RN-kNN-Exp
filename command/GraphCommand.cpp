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

#include "GraphCommand.h"

#include "../processing/Graph.h"
#include "../tuple/IndexTuple.h"
#include "../common.h"
#include "../utility/StopWatch.h"
#include "../utility/serialization.h"
#include "../processing/Junction.h"

#include <fstream>
#include <boost/archive/binary_oarchive.hpp>

void GraphCommand::execute(int argc, char* argv[]) {
    std::string coFilePath = "";
    std::string grFilePath = "";
    std::string indexOutputFilePath = "";
    std::string statsOutputFilePath = "";

    /*
     * Process Command Line Arguments
     */
    if (argc < 9) {
        // Arguments: -e <graph file> -n <coordinates file> -o <output file> -s <stats file>
        std::cerr << "Too few arguments!\n\n";
        this->showCommandUsage(argv[0]);
        return;
    }
    
    int opt;
    while ((opt = getopt (argc, argv, "e:n:o:s:")) != -1) {
        switch (opt) {
            case 'e':
                grFilePath = optarg;
                break;
            case 'n':
                coFilePath = optarg;
                break;
            case 'o':
                indexOutputFilePath = optarg;
                break;
            case 's':
                statsOutputFilePath = optarg;
                break;
            default:
                std::cerr << "Unknown option(s) provided!\n\n";
                showCommandUsage(argv[0]);
                exit(1);
        }
    }    
    
    if (grFilePath == "" || indexOutputFilePath == "" 
        || coFilePath == "" || statsOutputFilePath == "") {
        std::cerr << "Invalid argument(s)!\n\n";
        this->showCommandUsage(argv[0]);
        exit(1);
    }
    
    /*
     * Parse Graph from Text File and Serialize to Binary File
     */
    Graph graph;
    StopWatch sw;
    
    sw.start();
    graph.parseGraphFile(grFilePath,coFilePath);
    sw.stop();
    double processingTimeMs = sw.getTimeMs();    

    // Print Additional Graph Statistics
    graph.printGraphStats();

    /*
     * Serialize to Binary File
     */    
    serialization::outputIndexToBinaryFile<Graph>(graph,indexOutputFilePath);
    
    /*
     * Collect Stats and Output
     */
    int nodes = graph.getNumNodes();
    int edges = graph.getNumEdges();
//     double memoryUsage = graph.computeIndexSize();
    double memoryUsage = graph.computeINEIndexSize(); // For experiments we report INE size

    IndexTuple stats(graph.getNetworkName(),nodes,edges,constants::IDX_GRAPH_CMD,processingTimeMs,memoryUsage);

    std::cout << stats.getMultilineTupleString();

    this->outputCommandStats(statsOutputFilePath,stats.getTupleString());
    
    std::cout << "Binary graph index successfully created!" << std::endl;
    
    /*
     * Supplementary Junction Index
     */
    
    sw.reset();
    sw.start();
    Junction junc(graph.getNetworkName(),graph.getNumNodes(),graph.getNumEdges());
    junc.buildJunction(graph);
    sw.stop();

    // Print Additional Graph Statistics
    //junc.printGraphStats();

    /*
     * Serialize to Binary File
     */
    int extIndex = indexOutputFilePath.find_last_of(".");
    std::string idxFileName = indexOutputFilePath.substr(0,extIndex)+"."+constants::IDX_JUNC_CMD;
    serialization::outputIndexToBinaryFile<Junction>(junc,idxFileName);

    /*
     * Collect Stats and Output
     */
    processingTimeMs = sw.getTimeMs();
    memoryUsage = junc.computeIndexSize();
    
    IndexTuple stats2(junc.getNetworkName(),junc.getNumNodes(),junc.getNumEdges(),constants::IDX_JUNC_CMD,processingTimeMs,memoryUsage);

    //this->outputCommandStats(statsOutputFilePath,stats2.getTupleString()); // No need to output junction stats
    
    std::cout << "Junction index successfully created!" << std::endl;
 
    std::cout << stats2.getMultilineTupleString();
}

void GraphCommand::showCommandUsage(std::string programName) {
    std::cerr << "Usage: " << programName << " -c " + constants::IDX_GRAPH_CMD + " -e <text graph file>"
              << " -n <coordinates file>\n-o <binary graph output file> -s <stats file>\n";
}
