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

#include <iostream>
#include <unistd.h>
#include <getopt.h>
#include <vector>
#include <string>

#include "utility/utility.h"
#include "common.h"
#include "command/Command.h"
#include "command/GraphCommand.h"
#include "command/GenerateSampleSetCommand.h"
#include "command/KnnQueryCommand.h"
#include "command/ObjectIndexCommand.h"
#include "command/DynamicGraphCommand.h"
#include "command/TransformInputCommand.h"
#include "command/ExperimentsCommand.h"
#include "command/RoadNetworkIndexCommand.h"

void showUsage(std::string programName);
Command* getCommand(std::string);

int main(int argc, char* argv[])
{
    if (argc < 3) {
        // Arguments: -c <command>
        if (argc == 2) {
            std::cerr << "Too few arguments!\n\n";
        }
        showUsage(argv[0]);
        exit(1);
    }
    
    /*
     * Determine Command to Execute
     */ 
    Command* command;
    int opt = getopt(argc,argv, "c:");
    switch (opt) {
        case 'c':
            command = getCommand(optarg);
            break;
        default:
            std::cerr << "Unknown option provided!\n\n";
            showUsage(argv[0]);
            exit(1);
    }
    
    if (command != NULL) {
        if (argc == 3) {
            //  Commands are expected to have more that 2 options/arguments
            command->showCommandUsage(argv[0]);
        } else {
            command->execute(argc, argv);
        }
        delete command;
    } else {
        std::cerr << "Command not found!\n\n";
        showUsage(argv[0]);
        exit(1);
    }
    
}

void showUsage(std::string programName) {
    std::cerr << "Usage: " << programName << " -c <command>\n\n"
              << "Index Commands:\n"
              << utility::getFormattedUsageString(constants::IDX_GRAPH_CMD,"Fast in-memory graph") + "\n"
              << utility::getFormattedUsageString(constants::IDX_DYNAMICGRAPH_CMD,"Slower but updateable graph") + "\n"
              << utility::getFormattedUsageString(constants::NETWORK_IDX_CMD,"Other types of graph indexes") + "\n"
              << utility::getFormattedUsageString(constants::OBJ_IDX_CMD,"Object indexes for various graph indexes") + "\n\n"
              << "Query Commands:\n"
              << utility::getFormattedUsageString(constants::QUERY_KNN_CMD,"kNN query") + "\n"
              << "Other Commands:\n"
              << utility::getFormattedUsageString(constants::SAMPLE_SET_CMD,"Generate different types of object sets") + "\n"
              << utility::getFormattedUsageString(constants::TRANSFORM_INPUT_CMD,"Process graphs from known input data sources") + "\n"
              << utility::getFormattedUsageString(constants::EXPERIMENTS_CMD,"Run knn experiments") + "\n\n"
              << "Note: Execute commands without parameters to view usage and all options\n";
}

Command* getCommand(std::string commandName) {
    /*
     * Determine Command
     */ 
    Command* command = NULL;
    if (commandName == constants::IDX_GRAPH_CMD) {
        command = new GraphCommand();
    } else if (commandName == constants::IDX_DYNAMICGRAPH_CMD) {
        command = new DynamicGraphCommand();
    } else if (commandName == constants::QUERY_KNN_CMD) {
        command = new KnnQueryCommand();
    } else if (commandName == constants::SAMPLE_SET_CMD) {
        command = new GenerateSampleSetCommand();
    } else if (commandName == constants::OBJ_IDX_CMD) {
        command = new ObjectIndexCommand();
    } else if (commandName == constants::TRANSFORM_INPUT_CMD) {
        command = new TransformInputCommand();
    } else if (commandName == constants::EXPERIMENTS_CMD) {
        command = new ExperimentsCommand();
    } else if (commandName == constants::NETWORK_IDX_CMD) {
        command = new RoadNetworkIndexCommand();
    }
    return command;
}