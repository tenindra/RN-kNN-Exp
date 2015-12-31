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

#ifndef _KNNQUERYCOMMAND_H
#define _KNNQUERYCOMMAND_H

#include "Command.h"
#include "../tuple/KnnQueryTuple.h"
#include "../common.h"

#include <vector>

class KnnQueryCommand: public Command {

    public:
        void execute(int argc, char* argv[]);
        void showCommandUsage(std::string programName);
        void showMethodUsage(std::string method, std::string programName);
        
    private:
        std::vector<KnnQueryTuple> executeINEQuery(int k, std::vector<NodeID>& queryNodesIDs, 
                                                   std::string objSetFilePath, std::string bgrFilePath);
        std::vector<KnnQueryTuple> executeGtreeQuery(int k, std::vector<NodeID>& queryNodesIDs, 
                                                     std::string idxFilePath, std::string objIdxFilePath,
                                                     std::string bgrFilePath);
        std::vector<KnnQueryTuple> executeROADQuery(int k, std::vector<NodeID>& queryNodesIDs, 
                                                    std::string idxFilePath, std::string objIdxFilePath);
        std::vector<KnnQueryTuple> executeSILCQuery(std::string method, int k, std::vector<NodeID>& queryNodesIDs, 
                                                    std::string idxFilePath, std::string objIdxFilePath,
                                                    std::string bgrFilePath, std::string addIdxFile = "", 
                                                    std::string objSetFilePath = "");
        std::vector<KnnQueryTuple> executeIERQuery(std::string method, int k, std::vector<NodeID>& queryNodesIDs, 
                                                   std::string objIdxFilePath, std::string bgrFilePath, 
                                                   std::string idxFilePath, std::string addIdxFile = "");
};

#endif // _KNNQUERYCOMMAND_H
