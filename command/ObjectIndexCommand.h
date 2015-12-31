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

#ifndef _OBJECTINDEXCOMMAND_H
#define _OBJECTINDEXCOMMAND_H

#include "Command.h"

class ObjectIndexCommand: public Command {

    public:
        void execute(int argc, char* argv[]);
        void showCommandUsage(std::string programName);
        void showIndexTypeUsage(std::string index, std::string programName);
        void createGtreeOccurenceList(std::string objSetFilePath, 
                                      std::string binaryIndexFilePath, std::string outputFilePath,
                                      std::string statsFilePath);
        void createROADAssociationDirectory(std::string objSetFilePath, std::string binaryIndexFilePath, 
                                            std::string outputFilePath, std::string statsFilePath);
        void createQuadtreeObjectHierarchy(std::string objSetFilePath, std::string bgrFilePath, 
                                           std::string outputFilePath, std::string statsFilePath, 
                                           int maxQuadTreeLeafSize);
        void createRtree(std::string objSetFilePath, std::string bgrFilePath, 
                         std::string outputFilePath, std::string statsFilePath, 
                         int branchFactor);

};

#endif // _OBJECTINDEXCOMMAND_H