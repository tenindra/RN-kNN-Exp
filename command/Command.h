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

#ifndef _COMMAND_H
#define _COMMAND_H

#include <string>
#include <unistd.h>
#include <getopt.h>
#include <iostream>
#include <fstream>

class Command {

    public:
        virtual void execute(int argc, char* argv[]) = 0;
        virtual void showCommandUsage(std::string programName) = 0;
        virtual ~Command() {};
        
    protected:
        void outputCommandStats(std::string statsFilePath, std::string tupleString) {
            std::ofstream statsFS(statsFilePath, std::ios::out | std::ios::app);
            if (statsFS.is_open()) {
                statsFS << tupleString << std::endl;
            } else {
                std::cerr << "Cannot open stats file " << statsFilePath << std::endl;
            }
        }
};

#endif // _COMMAND_H