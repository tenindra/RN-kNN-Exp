# Introduction

Road Network kNN Experimental Evaluation


# Requirements


## Software

Before the executable can be compiled the following packages/libraries must be installed. Note that we were able to install METIS from the Ubuntu repositories if you are using that operating system.

- g++ v4.9 or higher
- Boost v1.57 or higher for serialization
- [METIS] (http://glaros.dtc.umn.edu/gkhome/metis/metis/download) version 5.1 or higher
- CMake v2.8 or higher
- gnuplot for figure generation
- Optional: [Google Sparsehash] (http://code.google.com/p/sparsehash/) for comparison of G-tree hashtables implementations

Note: We found Boost versions lower than 1.57 did not correctly support serialization of some STL data structures (such as std::unordered_map) thus they cannot be used here. This means you need to remove any existing Boost version, e.g. by uninstalling libboost-all-dev through your package manager.


## System & Hardware

We recommend a 64-bit OS due to the size of indexes (and we have not tested on 32-bit). Furthermore at least 36GB RAM is required to reproduce experiments on all datasets. Results can be reproduced using a smaller amount of RAM by omitting some of the larger datasets (this will be discussed later).


# Files

The C++ source code can be found in the top-level "rn_knn" directory, and in the subdirectories command, external, processing, queue, tuple and utility. All bash scripts for various tasks such as setting up 
the directories, running experiments, and generating figures are found in the scripts directory.


# Compilation

1. Change CMAKE_C_COMPILER and CMAKE_CXX_COMPILER in CMakeLists.txt to point to the correct g++ version (i.e. 4.9 or higher). CMakeLists.txt can be found in the "knn" directory.

2. Create a directory called "build" in the top-level (i.e. knn/) directory

```
mkdir build
```

3. Change into this directory and Generate makefiles using CMake

```
cd build
cmake -G "Unix Makefiles" ../CMakeLists.txt ..
```

4. Compile all executables using make

```
make
```


# Setup

1. Open the scripts/globalVariables file in your favourite editor

2. Change the "output_path" variable to the full-path of an *existing* location where you wish to store all data (e.g. indexes, objects sets, figures, etc...)

Note: This can be different to the path containing the code (recommended)

3. Change the "exe_path" to the full path of the build directory created above

4. Run `resetExperimentalSetup` in the scripts directory to create all necessary sub-directories

```
cd scripts
bash resetExperimentalSetup
```

5. Download the DIMACS distance edge-weight graph files (with extension .gr.gz and prefixed by USA-road-*d*) and coordinate files (with extension .co.gz) for DE, VT, ME and NW from https://github.com/tenindra/rn_knn_exp_data to the output_path/data/dimacs directory

6. Download the DIMACS distance edge-weight graph files and coordinate files for COL, CAL, E, W, CTR, USA from http://www.dis.uniroma1.it/challenge9/ to the output_path/data/dimacs directory

Note: You *must* rename the "CAL" in DIMACS graph and coordinate filenames to CALNV.

7. Download the node (with extension .cnode) and edge (with extension .cedge) files for North America (NA) from http://www.cs.utah.edu/~lifeifei/SpatialDataset to the output_path/data/tpq directory

Note: You *must* gzip each file so you have two files, NA.cnode.gz and NA.cedge.gz, in the output_path/data/tpq directory

8. Download the real_world_pois.tar.gz from https://github.com/tenindra/rn_knn_exp_data to the output_path/ directory

9. Unzip the real_world_pois.tar.gz archive (ensure this creates a directory at output_path/real_world_pois populated with subdirectories with POI sets for several road networks)


# Less Than 16GB RAM

Experiments can still be run with less than 16GB RAM, but for fewer datasets. To do this go to the globalVariables file and modify the road_networks array. The road networks are in size order, so simply remove the road networks from the end. E.g. to run experiments up to the Eastern US dataset change it to road_networks=("DE" "VT" "ME" "COL" "FLA" "CALNV" "E").

Note 1: All experiments on the US dataset will not be possible, of course. In this case it is best to comment out those experiments in the runPaperExperiments script.

Note 2: SILC (index used by Distance Browsing) can only be built upto the default Colorado dataset. SILC on the Colorado requires at least 8GB of memory. Less than this will cause all default experiments to be missing Distance Browsing comparisons.


# Running Distance Experiments

The follow instructions can be followed to re-create all figures from the paper. Assuming your are already in the scripts directory (otherwise cd into scripts):

1. Clean the DIMACS and TPQ datasets for errors and redundancy

```
bash transformInputData
```

2. Build the binary files for the basic graph representations

```
bash buildBinaryGraphs
```

3. Build all road network indexes (this may take a while - ~15 hours on our machine)

```
bash buildRoadNetworkIndexes
```

4. Generate query sets

```
bash generateQuerySets
```

5. Generate random object sets and build corresponding object indexes

```
bash buildObjectIndexes
```

Note: You can run `resetExperimentalSetup` again to remove all object indexes and query sets

6. Run all paper experiments and produce figures

```
bash runPaperExperiments
```

Note: The above command also creates figures, but if you wish to recreate figures without running experiments again (which can take sometime), you can use: 

```
bash createPaperFigures
```

Note: All above commands may be batched, just enter each separated by semi-colon ";"


# Running Travel Time Experiments

We also provided some additional experiments in the appendices on travel time graphs. These
can also be re-produced using the following procedure:

1. Create a new directory to store data for travel time experiments (e.g. figures etc...)

2. Change the "output_path" variable in globalVariables the full-path of this new location

3. Change the "edge_type" variable to "t"

Note: The "edge_type" variable must be changed back to "d" to run travel distance experiments

4. Run the `resetExperimentalSetup` to create all necessary sub-directories

```
bash resetExperimentalSetup
```

5. Download the DIMACS travel time edge-weight graph files (with extension .gr.gz and prefixed by USA-road-*t*) and coordinate files (with extension .co.gz) for DE, VT, ME and NW from https://github.com/tenindra/rn_knn_exp_data to the output_path/data/dimacs directory

6. Download the DIMACS travel time edge-weight graph files and coordinate files for COL, CAL, E, W, CTR, USA from http://www.dis.uniroma1.it/challenge9/ to the output_path/data/dimacs directory

Note: You *must* rename the "CAL" in DIMACS graph and coordinate filenames to CALNV.

7. Download the real_world_pois.tar.gz from https://github.com/tenindra/rn_knn_exp_data to the output_path/ directory

8. Unzip the real_world_pois.tar.gz archive (ensure this creates a directory at output_path/real_world_pois populated with subdirectories with POI sets for several road networks)

9. Clean the DIMACS datasets for errors and redundancy for travel time edge weights

```
bash transformInputData
```

10. Build the binary files for the basic graph representations for travel time edge weights

```
bash buildBinaryGraphs
```

11. As in steps 3-5 in "Running Experiments", execute:

```
bash buildRoadNetworkIndexes
bash generateQuerySets
bash buildObjectIndexes
```

12. Run all travel time experiments and produce figures

```
bash runTravelTimeExperiments
```

Note: The above command also creates figures, but if you wish to recreate figures without running experiments again (which can take sometime), you can use:

```
bash createTravelTimeFigures
```


# Licence

Road Network kNN Experimental Evaluation is free software; you can redistribute it and/or modify it under the terms of the GNU Affero General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

Road Network kNN Experimental Evaluation is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License along with Road Network kNN Experimental Evaluation; see LICENSE.txt; if not, see <http://www.gnu.org/licenses/>.
