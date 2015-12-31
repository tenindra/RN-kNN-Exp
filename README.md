# Road Network kNN Experimental Evaluation

This project consists of implementations of several kNN algorithms for road networks and the experimental framework to compare them. This has primarily been released to allow readers to reproduce results from a paper to appear at VLDB 2016 (details to be updated) and to use in future studies. If you use the code in a publication, and our work is relevant, please consider citing our paper. Please refer to the Requirements below and [FAQ](https://github.com/tenindra/RN-kNN-Exp/wiki/FAQ) if you have any issues. If you still have problems [contact us](http://users.monash.edu.au/~tenindra/). 


# Requirements


## Software

Before the executable can be compiled the following packages/libraries must be installed. Note that we were able to install METIS from the Ubuntu repositories if you are using that operating system.

- g++ version  4.9 or higher
- Boost version  1.57 or higher for serialization
- [METIS] (http://glaros.dtc.umn.edu/gkhome/metis/metis/download) version 5.1 or higher
- CMake version 2.8 or higher
- gnuplot for figure generation
- Optional: [Google Sparsehash] (http://code.google.com/p/sparsehash/) for comparison of G-tree hashtables implementations

Note: We found Boost versions lower than 1.57 did not correctly support serialization of some STL data structures (such as `std::unordered_map`) thus they cannot be used here. This means you need to remove any existing Boost version, e.g. by uninstalling `libboost-all-dev` through your package manager.


## System & Hardware

We recommend a 64-bit OS due to the size of indexes (and we have not tested on 32-bit). Furthermore 32GB RAM is required to reproduce all experiments as they are in the paper. Results can be reproduced using a smaller amount of RAM by omitting some of the larger datasets (see below). Also ensure that there is 200GB of disk space freely available on the hard disk storing indexes (for both travel distance and travel time experiments).


# Files

The C++ source code can be found in the top-level `RN-kNN-Exp` directory, and in the subdirectories `command`, `external`, `processing`, `queue`, `tuple` and `utility`. All bash scripts for various tasks such as setting up the directories, running experiments, and generating figures are found in the `scripts` subdirectory.


# Compilation

1. Change `CMAKE_C_COMPILER` and `CMAKE_CXX_COMPILER` in CMakeLists.txt to point to the correct `g++` version (i.e. 4.9 or higher). `CMakeLists.txt` can be found in the `RN-kNN-Exp` directory.

2. Create a directory called `build` in the top-level (e.g. RN-kNN-Exp/) directory

    ```
    mkdir build
    ```

3. Change into this directory and generate makefiles using `CMake`

    ```
    cd build
    cmake -G "Unix Makefiles" ../CMakeLists.txt ..
    ```

4. Compile all executables using make

    ```
    make
    ```


# Setup

1. Open the `globalVariables` bash script in `scripts` directory in your favourite editor

2. Change the `output_path` variable to the full-path of an **existing** (i.e. create it first) location where you wish to store all data (e.g. indexes, objects sets, figures, etc...)

    Note: This can be different to the path containing the code (recommended)

3. Change the `exe_path` variable to the full path of the build directory created above (e.g. `/home/user/Downloads/rn_knn/build`)

4. Run `resetExperimentalSetup` in the `scripts` directory to create all necessary sub-directories

    ```
    cd scripts
    bash resetExperimentalSetup
    ```

5. Download the DIMACS distance edge-weight graph files (with extension .gr.gz and prefixed by USA-road-*d*) and coordinate files (with extension .co.gz) for DE, VT, ME and NW from https://github.com/tenindra/RN-kNN-Exp-Data to the `output_path/data/dimacs` directory

6. Download the DIMACS distance edge-weight graph files and coordinate files for COL, CAL, E, W, CTR, USA from http://www.dis.uniroma1.it/challenge9/ to the `output_path/data/dimacs` directory

7. Download the node (with extension .cnode) and edge (with extension .cedge) files for North America (NA) from http://www.cs.utah.edu/~lifeifei/SpatialDataset to the `output_path/data/tpq` directory

    Note: You *must* gzip these files so that `output_path/data/tpq` directory contains NA.cnode.gz and NA.cedge.gz

8. Download the real_world_pois.tar.gz from https://github.com/tenindra/RN-kNN-Exp-Data to the `output_path` directory

9. Unzip the real_world_pois.tar.gz archive (this should create a `output_path/real_world_pois` directory populated with subdirectories with POI sets for several road networks)


# Less Than 32GB RAM

Experiments can still be run with less than 32GB RAM, but for fewer datasets. To do this go to the `globalVariables` file and modify the `road_networks` array. The road networks are in size order, so simply remove the road networks from the end. E.g. to run experiments up to the Eastern US dataset change it to `road_networks=("DE" "VT" "ME" "COL" "FLA" "CAL" "E")`.

Note 1: All in-depth experiments on the US dataset will not be possible, of course. In this case it is best to comment out those experiments in the `runPaperExperiments` script.

Note 2: SILC (index used by Distance Browsing) can only be built upto the default NW dataset. SILC on the NW requires at least 20GB of memory. Less than this will cause all default experiments to be missing Distance Browsing comparisons. In this case we suggest changing the default road network to COL (for which SILC only requires 8GB), by changing `default_network` and `default_parameters` to COL.


# Running Travel Distance Experiments

The follow instructions can be followed to re-create all figures from the paper. Assuming your are already in the `scripts` directory (otherwise cd into scripts):

1. Clean the DIMACS and TPQ datasets for errors and redundancy

    ```
    bash transformInputData
    ```

2. Build the binary files for the basic graph representations

    ```
    bash buildBinaryGraphs
    ```

3. Build all road network indexes (this may take a while... ~15 hours on our machine)

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

Note: All above commands may be batched in the shell terminal, just enter each separated by semi-colon ";"


# Running Travel Time Experiments

Travel times experiments must be reproduced indepedently, as they require different indexes. These can be re-produced using the following procedure:

1. Create a new directory to store data for travel time experiments (e.g. figures etc...) somewhere

2. Change the `output_path` variable in globalVariables the full-path of this new location

3. Change the `edge_type` variable to `edge_type=t`

    Note: The `edge_type` variable must be changed back to be "d" to run travel distance experiments again

4. Run the `resetExperimentalSetup` to create all necessary sub-directories

    ```
    bash resetExperimentalSetup
    ```

5. Download the DIMACS travel time edge-weight graph files (with extension .gr.gz and prefixed by USA-road-*t*) and coordinate files (with extension .co.gz) for DE, VT, ME and NW from https://github.com/tenindra/RN-kNN-Exp-Data to the `output_path/data/dimacs` directory

6. Download the DIMACS travel time edge-weight graph files and coordinate files for COL, CAL, E, W, CTR, USA from http://www.dis.uniroma1.it/challenge9/ to the `output_path/data/dimacs` directory

7. Download the real_world_pois.tar.gz from https://github.com/tenindra/RN-kNN-Exp-Data to the `output_path/ directory`

8. Unzip the real_world_pois.tar.gz archive (this should create a `output_path/real_world_pois` directory populated with subdirectories with POI sets for several road networks)

9. Clean the DIMACS datasets for errors and redundancy for travel time edge weights

    ```
    bash transformInputData
    ```

10. Build the binary files for the basic graph representations for travel time edge weights

    ```
    bash buildBinaryGraphs
    ```

11. As in steps 3-5 in "Running Travel Distance Experiments", execute:

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

