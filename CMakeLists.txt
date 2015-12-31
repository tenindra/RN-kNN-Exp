cmake_minimum_required(VERSION 2.8)

project(nd_knn)

set(CMAKE_C_COMPILER "/usr/bin/gcc")
set(CMAKE_CXX_COMPILER "/usr/bin/g++")

#set (CMAKE_CXX_FLAGS "-g -Wall -O1 -std=c++11")
set (CMAKE_CXX_FLAGS "-O3 -std=c++11")

file(GLOB MAIN_SRC
    "common.h"
    "main.cpp"
    "command/*.h"
    "command/*.cpp"
    "tuple/*.h"
    "tuple/*.cpp"
    "utility/*.h"
    "utility/*.cpp"
    "queue/*.h"
    "queue/*.cpp"
    "processing/*.h"
    "processing/*.cpp"
    "external/*.h"
    "external/*.cpp"
    "external/command/*.h"
    "external/command/*.cpp"
    "external/datastr/*.h"
    "external/datastr/*.cpp"
    "external/datastr/graph/*.h"
    "external/datastr/graph/*.cpp"
    "external/datastr/pqueue/*.h"
    "external/datastr/pqueue/*.cpp"
    "external/io/*.h"
    "external/io/*.cpp"
    "external/processing/*.h"
    "external/processing/*.cpp"
    "external/stats/*.h"
    "external/stats/*.cpp"
)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

add_executable(nd_knn ${MAIN_SRC})
target_link_libraries (nd_knn boost_serialization metis)

# Standard Optimisations
# By default all optimisations are already applied without pre-processor directives
set_target_properties(nd_knn PROPERTIES COMPILE_FLAGS "${DefaultDefinitions}")

# Additional Executables for Experiments

# Version of nd_knn that collects other query statistics (e.g. nodes visited in INE)
# Note: This is not to be used for query times due to collection overhead
add_executable(nd_knn_stats ${MAIN_SRC})
target_link_libraries (nd_knn_stats boost_serialization metis)
set_target_properties(nd_knn_stats PROPERTIES COMPILE_FLAGS "${DefaultDefinitions} -DCOLLECT_STATISTICS")

# Versions with INE using DynamicGraph is progressively improved
add_executable(nd_knn_ine_queue_opt ${MAIN_SRC})
target_link_libraries (nd_knn_ine_queue_opt boost_serialization metis)
set_target_properties(nd_knn_ine_queue_opt PROPERTIES COMPILE_FLAGS "${DefaultDefinitions} -DINE_QUEUE_OPTIMISED")
add_executable(nd_knn_ine_visited_opt ${MAIN_SRC})
target_link_libraries (nd_knn_ine_visited_opt boost_serialization metis)
set_target_properties(nd_knn_ine_visited_opt PROPERTIES COMPILE_FLAGS "${DefaultDefinitions} -DINE_QUEUE_OPTIMISED -DINE_VISITED_OPTIMISED")

# Version with G-tree leaf optimisation applied
add_executable(nd_knn_gtree_leaf_unoptimised ${MAIN_SRC})
target_link_libraries (nd_knn_gtree_leaf_unoptimised boost_serialization metis)
set_target_properties(nd_knn_gtree_leaf_unoptimised PROPERTIES COMPILE_FLAGS "${DefaultDefinitions} -DUNOPTIMISED_GTREE_LEAF_SEARCH")

# Version with G-tree distance matrix implemented using hash tables
add_executable(nd_knn_gtree_stl_hashtable_dist_matrix ${MAIN_SRC})
target_link_libraries (nd_knn_gtree_stl_hashtable_dist_matrix boost_serialization metis)
set_target_properties(nd_knn_gtree_stl_hashtable_dist_matrix PROPERTIES COMPILE_FLAGS "${DefaultDefinitions} -DGTREE_STL_HASH_TABLE_DIST_MATRIX")
add_executable(nd_knn_gtree_google_densehash_dist_matrix ${MAIN_SRC})
target_link_libraries (nd_knn_gtree_google_densehash_dist_matrix boost_serialization metis)
set_target_properties(nd_knn_gtree_google_densehash_dist_matrix PROPERTIES COMPILE_FLAGS "${DefaultDefinitions} -DGTREE_GOOGLE_DENSEHASH_DIST_MATRIX")

# Version with IER using unoptimsed shortest path methods (e.g. G-tree without materialization, SILC without skipping)
add_executable(nd_knn_ier_unoptimised ${MAIN_SRC})
target_link_libraries (nd_knn_ier_unoptimised boost_serialization metis)
set_target_properties(nd_knn_ier_unoptimised PROPERTIES COMPILE_FLAGS "${DefaultDefinitions} -DUNOPTIMISED_IER")

# Version with optimised SILC kNN and DistBrws without using chain skipping optimisation
add_executable(nd_knn_no_chain_skip ${MAIN_SRC})
target_link_libraries (nd_knn_no_chain_skip boost_serialization metis)
set_target_properties(nd_knn_no_chain_skip PROPERTIES COMPILE_FLAGS "${DefaultDefinitions} -DNO_CHAIN_SKIPPING")
