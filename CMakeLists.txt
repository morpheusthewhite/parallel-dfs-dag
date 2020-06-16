cmake_minimum_required(VERSION 3.16.5)
project(parallel-dfs-dag)

set(CMAKE_CXX_STANDARD 17)

add_executable(parallel-dfs-dag src/main.cpp src/GraphParallelDFS.cpp include/GraphParallelDFS.h)