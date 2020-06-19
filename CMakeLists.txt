cmake_minimum_required(VERSION 3.16.5)
project(parallel-dfs-dag)

set(CMAKE_CXX_STANDARD 17)

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

include_directories(include)

add_executable(parallel-dfs-dag src/main.cpp src/GraphParallelDFS.cpp)
target_link_libraries(parallel-dfs-dag Threads::Threads)
