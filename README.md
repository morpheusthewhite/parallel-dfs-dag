# Parallel DFS for Direct Acyclic Graphs

This is a C++ implementation of a parallel algorithm of the DFS traversal, according to
[this paper](https://research.nvidia.com/publication/parallel-depth-first-search-directed-acyclic-graphs).

More specifically, the algorithm is used to compute post order (pre order can
be as well computed with very minimal changes to the code) which are later used
to compute outer and inner rank, respectively defined as  
- post order + 1 (e_v or outer rank)
- equal to the outer rank if the node has no children, the minimum of the outer
  ranks of the children otherwise (s_v or inner rank).

More info (except regarding ranks) in the paper. 

## Running 

For building the project use the provided `CMakeLists.txt`, with the following
commands

```shell
$ mkdir build && cd build && cmake .. && make && cp parallel-dfs-dag .. && cd ..
```

The program receive a file containing the initial dag with the following format

```
<number of nodes>
0: <node1> <node2> ... #
...

<nodeId>: <nodeN> <nodeM> ... #

...
```

The initial line contains the number of nodes while all the next lines have the
same format, starting with the (numeric) node identifier (must be incremental)
followed by a colon and the the list of the nodes to which the current node
points to.

For example

```
7
0: 1 2 #
1: 3 4 #
2: 4 5 #
3: #
4: 6 #
5: 6 #
6: #
```

is a valid input file.

The third argument of the executable is the name of the file into which the
ranks are going to be saved.
