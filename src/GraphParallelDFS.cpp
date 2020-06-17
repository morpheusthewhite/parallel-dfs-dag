//
// Created by francesco on 6/16/20.
//

#include "../include/GraphParallelDFS.h"


int GraphParallelDFS::getNNodes() const {
    return n_nodes;
}

const vector<int> &GraphParallelDFS::getAp() const {
    return Ap;
}

const vector<int> &GraphParallelDFS::getAi() const {
    return Ai;
}

const vector<int> &GraphParallelDFS::getRoots() const {
    return roots;
}

const vector<int> &GraphParallelDFS::getEV() const {
    return e_v;
}

const vector<int> &GraphParallelDFS::getSV() const {
    return s_v;
}
