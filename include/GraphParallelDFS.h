#ifndef PARALLEL_DFS_DAG_GRAPHPARALLELDFS_H
#define PARALLEL_DFS_DAG_GRAPHPARALLELDFS_H

#include <string>
#include <vector>
#include <mutex>
#include <future>

using namespace std;

class GraphParallelDFS {
public:
    GraphParallelDFS() = delete;

    GraphParallelDFS(const string& filename);

    // copy constructor
    GraphParallelDFS(const GraphParallelDFS& other);

    // move constructor
    GraphParallelDFS(GraphParallelDFS&& other);

    // copy assignment
    GraphParallelDFS& operator=(const GraphParallelDFS& other);

    // move assignment
    GraphParallelDFS& operator=(GraphParallelDFS&& other);

    // this function will calculate the final labels, going through all the phases
    void computeLabels();

    // this function will convert the current dag to a dt and update its internal structures (phase 1)
    void convertToDT();

    friend ostream& operator<<(ostream& os, GraphParallelDFS& graphParallelDfs);

    // getters here
    int getNNodes() const;

    const vector<int> &getAp() const;

    const vector<int> &getAi() const;

    const vector<int> &getRoots() const;

    const vector<int> &getEV() const;

    const vector<int> &getSV() const;

private:
    // phase 2
    void computeSubGraphSize();

    // phase 3
    void computePostOrder();

    // phase 4
    void computeRanks();

    int n_nodes;
    vector<int> Ap;
    vector<int> Ai;

    // precomputed leaves and roots
    vector<int> leaves;
    vector<int> roots;

    // count of edges entering each node
    vector<int> incoming_edges;
    // count of edges outgoing each node
    vector<int> outgoing_edges;

    // actual labels
    vector<int> e_v;
    vector<int> s_v;

    vector<int> gamma;
    vector<int> gamma_tilde;
    vector<int> parents;
    vector<int> post_order;
};


#endif //PARALLEL_DFS_DAG_GRAPHPARALLELDFS_H
