#ifndef PARALLEL_DFS_DAG_GRAPHPARALLELDFS_H
#define PARALLEL_DFS_DAG_GRAPHPARALLELDFS_H

#include <string>
#include <vector>
#include <mutex>
#include <future>
#include "taskflow/taskflow.hpp"

using namespace std;
using namespace tf;

class GraphParallelDFS {
public:
    GraphParallelDFS() = delete;

    GraphParallelDFS(const string& filename);

    // this function will calculate the final labels, going through all the phases
    void computeLabels();

    // this function will convert the current dag to a dt and update its internal structures (phase 1)
    void convertToDT();

    friend ostream& operator<<(ostream& os, GraphParallelDFS& graphParallelDfs);

    // create output file
    void saveTo(const string& filename);

    // getters here
    int getNNodes() const;

    const vector<int> &getAp_dag() const;

    const vector<int> &getAi_dag() const;

    const vector<int> &getAp_dt() const;

    const vector<int> &getAi_dt() const;

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
    vector<int> Ap_dag;
    vector<int> Ai_dag;

    vector<int> Ap_dt;
    vector<int> Ai_dt;

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
    vector<int> post_order;

    // contains the parent of each node in the dt
    vector<int> parents;

    // contains all the parents of each node in the dag
    vector<vector<int>> parents_dag;

    Executor executor;
};


#endif //PARALLEL_DFS_DAG_GRAPHPARALLELDFS_H
