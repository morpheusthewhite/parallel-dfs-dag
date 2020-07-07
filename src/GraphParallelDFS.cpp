#include <fstream>
#include <sstream>
#include <cctype>
#include <algorithm>
#include <climits>
#include <chrono>
#include <iostream>
#include "InvalidGraphInputFile.h"
#include "OutputFileException.h"
#include "GraphParallelDFS.h"


using namespace std;
using namespace std::chrono;

GraphParallelDFS::GraphParallelDFS(const string &filename) : n_nodes(0) {

    ifstream inputFile(filename);
    string buffer;

    if(!inputFile.is_open()){
        throw InvalidGraphInputFile("Unable to open the file");
    }

#ifdef PRINT_TIME
    auto function_start = high_resolution_clock::now();
#endif

    // Get the number of nodes in the graph
    if(getline(inputFile, buffer)){
        istringstream stream(buffer);
        stream >> this->n_nodes;
    }

    //check the validity of the number of nodes
    if(this->n_nodes > 0){
        this->incoming_edges.resize(n_nodes);
        this->Ap_dag.resize(n_nodes);
    }else{
        throw InvalidGraphInputFile("Number of node in the graph not valid");
    }

    while(getline(inputFile, buffer)){
        int node;
        char c_buffer;
        istringstream stream(buffer);

        // read the first node number in the line
        stream >> node;

        if(node >= n_nodes || node < 0)
           throw InvalidGraphInputFile("Number of node in the graph not valid");

        // save node index (starting position in Ai) in Ap
        this->Ap_dag[node] = this->Ai_dag.size();

        // add the node itself in Ai
        this->Ai_dag.push_back(node);

        while(stream >> c_buffer && c_buffer != '#'){
            if(isdigit(c_buffer)){
                int child;
                // go one position back to read the entire number
                int pos = stream.tellg();
                stream.seekg(pos-1);

                // save the child
                stream >> child;

                // check the validity of the node value
                if(child < n_nodes && child >= 0){
                    this->Ai_dag.push_back(child);
                }else{
                    throw InvalidGraphInputFile("Out of bound node");
                }

                //update incoming nodes
                this->incoming_edges[child]++;
            }
        }
    }

    //check the saved number of nodes
    if(this->Ap_dag.size() != n_nodes){
        throw InvalidGraphInputFile("The number of read nodes doesn't match the declared one");
    }

    // mark end of Ai_dag array
    this->Ap_dag.push_back(Ai_dag.size());

    // work out roots of the graph
    // they will be necessary in the next phases
    for(unsigned int i = 0; i < this->n_nodes; i++){
        if(!this->incoming_edges[i]){
            this->roots.push_back(i);
        }
    }

#ifdef PRINT_TIME
    auto function_end = high_resolution_clock::now();
    auto function_interval = duration_cast<milliseconds>(function_end-function_start);
    cout << "File parsing ends in: " << function_interval.count() << " milliseconds"  << endl;
#endif
}

int GraphParallelDFS::getNNodes() const {
    return n_nodes;
}

const vector<int> &GraphParallelDFS::getAp_dag() const {
    return Ap_dag;
}

const vector<int> &GraphParallelDFS::getAi_dag() const {
    return Ai_dag;
}

const vector<int> &GraphParallelDFS::getAp_dt() const {
    return Ap_dt;
}

const vector<int> &GraphParallelDFS::getAi_dt() const {
    return Ai_dt;
}

const vector<int> &GraphParallelDFS::getEV() const {
    return e_v;
}

const vector<int> &GraphParallelDFS::getSV() const {
    return s_v;
}

void GraphParallelDFS::convertToDT() {
#ifdef PRINT_TIME
    auto function_start = high_resolution_clock::now();
#endif
    // initialize parent vector
    this->parents.resize(this->n_nodes, -1);

    // vector holding the best current path for each node
    vector<vector<int>> paths = vector<vector<int>>(this->n_nodes);

    // vector of mutexes, one for each node
    vector<mutex> node_mutexes = vector<mutex>(this->n_nodes);

    // copy roots since they will be used again later
    vector<int> Q = this->roots;

    vector<int> P;

    // mutex for protecting P
    mutex mP;

    // store the list of the parents in the dag for each node except the one in the dt
    this->parents_dag.resize(n_nodes);

    vector<future<void>> futures;

    while(!Q.empty()) {
        P = vector<int>();
        futures.clear();

        // get available parallelism in the current architecture to tune the number of generated async
        int nThreads = thread::hardware_concurrency();

// work out the number of iterations to give to the threads
        int totalSize = Q.size();

        // check the total number of iterations
        // if the number of iterations is lower than the number of threads, the number of threads started will be
        // equal to the number of iterations. Otherwise use all the available threads
        if (totalSize < nThreads)
            nThreads = totalSize;

        // the minimum number of iterations assigned to each thread
        int threadPortion = totalSize / nThreads;

        // the portion assigned to the last thread (will need to take care of the remainder)
        int lastPortion = totalSize - (threadPortion * nThreads);

        for (int thIndex = 0; thIndex < nThreads; thIndex++) {
            future<void> f = async([this, &paths, &node_mutexes, &mP, &P, &threadPortion, &lastPortion, &nThreads, &Q, thIndex]() {

                // calculate end index on Q: if this is the last thread then it will take care of the remaining portion
                int max = (thIndex == (nThreads - 1)) ? ((thIndex + 1) * threadPortion + lastPortion) : (threadPortion * (thIndex + 1));

                for (int i = threadPortion * thIndex; i < max; i++) {

                    int node = Q[i];

                    // iterate over the children of node
                    int first_child = this->Ap_dag[node] + 1;
                    int ending_child = this->Ap_dag[node + 1];

                    // path to the current node with node itself (used in child for comparison)
                    vector<int> Br = paths[node];

                    // check that node has more than 1 child before adding it to the path
                    // otherwise it will never be a decision point for the path selection
                    if (first_child != ending_child)
                        Br.push_back(node);

                    for (int j = first_child; j < ending_child; j++) {
                        int child = this->Ai_dag[j];

                        // lock mutex to access shared resources (paths and incoming_edges)
                        node_mutexes[child].lock();

                        // best path to child
                        vector<int> Qr = paths[child];

                        // update the path in the case in which
                        // - The path is empty, so this is the first time we meet this child
                        // - We found a path Br which is better than the current one
                        if (Qr.empty() || Br <= Qr) {
                            paths[child] = Br;

                            // the previous one, if existing, is a parent in the dag which is not present in the dt
                            if (parents[child] != -1)
                                this->parents_dag[child].push_back(parents[child]);

                            parents[child] = node;
                        } else {
                            // the current parent is not part of the dt but it is actually part of the dag
                            this->parents_dag[child].push_back(node);
                        }

                        // decrement the count of incoming edges which needs to be visited yet for the node child
                        int remaining = --(this->incoming_edges[child]);

                        node_mutexes[child].unlock();

                        // if all incoming edges into this child have been visited, add it to P
                        if (remaining == 0) {
                            mP.lock();
                            P.push_back(child);
                            mP.unlock();
                        }
                    }
                }
            });

            futures.push_back(move(f));
        }

        // wait for all the async
        for(auto& f : futures){
            f.get();
        }

        // move all content of P to Q
        Q = move(P);
    }

    // calculate new Ai_dt, new Ap_dt, number of outgoing edges, leaves and their gamma (obvious: it's 1)
    this->outgoing_edges.resize(n_nodes);
    this->Ap_dt.resize(n_nodes + 1, 0);

    this->gamma.resize(this->n_nodes, 0);

    for(int i=0; i<n_nodes; i++){
        Ap_dt[i] = Ai_dt.size();
        Ai_dt.push_back(i);

        int first_child = Ap_dag[i] + 1;
        int end_child = Ap_dag[i + 1];

        // iterate over child of current node i
        for(int j=first_child; j<end_child; j++){
            int child = Ai_dag[j];

            // verify that this child is a child in dt
            if(this->parents[child] == i)
                Ai_dt.push_back(child);
        }

        // skip first iteration since it operates also on the previously inserted element
        if(i > 0){
            this->outgoing_edges[i-1] = Ap_dt[i] - Ap_dt[i-1] - 1;

            if(!this->outgoing_edges[i-1]){
                this->leaves.push_back(i-1);
                this->gamma[i-1] = 1;
            }
        }
    }

    // Add Ai_dt dimension to dt
    Ap_dt[this->n_nodes] = Ai_dt.size();

    // update of the number of outgoing edges for the last node (since in the previous cycle it is
    // ignored (we skip the first iteration))
    this->outgoing_edges[this->n_nodes-1] = Ap_dt[n_nodes] - Ap_dt[n_nodes-1] - 1;
    if(!this->outgoing_edges[n_nodes-1]){
        this->leaves.push_back(n_nodes-1);
        this->gamma[n_nodes - 1] = 1;
    }
#ifdef PRINT_TIME
    auto function_end = high_resolution_clock::now();
    auto function_interval = duration_cast<milliseconds>(function_end-function_start);
    cout << "DT conversion ends in: " << function_interval.count() << " milliseconds"  << endl;
#endif
}

void GraphParallelDFS::computePostOrder() {
#ifdef PRINT_TIME
    auto function_start = high_resolution_clock::now();
#endif
    // Initialize the post order vector
    this->post_order.resize(this->n_nodes, 0);

    // initialize post order value for the roots
    // this is necessary in case the conversion from dag to dt generates many detached dts
    // without this, they will count from the same starting post-order and hence
    // share some post-orders, which is wrong, since they are unique by definition
    // this can be solved by initializing their initial post-order with their gamma_tilde,
    // as if all them are children of a common parent
    int gamma_tilde_root = 0;
    for(int root : roots){
        post_order[root] = gamma_tilde_root;
        gamma_tilde_root += gamma[root];
    }

    // Use the precomputed roots
    // Move them for performance reason since they won't be used anymore
    vector<int> Q = move(this->roots);

    // mutex which protect P vector modifications
    mutex mP;

    vector<future<void>> futures;

    while(!Q.empty()){

        vector<int> P;
        futures.clear();

        // work out iteration to give to the threads
        int totalSize = Q.size();
        int nThreads = thread::hardware_concurrency();

        // check the total number of iterations
        // if the number of iterations is lower than the number of threads, the number of threads started will be
        // equal to the number of iterations. Otherwise use all the available threads
        if (totalSize < nThreads) {
            nThreads = totalSize;
        }

        int threadPortion = totalSize / nThreads;
        int lastPortion = totalSize - (threadPortion * nThreads);

        for (int thIndex = 0; thIndex < nThreads; thIndex++) {

            future<void> f = async([this, &mP, &P, &threadPortion, &nThreads, &Q, &lastPortion, thIndex]() {

                // calculate end index on Q: if this is the last thread then it will take care of the remaining portion
                int max = (thIndex == (nThreads - 1)) ? ((thIndex + 1) * threadPortion + lastPortion) : (threadPortion * (thIndex + 1));
                for (int i = threadPortion * thIndex; i < max; i++) {
                    int node = Q[i];
                    int post = this->post_order[node];

                    int child_start = this->Ap_dt[node] + 1;
                    int child_end = this->Ap_dt[node + 1];

                    // iterate over the children of the current node
                    for (int j = child_start; j < child_end; j++) {
                        int child = this->Ai_dt[j];

                        // accumulate the post order from the roots to the current node
                        // this corresponds to calculating tau
                        post_order[child] = post + this->gamma_tilde[child];

                        // add child in P for the next iteration
                        mP.lock();
                        P.push_back(child);
                        mP.unlock();
                    }

                    // work out post-order of node
                    this->post_order[node] = post + this->gamma[node] - 1;
                }
            });
            futures.push_back(move(f));
        }

        for(auto& f : futures){
            f.get();
        }

        // move P in Q for the following iteration
        Q = move(P);
    }
#ifdef PRINT_TIME
    auto function_end = high_resolution_clock::now();
    auto function_interval = duration_cast<milliseconds>(function_end-function_start);
    cout << "Post order ends in: " << function_interval.count() << " milliseconds"  << endl;
#endif
}

void GraphParallelDFS::computeSubGraphSize(){
#ifdef PRINT_TIME
    auto function_start = high_resolution_clock::now();
#endif
    // use precomputed leaves
    // since they will not be used in the next phases they can be directly moved
    vector<int> Q = move(this->leaves);

    // initialize gamma_tilde vector
    // note: gamma vector does not need to be initialized since it was done
    // in the phase 1 (when also computing its (obvious) value for the leaves)
    this->gamma_tilde.resize(this->n_nodes, 0);

    vector<int> C;

    // mutex to protect C modifications
    mutex mC;

    // copy count of outgoing edges into an atomic vector
    // since there will be concurrency in the next cycle on these
    // variables
    vector<atomic_int> outgoing_edges_atomic(this->n_nodes);
    for(int i=0; i<n_nodes; i++)
        outgoing_edges_atomic[i].store(this->outgoing_edges[i]);

    vector<future<void>> futures;

    while(!Q.empty()){
        C = vector<int>();

        futures.clear();

        // work out iteration to give to the threads
        int totalSize = Q.size();
        int nThreads = thread::hardware_concurrency();

        // check the total number of iterations
        // if the number of iterations is lower than the number of threads, the number of threads started will be
        // equal to the number of iterations. Otherwise use all the available threads
        if (totalSize < nThreads)
            nThreads = totalSize;

        int threadPortion = totalSize / nThreads;
        int lastPortion = totalSize - (threadPortion * nThreads);

        for (int thIndex = 0; thIndex < nThreads; thIndex++) {
            future<void> f = async([this, &mC, &C, &outgoing_edges_atomic, &threadPortion, &lastPortion, &nThreads, &Q, thIndex]() {

                // calculate end index on Q: if this is the last thread then it will take care of the remaining portion
                int max = (thIndex == (nThreads - 1)) ? ((thIndex + 1) * threadPortion + lastPortion) : (threadPortion * (thIndex + 1));
                for (int i = threadPortion * thIndex; i < max; i++) {
                    int node = Q[i];

                    // differently from the paper here it is not present another loop
                    // this because in the previous phase the DAG was converted to a
                    // directed tree, hence each node will have at most one parent
                    int parent = this->parents[node];

                    // verify that it has actually a parent, i. e. it is not -1
                    if (parent != -1) {

                        // check if no more outgoing edges needs to be visited yet
                        // note: fetch_sub returns the previous value
                        if (outgoing_edges_atomic[parent].fetch_sub(1) == 1) {
                            mC.lock();
                            C.push_back(parent);
                            mC.unlock();
                        }
                    }
                }
            });
            futures.push_back(move(f));
        }

        // wait for the async
        for(auto& f : futures){
            f.get();
        }

        futures.clear();

        totalSize = C.size();

        // work out iteration to give to the threads

        if(totalSize > 0){
            nThreads = thread::hardware_concurrency();
            // check the dimension of the vector
            if (totalSize < nThreads) {
                nThreads = totalSize;
            }

            threadPortion = totalSize / nThreads;
            lastPortion = totalSize - (threadPortion * nThreads);
        } else
            nThreads = 0;

        for (int thIndex = 0; thIndex < nThreads; thIndex++) {
            future<void> f = async([this, &threadPortion, &lastPortion, &nThreads, &C, thIndex]() {

                // calculate end index on Q: if this is the last thread then it will take care of the remaining portion
                int max = (thIndex == (nThreads - 1)) ? ((thIndex + 1) * threadPortion + lastPortion) : (threadPortion * (thIndex + 1));
                for (int i = threadPortion * thIndex; i < max; i++) {
                    int p = C[i];

                    int first_child = Ap_dt[p] + 1;
                    int end_child = Ap_dt[p + 1];

                    // sort children of p in order to iterate on them in lexicographical order for computing gamma_tilde value
                    // no race condition since each async will operate on a different section of Ai_dt
                    sort(this->Ai_dt.begin() + first_child, this->Ai_dt.begin() + end_child);

                    // iterate over children of p
                    for (int j = first_child; j < end_child; j++) {
                        int child = Ai_dt[j];

                        // gamma_tilde is calculated accumulating the gamma of the previous children of this parent p,
                        // stored in gamma[p]
                        this->gamma_tilde[child] = gamma[p];
                        this->gamma[p] += gamma[child];
                    }

                    // count also p itself in the size of its subgraph
                    this->gamma[p]++;
                }
            });
            futures.push_back(move(f));
        }

        // wait for all launched async to terminate
        for(auto& f : futures){
            f.get();
        }

        // move all content of C to Q
        Q = move(C);
    }
#ifdef PRINT_TIME
    auto function_end = high_resolution_clock::now();
    auto function_interval = duration_cast<milliseconds>(function_end-function_start);
    cout << "Sub graph size calculation ends in: " << function_interval.count() << " milliseconds"  << endl;
#endif
}

void GraphParallelDFS::computeRanks(){
#ifdef PRINT_TIME
    auto function_start = high_resolution_clock::now();
#endif
    this->e_v.resize(this->n_nodes);
    this->s_v.resize(this->n_nodes);

    vector<int> P;

    // mutex to protect P modifications
    mutex mP;

    vector<int> Q;

    // initialize vector of atomic of outgoing edges through Ap of the dag
    // calculate also leaves of the dag
    vector<atomic_int> outgoing(this->n_nodes);
    for(int i=0; i<this->n_nodes; i++){
        int start_child = this->Ap_dag[i] + 1;
        int end_child = this->Ap_dag[i+1];
        
        int n_children = end_child - start_child;
        outgoing[i].store(n_children);

        if(n_children == 0)
            Q.push_back(i);
    }

    vector<future<void>> futures;

    while(!Q.empty()){
        P = vector<int>();
        futures.clear();

        // work out iteration to give to the threads
        int totalSize = Q.size();
        int nThreads = thread::hardware_concurrency();

        // check the total number of iterations
        // if the number of iterations is lower than the number of threads, the number of threads started will be
        // equal to the number of iterations. Otherwise use all the available threads
        if (totalSize < nThreads) {
            nThreads = totalSize;
        }

        int threadPortion = totalSize / nThreads;
        int lastPortion = totalSize - (threadPortion * nThreads);

        for (int thIndex = 0; thIndex < nThreads; thIndex++) {
            future<void> f = async([this, &outgoing, &mP, &P, &threadPortion, &lastPortion, &nThreads, &Q, thIndex]() {

            // calculate end index on Q: if this is the last thread then it will take care of the remaining portion
            int max = (thIndex == (nThreads - 1)) ? ((thIndex + 1) * threadPortion + lastPortion) : (threadPortion * (thIndex + 1));
            for (int i = threadPortion * thIndex; i < max; i++) {
                int node = Q[i];

                // e_v is, for definition, the corresponding post order + 1
                this->e_v[node] = this->post_order[node] + 1;

                int first_child = Ap_dag[node] + 1;
                int end_child = Ap_dag[node + 1];

                if (first_child == end_child) {
                    // this is a leaf so e_v and s_v will correspond
                    s_v[node] = e_v[node];
                } else {
                    // iterate over children of node
                    int min = INT_MAX;

                    // we compute the s_v in this case as the minimum of the s_v of the children
                    // since the s_v is:
                    // - equal to e_v in case the node is a leaf
                    // - equal to the minimum of the s_v of the children otherwise
                    for (int j = first_child; j < end_child; j++) {
                        int child = Ai_dag[j];

                        // once here the values s_v of the children has been already computed
                        if (this->s_v[child] < min) min = s_v[child];
                    }

                    s_v[node] = min;
                }

                // update the count of the visited outgoing edges for the parents of the current node
                // 1. first consider the parent in the dt

                int parent_dt = this->parents[node];
                // NOTE: root nodes will never enter this condition
                if (parent_dt != -1) {
                    int remaining = outgoing[parent_dt].fetch_sub(1);
                    // check that no more children (of the dt parent) need to be visited yet
                    if (remaining == 1) {
                        mP.lock();
                        P.push_back(parent_dt);
                        mP.unlock();
                    }

                    // 2. then consider the parent(s) in the dag except the one in the dt (already considered)
                    for (int parent : parents_dag[node]) {
                        remaining = outgoing[parent].fetch_sub(1);

                        // check that no more children (of this parent) need to be visited yet
                        if (remaining == 1) {
                            mP.lock();
                            P.push_back(parent);
                            mP.unlock();
                        }
                    }
                }
            }
          });
          futures.push_back(move(f));
        }

        // wait for all launched async to terminate
        for(auto& f : futures){
            f.get();
        }

        // move all content of P to Q
        Q = move(P);
    }
#ifdef PRINT_TIME
    auto function_end = high_resolution_clock::now();
    auto function_interval = duration_cast<milliseconds>(function_end-function_start);
    cout << "Rank computation ends in: " << function_interval.count() << " milliseconds"  << endl;
#endif
}

ostream& operator<<(ostream &os, GraphParallelDFS &graphParallelDfs) {

    for(int node = 0; node < graphParallelDfs.getNNodes(); node++){
        os << node << " " << graphParallelDfs.s_v[node] << " " << graphParallelDfs.e_v[node] << endl;
    }

    return os;
}

void GraphParallelDFS::saveTo(const string& filename){
    ofstream outputFile(filename);

    // check errors in file opening
    if(!outputFile.is_open()){
        throw OutputFileException();
    }

    outputFile << *this;
}

void GraphParallelDFS::computeLabels() {
    // phase 1
    this->convertToDT();

    // phase 2
    this->computeSubGraphSize();

    // phase 3
    this->computePostOrder();

    // phase 4
    this->computeRanks();
}