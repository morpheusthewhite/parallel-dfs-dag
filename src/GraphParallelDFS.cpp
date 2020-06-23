#include <fstream>
#include <sstream>
#include <cctype>
#include <algorithm>
#include <climits>
#include "InvalidGraphInputFile.h"
#include "OutputFileException.h"
#include "GraphParallelDFS.h"

using namespace std;

GraphParallelDFS::GraphParallelDFS(const string &filename) : n_nodes(0){
    ifstream inputFile(filename);
    string buffer;

    if(!inputFile.is_open()){
        throw InvalidGraphInputFile("Unable to open the file");
    }

    // Get the number of nodes in the graph
    if(getline(inputFile, buffer)){
        istringstream stream(buffer);
        stream >> this->n_nodes;
    }

    //check the validity of the number of nodes
    if(this->n_nodes > 0){
        this->incoming_edges.resize(n_nodes);
    }else{
        throw InvalidGraphInputFile("Number of node in the graph not valid");
    }

    while(getline(inputFile, buffer)){
        int node;
        char c_buffer;
        istringstream stream(buffer);

        // save node index in Ap_dt
        this->Ap_dag.push_back(this->Ai_dag.size());

        // read the first node
        stream >> node;
        this->Ai_dag.push_back(node);

        while(stream >> c_buffer && c_buffer != '#'){
            if(isdigit(c_buffer)){
                int child;
                //go one position back to read the number
                int pos = stream.tellg();
                stream.seekg(pos-1);

                //save the child
                stream >> child;

                //check the validity of the node value
                if(child < n_nodes){
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

    // mark end of Ai_dt array
    this->Ap_dag.push_back(Ai_dag.size());

    //work out roots of the graph
    for(unsigned int i = 0; i < this->n_nodes; i++){
        if(!this->incoming_edges[i]){
            this->roots.push_back(i);
        }
    }
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
    // initialize parent vector
    this->parents = vector<int>(this->n_nodes, -1);

    // vector holding the best current path for each node
    vector<vector<int>> paths = vector<vector<int>>(this->n_nodes);

    // vector of mutexes, one for each node
    vector<mutex> node_mutexes = vector<mutex>(this->n_nodes);

    vector<int> Q = this->roots;

    // mutex for protecting P
    mutex mP;

    // vector to collect all the node futures
    vector<future<void>> node_futures;

    vector<int> P;

    // store the list of the parents of the dag for each node
    parents_dag.resize(n_nodes);

    while(!Q.empty()){
        node_futures.clear();
        P = vector<int>();

        for(int node : Q){
            // create and launch a task for each node in Q
            packaged_task<void(int)> task([this, &paths, &node_mutexes, &mP, &P](int node) {
                // iterate over the children of node
                int first_child = this->Ap_dag[node] + 1;
                int ending_child = this->Ap_dag[node + 1];

                // path to the current node with node itself (used in child for comparison)
                vector<int> Br = paths[node];
                Br.push_back(node);

                // vector to collect children futures; needed to wait on them
                vector<future<void>> child_futures;

                // create and launch a task for each child of node
                for(int i=first_child; i<ending_child; i++){
                    packaged_task<void(int, int)> task_child([this, &paths, &node_mutexes, &mP, &P, &Br](int index, int current_parent){
                        int child = this->Ai_dag[index];

                        // existing path
                        vector<int> Qr = paths[child];

                        // lock mutex to access shared resources (paths and incoming_edges)
                        node_mutexes[child].lock();

                        // update the path in the case in which
                        // - The path is empty, so this is the first time we meet this child
                        // - We found a path Br which is better than the current one
                        if(Qr.empty() || Br <= Qr){
                            paths[child] = Br;
                            parents[child] = current_parent;
                        }

                        // decrement the count of incoming edges which needs to be visited yet for the node child
                        int remaining = --(this->incoming_edges[child]);

                        // store the current parent of the current child
                        this->parents_dag[child].push_back(current_parent);

                        node_mutexes[child].unlock();

                        // if all incoming edges into this child have been visited, add it to P
                        if(remaining == 0){
                            mP.lock();
                            P.push_back(child);
                            mP.unlock();
                        }
                    });

                    child_futures.push_back(move(task_child.get_future()));

                    // actually launch task
                    task_child(i, node);
                }

                // wait for all children task to terminate
                for(auto& child_future : child_futures)
                    child_future.get();

            });

            node_futures.push_back(move(task.get_future()));

            // actually launch task
            task(node);
        }

        // wait for all launched tasks to terminate
        for(auto& node_future : node_futures)
            node_future.get();

        // move all content of P to Q
        Q = move(P);
    }

    // calculate new Ai_dt, new Ap_dt, number of outgoing edges, leaves and their gamma (obvious: it's 1)
    this->outgoing_edges.resize(n_nodes);
    Ap_dt.resize(n_nodes + 1, 0);

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

        // skip first iteration
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
    // ignored (we skip the first element))
    this->outgoing_edges[this->n_nodes-1] = Ap_dt[n_nodes] - Ap_dt[n_nodes-1] - 1;
    if(!this->outgoing_edges[n_nodes-1]){
        this->leaves.push_back(n_nodes-1);
        this->gamma[n_nodes - 1] = 1;
    }

}

void GraphParallelDFS::computePostOrder() {
    // Initialize the post order vector
    this->post_order.resize(this->n_nodes, 0);

    // Use the precomputed roots
    // Move them for performance reason since they won't be used anymore
    vector<int> Q = move(this->roots);

    //mutex which protect P vector modifications
    mutex mP;

    // prepare vector of future for tasks
    vector<future<void>> node_futures;

    while(!Q.empty()){
        // clear the vector of future
        node_futures.clear();

        vector<int> P;

        for(int node : Q){
            // create and launch a task for each node in Q
            packaged_task<void(int)> task([this, &mP, &P](int node) {
                int post = this->post_order[node];

                int children_start = this->Ap_dt[node] + 1;
                int children_end = this->Ap_dt[node + 1];

                // vector which collect the children futures
                vector<future<void>> child_futures;

                //iterate over the children of the current node
                for(int i = children_start; i < children_end; i++){
                    // create and launch task for each child of the node
                    packaged_task<void(int)> task_child([this, &post, &mP, &P](int index){
                        int child = this->Ai_dt[index];

                        // pre-compute post-order
                        post_order[child] = post + this->gamma_tilde[child];

                        // add child in P for the next iteration
                        mP.lock();
                        P.push_back(child);
                        mP.unlock();
                    });

                    child_futures.push_back(move(task_child.get_future()));

                    // launch task
                    task_child(i);
                }
                // wait termination of child tasks
                for(auto& child_future : child_futures)
                    child_future.get();

                // work out post-order of node
                this->post_order[node] = post + this->gamma[node] -1;
            });

            node_futures.push_back(move(task.get_future()));

            // launch task
            task(node);
        }

        // wait termination of all the node tasks
        for(auto& node_future : node_futures)
            node_future.get();

        // move P in Q for the following iteration
        Q = move(P);
    }

}

void GraphParallelDFS::computeSubGraphSize(){
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

    // vector to collect all the node futures
    vector<future<void>> node_futures;

    while(!Q.empty()){
        node_futures.clear();
        C = vector<int>();

        for(int node : Q){
            // create and launch a task for each node in Q
            packaged_task<void(int)> task([this, &mC, &C, &outgoing_edges_atomic](int node) {

                // differently from the paper here it is not present another loop
                // this because in the previous phase the DAG was converted to a
                // directed tree, hence each node will have at most one parent
                int parent = this->parents[node];

                // verify that it has actually a parent, i. e. it is not -1
                if(parent != -1){

                    // check if no more outgoing edges needs to be visited yet
                    // note: fetch_sub returns the previous value
                    if(outgoing_edges_atomic[parent].fetch_sub(1) == 1){
                        mC.lock();
                        C.push_back(parent);
                        mC.unlock();
                    }
                }

            });

            node_futures.push_back(move(task.get_future()));

            // actually launch task
            task(node);
        }

        // wait for all launched tasks to terminate
        for(auto& node_future : node_futures)
            node_future.get();

        node_futures.clear();

        for(int p : C) {

            // create and launch a task for each node in C
            packaged_task<void(int)> task([this](int p) {
                int first_child = Ap_dt[p] + 1;
                int end_child = Ap_dt[p + 1];

                // order children of P
                // no race condition since each task will operate on a different section of Ai_dt
                sort(this->Ai_dt.begin() + first_child, this->Ai_dt.begin() + end_child);

                // iterate over children of P
                for (int i = first_child; i < end_child; i++) {
                    int child = Ai_dt[i];

                    this->gamma_tilde[child] = gamma[p];
                    this->gamma[p] += gamma[child];
                }

                // count also p itself in the size of its subgraph
                this->gamma[p]++;
            });

            node_futures.push_back(move(task.get_future()));

            // actually launch task
            task(p);
        }

        // wait for all launched tasks to terminate
        for(auto& node_future : node_futures)
            node_future.get();

        // move all content of C to Q
        Q = move(C);
    }
}

void GraphParallelDFS::computeRanks(){
    this->e_v.resize(this->n_nodes, 0);
    this->s_v.resize(this->n_nodes, 0);

    vector<int> P;
    
    // mutex to protect P modifications
    mutex mP;
    
    // move leaves since they won't be used anymore
    vector<int> Q;
    
    // list to collect children futures
    vector<future<void>> node_futures;

    // initialize vector of atomic variables with the precomputed outgoing_edges
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

    while(!Q.empty()){
        P = vector<int>();
        node_futures.clear();

        for(int node : Q){
            // create and launch a task for each node in Q
            packaged_task<void(int)> task([this, &outgoing, &mP, &P](int node) {
                // e_v is, for definition, the corresponding post order + 1
                this->e_v[node] = this->post_order[node] + 1;
                
                int first_child = Ap_dag[node] + 1;
                int end_child = Ap_dag[node + 1];

                int min = INT_MAX;

                if(first_child == end_child){
                    // this is a leaf so e_v and s_v will correspond
                    s_v[node] = e_v[node];
                } else {
                    // iterate over children of node

                    // we compute the s_v in this case as the minimum of the s_v of the children
                    // since they s_v is:
                    // - equal to e_v in case they are a leaf
                    // - equal to the mimimum of the s_v of the children otherwise
                    for (int i = first_child; i < end_child; i++) {
                        // no parallelization in this cycle: min should be protected in that case
                        // which would reduce the benefits of the parallelization, considering also that
                        // the actual parallelizable content of the cycle is almost none
                        int child = Ai_dag[i];

                        // once here we visited the children of the node node
                        if(this->s_v[child] < min) min = s_v[child];
                    }

                    s_v[node] = min;
                }

                // update the count of the visited outgoing edges for the parentS of the current node
                // NOTE: root nodes will never enter this cycle
                for(int parent : parents_dag[node]){
                    int remaining = outgoing[parent].fetch_sub(1);

                    // check that no more children (of this parent) needs to be visited yet
                    if(remaining == 1){
                        mP.lock();
                        P.push_back(parent);
                        mP.unlock();
                    }
                }
            });

            node_futures.push_back(move(task.get_future()));
            
            // actually launch task
            task(node);

        }

        // wait for all launched tasks to terminate
        for(auto& node_future : node_futures)
            node_future.get();

        // move all content of P to Q
        Q = move(P);
    }
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