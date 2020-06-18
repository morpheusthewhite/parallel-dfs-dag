//
// Created by francesco on 6/16/20.
//

#include "../include/GraphParallelDFS.h"
#include <fstream>
#include <sstream>
#include "../include/InvalidGraphInputFile.h"
#include <cctype>

using namespace std;

GraphParallelDFS::GraphParallelDFS(const string &filename) {
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

        // save node index in Ap
        this->Ap.push_back(this->Ai.size());

        // read the first node
        stream >> node;
        this->Ai.push_back(node);

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
                    this->Ai.push_back(child);
                }else{
                    throw InvalidGraphInputFile("Out of bound node");
                }

                //update incoming nodes
                this->incoming_edges[child]++;
            }
        }
    }

    //check the saved number of nodes
    if(this->Ap.size() != n_nodes){
        throw InvalidGraphInputFile("The number of read nodes doesn't match the declared one");
    }

    // mark end of Ai array
    this->Ap.push_back(Ai.size());

    //work out roots of the graph
    for(int i = 0; i < this->n_nodes; i++){
        if(!this->incoming_edges[i]){
            this->roots.push_back(i);
        }
    }
}

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

void GraphParallelDFS::convertToDT() {
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

    while(!Q.empty()){
        node_futures.clear();
        P.clear();

        for(int node : Q){
            // create and launch a task for each node in Q
            packaged_task<void(int)> task([this, &paths, &node_mutexes, &mP, &P](int node) {
                // iterate over the children of node
                int first_child = this->Ap[node] + 1;
                int ending_child = this->Ap[node + 1];

                // path to the current node with node itself (used in child for comparison)
                vector<int> Br = paths[node];
                Br.push_back(node);

                // vector to collect children futures; needed to wait on them
                vector<future<void>> child_futures;

                // create and launch a task for each child of node
                for(int i=first_child; i<ending_child; i++){
                    packaged_task<void(int, int)> task_child([this, &paths, &node_mutexes, &mP, &P, &Br](int index, int current_parent){
                        int child = this->Ai[index];

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

}
