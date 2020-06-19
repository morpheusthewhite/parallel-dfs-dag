#ifndef PARALLEL_DFS_DAG_INVALIDGRAPHINPUTFILE_H
#define PARALLEL_DFS_DAG_INVALIDGRAPHINPUTFILE_H

using namespace std;

class InvalidGraphInputFile : exception {
public:
    InvalidGraphInputFile() : message("The graph file format is not valid") {}
    explicit InvalidGraphInputFile(const char * message) : message(message) {}

    const char * what() {
        return message;
    }
private:
    const char * message;
};

#endif //PARALLEL_DFS_DAG_INVALIDGRAPHINPUTFILE_H
