//
// Created by claud on 18/06/2020.
//

#ifndef PARALLEL_DFS_DAG_INVALIDGRAPHINPUTFILE_H
#define PARALLEL_DFS_DAG_INVALIDGRAPHINPUTFILE_H

class InvalidGraphInputFile {
public:
    InvalidGraphInputFile() : message("The graph file format is not valid") {}
    InvalidGraphInputFile(char * message) : message(message) {}

    const char * what() {
        return message;
    }
private:
    const char * message;
};

#endif //PARALLEL_DFS_DAG_INVALIDGRAPHINPUTFILE_H
