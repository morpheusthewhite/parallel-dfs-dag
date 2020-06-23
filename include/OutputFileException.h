//
// Created by claud on 23/06/2020.
//

#ifndef PARALLEL_DFS_DAG_OUTPUTFILEEXCEPTION_H
#define PARALLEL_DFS_DAG_OUTPUTFILEEXCEPTION_H

using namespace std;

class OutputFileException : exception {
public:
    OutputFileException() : message("Unable to open the output file") {}
    explicit OutputFileException(const char * message) : message(message) {}

    const char * what() {
        return message;
    }
private:
    const char * message;
};

#endif //PARALLEL_DFS_DAG_OUTPUTFILEEXCEPTION_H
