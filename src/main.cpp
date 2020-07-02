#include <iostream>
#include <InvalidGraphInputFile.h>
#include <OutputFileException.h>
#include "GraphParallelDFS.h"

using namespace std;

int main(int argc, char * argv[]){
    if(argc < 3){
        cout << "Usage: input_file output_file" << endl;
        return 1;
    }

    try{
        GraphParallelDFS graph(argv[1]);
        graph.computeLabels();
        graph.saveTo(argv[2]);
    } catch (InvalidGraphInputFile &e) {
        cout << e.what();
        return 1;
    } catch (OutputFileException &e) {
        cout << e.what();
        return 1;
    }

    return 0;
}
