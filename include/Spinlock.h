#ifndef PARALLEL_DFS_DAG_SPINLOCK_H
#define PARALLEL_DFS_DAG_SPINLOCK_H
#include <atomic>

using namespace std;

class Spinlock {
public:
    Spinlock();

    void lock();

    void unlock();

private:
    atomic_bool value;
};


#endif //PARALLEL_DFS_DAG_SPINLOCK_H
