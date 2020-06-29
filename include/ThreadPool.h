#ifndef THREAD_POOL_H
#define THREAD_POOL_H

#include <vector>
#include <queue>
#include <memory>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <future>
#include <functional>
#include <stdexcept>

class ThreadPool {
public:
    ThreadPool(size_t);
    template<class F, class... Args>
    std::future<typename std::result_of<F(Args...)>::type> enqueue(F&& f, Args&&... args);
    ~ThreadPool();
private:
    // need to keep track of threads so we can join them
    std::vector< std::thread > workers;
    // the task queue
    std::queue< std::function<void()> > tasks;

    // synchronization; needed to handle thread start/end
    std::mutex queue_mutex;
    std::condition_variable condition;
    bool stop;
};

// the constructor launches some amount of workers
// create the number of thread passed as only argument
inline ThreadPool::ThreadPool(size_t threads) : stop(false) {
    for(size_t i = 0; i < threads; ++i) {

        // save this thread into the workers vector
        workers.emplace_back( [this]() {

            // each one of this tasks loop indefinitely until the stop variable is true and all the work was done
            for (;;) {
                std::function<void()> task;

                // scope this part to automatically release mutex
                {
                    std::unique_lock<std::mutex> lock(this->queue_mutex);

                    // wait until a new task is available or the stop is reached
                    this->condition.wait(lock,
                                         [this] { return this->stop || !this->tasks.empty(); });
                    if (this->stop && this->tasks.empty())
                        return;

                    // get a task from the queue
                    task = std::move(this->tasks.front());
                    this->tasks.pop();
                }

                // and execute it
                task();
            }
        });
    }
}

// add new work item to the pool
template<class F, class... Args>
std::future<typename std::result_of<F(Args...)>::type> ThreadPool::enqueue(F&& f, Args&&... args) {
    using return_type = typename std::result_of<F(Args...)>::type;

    // bind is used to create a wrapper to call f with the given argument(s), without the need to pass them again when
    // actually executing the task
    // forward is needed since it is not known the type (and number) of the parameters which will be passed to f
    auto task = std::make_shared< std::packaged_task<return_type()> >(
        std::bind(std::forward<F>(f), std::forward<Args>(args)...)
    );

    // future to wait for the task completion
    std::future<return_type> res = task->get_future();

    // scope to release mutex
    {
        std::unique_lock<std::mutex> lock(queue_mutex);

        // don't allow enqueueing after stopping the pool
        if(stop)
            throw std::runtime_error("enqueue on stopped ThreadPool");

        // add a task to the queue
        tasks.emplace([task](){ (*task)(); });
    }

    // wake up a thread if someone is waiting
    condition.notify_one();
    return res;
}

// the destructor joins all threads
inline ThreadPool::~ThreadPool() {
    {
        std::unique_lock<std::mutex> lock(queue_mutex);
        stop = true;
    }

    // wake all threads notifying the stop of the thread pool
    condition.notify_all();

    // and wait fot their termination before continuing
    for(std::thread &worker: workers)
        worker.join();
}

#endif
