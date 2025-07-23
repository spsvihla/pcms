/**
 * @file thread-pool.hpp
 * @brief The thread pool class.
 * @author Sean Svihla
 */
#ifndef THREAD_POOL_HPP
#define THREAD_POOL_HPP

#include <vector>
#include <thread>
#include <queue>
#include <future>
#include <functional>
#include <mutex>
#include <condition_variable>
#include <atomic>

/**
 * @class ThreadPool
 * @brief A simple fixed-size thread pool for executing asynchronous tasks.
 *
 * This class creates a pool of worker threads which continuously
 * fetch and execute tasks from an internal task queue.
 * Users can submit any callable task with arguments and receive a
 * std::future to get the result asynchronously.
 *
 * The pool manages thread lifecycle internally: threads are started
 * on construction and stopped/joined on destruction.
 *
 * Thread-safety:
 * - submit() is thread-safe.
 * - Destruction waits for all running tasks to complete.
 *
 * Usage example:
 * @code
 * ThreadPool pool(4); // create a pool with 4 worker threads
 * auto future = pool.submit([](int x){ return x*x; }, 5);
 * int result = future.get(); // result == 25
 * @endcode
 */
class ThreadPool {
public:
    /**
     * @brief Construct a ThreadPool with a fixed number of threads.
     * @param num_threads Number of worker threads to create.
     *
     * Creates and launches `num_threads` worker threads.
     */
    explicit ThreadPool(size_t num_threads)
        : done(false)
    {
        start(num_threads);
    }

    /**
     * @brief Destructor.
     *
     * Signals all worker threads to stop,
     * joins them, and cleans up resources.
     * Waits for all currently queued tasks to finish.
     */
    ~ThreadPool() 
    {
        stop();
    }

    /**
     * @brief Submit a task to be executed asynchronously by the thread pool.
     * @tparam Func Callable type.
     * @tparam Args Argument types for the callable.
     * @param func Callable object (function, lambda, etc).
     * @param args Arguments to pass to the callable.
     * @return std::future holding the result of the callable.
     */
    template<typename Func, typename... Args>
    auto submit(Func&& func, Args&&... args)
        -> std::future<typename std::invoke_result<Func, Args...>::type>
    {
        using ReturnType = typename std::invoke_result<Func, Args...>::type;

        // Package the task with bound arguments into a packaged_task
        auto task = std::make_shared<std::packaged_task<ReturnType()>>(
            std::bind(std::forward<Func>(func), std::forward<Args>(args)...)
        );

        std::future<ReturnType> res = task->get_future();

        {
            // Lock queue and push new task
            std::unique_lock<std::mutex> lock(queue_mutex);
            if(done)
            {
                throw std::runtime_error("ThreadPool is stopped");
            }

            tasks.emplace([task]() { (*task)(); });
        }

        // Notify one worker thread
        cond_var.notify_one();

        return res;
    }

private:
    std::vector<std::thread> workers;                   ///< Vector of worker threads.
    std::queue<std::function<void()>> tasks;            ///< Task queue holding pending tasks.

    std::mutex queue_mutex;                             ///< Mutex for protecting access to the task queue.
    std::condition_variable cond_var;                   ///< Condition variable to notify worker threads.
    std::atomic<bool> done;                             ///< Atomic flag to signal shutdown.

    /**
     * @brief Launches worker threads.
     * @param num_threads Number of threads to create.
     *
     * Each worker thread loops waiting for tasks in the queue.
     * When a task is available, it is executed.
     * Threads exit when the `done` flag is set and task queue is empty.
     */
    void start(size_t num_threads) 
    {
        for(size_t i = 0; i < num_threads; ++i) 
        {
            workers.emplace_back([this] {
                while (true) {
                    std::function<void()> task;

                    {
                        std::unique_lock<std::mutex> lock(this->queue_mutex);
                        this->cond_var.wait(lock, [this] {
                            return done || !this->tasks.empty();
                        });
                        if(done && this->tasks.empty())
                        {
                            return; // Exit thread
                        }
                        task = std::move(this->tasks.front());
                        this->tasks.pop();
                    }

                    // Execute the retrieved task outside the lock
                    task();
                }
            });
        }
    }

    /**
     * @brief Signals all threads to stop and joins them.
     *
     * Notifies all workers to exit after completing current tasks,
     * then waits for threads to join.
     */
    void stop() 
    {
        {
            std::unique_lock<std::mutex> lock(queue_mutex);
            done = true;
        }
        cond_var.notify_all();

        for(auto& thread : workers) 
        {
            if(thread.joinable())
            {
                thread.join();
            }
        }
    }
};

#endif // THREAD_POOL_HPP