//
// Created by jakob on 9/17/24.
//

#ifndef CPP_MS_PARALLEL_Z_H
#define CPP_MS_PARALLEL_Z_H

// Local includes
#include "./ms_base.h"
#include "./ms_no_optimization.h"
#include "./ms_split_z.h"
#include "../mombf/modselIntegrals.h"
#include "../SafeQueue.h"

// Standard library includes
#include <thread>
#include <atomic>
#include <condition_variable>
#include <future>


#define NUM_THREADS 16

Eigen::VectorXi model_selection_parallel_z(
        Eigen::VectorXd y,
        Eigen::MatrixXi x,
        int n_iter,
        msPriorSpec prior_coef,
        msPriorSpec prior_delta,
        double phi,
        Eigen::VectorXi w_i,
        int n_observations,
        int n_timeperiods
);

struct Task {
    int index;
    Eigen::MatrixXi x;
    Eigen::VectorXd y;
    Eigen::VectorXi wi;
    int n_iter;
    const msPriorSpec *prior_coef;
    const msPriorSpec *prior_delta;
    double phi;
    int n_observations;
    int n_timeperiods;
    std::promise<Eigen::MatrixXi> result_promise;
};

void worker(SafeQueue<Task> &task_queue,
            std::vector<Eigen::MatrixXi> &results,
            int n_iter,
            const msPriorSpec &prior_coef,
            const msPriorSpec &prior_delta,
            double phi,
            int n_observations,
            int n_timeperiods,
            std::atomic<int> &completed_tasks);

class GlobalThreadPool {
private:
    SafeQueue<std::shared_ptr<Task>> task_queue;
    std::vector<std::thread> threads;
    std::atomic<bool> stop_flag{false};
    std::condition_variable cv;
    std::mutex cv_mutex;

    void worker_thread() {
        while (!stop_flag) {
            std::shared_ptr<Task> task;
            {
                std::unique_lock<std::mutex> lock(cv_mutex);
                cv.wait(lock, [this] { return !task_queue.Empty() || stop_flag; });
                if (stop_flag) break;
                if (!task_queue.Consume(task)) continue;
            }

            Eigen::MatrixXi result = model_selection_no_optimization(
                    task->y, task->x, task->n_iter,
                    *(task->prior_coef), *(task->prior_delta), task->phi, task->wi,
                    task->n_observations, task->n_timeperiods);
            task->result_promise.set_value(std::move(result));
        }
    }

public:
    GlobalThreadPool() {
        for (int i = 0; i < NUM_THREADS; ++i) {
            threads.emplace_back(&GlobalThreadPool::worker_thread, this);
        }
    }

    ~GlobalThreadPool() {
        {
            std::lock_guard<std::mutex> lock(cv_mutex);
            stop_flag = true;
        }
        cv.notify_all();
        for (auto &thread: threads) {
            thread.join();
        }
    }

    std::future<Eigen::MatrixXi> add_task(std::shared_ptr<Task> task) {
        std::future<Eigen::MatrixXi> future = task->result_promise.get_future();
        {
            std::lock_guard<std::mutex> lock(cv_mutex);
            task_queue.Produce(std::move(task));
        }
        cv.notify_one();
        return future;
    }
};

#endif //CPP_MS_PARALLEL_Z_H
