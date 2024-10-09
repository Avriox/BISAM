//
// Created by jakob on 9/17/24.
//

#include "./ms_parallel_z.h"


// Worker function for thread pool
void worker(SafeQueue<Task> &task_queue,
            std::vector<Eigen::MatrixXi> &results,
            int n_iter,
            const msPriorSpec &prior_coef,
            const msPriorSpec &prior_delta,
            double phi,
            int n_observations,
            int n_timeperiods,
            std::atomic<int> &completed_tasks) {
    Task task;
    while (task_queue.Consume(task)) {
        results[task.index] = model_selection_no_optimization(
                task.y, task.x, n_iter,
                prior_coef, prior_delta, phi, task.wi,
                n_observations, n_timeperiods, task.standardize);
        completed_tasks++;
    }
}

// Global thread pool instance
GlobalThreadPool global_thread_pool;

Eigen::VectorXi model_selection_parallel_z(
        Eigen::VectorXd y,
        Eigen::MatrixXi x,
        int n_iter,
        msPriorSpec prior_coef,
        msPriorSpec prior_delta,
        double phi,
        Eigen::VectorXi w_i,
        int n_observations,
        int n_timeperiods,
        bool standardize) {

    std::vector<Eigen::MatrixXi> split_xs = splitMatrix(x, n_timeperiods, n_observations);
    std::vector<Eigen::VectorXd> split_ys = splitVector(y, n_timeperiods);
    std::vector<Eigen::VectorXi> split_wis = splitVector(w_i, n_observations);

    int num_segments = split_xs.size();
    std::vector<std::future<Eigen::MatrixXi>> futures(num_segments);

    // Submit tasks to the global thread pool
    for (int i = 0; i < num_segments; ++i) {
        auto task = std::make_shared<Task>(Task{
                i, std::move(split_xs[i]), std::move(split_ys[i]), std::move(split_wis[i]),
                n_iter, &prior_coef, &prior_delta, phi, n_observations, n_timeperiods, standardize
        });
        futures[i] = global_thread_pool.add_task(std::move(task));
    }

    // Collect results
    std::vector<Eigen::MatrixXi> results(num_segments);
    for (int i = 0; i < num_segments; ++i) {
        results[i] = futures[i].get();
    }

    // Combine the results
    Eigen::VectorXi post_samples = Eigen::VectorXi(x.cols());
    for (int i = 0; i < num_segments; i++) {
        post_samples.segment(n_observations * i, n_observations) = results[i];
    }

    return post_samples;
}