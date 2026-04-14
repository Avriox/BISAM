//
// Created by jakob on 4/13/25.
//
#include "modelselection_strategy.h"
#include <stdexcept>
#include <algorithm>
#include <iostream>

namespace bisam {
    // Global instance of the parallel executor to be reused across calls
    static ModelSelectionParallelExecutor g_parallel_executor;

    // Flag to track if we've initialized the persistent parallel region
    static bool g_omp_pool_initialized = false;

    // ModelSelectionParallelExecutor implementation
    ModelSelectionParallelExecutor::ModelSelectionParallelExecutor(int num_threads) : initialized(false) {
        if (num_threads <= 0) {
            // Use the maximum number of threads available
#ifdef _OPENMP
            this->num_threads = omp_get_max_threads();
#else
            this->num_threads = 1;
#endif
        } else {
            this->num_threads = num_threads;
        }
    }

    void ModelSelectionParallelExecutor::set_max_threads(int max_threads) {
        // Handle special cases: 0 or negative means use all available
        if (max_threads <= 0) {
#ifdef _OPENMP
            this->num_threads = omp_get_max_threads();
#else
            this->num_threads = 1;
#endif
        } else {
            // Use exactly the specified number of threads
            this->num_threads = max_threads;
        }

#ifdef _OPENMP
        // Force re-initialization with new thread count
        initialized            = false;
        g_omp_pool_initialized = false;
#endif
    }

    void ModelSelectionParallelExecutor::initialize() {
        if (!initialized) {
#ifdef _OPENMP
            // Control nested parallelism - disable it to prevent excessive thread creation
            omp_set_nested(0);

            // Set dynamic adjustment off to enforce exact thread count
            omp_set_dynamic(0);

            // Set the OpenMP thread pool size
            omp_set_num_threads(num_threads);

            // Force thread pool creation with exact thread count
            if (!g_omp_pool_initialized) {
                // Create thread pool with a dummy parallel region
#pragma omp parallel num_threads(num_threads)
                {
#pragma omp single
                    {
                        g_omp_pool_initialized = true;
                        // Debug output to verify thread count
                        // std::cout << "OpenMP initialized with " << omp_get_num_threads() << " threads" << std::endl;
                    }
                }
            }
#endif
            initialized = true;
        }
    }

    // The main execution method that handles parallel processing
    ModelSelectionOutput ModelSelectionParallelExecutor::execute_parallel( // CHANGED
        const arma::vec &y,
        const arma::mat &x,
        int niter,
        int thinning,
        int burnin,
        arma::Col<int> &deltaini_input,
        bool center,
        bool scale,
        bool XtXprecomp,
        double phi,
        double tau,
        double priorSkew,
        int prDelta,
        double prDeltap,
        std::vector<double> parprDeltap,
        arma::vec thinit,
        InitType initpar_type,
        // arma::Col<int> &include_vars,
        int method,
        int hesstype,
        int optimMethod,
        int optim_maxit,
        int B,
        int knownphi,
        int r,
        double alpha,
        double lambda,
        int n,
        int max_threads // NEW: Use this parameter
    ) {
        // Set thread count based on max_threads parameter
        int desired_threads;
        if (max_threads <= 0) {
            // Use all available threads
#ifdef _OPENMP
            desired_threads = omp_get_max_threads();
#else
            desired_threads = 1;
#endif
        } else {
            // Use exactly the specified number
            desired_threads = max_threads;
        }

        // Only update if different from current setting
        if (desired_threads != num_threads || !initialized) {
            set_max_threads(desired_threads);
            initialize();
        }

        // Prepare the split data
        DataPartition split_data = partition_data(y,
                                                  x,
                                                  deltaini_input,
                                                  thinit,
                                                  n);

        // Initialize vector for results
        std::vector<ModelSelectionOutput> results(n); // CHANGED

#ifdef _OPENMP
        // CRITICAL: Set threads before any parallel region
        int original_threads = omp_get_max_threads();
        omp_set_num_threads(num_threads);

        // Use parallel for with explicit thread count instead of tasks
        // This gives better control over thread usage
#pragma omp parallel for num_threads(num_threads) schedule(dynamic)
        for (int part = 0; part < n; part++) {
            // Set thread count for nested regions within modelSelection
            omp_set_num_threads(1); // Force serial execution in nested regions

            results[part] = modelSelection(
                split_data.y_parts[part],
                split_data.common_x,
                niter,
                thinning,
                burnin,
                split_data.delta_init_parts[part],
                center,
                scale,
                XtXprecomp,
                phi,
                tau,
                priorSkew,
                prDelta,
                prDeltap,
                parprDeltap,
                split_data.theta_init_parts[part],
                initpar_type,
                method,
                hesstype,
                optimMethod,
                optim_maxit,
                B,
                knownphi,
                r,
                alpha,
                lambda
            );
        }

        // Restore original thread count
        omp_set_num_threads(original_threads);
#else
        // Fallback to sequential execution if OpenMP is not available
        for (int part = 0; part < n; part++) {
            results[part] = modelSelection(
                split_data.y_parts[part],
                split_data.common_x,
                niter,
                thinning,
                burnin,
                split_data.delta_init_parts[part],
                center,
                scale,
                XtXprecomp,
                phi,
                tau,
                priorSkew,
                prDelta,
                prDeltap,
                parprDeltap,
                split_data.theta_init_parts[part],
                initpar_type,
                method,
                hesstype,
                optimMethod,
                optim_maxit,
                B,
                knownphi,
                r,
                alpha,
                lambda
            );
        }
#endif

        // Combine and return results
        return combine_partition_results(results, split_data.start_columns, split_data.end_columns, x.n_cols);
        // CHANGED (results type)
    }

    ModelSelectionOutput model_selection_with_strategy( // CHANGED
        const arma::vec &y,
        const arma::mat &x,
        int niter,
        int thinning,
        int burnin,
        arma::Col<int> &deltaini_input,
        bool center,
        bool scale,
        bool XtXprecomp,
        double phi,
        double tau,
        double priorSkew,
        arma::vec thinit,
        InitType initpar_type,
        // arma::Col<int> &include_vars,
        int method,
        int hesstype,
        int optimMethod,
        int optim_maxit,
        int B,
        int knownphi,
        int r,
        double alpha,
        double lambda,
        ComputationStrategy strategy,
        int n,
        int prDelta,
        double prDeltap,
        std::vector<double> parprDeltap,
        int max_threads // NEW: Thread control parameter
    ) {
        // Only apply thread limiting for SPLIT_PARALLEL strategy
        if (strategy == ComputationStrategy::SPLIT_PARALLEL) {
#ifdef _OPENMP
            // Configure thread count based on max_threads parameter
            int desired_threads;
            if (max_threads <= 0) {
                // 0 or negative means use all available
                desired_threads = omp_get_max_threads();
            } else {
                // Use exactly the specified number
                desired_threads = max_threads;
            }

            // Update the global executor with the desired thread count
            g_parallel_executor.set_max_threads(desired_threads);
            g_parallel_executor.initialize();

            // Debug output
            // std::cout << "Using " << g_parallel_executor.get_num_threads() << " threads for parallel execution" << std::endl;
#endif
        }

        switch (strategy) {
            case ComputationStrategy::STANDARD:
                // Simple passthrough to modelSelection
                return modelSelection(y,
                                      x,
                                      niter,
                                      thinning,
                                      burnin,
                                      deltaini_input,
                                      center,
                                      scale,
                                      XtXprecomp,
                                      phi,
                                      tau,
                                      priorSkew,
                                      prDelta,
                                      prDeltap,
                                      parprDeltap,
                                      thinit,
                                      initpar_type,
                                      method,
                                      hesstype,
                                      optimMethod,
                                      optim_maxit,
                                      B,
                                      knownphi,
                                      r,
                                      alpha,
                                      lambda);

            case ComputationStrategy::SPLIT_SEQUENTIAL: {
                // Prepare the split data
                DataPartition split_data = partition_data(y,
                                                          x,
                                                          deltaini_input,
                                                          thinit,
                                                          n);

                // Process each part sequentially
                std::vector<ModelSelectionOutput> results(n); // CHANGED
                for (int part = 0; part < n; part++) {
                    results[part] = modelSelection(
                        split_data.y_parts[part],
                        split_data.common_x,
                        niter,
                        thinning,
                        burnin,
                        split_data.delta_init_parts[part],
                        center,
                        scale,
                        XtXprecomp,
                        phi,
                        tau,
                        priorSkew,
                        prDelta,
                        prDeltap,
                        parprDeltap,
                        split_data.theta_init_parts[part],
                        initpar_type,
                        method,
                        hesstype,
                        optimMethod,
                        optim_maxit,
                        B,
                        knownphi,
                        r,
                        alpha,
                        lambda);
                }

                // Combine and return results
                return combine_partition_results(results,
                                                 split_data.start_columns,
                                                 split_data.end_columns,
                                                 x.n_cols); // CHANGED (results type)
            }

            case ComputationStrategy::SPLIT_PARALLEL: {
                // Use the global parallel executor with specified thread count
                return g_parallel_executor.execute_parallel(
                    y,
                    x,
                    niter,
                    thinning,
                    burnin,
                    deltaini_input,
                    center,
                    scale,
                    XtXprecomp,
                    phi,
                    tau,
                    priorSkew,
                    prDelta,
                    prDeltap,
                    parprDeltap,
                    thinit,
                    initpar_type,
                    method,
                    hesstype,
                    optimMethod,
                    optim_maxit,
                    B,
                    knownphi,
                    r,
                    alpha,
                    lambda,
                    n,
                    max_threads // Pass the thread control parameter
                );
            }

            default:
                throw std::runtime_error("Unknown model selection strategy");
        }
    }

    // Updated to handle include_vars splitting
    DataPartition partition_data(
        const arma::vec &y,
        const arma::mat &x,
        arma::Col<int> &delta_initial,
        arma::vec &theta_init,
        int num_partitions
    ) {
        DataPartition data;

        size_t n_rows = y.size();
        size_t n_cols = x.n_cols;

        size_t thinit_size = theta_init.size();

        // Calculate sizes for each part
        size_t rows_per_part  = n_rows / num_partitions;
        size_t cols_per_part  = n_cols / num_partitions;
        size_t rows_remainder = n_rows % num_partitions;
        size_t cols_remainder = n_cols % num_partitions;

        // Extract the same submatrix of x for all parts
        size_t x_start_row = 0;
        size_t x_end_row   = rows_per_part + (rows_remainder > 0 ? 1 : 0) - 1;
        size_t x_start_col = 0;
        size_t x_end_col   = cols_per_part + (cols_remainder > 0 ? 1 : 0) - 1;

        // Extract the submatrix for x (same for all parts)
        data.common_x = x.submat(x_start_row, x_start_col, x_end_row, x_end_col);

        // Prepare vectors for all parts
        data.y_parts.resize(num_partitions);
        data.delta_init_parts.resize(num_partitions);
        data.theta_init_parts.resize(num_partitions);
        data.start_columns.resize(num_partitions);
        data.end_columns.resize(num_partitions);

        // Check if thinit needs to be split
        bool split_thinit = (thinit_size == n_cols);

        // Split y and deltaini_input
        for (int part = 0; part < num_partitions; part++) {
            // Calculate row range for this part
            size_t start_row = part * rows_per_part + std::min(static_cast<size_t>(part), rows_remainder);
            size_t end_row   = (part + 1) * rows_per_part + std::min(static_cast<size_t>(part + 1), rows_remainder) - 1;

            // Calculate column range for this part
            size_t start_col = part * cols_per_part + std::min(static_cast<size_t>(part), cols_remainder);
            size_t end_col   = (part + 1) * cols_per_part + std::min(static_cast<size_t>(part + 1), cols_remainder) - 1;

            // Extract the corresponding part of y
            data.y_parts[part] = y.subvec(start_row, end_row);

            // Extract the corresponding part of deltaini_input
            data.delta_init_parts[part] = delta_initial.subvec(start_col, end_col);

            // Handle thinit appropriately based on its size
            if (split_thinit) {
                // If thinit has the same length as x.n_cols, split it accordingly
                data.theta_init_parts[part] = theta_init.subvec(start_col, end_col);
            } else if (thinit_size > 0) {
                // If thinit is not empty but doesn't match x.n_cols, use the full vector
                data.theta_init_parts[part] = theta_init;
            } else {
                // If thinit is empty, create an empty vector
                data.theta_init_parts[part] = arma::vec();
            }

            // Store the column indices for later reconstruction
            data.start_columns[part] = start_col;
            data.end_columns[part]   = end_col;
        }

        return data;
    }

    ModelSelectionOutput combine_partition_results(       // CHANGED
        const std::vector<ModelSelectionOutput> &results, // CHANGED
        const std::vector<size_t> &start_columns,
        const std::vector<size_t> &end_columns,
        size_t total_columns
    ) {
        ModelSelectionOutput combined_result; // CHANGED
        combined_result.post_sample.set_size(total_columns);
        combined_result.post_sample.fill(0);
        combined_result.margpp.set_size(total_columns);
        combined_result.margpp.fill(0.0);

        for (size_t part = 0; part < results.size(); part++) {
            for (size_t j = 0; j <= (end_columns[part] - start_columns[part]); j++) {
                combined_result.post_sample(start_columns[part] + j) = results[part].post_sample(j); // CHANGED
                combined_result.margpp(start_columns[part] + j)      = results[part].margpp(j);      // CHANGED
            }
        }

        return combined_result;
    }
}
