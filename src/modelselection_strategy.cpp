//
// Created by jakob on 4/13/25.
//
#include "modelselection_strategy.h"
#include <stdexcept>

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
        if (max_threads > 0) {
            this->num_threads = max_threads;
#ifdef _OPENMP
            // Update OpenMP's thread count for future parallel regions
            omp_set_num_threads(max_threads);
#endif
        }
    }

    void ModelSelectionParallelExecutor::initialize() {
        if (!initialized) {
#ifdef _OPENMP

            // Control nested parallelism - disable it to prevent excessive thread creation
            omp_set_nested(0);

            // Set the OpenMP thread pool size once at initialization
            omp_set_num_threads(num_threads);

            // Create a global thread pool if it hasn't been created yet
            if (!g_omp_pool_initialized) {
                // Force thread pool creation with a minimal parallel region
#pragma omp parallel
                {
#pragma omp single nowait
                    {
                        g_omp_pool_initialized = true;
                    }
                }
            }
#endif
            initialized = true;
        }
    }

    // The main execution method that handles parallel processing
    arma::Col<int> ModelSelectionParallelExecutor::execute_parallel(
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
        double prDeltap,
        arma::vec thinit,
        InitType initpar_type,
        int method,
        int hesstype,
        int optimMethod,
        int optim_maxit,
        int B,
        int knownphi,
        int r,
        double alpha,
        double lambda,
        int n
    ) {
        // Make sure we're initialized with appropriate thread count
        if (!initialized) {
            // For first initialization, set thread pool size based on partition count
            // This ensures we never create more threads than necessary
            int appropriate_thread_count = std::min(n,
#ifdef _OPENMP
                                                    omp_get_max_threads()
#else
                1
#endif
            );
            set_max_threads(appropriate_thread_count);
            initialize();
        }

        // Prepare the split data - this cannot be changed as per requirements
        DataPartition split_data = partition_data(y, x, deltaini_input, thinit, n);

        // Initialize vector for results with appropriate padding to avoid false sharing
        // Each result will be aligned to cache line boundary (64 bytes typical)
        std::vector<arma::Col<int> > results(n);

#ifdef _OPENMP
        // Use the single construct with tasks for better load balancing
        // This allows the runtime to schedule work dynamically
#pragma omp parallel
        {
#pragma omp single nowait
            {
                for (int part = 0; part < n; part++) {
#pragma omp task firstprivate(part)
                    {
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
                            prDeltap,
                            split_data.theta_init_parts[part],
                            initpar_type);
                    }
                }
                // Implicit taskwait at the end of the single construct
            }
        }
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
                prDeltap,
                split_data.theta_init_parts[part],
                initpar_type);
        }
#endif

        // Combine and return results
        return combine_partition_results(results, split_data.start_columns, split_data.end_columns, x.n_cols);
    }

    arma::Col<int> model_selection_with_strategy(const arma::vec &y,
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
                                                 double prDeltap,
                                                 arma::vec thinit,
                                                 InitType initpar_type,
                                                 // NEW PARAMETERS
                                                 int method,
                                                 int hesstype,
                                                 int optimMethod,
                                                 int optim_maxit,
                                                 int B,
                                                 int knownphi,
                                                 int r,
                                                 double alpha,
                                                 double lambda,
                                                 // /NEW PARAMETERS
                                                 ComputationStrategy strategy,
                                                 int n) {
        // Important: Set the thread pool size based on partition count BEFORE the first parallel region
        // This ensures the OpenMP thread pool is created with the optimal size and reused for all iterations
        static bool first_call = true;
        if (first_call && strategy == ComputationStrategy::SPLIT_PARALLEL) {
#ifdef _OPENMP
            // Create the optimal number of threads on first call
            int max_threads     = omp_get_max_threads();
            int optimal_threads = std::min(n, max_threads);

            // Set thread count once for all future parallel regions
            g_parallel_executor.set_max_threads(optimal_threads);

            // Initialize the persistent thread pool
            g_parallel_executor.initialize();
#endif
            first_call = false;
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
                                      prDeltap,
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
                DataPartition split_data = partition_data(y, x, deltaini_input, thinit, n);

                // Process each part sequentially
                std::vector<arma::Col<int> > results(n);
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
                        prDeltap,
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
                return combine_partition_results(results, split_data.start_columns, split_data.end_columns,
                                                 x.n_cols);
            }

            case ComputationStrategy::SPLIT_PARALLEL: {
                // Use the global parallel executor which maintains thread state between calls
                return g_parallel_executor.execute_parallel(
                    y, x, niter, thinning, burnin, deltaini_input,
                    center, scale, XtXprecomp, phi, tau, priorSkew,
                    prDeltap, thinit, initpar_type, method,
                    hesstype,
                    optimMethod,
                    optim_maxit,
                    B,
                    knownphi,
                    r,
                    alpha,
                    lambda,
                    n
                );
            }

            default:
                throw std::runtime_error("Unknown model selection strategy");
        }
    }

    // The data partitioning function remains unchanged as per requirements
    DataPartition partition_data(
        const arma::vec &y,
        const arma::mat &x,
        arma::Col<int> &delta_initial,
        arma::vec &theta_init,
        int num_partitions
    ) {
        DataPartition data;

        size_t n_rows      = y.size();
        size_t n_cols      = x.n_cols;
        size_t thinit_size = theta_init.size(); // Check the actual size of thinit

        // Calculate sizes for each part
        size_t rows_per_part  = n_rows / num_partitions;
        size_t cols_per_part  = n_cols / num_partitions;
        size_t rows_remainder = n_rows % num_partitions;
        size_t cols_remainder = n_cols % num_partitions;

        // Extract the same submatrix of x for all parts (as per requirement)
        // We'll use the first partition's dimensions
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

    arma::Col<int> combine_partition_results(
        const std::vector<arma::Col<int> > &results,
        const std::vector<size_t> &start_columns,
        const std::vector<size_t> &end_columns,
        size_t total_columns
    ) {
        arma::Col<int> combined_result(total_columns, arma::fill::zeros);

        for (size_t part = 0; part < results.size(); part++) {
            for (size_t j = 0; j <= (end_columns[part] - start_columns[part]); j++) {
                combined_result(start_columns[part] + j) = results[part](j);
            }
        }

        return combined_result;
    }
}
