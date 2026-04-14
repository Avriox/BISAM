//
// Created by jakob on 4/13/25.
//

#ifndef MODELSELECTION_STRATEGY_H
#define MODELSELECTION_STRATEGY_H
#include <vector>
#include <RcppArmadillo.h>
#include "bisam_types.h"
#include "mombf_bridge.h"

// Check for OpenMP availability
#ifdef _OPENMP
#include <omp.h>
#endif

namespace bisam {
    ModelSelectionOutput model_selection_with_strategy(const arma::vec &y,
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
                                                       // NEW PARAMETERS
                                                       // arma::Col<int> &include_vars,
                                                       int method      = 0,
                                                       int hesstype    = 1,
                                                       int optimMethod = 2,
                                                       int optim_maxit = 0,
                                                       int B           = 100000,
                                                       int knownphi    = 1,
                                                       int r           = 1,
                                                       double alpha    = 0.01,
                                                       double lambda   = 0.01,

                                                       // /NEW PARAMETERS
                                                       ComputationStrategy strategy =
                                                               ComputationStrategy::SPLIT_SEQUENTIAL,
                                                       int n                           = 3,
                                                       int prDelta                     = 1,
                                                       double prDeltap                 = 0.5,
                                                       std::vector<double> parprDeltap = {1, 1},
                                                       int max_threads                 = 0);

    // NEW: Thread control parameter (0 = use all)

    DataPartition partition_data(
        const arma::vec &y,
        const arma::mat &x,
        arma::Col<int> &delta_initial,
        arma::vec &theta_init,
        // arma::Col<int> &include_vars, // Added parameter to match deltaini_input handling
        int num_partitions
    );

    ModelSelectionOutput combine_partition_results(
        const std::vector<ModelSelectionOutput> &results,
        const std::vector<size_t> &start_columns,
        const std::vector<size_t> &end_columns,
        size_t total_columns
    );

    // Persistent executor class to minimize thread creation overhead
    class ModelSelectionParallelExecutor {
    private:
        int num_threads;
        bool initialized;
        int saved_num_threads; // Store original thread count for restoration

    public:
        // Constructor - default to 0 which means use system max
        ModelSelectionParallelExecutor(int num_threads = 0);

        // Initialize the executor and configure the OpenMP environment
        void initialize();

        // Set the maximum number of threads for future operations
        void set_max_threads(int max_threads);

        // Get the current thread count
        int get_num_threads() const { return num_threads; }

        // Cleanup and restore original OpenMP state
        void cleanup();

        // The main execution method that handles parallel processing
        ModelSelectionOutput execute_parallel(
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
            int max_threads // NEW: Thread control parameter
        );
    };
}

#endif //MODELSELECTION_STRATEGY_H
