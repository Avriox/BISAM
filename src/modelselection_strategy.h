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
                                                 ComputationStrategy strategy = ComputationStrategy::SPLIT_SEQUENTIAL,
                                                 int n                        = 3);

    DataPartition partition_data(
        const arma::vec &y,
        const arma::mat &x,
        arma::Col<int> &delta_initial,
        arma::vec &theta_init,
        int num_partitions
    );

    arma::Col<int> combine_partition_results(
        const std::vector<arma::Col<int> > &results,
        const std::vector<size_t> &start_columns,
        const std::vector<size_t> &end_columns,
        size_t total_columns
    );

    // Persistent executor class to minimize thread creation overhead
    class ModelSelectionParallelExecutor {
    private:
        int num_threads;
        bool initialized;

    public:
        // Constructor - default to 0 which means use system max
        ModelSelectionParallelExecutor(int num_threads = 0);

        // Initialize the executor and configure the OpenMP environment
        void initialize();

        // Set the maximum number of threads for future operations
        void set_max_threads(int max_threads);

        // The main execution method that handles parallel processing
        arma::Col<int> execute_parallel(
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
        );
    };
}

#endif //MODELSELECTION_STRATEGY_H
