//
// Created by jakob on 4/13/25.
//

#ifndef BISAM_TYPES_H
#define BISAM_TYPES_H

#include <RcppArmadillo.h>

namespace bisam {
    // Output structure (renamed from BIsmOutput)
    struct BisamResult {
        arma::mat beta_samples;           // Previously b_store
        arma::mat gamma_samples;          // Previously g_store
        arma::mat sigma2_samples;         // Previously s2_store
        arma::Mat<int> indicator_samples; // Previously w_store
        arma::rowvec indicator_means;     // Previously w_store_means
        arma::rowvec beta_means;          // Previously b_store_means
        arma::rowvec sigma2_means;        // Previously s2_store_means
    };

    // Prior type (replacing string codes)
    enum class PriorType {
        GAUSSIAN, // Previously "g"
        FIXED,    // Previously "f"
        HORSESHOE // Previously "hs"
    };

    // Selection strategy (replacing ModelSelectionStrategy)
    enum class SelectionStrategy {
        STANDARD, // Previously NO_SPLIT
        SPLIT_SEQUENTIAL,
        SPLIT_PARALLEL
    };

    // Initialization type (replacing InitparType)
    enum class InitType {
        NONE,
        MLE,
        MLE_AISG,
        L1,
        L2_AISGD,
        AUTO
    };

    enum class ComputationStrategy {
        // Previously SelectionStrategy
        STANDARD,         // Single computation on full dataset
        SPLIT_SEQUENTIAL, // Split data, process sequentially
        SPLIT_PARALLEL    // Split data, process in parallel
    };

    // In computation_strategy.h
    struct DataPartition {
        // Previously DataPartition (kept, but in better-named file)
        arma::mat common_x;
        std::vector<arma::vec> y_parts;
        std::vector<arma::Col<int> > delta_init_parts;
        std::vector<arma::vec> theta_init_parts;
        std::vector<size_t> start_columns;
        std::vector<size_t> end_columns;
    };
}
#endif //BISAM_TYPES_H
