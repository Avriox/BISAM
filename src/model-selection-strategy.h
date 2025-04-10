//
// Created by Jakob Goldmann on 02.04.25.
//

#ifndef MODEL_SELECTION_STRATEGY_H
#define MODEL_SELECTION_STRATEGY_H
#include "mombf-bridge.h"

/* -------------------------------------------------------------------------- */
/*                            Split & Combine Data                            */
/* -------------------------------------------------------------------------- */

// Helper struct to hold split data
struct SplitData {
    arma::mat x_sub;
    std::vector<arma::vec> y_parts;
    std::vector<arma::Col<int> > deltaini_parts;
    std::vector<arma::Col<double> > thinit_parts;
    std::vector<size_t> start_cols;
    std::vector<size_t> end_cols;
};

// Helper function to prepare split data
inline SplitData prepare_split_data(const arma::vec &y, const arma::mat &x, arma::Col<int> &deltaini_input,
                                    arma::vec &thinit, int n) {
    SplitData data;

    size_t n_rows      = y.size();
    size_t n_cols      = x.n_cols;
    size_t thinit_size = thinit.size(); // Check the actual size of thinit

    // Calculate sizes for each part
    size_t rows_per_part  = n_rows / n;
    size_t cols_per_part  = n_cols / n;
    size_t rows_remainder = n_rows % n;
    size_t cols_remainder = n_cols % n;

    // Extract the same submatrix of x for all parts (as per requirement)
    // We'll use the first partition's dimensions
    size_t x_start_row = 0;
    size_t x_end_row   = rows_per_part + (rows_remainder > 0 ? 1 : 0) - 1;
    size_t x_start_col = 0;
    size_t x_end_col   = cols_per_part + (cols_remainder > 0 ? 1 : 0) - 1;

    // Extract the submatrix for x (same for all parts)
    data.x_sub = x.submat(x_start_row, x_start_col, x_end_row, x_end_col);

    // Prepare vectors for all parts
    data.y_parts.resize(n);
    data.deltaini_parts.resize(n);
    data.thinit_parts.resize(n);
    data.start_cols.resize(n);
    data.end_cols.resize(n);

    // Check if thinit needs to be split
    bool split_thinit = (thinit_size == n_cols);

    // Split y and deltaini_input
    for (int part = 0; part < n; part++) {
        // Calculate row range for this part
        size_t start_row = part * rows_per_part + std::min(static_cast<size_t>(part), rows_remainder);
        size_t end_row   = (part + 1) * rows_per_part + std::min(static_cast<size_t>(part + 1), rows_remainder) - 1;

        // Calculate column range for this part
        size_t start_col = part * cols_per_part + std::min(static_cast<size_t>(part), cols_remainder);
        size_t end_col   = (part + 1) * cols_per_part + std::min(static_cast<size_t>(part + 1), cols_remainder) - 1;

        // Extract the corresponding part of y
        data.y_parts[part] = y.subvec(start_row, end_row);

        // Extract the corresponding part of deltaini_input
        data.deltaini_parts[part] = deltaini_input.subvec(start_col, end_col);

        // Handle thinit appropriately based on its size
        if (split_thinit) {
            // If thinit has the same length as x.n_cols, split it accordingly
            data.thinit_parts[part] = thinit.subvec(start_col, end_col);
        } else if (thinit_size > 0) {
            // If thinit is not empty but doesn't match x.n_cols, use the full vector
            data.thinit_parts[part] = thinit;
        } else {
            // If thinit is empty, create an empty vector
            data.thinit_parts[part] = arma::vec();
        }

        // Store the column indices for later reconstruction
        data.start_cols[part] = start_col;
        data.end_cols[part]   = end_col;
    }

    return data;
}


// Helper function to combine results
inline arma::Col<int> combine_results(const std::vector<arma::Col<int> > &results,
                                      const std::vector<size_t> &start_cols,
                                      const std::vector<size_t> &end_cols,
                                      size_t total_cols) {
    arma::Col<int> combined_result(total_cols, arma::fill::zeros);

    for (size_t part = 0; part < results.size(); part++) {
        for (size_t j = 0; j <= (end_cols[part] - start_cols[part]); j++) {
            combined_result(start_cols[part] + j) = results[part](j);
        }
    }

    return combined_result;
}

/* -------------------------------------------------------------------------- */
/*                          Model Selection Strategy                          */
/* -------------------------------------------------------------------------- */

// This is more or less a wrapper for modelSelection(). Which is used to provide a unified
// method interface to the main b_ism function so we can easily switch between Z matrix
// handling strategies.

enum ModelSelectionStrategy {
    NO_SPLIT,
    SPLIT_SEQUENTIAL,
    SPLIT_PARALLEL
};


// Main function with strategy selection
inline arma::Col<int> model_selection_strategy(const arma::vec &y, const arma::mat &x, int niter, int thinning,
                                               int burnin,
                                               arma::Col<int> &deltaini_input, bool center, bool scale,
                                               bool XtXprecomp, double phi, double tau, double priorSkew,
                                               double prDeltap, arma::vec thinit, MombfBridge::InitparType initpar_type,
                                               ModelSelectionStrategy strategy, int n = 3) {
    switch (strategy) {
        case NO_SPLIT:
            // Simple passthrough to modelSelection
            return MombfBridge::modelSelection(y,
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
                                               initpar_type);

        case SPLIT_SEQUENTIAL: {
            // Prepare the split data
            SplitData split_data = prepare_split_data(y, x, deltaini_input, thinit, n);

            // Process each part sequentially
            std::vector<arma::Col<int> > results(n);
            for (int part = 0; part < n; part++) {
                results[part] = MombfBridge::modelSelection(
                    split_data.y_parts[part],
                    split_data.x_sub,
                    niter,
                    thinning,
                    burnin,
                    split_data.deltaini_parts[part],
                    center,
                    scale,
                    XtXprecomp,
                    phi,
                    tau,
                    priorSkew,
                    prDeltap,
                    split_data.thinit_parts[part],
                    initpar_type);
            }

            // Combine and return results
            return combine_results(results, split_data.start_cols, split_data.end_cols, x.n_cols);
        }

        case SPLIT_PARALLEL: {
            // Prepare the split data (same as sequential)
            SplitData split_data = prepare_split_data(y, x, deltaini_input, thinit, n);

            // Process each part in parallel (placeholder for now)
            std::vector<arma::Col<int> > results(n);

            // Here you would add parallel execution code
            // For example using OpenMP:
            // #pragma omp parallel for
            for (int part = 0; part < n; part++) {
                results[part] = MombfBridge::modelSelection(
                    split_data.y_parts[part],
                    split_data.x_sub,
                    niter,
                    thinning,
                    burnin,
                    split_data.deltaini_parts[part],
                    center,
                    scale,
                    XtXprecomp,
                    phi,
                    tau,
                    priorSkew,
                    prDeltap,
                    split_data.thinit_parts[part],
                    initpar_type);
            }

            // Combine and return results (same as sequential)
            return combine_results(results, split_data.start_cols, split_data.end_cols, x.n_cols);
        }

        default:
            throw std::runtime_error("Unknown model selection strategy");
    }
}


#endif //MODEL_SELECTION_STRATEGY_H
