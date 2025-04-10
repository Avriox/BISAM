//
// Created by jakob on 3/19/25.
//

#ifndef CONTROLLED_SIMULATION_H
#define CONTROLLED_SIMULATION_H

#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <armadillo>
#include <string>
#include <iomanip>
#include <optional>

/* -------------------------------------------------------------------------- */
/*                           Controlled Simulation                            */
/* -------------------------------------------------------------------------- */

// This should be a one to one translation of the original file `contr_sim_breaks_fun.R`. It has been auto
// translated using gemini.

// Helper function to create a sequence of repeated values
inline arma::ivec rep_seq(int value, int times) {
    arma::ivec result(times);
    result.fill(value);
    return result;
}

// Helper function for Kronecker product of two matrices
inline arma::mat kronecker_product(const arma::mat &A, const arma::mat &B) {
    arma::mat result(A.n_rows * B.n_rows, A.n_cols * B.n_cols);
    for (arma::uword i = 0; i < A.n_rows; ++i) {
        for (arma::uword j = 0; j < A.n_cols; ++j) {
            result.submat(i * B.n_rows, j * B.n_cols, (i + 1) * B.n_rows - 1, (j + 1) * B.n_cols - 1) = A(i, j) * B;
        }
    }
    return result;
}

struct SimulationOutput {
    arma::mat data;
    arma::vec true_b;
    arma::vec errors;
    std::optional<double> true_const;
    std::optional<arma::vec> true_ife;
    std::optional<arma::vec> true_tfe;
    arma::mat tr_idx; // Treated index information
};

inline SimulationOutput contr_sim_breaks(
    int n,
    int t,
    int nx,
    bool iis                   = true,
    bool sis                   = true,
    int pos_outl               = 0,  // count in total position of panel setup
    const arma::ivec &pos_step = {}, // count in total position of panel setup
    bool const_val             = false,
    bool ife                   = false,
    bool tfe                   = false,
    int outl_mean              = 0,
    int step_mean              = 5,
    double error_sd            = 1.0
) {
    std::mt19937 rng(std::rand()); // Use Mersenne Twister engine
    std::normal_distribution<> normal_dist(0.0, error_sd);
    std::uniform_real_distribution<> uniform_dist(-10.0, 10.0);
    std::uniform_int_distribution<> uniform_int_dist(-10, 10);

    int n_obs = n * t;

    // Create individual and time indicators
    arma::ivec n_(n_obs);
    arma::ivec t_(n_obs);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < t; ++j) {
            n_(i * t + j) = i + 1;
            t_(i * t + j) = j + 1;
        }
    }

    arma::mat X_;
    arma::vec b_;
    arma::mat X;
    arma::vec b;

    // Create obs. matrices and corresponding betas
    if (nx > 0) {
        X_ = arma::randn(n_obs, nx);
        X  = X_;
        b_ = arma::vec(nx);
        for (int i = 0; i < nx; ++i) {
            b_(i) = uniform_int_dist(rng);
        }
        b = b_;
    } else {
        X_ = arma::mat();
        X  = arma::mat();
        b_ = arma::vec();
        b  = arma::vec();
    }

    std::optional<double> b0;
    std::optional<arma::vec> bife;
    std::optional<arma::vec> btfe;

    if (const_val) {
        // Intercept
        arma::vec const_col(n_obs, arma::fill::ones);
        X.insert_cols(X.n_cols, const_col);
        b0 = uniform_int_dist(rng);
        b.resize(b.n_elem + 1);
        b.tail(1)(0) = *b0;
    }

    if (ife) {
        // Individual fixed effects
        arma::mat IFE = kronecker_product(arma::eye(n, n), arma::ones(t, 1));
        if (const_val) {
            IFE = IFE.cols(1, IFE.n_cols - 1);
        }
        X.insert_cols(X.n_cols, IFE);
        int n_ife = IFE.n_cols;
        bife      = arma::vec(n_ife);
        for (int i = 0; i < n_ife; ++i) {
            (*bife)(i) = uniform_int_dist(rng);
        }
        b.resize(b.n_elem + n_ife);
        b.tail(n_ife) = *bife;
    }

    if (tfe) {
        // Time fixed effects
        arma::mat TFE = kronecker_product(arma::ones(n, 1), arma::eye(t, t));
        if (const_val) {
            TFE = TFE.cols(1, TFE.n_cols - 1);
        }
        X.insert_cols(X.n_cols, TFE);
        int n_tfe = TFE.n_cols;
        btfe      = arma::vec(n_tfe);
        for (int i = 0; i < n_tfe; ++i) {
            (*btfe)(i) = uniform_int_dist(rng);
        }
        b.resize(b.n_elem + n_tfe);
        b.tail(n_tfe) = *btfe;
    }

    if (tfe && ife && !const_val) {
        std::cerr <<
                "Warning: Both time and unit fixed effects used.\n Dropping first indiv. FE dummy to avoid perfect colinearity"
                << std::endl;
        if (nx > 0) {
            X.shed_col(nx);
        } else if (const_val) {
            X.shed_col(0);
        } else if (ife) {
            X.shed_col(0);
        }
        if (!b.empty()) {
            b.shed_row(b.n_elem - 1);
        }
    }

    // IIS
    arma::mat I = arma::eye(n_obs, n_obs);
    arma::vec a = arma::zeros(n_obs);
    if (pos_outl > 0 && pos_outl <= n_obs) {
        a(pos_outl - 1) = outl_mean;
    }

    // SIS
    int n_step_cols = n * (t - 1 - (iis ? 1 : 0));
    arma::mat Z     = arma::zeros(n_obs, n_step_cols);
    if (t > 1 + (iis ? 1 : 0)) {
        arma::mat lower_tri = arma::trimatu(arma::ones(t, t), 0);
        lower_tri.shed_col(0);
        if (iis) {
            lower_tri.shed_col(lower_tri.n_cols - 1);
        }
        Z = kronecker_product(arma::eye(n, n), lower_tri);
    }
    arma::vec g = arma::zeros(n_step_cols);
    for (const auto &pos: pos_step) {
        if (pos > 0 && pos <= n_obs) {
            int unit_index = (pos - 1) / t;
            int time_index = (pos - 1) % t;
            if (time_index > 0 && time_index < t - (iis ? 1 : 0)) {
                int col_index = unit_index * (t - 1 - (iis ? 1 : 0)) + time_index - 1;
                if (col_index >= 0 && col_index < g.n_elem) {
                    g(col_index) = step_mean;
                }
            }
        }
    }

    arma::vec e(n_obs);
    for (int i = 0; i < n_obs; ++i) {
        e(i) = normal_dist(rng);
    }

    arma::vec y;
    if (X.is_empty()) {
        y = I * a + Z * g + e;
    } else if (b.is_empty()) {
        y = I * a + Z * g + e;
    } else {
        y = X * b + I * a + Z * g + e;
    }

    arma::mat tr_ind_arma;
    arma::mat tr_stp_arma;

    arma::umat a_nonzero_idx = arma::find(a != 0);
    if (!a_nonzero_idx.empty()) {
        arma::mat a_reshaped = arma::reshape(a, t, n).t();
        std::vector<arma::uvec> rows_cols;
        std::vector<double> values;
        for (arma::uword i = 0; i < a_reshaped.n_rows; ++i) {
            for (arma::uword j = 0; j < a_reshaped.n_cols; ++j) {
                if (a_reshaped(i, j) != 0) {
                    rows_cols.push_back({j + 1, i + 1});
                    values.push_back(a_reshaped(i, j));
                }
            }
        }
        if (!rows_cols.empty()) {
            tr_ind_arma.set_size(rows_cols.size(), 3);
            for (size_t i = 0; i < rows_cols.size(); ++i) {
                tr_ind_arma(i, 0) = rows_cols[i](0); // column (time)
                tr_ind_arma(i, 1) = rows_cols[i](1); // row (unit)
                tr_ind_arma(i, 2) = values[i];       // value
            }
        }
    }

    arma::umat g_nonzero_idx = arma::find(g != 0);
    if (!g_nonzero_idx.empty()) {
        arma::mat g_padded = arma::join_cols(arma::zeros(n, 1), arma::reshape(g, (t - 1 - (iis ? 1 : 0)), n).t());
        std::vector<arma::uvec> rows_cols_step;
        std::vector<double> values_step;
        for (arma::uword i = 0; i < g_padded.n_rows; ++i) {
            for (arma::uword j = 0; j < g_padded.n_cols; ++j) {
                if (g_padded(i, j) != 0) {
                    rows_cols_step.push_back({j + 1, i + 1});
                    values_step.push_back(g_padded(i, j));
                }
            }
        }
        if (!rows_cols_step.empty()) {
            tr_stp_arma.set_size(rows_cols_step.size(), 3);
            for (size_t i = 0; i < rows_cols_step.size(); ++i) {
                tr_stp_arma(i, 0) = rows_cols_step[i](0); // column (time)
                tr_stp_arma(i, 1) = rows_cols_step[i](1); // row (unit)
                tr_stp_arma(i, 2) = values_step[i];       // value
            }
        }
    }

    arma::mat treated;
    arma::vec treated_indices;
    arma::vec treated_sizes;

    if (!tr_ind_arma.empty()) {
        arma::mat temp_ind(tr_ind_arma.n_rows, 2);
        temp_ind.col(0) = tr_ind_arma.col(2);                                // size
        temp_ind.col(1) = (tr_ind_arma.col(1) - 1) * t + tr_ind_arma.col(0); // index
        if (treated.is_empty()) {
            treated = temp_ind;
        } else {
            treated = arma::join_cols(treated, temp_ind);
        }
    }
    if (!tr_stp_arma.empty()) {
        arma::mat temp_stp(tr_stp_arma.n_rows, 2);
        temp_stp.col(0) = tr_stp_arma.col(2);                                // size
        temp_stp.col(1) = (tr_stp_arma.col(1) - 1) * t + tr_stp_arma.col(0); // index
        if (treated.is_empty()) {
            treated = temp_stp;
        } else {
            treated = arma::join_cols(treated, temp_stp);
        }
    }

    arma::mat treated_with_effects;
    if (!treated.is_empty()) {
        treated_with_effects.set_size(treated.n_rows, 3);
        treated_with_effects.col(0) = treated.col(0); // size
        treated_with_effects.col(1) = treated.col(1); // index
        arma::vec e_subset(treated.n_rows);
        for (arma::uword i = 0; i < treated.n_rows; ++i) {
            e_subset(i) = e(treated(i, 1) - 1);
        }
        treated_with_effects.col(2) = (treated.col(0) + e_subset) / error_sd; // rel_net_eff
    }

    arma::vec errors(2 * n_obs);
    arma::vec iis_errors                = e;
    arma::vec sis_errors                = e;
    errors.subvec(0, n_obs - 1)         = iis_errors;
    errors.subvec(n_obs, 2 * n_obs - 1) = sis_errors;

    arma::mat data(n_obs, 3 + (nx > 0 ? nx : 0));
    data.col(0) = arma::conv_to<arma::colvec>::from(n_);
    data.col(1) = arma::conv_to<arma::colvec>::from(t_);
    data.col(2) = y;
    if (nx > 0) {
        data.cols(3, 3 + nx - 1) = X_;
    }

    SimulationOutput sim_output;
    sim_output.data       = data;
    sim_output.true_b     = b_;
    sim_output.errors     = errors;
    sim_output.true_const = b0;
    sim_output.true_ife   = bife;
    sim_output.true_tfe   = btfe;
    if (!treated_with_effects.is_empty()) {
        sim_output.tr_idx = treated_with_effects;
    } else {
        sim_output.tr_idx = arma::mat();
    }

    return sim_output;
}

#endif //CONTROLLED_SIMULATION_H
