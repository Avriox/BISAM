//
// Created by jakob on 9/17/24.
//

#ifndef CPP_MS_SPLIT_Z_H
#define CPP_MS_SPLIT_Z_H

#include "./ms_base.h"
#include "../mombf/modselIntegrals.h"
#include "./ms_no_optimization.h"

std::vector<Eigen::MatrixXi> splitMatrix(const Eigen::MatrixXi &x, int sub_rows, int sub_cols);

// Function to split VectorXi
std::vector<Eigen::VectorXi> splitVector(const Eigen::VectorXi &x, int sub_size);

// Function to split VectorXd
std::vector<Eigen::VectorXd> splitVector(const Eigen::VectorXd &x, int sub_size);

Eigen::VectorXi model_selection_split_z(
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


#endif //CPP_MS_SPLIT_Z_H
