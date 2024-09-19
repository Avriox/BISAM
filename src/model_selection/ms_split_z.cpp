//
// Created by jakob on 9/17/24.
//

#include <iostream>
#include "ms_split_z.h"

std::vector<Eigen::MatrixXi> splitMatrix(const Eigen::MatrixXi &x, int sub_rows, int sub_cols) {
    int n = x.cols() / sub_cols;
    std::vector<Eigen::MatrixXi> result;

    for (int i = 0; i < n; ++i) {
        int start_col = i * sub_cols;
        int start_row = (i == 0) ? 0 : (i * sub_rows);

        Eigen::MatrixXi block = x.block(start_row, start_col, sub_rows, sub_cols);
        result.push_back(block);
    }

    return result;
}

// Function to split VectorXi
std::vector<Eigen::VectorXi> splitVector(const Eigen::VectorXi &x, int sub_size) {
    int n = x.size() / sub_size;
    std::vector<Eigen::VectorXi> result;

    for (int i = 0; i < n; ++i) {
        int start = i * sub_size;
        Eigen::VectorXi block = x.segment(start, sub_size);
        result.push_back(block);
    }

    return result;
}

// Function to split VectorXd
std::vector<Eigen::VectorXd> splitVector(const Eigen::VectorXd &x, int sub_size) {
    int n = x.size() / sub_size;
    std::vector<Eigen::VectorXd> result;

    for (int i = 0; i < n; ++i) {
        int start = i * sub_size;
        Eigen::VectorXd block = x.segment(start, sub_size);
        result.push_back(block);
    }

    return result;
}


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
) {



//    std::cout << "x:" << std::endl << x << std::endl << std::endl;
//
//    std::cout << "y:" << std::endl << y << std::endl << std::endl;
//
//    std::cout << "w_i:" << std::endl << w_i << std::endl << std::endl;

    std::vector<Eigen::MatrixXi> split_xs = splitMatrix(x, n_timeperiods, n_observations);
    std::vector<Eigen::VectorXd> split_ys = splitVector(y, n_timeperiods);
    std::vector<Eigen::VectorXi> split_wis = splitVector(w_i, n_observations);

    Eigen::VectorXi post_samples = Eigen::VectorXi(x.cols());

    for (int i = 0; i < split_xs.size(); i++) {
        Eigen::MatrixXi split_x = split_xs[i];
        Eigen::VectorXd split_y = split_ys[i];
        Eigen::VectorXi split_wi = split_wis[i];

        Eigen::MatrixXi post_sample = model_selection_no_optimization(split_y, split_x, n_iter, prior_coef, prior_delta,
                                                                      phi, split_wi, n_observations, n_timeperiods);

        post_samples.segment(n_observations * i, n_observations) = post_sample;

    }

    return post_samples;


}
