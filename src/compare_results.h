//
// Created by jakob on 8/26/24.
//

#ifndef CPP_COMPARE_RESULTS_H
#define CPP_COMPARE_RESULTS_H

#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <cmath>

#ifdef RCPP_EIGEN
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
#else

#include <Eigen/Dense>

#endif
// Assume BismResults struct is defined here

void print_colored(double value, const std::string &color) {
    std::cout << "\033[" << color << "m" << std::setw(15) << std::fixed << std::setprecision(6) << value << "\033[0m";
}

void compare_bism_results(const BismResults &result1, const BismResults &result2,
                          const std::string &id1, const std::string &id2) {
    if (result1.w_store.cols() != result2.w_store.cols()) {
        std::cerr << "Error: w_store matrices have different number of columns." << std::endl;
        return;
    }

    std::cout << std::endl;

    int num_cols = result1.w_store.cols();
    Eigen::VectorXd means1 = result1.w_store.cast<double>().colwise().mean();
    Eigen::VectorXd means2 = result2.w_store.cast<double>().colwise().mean();
    Eigen::VectorXd abs_diff = (means1 - means2).cwiseAbs();

    // Find min and max absolute differences
    int min_abs_idx = 0, max_abs_idx = 0;
    for (int i = 1; i < num_cols; ++i) {
        if (abs_diff[i] < abs_diff[min_abs_idx]) min_abs_idx = i;
        if (abs_diff[i] > abs_diff[max_abs_idx]) max_abs_idx = i;
    }

    // Define threshold for "significantly different from 0"
    const double threshold = 0.09;

    // Print header
    std::cout << std::setw(10) << "Column" << " | "
              << std::setw(15) << id1 << " | "
              << std::setw(15) << id2 << " | "
              << std::setw(15) << "Abs Diff" << std::endl;
    std::cout << std::string(60, '-') << std::endl;

    // Print data
    for (int i = 0; i < num_cols; ++i) {
        std::cout << std::setw(10) << i << " | ";

        // Print Version1 value
        if (std::abs(means1[i]) > threshold) {
            print_colored(means1[i], "32"); // Green for significant
        } else {
            std::cout << std::setw(15) << std::fixed << std::setprecision(6) << means1[i];
        }
        std::cout << " | ";

        // Print Version2 value
        if (std::abs(means2[i]) > threshold) {
            print_colored(means2[i], "32"); // Green for significant
        } else {
            std::cout << std::setw(15) << std::fixed << std::setprecision(6) << means2[i];
        }
        std::cout << " | ";

        // Print Abs Diff
        if (i == min_abs_idx) {
            print_colored(abs_diff[i], "32"); // Green for minimum
        } else if (i == max_abs_idx) {
            print_colored(abs_diff[i], "31"); // Red for maximum
        } else {
            std::cout << std::setw(15) << std::fixed << std::setprecision(6) << abs_diff[i];
        }

        std::cout << std::endl;
    }
}

#endif //CPP_COMPARE_RESULTS_H
