//
// Created by jakob on 6/25/24.
//

#ifndef CPP_B_ISM_H
#define CPP_B_ISM_H

#include <Eigen/src/Cholesky/LDLT.h>
#include <Rcpp.h>
#include <RcppEigen.h>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <algorithm>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <random>
#include <set>

#include "./model_selection/ms_base.h"
#include "model_selection/ms_no_optimization.h"
#include "model_selection/ms_parallel_z.h"
#include "model_selection/ms_split_z.h"
#include "mombf/cstat.h"
#include "progressbar.h"

// [[Rcpp::depends(RcppEigen)]]

enum ModelSelectionVersion {
    NO_OPTIMIZATION,
    SPLIT_MATRIX,
    SPLIT_MATRIX_PARALLEL
};

struct BismResults {
    Eigen::MatrixXd b_store;
    Eigen::VectorXd s2_store;
    Eigen::MatrixXd g_store;
    Eigen::MatrixXi w_store;
    Eigen::VectorXd v_store;
    Eigen::MatrixXd w_pip_store;
};

void b_ism(
    const Eigen::MatrixXd &data,
    bool include_constant,
    bool tfe,
    bool ife,
    bool iis,
    bool sis,
    int y_index,
    int i_index,
    int t_index,
    long Ndraw,
    long Nburn,
    std::string b_prior,
    double lambda_b,
    double c0,
    double C0,
    bool geweke,
    BismResults &results,
    ModelSelectionVersion model_selection_version);

Rcpp::List b_ism_wrapper(
    Rcpp::NumericMatrix data,
    bool include_constant,
    bool tfe,
    bool ife,
    bool iis,
    bool sis,
    int y_index,
    int i_index,
    int t_index,
    long Ndraw,
    long Nburn,
    std::string b_prior,
    double lambda_a,
    double lambda_b,
    double lambda_g,
    double c0,
    double C0,
    double oa,
    double ob,
    double va,
    double vb,
    double pip_thr,
    bool w_sis_rand,
    bool geweke);

#endif  // CPP_B_ISM_H
