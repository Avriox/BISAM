//
// Created by Jakob Goldmann on 08.04.25.
//

#ifndef B_ISM_H
#define B_ISM_H
#include <RcppArmadillo.h>

#include "model-selection-strategy.h"

// Define struct for b_ism output
struct BIsmOutput {
    arma::mat b_store;           // Store for beta coefficients
    arma::mat g_store;           // Store for gamma coefficients
    arma::mat s2_store;          // Store for variance
    arma::Mat<int> w_store;      // Store for model indicators
    arma::rowvec w_store_means;  // Mean of model indicators
    arma::rowvec b_store_means;  // Mean of beta coefficients
    arma::rowvec s2_store_means; // Mean of variance
};


BIsmOutput b_ism(
    arma::mat &data,
    int i_index,
    int t_index,
    int y_index,
    long long Ndraw,
    long long Nburn,
    std::string b_prior,
    double lambda_b,
    double c0,
    double C0,
    double va,
    double vb,
    double tau,
    bool geweke,
    bool use_phiinit,
    bool const_val,
    bool ife,
    bool tfe,
    bool iis,
    bool sis
);


#endif //B_ISM_H
