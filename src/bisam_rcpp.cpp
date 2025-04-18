//
// Created by Jakob Goldmann on 08.04.25.
//
// #include <Rcpp.h>
#include <RcppArmadillo.h>
#include "biasm_model.h"
#include "bisam_types.h"

// TODO HARD CODED VALUES!!!

// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List rcpp_estimate_model(
    arma::mat data,
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
    bool use_phiinit,
    bool const_val,
    bool ife,
    bool tfe,
    bool iis,
    bool sis,
    int new_par_method,
    int new_par_hesstype,
    int new_par_optim_method,
    int new_par_optim_maxit,
    int new_par_B,
    int new_par_knownphi,
    int new_par_r,
    double new_par_alpha,
    double new_par_lambda,
    int computation_strategy
) {
    // Call the C++ function

    bisam::BisamResult result = bisam::estimate_model(
        data,
        i_index,
        t_index,
        y_index,
        Ndraw,
        Nburn,
        b_prior,
        lambda_b,
        c0,
        C0,
        va,
        vb,
        tau,
        use_phiinit,
        const_val,
        ife,
        tfe,
        iis,
        sis,
        new_par_method,
        new_par_hesstype,
        new_par_optim_method,
        new_par_optim_maxit,
        new_par_B,
        new_par_knownphi,
        new_par_r,
        new_par_alpha,
        new_par_lambda,
        static_cast<bisam::ComputationStrategy>(computation_strategy)
    );

    // Convert the output struct to an R list
    Rcpp::List ret;
    ret["b_store"]        = result.beta_samples;
    ret["g_store"]        = result.gamma_samples;
    ret["s2_store"]       = result.sigma2_samples;
    ret["w_store"]        = result.indicator_samples;
    ret["w_store_means"]  = result.indicator_means;
    ret["b_store_means"]  = result.beta_means;
    ret["s2_store_means"] = result.sigma2_means;

    return ret;
}

// Export enum values as integers
// [[Rcpp::export]]
int comp_strategy_standard() {
    return static_cast<int>(bisam::ComputationStrategy::STANDARD);
}

// [[Rcpp::export]]
int comp_strategy_split_sequential() {
    return static_cast<int>(bisam::ComputationStrategy::SPLIT_SEQUENTIAL);
}

// [[Rcpp::export]]
int comp_strategy_split_parallel() {
    return static_cast<int>(bisam::ComputationStrategy::SPLIT_PARALLEL);
}
