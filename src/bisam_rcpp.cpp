//
// Created by Jakob Goldmann on 08.04.25.
//
// #include <Rcpp.h>
#include <RcppArmadillo.h>
#include "b_ism.h"


// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List b_ism_rcpp(
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
    bool geweke,
    bool use_phiinit,
    bool const_val,
    bool ife,
    bool tfe,
    bool iis,
    bool sis
) {
    // Call the C++ function
    BIsmOutput result = b_ism(
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
        geweke,
        use_phiinit,
        const_val,
        ife,
        tfe,
        iis,
        sis
    );

    // Convert the output struct to an R list
    Rcpp::List ret;
    ret["b_store"]        = result.b_store;
    ret["g_store"]        = result.g_store;
    ret["s2_store"]       = result.s2_store;
    ret["w_store"]        = result.w_store;
    ret["w_store_means"]  = result.w_store_means;
    ret["b_store_means"]  = result.b_store_means;
    ret["s2_store_means"] = result.s2_store_means;

    return ret;
}
