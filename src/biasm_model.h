//
// Created by jakob on 4/13/25.
//
#ifndef BIASM_MODEL_H
#define BIASM_MODEL_H
#include <string>

#include "bisam_types.h"
#include <RcppArmadillo.h>
#include "utils.h"

namespace bisam {
    BisamResult estimate_model(
        arma::mat &data,
        int i_index,
        int t_index,
        int y_index,
        long long num_draws,
        long long num_burnin,
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
        ComputationStrategy strategy
    );
}


#endif //BIASM_MODEL_H
