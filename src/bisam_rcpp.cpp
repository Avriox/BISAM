//
// Created by Jakob Goldmann on 08.04.25.
//
// #include <Rcpp.h>
#include <RcppArmadillo.h>
#include "biasm_model.h"
#include "bisam_types.h"
#include "modelselection_strategy.h"

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

    int prDelta,
    double p,
    double a,
    double b,

    int computation_strategy
) {
    // Call the C++ function

    // SInce include vars does currently not work, i hard code it here
    // arma::Col<int> new_par_include_vars,
    arma::Col<int> new_par_include_vars;

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
        // new_par_include_vars,
        new_par_method,
        new_par_hesstype,
        new_par_optim_method,
        new_par_optim_maxit,
        new_par_B,
        new_par_knownphi,
        new_par_r,
        new_par_alpha,
        new_par_lambda,
        prDelta,
        p,
        {a, b},
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

// [[Rcpp::export]]
Rcpp::List rcpp_fast_model_selection(
    const arma::vec &y,
    const arma::mat &x,
    int niter                = 5000,
    int thinning             = 1,
    int burnin               = -1,         // -1 => round(niter/10)
    SEXP deltaini            = R_NilValue, // NULL or int/logical length p
    bool center              = true,
    bool scale               = true,
    bool XtXprecomp          = true,
    double phi               = 1.0,
    double tau               = 0.348,
    double priorSkew         = 0.348,
    SEXP thinit              = R_NilValue, // NULL or numeric length p (optional)
    int initpar_type         = 0,          // cast to bisam::InitType
    int method               = 1,
    int hesstype             = 1,
    int optimMethod          = 1,
    int optim_maxit          = 10,
    int B                    = 100000,
    int knownphi             = 0,
    int r                    = 1,
    double alpha             = 0.05,
    double lambda            = 1.0,
    int prDelta              = 2,   // 1=binomial, 2=beta-binomial (your convention)
    double prDeltap          = 0.5, // used if prDelta==1
    double prDelta_a         = 1.0, // used if prDelta==2
    double prDelta_b         = 1.0, // used if prDelta==2
    int computation_strategy = 0,   // STANDARD / SPLIT_* as int
    int n_units              = 1,   // used for split strategies
    int max_threads          = 0    // used for SPLIT_PARALLEL
) {
    if (x.n_rows != y.n_elem) {
        Rcpp::stop("Dimension mismatch: nrow(x) must equal length(y).");
    }
    if (niter <= 0) Rcpp::stop("niter must be > 0.");
    if (thinning <= 0) Rcpp::stop("thinning must be > 0.");
    if (burnin < -1) Rcpp::stop("burnin must be >= 0, or -1 for default.");
    if (n_units <= 0) Rcpp::stop("n_units must be > 0.");

    if (burnin == -1) burnin = static_cast<int>(std::llround(niter / 10.0));

    const arma::uword p = x.n_cols;

    // ---- deltaini: manual conversion R -> arma::Col<int> (0/1)
    arma::Col<int> delta_init(p, arma::fill::zeros);
    if (!Rf_isNull(deltaini)) {
        // This coerces logical -> integer if needed
        Rcpp::IntegerVector d = Rcpp::as<Rcpp::IntegerVector>(deltaini);
        if (static_cast<arma::uword>(d.size()) != p) {
            Rcpp::stop("deltaini must have length ncol(x).");
        }
        for (arma::uword j = 0; j < p; ++j) {
            delta_init(j) = d[j];
        }
    }

    // ---- thinit: optional numeric vector length p (or empty)
    arma::vec thinit_arma;
    if (!Rf_isNull(thinit)) {
        Rcpp::NumericVector tv = Rcpp::as<Rcpp::NumericVector>(thinit);
        if (tv.size() != 0 && static_cast<arma::uword>(tv.size()) != p) {
            Rcpp::stop("thinit must be NULL/empty or have length ncol(x).");
        }
        thinit_arma = arma::vec(tv.begin(), tv.size(), /*copy_aux_mem=*/true);
    } else {
        thinit_arma.reset(); // empty
    }

    std::vector<double> parprDeltap = {prDelta_a, prDelta_b};

    bisam::ModelSelectionOutput out = bisam::model_selection_with_strategy(
        y,
        x,
        niter,
        thinning,
        burnin,
        delta_init,
        center,
        scale,
        XtXprecomp,
        phi,
        tau,
        priorSkew,
        thinit_arma,
        static_cast<bisam::InitType>(initpar_type),
        method,
        hesstype,
        optimMethod,
        optim_maxit,
        B,
        knownphi,
        r,
        alpha,
        lambda,
        static_cast<bisam::ComputationStrategy>(computation_strategy),
        n_units,
        prDelta,
        prDeltap,
        parprDeltap,
        max_threads
    );

    return Rcpp::List::create( // CHANGED: return both fields
        Rcpp::Named("post_sample") =
        Rcpp::IntegerVector(out.post_sample.begin(), out.post_sample.end()),
        Rcpp::Named("margpp") =
        Rcpp::NumericVector(out.margpp.begin(), out.margpp.end())
    );
}
