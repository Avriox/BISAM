//
// Created by jakob on 4/13/25.
//

#ifndef MOMBF_BRIDGE_H
#define MOMBF_BRIDGE_H
#include <RcppArmadillo.h>

#include "bisam_types.h"
#include "modelSel_regression.h"
#include "LassoRegression.h"

namespace bisam {
    struct ModelSelectionOutput {
        arma::Col<int> post_sample;
        arma::Col<double> margpp;
    };

    ModelSelectionOutput modelSelection(const arma::vec &y,
                                        const arma::mat &x,
                                        int niter,
                                        int thinning,
                                        int burnin,
                                        arma::Col<int> &deltaini_input,
                                        bool center,
                                        bool scale,
                                        bool XtXprecomp,
                                        double phi,
                                        double tau,
                                        double priorSkew,
                                        int prDelta,
                                        double prDeltap,
                                        std::vector<double> parprDeltap,
                                        arma::vec thinit,
                                        InitType initpar_type,
                                        // arma::Col<int> include_vars,
                                        int method      = 0,
                                        int hesstype    = 1,
                                        int optimMethod = 2,
                                        int optim_maxit = 0,
                                        int B           = 100000,
                                        int knownphi    = 1,
                                        int r           = 1,
                                        double alpha    = 0.01,
                                        double lambda   = 0.01
    );

    ModelSelectionOutput modelSelectionGibbsCI(const arma::vec &SpostModeini,
                                               double SpostModeiniProb,
                                               int Sknownphi,
                                               int Sfamily,
                                               int SpriorCoef,
                                               int SpriorGroup,
                                               int Sniter,
                                               int Sthinning,
                                               int Sburnin,
                                               int Sndeltaini,
                                               arma::Col<int> &Sdeltaini,
                                               arma::Col<int> &Sincludevars,
                                               int Sn,
                                               int Sp,
                                               arma::vec &Sy,
                                               int Suncens,
                                               double Ssumy2,
                                               double Ssumy,
                                               double Ssumlogyfact,
                                               arma::mat &Sx,
                                               arma::vec &Scolsumsx,
                                               bool ShasXtX,
                                               arma::mat &SXtX,
                                               arma::rowvec &SytX,
                                               int Smethod,
                                               int Sadjoverdisp,
                                               int Shesstype,
                                               int SoptimMethod,
                                               int Soptim_maxit,
                                               arma::vec Sthinit,
                                               int Susethinit,
                                               int SB,
                                               double Salpha,
                                               double Slambda,
                                               double Sphi,
                                               double Stau,
                                               double Staugroup,
                                               double Staualpha,
                                               double Sfixatanhalpha,
                                               int Sr,
                                               int SpriorDelta,
                                               double SprDeltap,
                                               std::vector<double> SparprDeltap,
                                               int SpriorConstr,
                                               double SprConstrp,
                                               std::vector<double> SparprConstrp,
                                               int *Sgroups,
                                               int Sngroups,
                                               arma::Col<int> &Snvaringroup,
                                               arma::ivec &Sconstraints,
                                               arma::ivec &Sinvconstraints,
                                               int Sverbose);

    arma::vec initParameters(const arma::vec &y, const arma::mat &x, int family,
                             InitType initpar);

    std::pair<arma::vec, int> getthinit(const arma::vec &y, const arma::mat &x, int family,
                                        const arma::vec &initpar, bool enumerate, InitType initpar_type);

    void countConstraints(int *nconstraints,
                          intptrvec *constraints,
                          int *ninvconstraints,
                          intptrvec *invconstraints,
                          int *ngroupsconstr,
                          int *isgroup,
                          int *ngroups,
                          int *nvaringroup,
                          arma::ivec &Sconstraints,
                          arma::ivec &Sinvconstraints);

    void rnlpPost_lm(double *ans, int niter, int burnin, int thinning, double *y, double *x, int n, int p, int r,
                     double tau, double a_phi, double b_phi, int prior, arma::dvec &thinit, bool use_thinit,
                     double phiinit, bool use_phiinit);
}

#endif //MOMBF_BRIDGE_H
