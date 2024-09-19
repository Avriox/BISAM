//
// Created by jakob on 8/17/24.
//

#ifndef CPP_MS_BASE_H
#define CPP_MS_BASE_H

#ifdef RCPP_EIGEN
// R environment
#include <RcppEigen.h>
  // [[Rcpp::depends(RcppEigen)]]
#else
// Pure C++ environment
#include <Eigen/Core>
#include <Eigen/Dense>

#endif

// Standard library includes (common to both environments)
#include <map>
#include <set>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <cstdint>
#include <iterator>
#include <memory>

// Local includes
#include "./sel2selnew_optimization.h"

struct ModelSelectionResult {
    Eigen::VectorXi post_sample;
    Eigen::VectorXd marg_pp;
    Eigen::VectorXd post_mode;
    double post_mode_prob;
    Eigen::VectorXd post_prob;
    std::string model_id;
    Eigen::VectorXd *post_mean;
    Eigen::VectorXd *post_var;
    std::string family;
    int p;
    bool enumerate;
    // Priors struct definition omitted for brevity
    Eigen::VectorXd y_std;
    Eigen::MatrixXd x_std;
    Eigen::VectorXi groups;
    std::vector<std::vector<int>> constraints;
    Eigen::MatrixXd std_constants;
    std::string type;
};

struct FormatMsMethodResult {
    int method;
    int optimMethod;
    int optim_maxit;
    int adj_overdisp;
    int hesstype;
};

class msPriorSpec {
public:
    msPriorSpec(std::string priorType, std::string priorDistr, std::map<std::string, double> priorPars);

    std::string priorType;
    std::string priorDistr;
    std::map<std::string, double> priorPars;
};

msPriorSpec imomprior(double tau, double tau_adj = 1e6);

msPriorSpec modelbbprior(double alpha_p = 1.0, double beta_p = 1.0);

msPriorSpec igprior(double alpha = 1.0, double lambda = 1.0);

msPriorSpec momprior(double tau = 1.0, double tau_adj = 1e6, double r = 1);

typedef std::vector<int *> intptrvec;


struct FormatMsPriorsMargResult {
    int r;
    int prior;
    int priorgr;
    double tau;
    double taugroup;
    double alpha;
    double lambda;
    double taualpha;
    double fixatanhalpha;
    msPriorSpec priorCoef;
    msPriorSpec priorGroup;
};

struct Priors {
    msPriorSpec priorCoef;
    msPriorSpec priorGroup;
    msPriorSpec priorDelta;
    msPriorSpec priorConstraints;
    msPriorSpec priorVar;
    msPriorSpec priorSkew;
};

struct GroupConstraintInfo {
    int ngroups;
    std::vector<std::vector<int>> constraints;
    std::vector<std::vector<int>> invconstraints;
    Eigen::VectorXi nvaringroup;
    Eigen::VectorXi groups;
};

GroupConstraintInfo
code_groups_and_constraints(int p, const Eigen::VectorXi &groups,
                            const std::vector<std::vector<int>> &constraints = {});

std::pair<int, int> get_family_info(const std::string &family, bool issurvival);


FormatMsPriorsMargResult format_ms_priors_marg(const msPriorSpec &priorCoef, const msPriorSpec &priorGroup,
                                               const msPriorSpec &priorVar, const msPriorSpec &priorSkew, int n);

void count_constraints(Eigen::VectorXi &n_constraints, intptrvec *constraints,
                       Eigen::VectorXi &n_inv_constraints, intptrvec *inv_constraints,
                       int *n_groups_constr, Eigen::VectorXi &is_group, int n_groups,
                       const Eigen::VectorXi &n_var_in_group,
                       const std::vector<std::vector<int>> &s_constraints,
                       const std::vector<std::vector<int>> &s_inv_constraints);


#endif //CPP_MODELSELECTION_H
