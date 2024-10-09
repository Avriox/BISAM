#include "./ms_base.h"
#include "mombf/modselIntegrals.h"

using namespace Rcpp;

// [[Rcpp::export]]
msPriorSpec::msPriorSpec(std::string priorType, std::string priorDistr, std::map<std::string, double> priorPars)
        : priorType(priorType), priorDistr(priorDistr), priorPars(priorPars) {}

// Convert msPriorSpec to R list
List msPriorSpec_to_list(const msPriorSpec &obj) {
    return List::create(
            _["priorType"] = obj.priorType,
            _["priorDistr"] = obj.priorDistr,
            _["priorPars"] = wrap(obj.priorPars)
    );
}

// Convert R list to msPriorSpec
msPriorSpec list_to_msPriorSpec(List lst) {
    std::string priorType = as<std::string>(lst["priorType"]);
    std::string priorDistr = as<std::string>(lst["priorDistr"]);
    std::map<std::string, double> priorPars;
    List pars = lst["priorPars"];
    CharacterVector names = pars.names();
    for (int i = 0; i < pars.size(); i++) {
        std::string name = as<std::string>(names[i]);
        double value = as<double>(pars[i]);
        priorPars[name] = value;
    }
    return msPriorSpec(priorType, priorDistr, priorPars);
}

msPriorSpec imomprior_cpp(double tau, double tau_adj) {
    std::map<std::string, double> priorPars = {
            {"tau",     tau},
            {"tau.adj", tau_adj}
    };

    return msPriorSpec("coefficients", "piMOM", priorPars);
}


msPriorSpec modelbbprior_cpp(double alpha_p, double beta_p) {
    std::map<std::string, double> priorPars = {
            {"alpha.p", alpha_p},
            {"beta.p",  beta_p}
    };

    return msPriorSpec("modelIndicator", "binomial", priorPars);
}

msPriorSpec igprior_cpp(double alpha, double lambda) {
    std::map<std::string, double> priorPars = {
            {"alpha",  alpha},
            {"lambda", lambda}
    };

    return msPriorSpec("nuisancePars", "invgamma", priorPars);
}

msPriorSpec momprior_cpp(double tau, double tau_adj, double r) {
    std::map<std::string, double> priorPars = {
            {"tau",     tau},
            {"tau_adj", tau_adj},
            {"r",       r}
    };

    return msPriorSpec("nuisancePars", "invgamma", priorPars);
}

// R-exportable versions
// [[Rcpp::export]]
List imomprior(double tau, double tau_adj = 1e6) {
    return msPriorSpec_to_list(imomprior_cpp(tau, tau_adj));
}

// [[Rcpp::export]]
List modelbbprior(double alpha_p = 1.0, double beta_p = 1.0) {
    return msPriorSpec_to_list(modelbbprior_cpp(alpha_p, beta_p));
}

// [[Rcpp::export]]
List igprior(double alpha = 1.0, double lambda = 1.0) {
    return msPriorSpec_to_list(igprior_cpp(alpha, lambda));
}

// [[Rcpp::export]]
List momprior(double tau = 1.0, double tau_adj = 1e6, double r = 1) {
    return msPriorSpec_to_list(momprior_cpp(tau, tau_adj, r));
}


void count_constraints(Eigen::VectorXi &n_constraints, intptrvec *constraints,
                       Eigen::VectorXi &n_inv_constraints, intptrvec *inv_constraints,
                       int *n_groups_constr, Eigen::VectorXi &is_group, int n_groups,
                       const Eigen::VectorXi &n_var_in_group,
                       const std::vector<std::vector<int>> &s_constraints,
                       const std::vector<std::vector<int>> &s_inv_constraints) {
    int jj = 0;
    *n_groups_constr = 0;
    for (int j = 0; j < n_groups; j++) {
        n_constraints(j) = s_constraints[j].size();
        n_inv_constraints(j) = s_inv_constraints[j].size();
        constraints->push_back(const_cast<int *>(s_constraints[j].data()));
        inv_constraints->push_back(const_cast<int *>(s_inv_constraints[j].data()));
        if (n_constraints(j) > 0) (*n_groups_constr)++;
        is_group(jj) = (n_var_in_group(j) > 1);
        jj++;
        for (int i = 1; i < n_var_in_group(j); i++, jj++) is_group(jj) = is_group(jj - 1);
    }
}

GroupConstraintInfo
code_groups_and_constraints(int p, const Eigen::VectorXi &groups,
                            const std::vector<std::vector<int>> &constraints) {
    Eigen::VectorXi groupsnum = groups.array() + 1;
    Eigen::VectorXi unique_groups = groupsnum.segment(0, 1);
    for (int i = 1; i < groupsnum.size(); ++i) {
        if (groupsnum[i] != groupsnum[i - 1]) {
            unique_groups.conservativeResize(unique_groups.size() + 1);
            unique_groups[unique_groups.size() - 1] = groupsnum[i];
        }
    }
    Eigen::VectorXi groupsnum_recoded = Eigen::VectorXi::Zero(groupsnum.size());
    for (int i = 0; i < groupsnum.size(); ++i) {
        groupsnum_recoded[i] =
                (std::find(unique_groups.data(), unique_groups.data() + unique_groups.size(), groupsnum[i]) -
                 unique_groups.data()) + 1;
    }

    int ngroups = groupsnum_recoded.maxCoeff();
    if (ngroups > p) {
        throw std::runtime_error("There cannot be more groups than variables (columns in x)");
    }

    std::vector<std::vector<int>> constraints_processed;
    if (constraints.empty()) {
        constraints_processed.resize(ngroups);
    } else {
        if (constraints.size() != ngroups) {
            throw std::runtime_error("constraints.size() must be equal to number of variable groups");
        }
        constraints_processed = constraints;
        for (auto &constraint: constraints_processed) {
            for (auto &c: constraint) {
                c = groupsnum_recoded[c] - 1;
            }
        }
    }

    Eigen::VectorXi nvaringroup;
    Eigen::VectorXi groups_processed;
    if (ngroups == p) {
        nvaringroup = Eigen::VectorXi::Ones(p);
        groups_processed = Eigen::VectorXi::LinSpaced(p, 0, p - 1);
    } else {
        std::map<int, int> group_counts;
        for (int i = 0; i < groupsnum_recoded.size(); ++i) {
            group_counts[groupsnum_recoded[i]]++;
        }
        nvaringroup.resize(ngroups);
        for (int i = 0; i < ngroups; ++i) {
            nvaringroup[i] = group_counts[i + 1];
        }
        groups_processed = groupsnum_recoded.array() - 1;
    }

    std::vector<std::vector<int>> invconstraints(ngroups);
    std::vector<std::pair<int, int>> tabconstr;
    for (int i = 0; i < constraints_processed.size(); ++i) {
        for (int j: constraints_processed[i]) {
            tabconstr.push_back({i, j});
        }
    }
    for (int i = 0; i < ngroups; ++i) {
        for (const auto &pair: tabconstr) {
            if (pair.second == i) {
                invconstraints[i].push_back(pair.first);
            }
        }
    }

    return {ngroups, constraints_processed, invconstraints, nvaringroup, groups_processed};
}


std::pair<int, int> get_family_info(const std::string &family, bool issurvival) {
    int familyint, familygreedy;

    if (family == "auto") {
        familyint = 0;
        familygreedy = 1;
    } else if (family == "normal") {
        familyint = familygreedy = (!issurvival) ? 1 : 11;
    } else if (family == "twopiecenormal") {
        familyint = 2;
        familygreedy = 1;
    } else if (family == "laplace") {
        familyint = 3;
        familygreedy = 1;
    } else if (family == "twopiecelaplace") {
        familyint = 4;
        familygreedy = 1;
    } else if (family == "binomial" || family == "binomial logit") {
        familyint = familygreedy = 21;
    } else if (family == "poisson" || family == "poisson log") {
        familyint = familygreedy = 22;
    } else {
        throw std::runtime_error("family not available");
    }

    return std::make_pair(familyint, familygreedy);
}


FormatMsPriorsMargResult format_ms_priors_marg(const msPriorSpec &priorCoef, const msPriorSpec &priorGroup,
                                               int n) { // const msPriorSpec &priorVar, const msPriorSpec &priorSkew,

    int r = 1;

    double tau = priorCoef.priorPars.at("tau");
    double taugroup = priorGroup.priorPars.at("tau");

    int prior = 1;
    int priorgr = 1;

//    double alpha = priorVar.priorPars.at("alpha");
//    double lambda = priorVar.priorPars.at("lambda");
//
//    double taualpha = priorSkew.priorPars.at("tau");
    double fixatanhalpha = -10000.0;

    FormatMsPriorsMargResult result{
            r,
            prior,
            priorgr,
            tau,
            taugroup,
//            alpha,
//            lambda,
//            taualpha,
//            fixatanhalpha,
            priorCoef,
            priorGroup
    };

    return result;
}

