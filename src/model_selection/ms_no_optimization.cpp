//
// Created by jakob on 9/13/24.
//


#include "ms_no_optimization.h"


Eigen::VectorXi model_selection_no_optimization(
        Eigen::VectorXd y,
        Eigen::MatrixXi x,
        int n_iter,
        msPriorSpec prior_coef,
        msPriorSpec prior_delta,
        double phi,
        Eigen::VectorXi w_i,
        int n_observations,
        int n_timeperiods
) {


    int thinning = 1;
    int burnin = 1;
    std::string family = "normal";

    int b = 10e5;

    int p = x.cols();
    int n = y.size();

    // Create groups vector
    Eigen::VectorXi groups = Eigen::VectorXi::LinSpaced(p, 1, p);

    // Create include_vars vector
//    Eigen::Array<bool, Eigen::Dynamic, 1> include_vars = Eigen::Array<bool, Eigen::Dynamic, 1>::Constant(p, false);
    Eigen::VectorXi include_vars = Eigen::VectorXi::Zero(p);

    // Code groups and constraints
    GroupConstraintInfo group_info = code_groups_and_constraints(p, groups);
    int n_groups = group_info.ngroups;
    std::vector<std::vector<int>> constraints = group_info.constraints;
    std::vector<std::vector<int>> inv_constraints = group_info.invconstraints;
    Eigen::VectorXi n_var_in_group = group_info.nvaringroup;
    groups = group_info.groups;

    // Calculate column means and standard deviations
    Eigen::VectorXd mx = x.cast<double>().colwise().mean();
    Eigen::VectorXd mean_square = x.cast<double>().array().square().colwise().mean();
    Eigen::VectorXd mx_squared = mx.array().square();
    Eigen::VectorXd difference = mean_square - mx_squared;
    Eigen::VectorXd sqrt_difference = difference.array().sqrt();
    double correction_factor = std::sqrt(static_cast<double>(n) / (n - 1));
    Eigen::VectorXd sx = sqrt_difference * correction_factor;

    // Check for constant columns
    Eigen::Array<bool, Eigen::Dynamic, 1> ct = (sx.array() == 0);
    if (ct.cast<int>().sum() > 1) {
        throw std::runtime_error("There are >1 constant columns in x (e.g. two intercepts)");
    }

    double my = 0.0;
    mx = Eigen::VectorXd::Zero(p);
    double sy = 1.0;
    sx = Eigen::VectorXd::Ones(p);

    // Standardize y
    Eigen::VectorXd y_std = (y.array() - my) / sy;

    // Initialize xstd as a copy of x
    Eigen::MatrixXd x_std = x.cast<double>();

    // Create a vector of indices where !ct is true
    Eigen::VectorXi non_constant_cols = (ct.array() == false).cast<int>();

    // Perform the normalization for non-constant columns
    for (int i = 0; i < p; ++i) {
        if (non_constant_cols(i)) {
            x_std.col(i) = (x_std.col(i).array() - mx(i)) / sx(i);
        }
    }

    int known_phi = 1;

    // Create the stdconstants matrix
    Eigen::MatrixXd std_constants(p + 1, 2);
    std_constants(0, 0) = my;
    std_constants(0, 1) = sy;
    std_constants.block(1, 0, p, 1) = mx;
    std_constants.block(1, 1, p, 1) = sx;

    // Initialize delta_ini
//    int delta_ini = 0;



    // Calculate sums and products
    double sum_y2 = y_std.array().square().sum();
    double sum_y = y_std.sum();
    Eigen::VectorXd yt_x = y_std.transpose() * x_std;
    Eigen::VectorXd col_sums_x = x_std.colwise().sum();
    Eigen::MatrixXd xt_x = x_std.transpose() * x_std;
    bool has_xt_x = true;

    // Get family info
    // TODO we should be able to get rid of the whole family thing by setting fixed values and removing the corresponding
    //  if statements later on
    std::pair<int, int> f_family = get_family_info(family, false);
    int family_int = f_family.first;
    int family_greedy = f_family.second;

    // Initialize priors
    msPriorSpec prior_var = igprior(0.01, 0.01);
    msPriorSpec prior_skew = momprior(0.348);
    msPriorSpec prior_group = prior_coef;
    msPriorSpec prior_constraints = prior_delta;

    // Format priors
    FormatMsPriorsMargResult tmp_pm = format_ms_priors_marg(prior_coef, prior_group, prior_var, prior_skew, n);

    int r = tmp_pm.r;
    int prior = tmp_pm.prior;
    int prior_gr = tmp_pm.priorgr;
    double tau = tmp_pm.tau;
    double tau_group = tmp_pm.taugroup;
    double alpha = tmp_pm.alpha;
    double lambda = tmp_pm.lambda;
    double tau_alpha = tmp_pm.taualpha;
    double fix_atanh_alpha = tmp_pm.fixatanhalpha;
    prior_coef = tmp_pm.priorCoef;
    prior_group = tmp_pm.priorGroup;

    int pr_delta = 2;
    double pr_delta_p = 0.5;
    std::vector<double> par_pr_delta_p = {prior_delta.priorPars.at("alpha.p"), prior_delta.priorPars.at("beta.p")};

    int pr_constr = 2;
    double pr_constr_p = 0.5;
    std::vector<double> par_pr_constr_p = {prior_constraints.priorPars.at("alpha.p"),
                                           prior_constraints.priorPars.at("beta.p")};

    // Convert include_vars to int
    Eigen::VectorXi include_vars_int = include_vars.cast<int>();

//    int deltaini;
    int ndeltaini = w_i.sum();

    Eigen::VectorXi deltaini;  // This will be our output

    // Perform the element-wise OR operation
    Eigen::Array<bool, Eigen::Dynamic, 1> mask = (w_i.array() != 0 || include_vars.array() != 0);

// Find the indices where mask is true and subtract 1
    std::vector<int> indices;
    for (Eigen::Index i = 0; i < mask.size(); ++i) {
        if (mask(i)) {
            indices.push_back(i);  // Subtracting 1 to match R's 1-based indexing
        }
    }

// Convert the vector of indices to Eigen::VectorXi
    deltaini = Eigen::Map<Eigen::VectorXi>(indices.data(), indices.size());

    // Initialize post_mode and post_mode_prob
    Eigen::VectorXi post_mode = Eigen::VectorXi::Zero(p);
    //    double post_mode_prob = 1.0;

    // Get method parameters
    int method_int = 0; // Laplace

    // Initialize th_init
    Eigen::VectorXd th_init = Eigen::VectorXd::Zero(p);

    // Prepare for model selection
    int mcmc_2_save = std::floor((n_iter - burnin + 0.0) / (thinning + 0.0));
    // Removed the need for family as we will always have family 1
    int my_cols = p;
    int my_cols2 = p;

    // Allocate memory for results
    Eigen::VectorXi post_sample = Eigen::VectorXi::Zero(mcmc_2_save * my_cols);
    Eigen::VectorXd marg_pp = Eigen::VectorXd::Zero(my_cols2);
    Eigen::VectorXi post_mode_vec = Eigen::VectorXi::Zero(my_cols);
    Eigen::VectorXd post_prob = Eigen::VectorXd::Zero(mcmc_2_save);
    Eigen::VectorXi is_group = Eigen::VectorXi::Zero(p);
    Eigen::VectorXi n_constraints = Eigen::VectorXi::Zero(n_groups);
    Eigen::VectorXi n_inv_constraints = Eigen::VectorXi::Zero(n_groups);


    // Count constraints
    int n_groups_constr = 0;
    intptrvec constraints_ptr, inv_constraints_ptr;
    count_constraints(n_constraints, &constraints_ptr, n_inv_constraints, &inv_constraints_ptr,
                      &n_groups_constr, is_group, n_groups, n_var_in_group,
                      constraints, inv_constraints);


    // Prepare marginal_pars structure
    marginalPars pars;
    pars.n = &n;
    pars.p = &p;
    pars.y = y_std.data();
    pars.sumy2 = &sum_y2;
    pars.XtX = new crossprodmat(xt_x.data(), n, p, has_xt_x);
    pars.ytX = yt_x.data();
    pars.method = &method_int;
    pars.B = &b;
    pars.phi = &phi;
    pars.tau = &tau;
    pars.r = &r;
    pars.parprDeltap = par_pr_delta_p.data();
    pars.parprConstrp = par_pr_constr_p.data();
    int log_scale = 1;
    pars.logscale = &log_scale;
    pars.groups = groups.data();
    pars.ngroups = &n_groups;
    pars.ngroupsconstr = &n_groups_constr;
    pars.nvaringroup = n_var_in_group.data();
    pars.nconstraints = n_constraints.data();


    int verbose = 0;

    // Start of modelSelectionGibbsMinimal logic
    bool copylast = false;
    bool validmodel = false;

    int jgroup;
    int niter10 = 1;
    int niterthin = (int) floor((n_iter - burnin + .0) / (thinning + .0));
    int nbvars = p;
    int *firstingroup = ivector(0, n_groups);
    int *addgroups = ivector(0, 1);
    int *dropgroups = ivector(0, 1);
    int naddgroups;
    int ndropgroups;
    int *sample = ivector(0, n_groups - 1);

    int nsel = ndeltaini;
    int nselnew;
    int *sel = ivector(0, nbvars);
    int *selnew = ivector(0, nbvars);
    int *selaux;

    double *newJ = dvector(0, 2);
    double *ppnew = dvector(0, 3);
    double ppnewsum;
    double u;

    for (int j = 0; j < ndeltaini; j++) {
        sel[j] = deltaini[j];
        post_mode[deltaini[j]] = 1;
//        postMode[deltaini[j]] = 1;
    }

    modselIntegrals *integrals = new modselIntegrals(pimomMarginalKC, betabinPrior, p);

    firstingroup[0] = 0;
    for (int j = 1; j < n_groups; j++) {
        firstingroup[j] = firstingroup[j - 1] + n_var_in_group[j - 1];
    }

    for (int j = 0; j < n_groups; j++) {
        sample[j] = 1;
    }

    for (int j = 0; j < p; j++) {
        marg_pp[j] = 0;
    }


    double currentJ = integrals->getJoint(sel, &nsel, &pars);
//    post_prob[0] = post_prob = currentJ;

    int ilow = -burnin;
    int savecnt = 0;
    int iupper = n_iter - burnin + 1;


    // TODO n_constraints always 0? n_var_in_grop always ones? perhaps not needed?

    double post_mode_prob = -std::numeric_limits<double>::infinity();

    // Iterate
    for (int i = ilow, itercount = 0; i < iupper; i++, itercount++) {
        int j = jgroup = 0;
        while (j < p) {
            if (sample[jgroup] > 0) { // if jgroup should be sampled
                sel2selnew(jgroup, sel, &nsel, selnew, &nselnew, copylast, &n_groups, n_var_in_group.data(),
                           firstingroup); // copy sel into selnew, adding/removing jth group
                if (nsel > nselnew) {
                    naddgroups = 0;
                    ndropgroups = 1;
                    dropgroups[0] = jgroup;
                } else {
                    naddgroups = 1;
                    ndropgroups = 0;
                    addgroups[0] = jgroup;
                }
                validmodel = checkConstraints(addgroups, &naddgroups, dropgroups, &ndropgroups, &constraints_ptr,
                                              n_constraints.data(), &inv_constraints_ptr, n_inv_constraints.data(),
                                              groups.data(),
                                              n_var_in_group.data(), sel, &nsel);
                if (nselnew > n) validmodel = false;

                ppnew[0] = ppnewsum = 1;
                if (include_vars_int[j] == 0 && validmodel) { // if proposed model is valid
//                    if (
                    newJ[0] = integrals->getJoint(selnew, &nselnew, &pars);

                    if (newJ[0] > post_mode_prob) {
                        post_mode_prob = newJ[0];
                    }

                    ppnew[1] = exp(newJ[0] - currentJ);
                    ppnewsum += ppnew[1];
                }

                if (include_vars_int[j] == 0 && validmodel) { // if proposed model is valid
                    ppnew[1] /= ppnewsum;
                    if (i >= 0) {
                        if (nselnew > nsel) {
                            marg_pp[j] += ppnew[1];
                        } else {
                            marg_pp[j] += (1 - ppnew[1]);
                        }
                    } // update Rao-Blackwellized inclusion probabilities
                    u = runif();
                    if (u < ppnew[1]) {
                        selaux = sel;
                        sel = selnew;
                        selnew = selaux;
                        nsel = nselnew;
                        currentJ = newJ[0];
                    } // update model
                }
            } // end if jgroup should be sampled

            j += n_var_in_group[jgroup];
            jgroup++;

        }  // end j for

        if ((i > 0) && ((i % thinning) == 0)) {
            for (j = 0; j < nsel; j++) {
                post_sample[sel[j] * niterthin + savecnt] = 1;
            }
            post_prob[savecnt] = currentJ;
            savecnt++;
        }
    }
    if (iupper > ilow) {
        for (int j = 0; j < p; j++) {
            marg_pp[j] /= (iupper - std::max(0, ilow) + .0);
        }
    } // from sum to average

    // Clean up
    delete pars.XtX;
    delete integrals;
    free_ivector(addgroups, 0, 1);
    free_ivector(dropgroups, 0, 1);
    free_dvector(newJ, 0, 2);
    free_dvector(ppnew, 0, 3);
    free_ivector(firstingroup, 0, n_groups);
    free_ivector(sel, 0, nbvars);
    free_ivector(selnew, 0, nbvars);
    free_ivector(sample, 0, n_groups - 1);

    // Prepare the result
    Priors priors = {
            prior_coef,
            prior_group,
            prior_delta,
            prior_constraints,
            prior_var,
            prior_skew
    };

    return post_sample;


}
