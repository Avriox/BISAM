// we only include RcppArmadillo.h which pulls Rcpp.h in for us
// #include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]


//Include other headers

#include "modelSel_regression.h"


//Global variables defined for minimization/integration routines
struct marginalPars f2opt_pars;





//********************************************************************************************
// GENERAL ALGEBRA
//********************************************************************************************

//multiply symmetric A[1..fi][1..fi] * x[sel[0]..sel[fi-1]]
//Note: A is indexed at 1. x and sel are indexed at 0. ans is indexed at 1.
void Asym_xsel(double **A, int fi, double *x, int *sel, double *ans) {
    int _i, _j;
    for (_i = 1; _i <= fi; _i++) {
        for (_j = _i, ans[_i] = 0; _j <= fi; _j++) { ans[_i] += A[_i][_j] * x[sel[_j - 1]]; }
        for (_j = 1; _j < _i; _j++) { ans[_i] += A[_j][_i] * x[sel[_j - 1]]; }
    }
}


//Add constant ct to diagonal elements in XtX[sel,sel]. XtX[0..p-1][0..p-1] is formatted as a vector indexed at 0, V[sel[0]..sel[nsel-1]][sel[0]..sel[nsel-1]] as a matrix indexed at 1, sel is indexed at 0
//Note: Only diagonal & upper-diagonal elements in V are set.
void addct2XtX(double *ct, crossprodmat *XtX, int *sel, int *nsel, int *p, double **V) {
    int i, j;
    for (i = 1; i <= (*nsel); i++) { V[i][i] = XtX->at(sel[i - 1] * (*p) + sel[i - 1]) + (*ct); }
    for (i = 1; i <= (*nsel); i++) {
        for (j = i + 1; j <= (*nsel); j++) {
            V[i][j] = XtX->at(sel[j - 1] * (*p) + sel[i - 1]);
        }
    }
}


void set_f2opt_pars(double *m, double **S, double *sumy2, crossprodmat *XtX, double *ytX, double *alpha, double *lambda,
                    double *phi, double *tau, int *r, int *n, int *p, int *sel, int *nsel) {
//    f2opt_pars.m = m;
//    f2opt_pars.S = S;
    f2opt_pars.sumy2 = sumy2;
    f2opt_pars.XtX = XtX;
    f2opt_pars.ytX = ytX;
//    f2opt_pars.alpha = alpha;
//    f2opt_pars.lambda = lambda;
    f2opt_pars.phi = phi;
    f2opt_pars.tau = tau;
    f2opt_pars.r = r;
//    f2opt_pars.n = n;
    f2opt_pars.p = p;
    f2opt_pars.sel = sel;
    f2opt_pars.nsel = nsel;
}


bool checkConstraints(int *addgroups, int *naddgroups, int *dropgroups, int *ndropgroups, intptrvec *constraints,
                      int *nconstraints, intptrvec *invconstraints, int *ninvconstraints, int *groups, int *nvaringroup,
                      int *sel, int *nsel) {
    /* Check if adding variables in addgroups and dropping variables from dropgroups to current model (sel) gives a valid model satisfying all constraints
   Input

   - addgroups: id of groups to be added to sel (vector of length naddgroups)
   - naddgroups:  number of groups to be added
   - dropgroups: id of groups to be dropped from sel (vector of length ndropgroups)
   - ndropgroups: number of groups to be dropped
   - constraints: constraints[j] is a vector indicating all groups required by group j. Assumed to be ordered increasingly
   - nconstraints: nconstraints[j] is the length of constraints[j] (number of constraints required by group j)
   - invconstraints: invconstraints[j] is a vector indicating all inverse constraints of group j (ids of other groups requiring group j)
   - ninvconstraints: ninvconstraints[j] is the length of invconstraints[j] (number of inverse constraints of group j)
   - firstingroup: index of the first variable in all groups
   - sel: index of variables currently in the model. Assumed to be ordered increasingly
   - nsel: number of variables currently in the model
  */

    bool valid = true;
    int j, k, l, *curconstraints, curgroup, nvalid;

    //For any group we want to add to sel, check that its constraints are also in sel
    for (k = 0; (k < *naddgroups) && valid; k++) {
        nvalid = j = l = 0;
        curgroup = addgroups[k];
        curconstraints = (*constraints)[curgroup];
        while ((l < nconstraints[curgroup]) && (j < *nsel) && valid) {
            if (groups[sel[j]] > curconstraints[l]) {
                valid = false;
            } else if (groups[sel[j]] == curconstraints[l]) {  //add groups and nvaringroup as parameters
                l++;
                nvalid++;
            } else {
                j += nvaringroup[sel[j]];
            }
        }
        if (nvalid < nconstraints[curgroup]) { valid = false; }
    }

    //For any group we want to drop from sel, check that its inverse constraints are not in sel
    for (k = 0; (k < *ndropgroups) && valid; k++) {
        j = l = 0;
        curgroup = dropgroups[k];
        curconstraints = (*invconstraints)[curgroup];
        while ((l < ninvconstraints[curgroup]) && (j < *nsel) && valid) {
            if (groups[sel[j]] > curconstraints[l]) {
                l++;
            } else if (groups[sel[j]] == curconstraints[l]) {
                valid = false;
            } else {
                j += nvaringroup[groups[sel[j]]];
            }
        }
    }

    return valid;
}


void
nselConstraints(int *ngroups0, int *ngroups1, int *sel, int *nsel, int *group, int *nconstraints, int *nvaringroup) {
    /* Return number of selected groups that have hierarchical constraints

     Input
     - sel: indexes of selected variables
     - nsel: number of selected variables (length of sel)
     - group: group[sel[j]] indicates the group that variable j belongs to. Variables are assumed to be ordered by groups (group 1, group 2 etc)
     - nconstraints: nconstraints[l] is the number of constraints for group l (nconstraints[l]==0 indicates no constraints)

     Output:
     - ngroups0: number of groups in sel that do not have hierarchical constraints
     - ngroups1: number of groups in sel that have hierarchical constraints
  */
    int j = 0, g;
    (*ngroups0) = (*ngroups1) = 0;
    while (j < *nsel) {
        g = group[sel[j]];
        if (nconstraints[g] == 0) { (*ngroups0)++; } else { (*ngroups1)++; }
        j += nvaringroup[g];
    }
}


//nsel ~ Beta-Binomial(prModelpar[0],prModelPar[1])
double betabinPrior(int *sel, int *nsel, struct marginalPars *pars) {
    int ngroups0, ngroups1;
    double ans;
    nselConstraints(&ngroups0, &ngroups1, sel, nsel, (*pars).groups, (*pars).nconstraints, (*pars).nvaringroup);
    ans = bbPrior(ngroups0, *(*pars).ngroups - *(*pars).ngroupsconstr, (*pars).parprDeltap[0], (*pars).parprDeltap[1],
                  1);
    if ((*(*pars).ngroupsconstr) > 0)
        ans += bbPrior(ngroups1, *(*pars).ngroupsconstr, (*pars).parprConstrp[0], (*pars).parprConstrp[1], 1);
    return ans;
    //return bbPrior(*nsel, *(*pars).p, (*pars).parprDeltap[0], (*pars).parprDeltap[1],1);
}

double vectBinom(int *sel, int *nsel, int len_prDeltap, int len_prConstrp, struct marginalPars *pars) {
    int i, sel_i = 0, delta_i = 0, constr_i = 0, ngroups = *(*pars).ngroups, *groups = (*pars).groups, *nconstraints = (*pars).nconstraints, *nvaringroup = (*pars).nvaringroup;
    double ans = 0, *prDeltap = (*pars).prDeltap, *prConstrp = (*pars).prConstrp;

    if (*nsel == 0) {
        for (i = 0; i < len_prDeltap; i++) ans += log(1 - prDeltap[(len_prDeltap > 1) ? i : 0]);
        if (*(*pars).ngroupsconstr > 0) {
            for (i = 0; i < len_prConstrp; i++) ans += log(1 - prConstrp[(len_prConstrp > 1) ? i : 0]);
        }
    } else {
        for (i = 0; i < ngroups; i++) {
            if (nconstraints[i] == 0) {
                if (i == groups[sel[sel_i]]) {
                    ans += log(prDeltap[delta_i]);
                    if (sel_i < *nsel - 1) sel_i += nvaringroup[groups[i]];
                } else {
                    ans += log(1 - prDeltap[delta_i]);
                }
                if (len_prDeltap > 1) delta_i++;
            } else {  // constrained, use prConstrp
                if (i == groups[sel[sel_i]]) {
                    ans += log(prConstrp[constr_i]);
                    if (sel_i < *nsel - 1) sel_i += nvaringroup[groups[i]];
                } else {
                    ans += log(1 - prConstrp[constr_i]);
                }
                if (len_prConstrp > 1) constr_i++;
            }
        }
    }
    return ans;
}

double binomPrior(int *sel, int *nsel, struct marginalPars *pars) {
    int ngroups0, ngroups1, n_notconstr, n_constr = *(*pars).ngroupsconstr,
            len_prDeltap = (int) *(*pars).parprDeltap, len_prConstrp = (int) *(*pars).parprConstrp;
    double ans, *prDeltap = (*pars).prDeltap, *prConstrp = (*pars).prConstrp;
    nselConstraints(&ngroups0, &ngroups1, sel, nsel, (*pars).groups, (*pars).nconstraints, (*pars).nvaringroup);
    n_notconstr = *(*pars).ngroups - n_constr;
    if ((len_prDeltap == 1) & (len_prConstrp == 1)) {
        ans = (ngroups0 + .0) * log(*prDeltap) + (n_notconstr - ngroups0 + .0) * log(1 - *prDeltap);
        if ((n_constr) > 0) {
            ans += (ngroups1 + .0) * log(*prConstrp) + (n_constr - ngroups1 + .0) * log(1 - *prConstrp);
        }
    } else {
        ans = vectBinom(sel, nsel, len_prDeltap, len_prConstrp, pars);
    }
    return ans;
    // return dbinomial(*nsel,*(*pars).p,*(*pars).prDeltap,1);
}

double f2opt_imom(double *th) {
    double ans;
    ans = fimomNegC_non0(th + 1, f2opt_pars.XtX, f2opt_pars.ytX, f2opt_pars.phi, f2opt_pars.tau, f2opt_pars.n,
                         f2opt_pars.p, f2opt_pars.sel, f2opt_pars.nsel);
    return (ans);
}

double fimomNegC_non0(double *th, crossprodmat *XtX, double *ytX, double *phi, double *tau, int *n, int *p, int *sel,
                      int *nsel) {
//same as fimomNegC but loops over all elem in th (i.e. th has length *nsel and contains non-zero elements only). th is indexed at 0.
    int i;
    double ans, ytXth, sumlogth, suminvth, th2;
    for (i = 0, ytXth = 0, sumlogth = 0, suminvth = 0; i < (*nsel); i++) {
        ytXth += ytX[sel[i]] * th[i];
        th2 = th[i] * th[i];
        suminvth += 1 / th2;
        sumlogth += log(th2);
    }
    ans = .5 * (quadratic_xtAselx(th, XtX, p, nsel, sel) - 2 * ytXth) / (*phi) + (*tau) * (*phi) * suminvth + sumlogth;
    return ans;
}

double fimomNegC_non0_eigen(const Eigen::VectorXd &th, const Eigen::MatrixXd &XtX, const Eigen::VectorXd &ytX,
                            double phi, double tau, int n, int p, const Eigen::VectorXi &sel) {
    int nsel = sel.size();

    // Calculate ytXth using Eigen operations
    double ytXth = (ytX(sel).array() * th.array()).sum();

    // Calculate sumlogth and suminvth
    double sumlogth = (th.array().square().log()).sum();
    double suminvth = (th.array().square().inverse()).sum();

    // Calculate the quadratic form
    double quadratic = quadratic_xtAselx_eigen(th, XtX, sel);

    // Calculate the final result
    double ans = 0.5 * (quadratic - 2 * ytXth) / phi + tau * phi * suminvth + sumlogth;

    return ans;
}

//Hessian of fimomNegC
// - ans: hessian matrix evaluated at th (indexed at 1, i.e. ans[1:(*nsel)][1:(*nsel)])
// - th: th[1:(*nsel)] indicates point at which to evaluate the hessian.
// - Other arguments as in fimomNegC_non0
void
fppimomNegC_non0(double **ans, double *th, crossprodmat *XtX, double *ytX, double *phi, double *tau, int *n, int *p,
                 int *sel, int *nsel) {
    int i, j;
    double th2;

    for (i = 1; i <= (*nsel); i++) {
        th2 = th[i] * th[i];
        ans[i][i] =
                (XtX->at(sel[i - 1] * (*p) + sel[i - 1])) / (*phi) + 6.0 * (*tau) * (*phi) / (th2 * th2) - 2.0 / th2;
    }
    for (i = 1; i <= (*nsel); i++) {
        for (j = i + 1; j <= (*nsel); j++) {
            ans[i][j] = ans[j][i] = (XtX->at(sel[i - 1] * (*p) + sel[j - 1])) / (*phi);
        }
    }
}


void imomModeK(double *th, PolynomialRootFinder::RootStatus_T *status, crossprodmat *XtX, double *ytX, double *phi,
               double *tau, int *sel, int *nsel, int *p) {
    //piMOM mode when phi is known using gradient algorithm
    // - th: contains initial estimate at input and mode at output
    // - status: indicates if root finding has been successful
    bool found = false;
    int i, j, niter = 0, root_count;
    double err = 1.0, *coef, *real_vector, *imag_vector;
    Polynomial poly;

    coef = dvector(0, 4);
    real_vector = dvector(0, 4);
    imag_vector = dvector(0, 4);

    coef[0] = 2.0 * (*tau) * (*phi);
    coef[1] = 0.0;
    coef[2] = -2;
    while ((err > 1.0e-5) & (niter < 50)) {
        err = 0;
        for (i = 1; i <= (*nsel); i++) {
            coef[3] = ytX[sel[i - 1]];
            for (j = 1; j < i; j++) { coef[3] -= (XtX->at(sel[i - 1] * (*p) + sel[j - 1])) * th[j]; }
            for (j = i + 1; j <= (*nsel); j++) { coef[3] -= (XtX->at(sel[i - 1] * (*p) + sel[j - 1])) * th[j]; }
            coef[3] = coef[3] / (*phi);
            coef[4] = -(XtX->at(sel[i - 1] * (*p) + sel[i - 1])) / (*phi);
            poly.SetCoefficients(coef, 4);
            (*status) = poly.FindRoots(real_vector, imag_vector, &root_count);

            j = 0;
            found = false;
            while ((!found) & (j <= 4)) {
                if (fabs(imag_vector[j]) < 1.0e-5) {
                    if (((real_vector[j] > 0) & (th[i] > 0)) | ((real_vector[j] < 0) & (th[i] < 0))) {
                        err = max_xy(err, fabs(th[i] - real_vector[j]));
                        th[i] = real_vector[j];
                        found = true;
                    }
                }
                j++;
            }

        }
        niter++;
    }

    free_dvector(coef, 0, 4);
    free_dvector(real_vector, 0, 4);
    free_dvector(imag_vector, 0, 4);
}


void imomIntegralApproxC(double *ILaplace, double *thopt, double **Voptinv, double *fopt, int *sel, int *nsel, int *n,
                         int *p, crossprodmat *XtX, double *ytX, double *phi, double *tau, int *logscale,
                         int *hessian) {
    bool posdef;
    int iter, maxit = 100, emptyint;
    double **V, **Vinv, ftol = 1.0e-5, **dirth, **Vopt, detVopt, emptydouble = 0, **emptymatrix;
    PolynomialRootFinder::RootStatus_T status;

    V = dmatrix(1, *nsel, 1, *nsel);
    Vinv = dmatrix(1, *nsel, 1, *nsel);
    Vopt = dmatrix(1, *nsel, 1, *nsel);
    dirth = dmatrix(1, *nsel, 1, *nsel);
    emptymatrix = dmatrix(1, 1, 1, 1);
    //Initialize
    addct2XtX(tau, XtX, sel, nsel, p, V); //add tau to XtX diagonal, store in V
    inv_posdef_upper(V, *nsel, Vinv, &posdef);
    Asym_xsel(Vinv, *nsel, ytX, sel, thopt);  //product Vinv * selected elements in ytX
    //Minimization
    imomModeK(thopt, &status, XtX, ytX, phi, tau, sel, nsel, p);
    set_f2opt_pars(&emptydouble, emptymatrix, &emptydouble, XtX, ytX, &emptydouble, &emptydouble, phi, tau, &emptyint,
                   n, p, sel, nsel);
    if (status == PolynomialRootFinder::SUCCESS) {
        (*fopt) = f2opt_imom(thopt);
    } else {
        // TODO Never hit in testing / coverage - but should potentially be included again!
//        ddiag(dirth, 1, *nsel);
//        minimize(thopt ddia, dirth, *nsel, ftol, &iter, fopt, f2opt_imom, maxit);
    }

    if (*hessian == 1) {
        //Laplace approx
        fppimomNegC_non0(Vopt, thopt, XtX, ytX, phi, tau, n, p, sel, nsel);
        invdet_posdef(Vopt, *nsel, Voptinv, &detVopt);
        (*ILaplace) = -(*fopt) - 0.5 * log(detVopt);
    } else {
        (*ILaplace) = -(*fopt) - 0.5 * (*nsel) * log(*n + .0);  //BIC-type approximation
    }

    free_dmatrix(V, 1, *nsel, 1, *nsel);
    free_dmatrix(Vinv, 1, *nsel, 1, *nsel);
    free_dmatrix(Vopt, 1, *nsel, 1, *nsel);
    free_dmatrix(dirth, 1, *nsel, 1, *nsel);
    free_dmatrix(emptymatrix, 1, 1, 1, 1);
    if ((*logscale) != 1) { (*ILaplace) = exp(*ILaplace); }
}


double pimomMarginalKC(int *sel, int *nsel, struct marginalPars *pars) {
    int one = 1, hessian;
    double k, ans, m, s, ILaplace, *thopt, **Voptinv, fopt;
    thopt = dvector(1, *nsel);
    Voptinv = dmatrix(1, *nsel, 1, *nsel);
    if ((*nsel) == 0) {
        m = 0;
        s = sqrt(*(*pars).phi);
        ans = dnormC_jvec((*pars).y, *(*pars).n, m, s, 1);
    } else {
        if (*(*pars).method == 2) { hessian = 0; } else { hessian = 1; }
        imomIntegralApproxC(&ILaplace, thopt, Voptinv, &fopt, sel, nsel, (*pars).n, (*pars).p, (*pars).XtX, (*pars).ytX,
                            (*pars).phi, (*pars).tau, &one, &hessian);
        k = .5 * ((*nsel) * log(*(*pars).tau) - (*(*pars).sumy2) / (*(*pars).phi) - (*(*pars).n + .0) * LOG_M_2PI -
                  (*(*pars).n - *nsel) * log(*(*pars).phi) - (*nsel) * LOG_M_PI);
        if (((*(*pars).method) == 0) || ((*(*pars).method) == 2)) {
            ans = k + ILaplace;
        } else {
//            ans = k + IS_imom(thopt, Voptinv, sel, nsel, (*pars).n, (*pars).p, (*pars).XtX, (*pars).ytX, (*pars).phi,
//                              (*pars).tau, (*pars).B);


        }
    }
    if ((*(*pars).logscale) != 1) { ans = exp(ans); }
    free_dvector(thopt, 1, *nsel);
    free_dmatrix(Voptinv, 1, *nsel, 1, *nsel);
    return (ans);
}


//Evaluation of iMOM integral via Importance Sampling (result is returned in log-scale)
//double IS_imom(double *thopt, double **Voptinv, int *sel, int *nsel, int *n, int *p, crossprodmat *XtX, double *ytX,
//               double *phi, double *tau, int *B) {
//    bool posdef;
//    int i, j;
//    double *sdprop, **Vprop, *sopt, **cholVprop, **cholVpropinv, detVpropinv, *mprop, *thsim, *logr, maxlogr, ans;
//
//    sdprop = dvector(1, *nsel);
//    sopt = dvector(1, *nsel);
//    mprop = dvector(1, *nsel);
//    thsim = dvector(1, *nsel);
//    logr = dvector(0, 999);
//    Vprop = dmatrix(1, *nsel, 1, *nsel);
//    cholVprop = dmatrix(1, *nsel, 1, *nsel);
//    cholVpropinv = dmatrix(1, *nsel, 1, *nsel);
//
//    for (i = 1; i <= (*nsel); i++) {
//        mprop[i] = 0;
//        sopt[i] = sqrt(Voptinv[i][i]);
//        sdprop[i] = .5 * fabs(thopt[i] + 2 * dsign(thopt[i]) * sopt[i]);
//    }
//    for (i = 1; i <= (*nsel); i++) {
//        for (j = i; j <= (*nsel); j++) {
//            Vprop[i][j] = Vprop[j][i] = sdprop[i] * sdprop[j] * Voptinv[i][j] / (sopt[i] * sopt[j]);
//        }
//    }
//    choldc(Vprop, *nsel, cholVprop, &posdef);
//    choldc_inv(Vprop, *nsel, cholVpropinv, &posdef);
//    detVpropinv = choldc_det(cholVpropinv, *nsel);
//    rmvtC(thsim, *nsel, mprop, cholVprop, 1);
//    maxlogr = logr[0] = -fimomNegC_non0(thsim + 1, XtX, ytX, phi, tau, n, p, sel, nsel) -
//                        dmvtC(thsim, *nsel, mprop, cholVpropinv, detVpropinv, 1, 1);
//    for (i = 1; i < 1000; i++) {
//        rmvtC(thsim, *nsel, mprop, cholVprop, 1);
//        logr[i] = -fimomNegC_non0(thsim + 1, XtX, ytX, phi, tau, n, p, sel, nsel) -
//                  dmvtC(thsim, *nsel, mprop, cholVpropinv, detVpropinv, 1, 1);
//        if (logr[i] > maxlogr) { maxlogr = logr[i]; }
//    }
//    for (i = 0, ans = 0; i < 1000; i++) { ans += exp(logr[i] - maxlogr + 500); }
//    for (i = 1000; i < (*B); i++) {
//        rmvtC(thsim, *nsel, mprop, cholVprop, 1);
//        ans += exp(-fimomNegC_non0(thsim + 1, XtX, ytX, phi, tau, n, p, sel, nsel) -
//                   dmvtC(thsim, *nsel, mprop, cholVpropinv, detVpropinv, 1, 1) - maxlogr + 500);
//    }
//    ans = log(ans / (.0 + (*B))) + maxlogr - 500;
//
//    free_dvector(sdprop, 1, *nsel);
//    free_dvector(sopt, 1, *nsel);
//    free_dvector(mprop, 1, *nsel);
//    free_dvector(thsim, 1, *nsel);
//    free_dvector(logr, 0, 999);
//    free_dmatrix(Vprop, 1, *nsel, 1, *nsel);
//    free_dmatrix(cholVprop, 1, *nsel, 1, *nsel);
//    free_dmatrix(cholVpropinv, 1, *nsel, 1, *nsel);
//    return (ans);
//}


int mspriorCode(int *prCoef, int *prGroup, struct marginalPars *pars) {
    // Returns a two-digit code indicating the prior on regression coefficients. The 1st digit is the prior on individual coef; The 2nd digit the prior on groups of coefficients
    //  Input
    //  - prCoef: 0 for pMOM; 1 for piMOM; 2 for peMOM; 3 for Zellner; 4 for normalid; 10 for group pMOM; 13 for group Zellner
    //  - prGroup: 0 for pMOM; 1 for piMOM; 2 for peMOM; 3 for Zellner; 4 for normalid; 10 for group pMOM; 11 for group iMOM; 12 for group eMOM; 13 for group Zellner
    //  Output
    //    0: pMOM on all coef
    //    1: peMOM on all coef
    //    2: piMOM on all coef
    //    3: Zellner on all coef
    //    4: normalid on all coef
    //    5: group pMOM (same as pMOM, standardized by n / X'X)
    //    9: group Zellner on all coef
    //   10: pMOM + group MOM
    //   13: pMOM + group Zellner
    //   32: peMOM + group eMOM
    //   33: peMOM + group Zellner
    //   43: Zellner + group Zellner
    //   50: group MOM + group MOM
    //   53: group MOM + group Zellner
    //   63: group Zellner + group Zellner
    //   73: normalid + group Zellner
    //  100: BIC (no prior, tells marginal likelihood routines to return -0.5 BIC, the BIC approx to the marginal likelihood)
    bool hasgroups = (*((*pars).ngroups)) < (*((*pars).p));
    int ans;
    if (*prCoef == 100) {
        ans = 100;  // BIC
    } else {
        if (!hasgroups) {
            if (*prCoef == 0) {  // pMOM on all coef
                ans = 0;
            } else if (*prCoef == 1) {  // piMOM on all coef
                ans = 1;
            } else if (*prCoef == 2) {  // peMOM on all coef
                ans = 2;
            } else if (*prCoef == 3) {  // Zellner on all coef
                ans = 3;
            } else if (*prCoef == 4) {  // normalid on all coef
                ans = 4;
            } else if (*prCoef == 10) {
                ans = 5;                 // group pMOM
            } else if (*prCoef == 13) {  // block Zellner on all coef
                ans = 9;
            } else {
                Rf_error("Prior specified by priorCoef not currently implemented\n");
            }
        } else {
            if ((*prCoef == 0) & (*prGroup == 0)) {  // pMOM on all coef
                ans = 0;
            } else if ((*prCoef == 1) & (*prGroup == 1)) {  // piMOM on all coef
                ans = 1;
            } else if ((*prCoef == 2) & (*prGroup == 2)) {  // peMOM on all coef
                ans = 2;
            } else if ((*prCoef == 3) & (*prGroup == 3)) {  // Zellner on all coef
                ans = 3;
            } else if ((*prCoef == 4) & (*prGroup == 4)) {  // normalid on all coef
                ans = 4;
            } else if ((*prCoef == 0) & (*prGroup == 10)) {  // pMOM + group MOM
                ans = 10;
            } else if ((*prCoef == 0) & (*prGroup == 13)) {  // pMOM + group Zellner
                ans = 13;
            } else if ((*prCoef == 2) & (*prGroup == 12)) {  // peMOM + group eMOM
                ans = 32;
            } else if ((*prCoef == 2) & (*prGroup == 13)) {  // peMOM + group Zellner
                ans = 33;
            } else if ((*prCoef == 3) & (*prGroup == 13)) {  // Zellner + group Zellner
                ans = 43;
            } else if ((*prCoef == 10) & (*prGroup == 10)) {  // group pMOM + group pMOM
                ans = 50;
            } else if ((*prCoef == 10) & (*prGroup == 13)) {  // group pMOM + group Zellner
                ans = 53;
            } else if ((*prCoef == 13) & (*prGroup == 13)) {  // group Zellner + group Zellner
                ans = 63;
            } else if ((*prCoef == 4) & (*prGroup == 13)) {  // normalid + group Zellner
                ans = 73;
            } else {
                Rf_error("Prior specified by priorCoef and priorGroup not currently implemented\n");
            }
        }
    }
    return ans;
}