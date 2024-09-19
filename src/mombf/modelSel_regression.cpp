// we only include RcppArmadillo.h which pulls Rcpp.h in for us
// #include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]


//Include other headers
#include "cstat.h"
#include "crossprodmat.h"
#include "modelSel_regression.h"
#include "modselIntegrals.h"
#include "Polynomial.h"

#include <csignal>
#include <map>
#include <string>

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
            ans = k + IS_imom(thopt, Voptinv, sel, nsel, (*pars).n, (*pars).p, (*pars).XtX, (*pars).ytX, (*pars).phi,
                              (*pars).tau, (*pars).B);
        }
    }
    if ((*(*pars).logscale) != 1) { ans = exp(ans); }
    free_dvector(thopt, 1, *nsel);
    free_dmatrix(Voptinv, 1, *nsel, 1, *nsel);
    return (ans);
}


//Evaluation of iMOM integral via Importance Sampling (result is returned in log-scale)
double IS_imom(double *thopt, double **Voptinv, int *sel, int *nsel, int *n, int *p, crossprodmat *XtX, double *ytX,
               double *phi, double *tau, int *B) {
    bool posdef;
    int i, j;
    double *sdprop, **Vprop, *sopt, **cholVprop, **cholVpropinv, detVpropinv, *mprop, *thsim, *logr, maxlogr, ans;

    sdprop = dvector(1, *nsel);
    sopt = dvector(1, *nsel);
    mprop = dvector(1, *nsel);
    thsim = dvector(1, *nsel);
    logr = dvector(0, 999);
    Vprop = dmatrix(1, *nsel, 1, *nsel);
    cholVprop = dmatrix(1, *nsel, 1, *nsel);
    cholVpropinv = dmatrix(1, *nsel, 1, *nsel);

    for (i = 1; i <= (*nsel); i++) {
        mprop[i] = 0;
        sopt[i] = sqrt(Voptinv[i][i]);
        sdprop[i] = .5 * fabs(thopt[i] + 2 * dsign(thopt[i]) * sopt[i]);
    }
    for (i = 1; i <= (*nsel); i++) {
        for (j = i; j <= (*nsel); j++) {
            Vprop[i][j] = Vprop[j][i] = sdprop[i] * sdprop[j] * Voptinv[i][j] / (sopt[i] * sopt[j]);
        }
    }
    choldc(Vprop, *nsel, cholVprop, &posdef);
    choldc_inv(Vprop, *nsel, cholVpropinv, &posdef);
    detVpropinv = choldc_det(cholVpropinv, *nsel);
    rmvtC(thsim, *nsel, mprop, cholVprop, 1);
    maxlogr = logr[0] = -fimomNegC_non0(thsim + 1, XtX, ytX, phi, tau, n, p, sel, nsel) -
                        dmvtC(thsim, *nsel, mprop, cholVpropinv, detVpropinv, 1, 1);
    for (i = 1; i < 1000; i++) {
        rmvtC(thsim, *nsel, mprop, cholVprop, 1);
        logr[i] = -fimomNegC_non0(thsim + 1, XtX, ytX, phi, tau, n, p, sel, nsel) -
                  dmvtC(thsim, *nsel, mprop, cholVpropinv, detVpropinv, 1, 1);
        if (logr[i] > maxlogr) { maxlogr = logr[i]; }
    }
    for (i = 0, ans = 0; i < 1000; i++) { ans += exp(logr[i] - maxlogr + 500); }
    for (i = 1000; i < (*B); i++) {
        rmvtC(thsim, *nsel, mprop, cholVprop, 1);
        ans += exp(-fimomNegC_non0(thsim + 1, XtX, ytX, phi, tau, n, p, sel, nsel) -
                   dmvtC(thsim, *nsel, mprop, cholVpropinv, detVpropinv, 1, 1) - maxlogr + 500);
    }
    ans = log(ans / (.0 + (*B))) + maxlogr - 500;

    free_dvector(sdprop, 1, *nsel);
    free_dvector(sopt, 1, *nsel);
    free_dvector(mprop, 1, *nsel);
    free_dvector(thsim, 1, *nsel);
    free_dvector(logr, 0, 999);
    free_dmatrix(Vprop, 1, *nsel, 1, *nsel);
    free_dmatrix(cholVprop, 1, *nsel, 1, *nsel);
    free_dmatrix(cholVpropinv, 1, *nsel, 1, *nsel);
    return (ans);
}