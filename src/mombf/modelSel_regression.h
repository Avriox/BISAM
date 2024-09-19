#ifndef MODELSEL_H
#define MODELSEL_H 1

//#include <R.h>
//#include <Rinternals.h>
#include <list>
#include <math.h>
#include <vector>
#include <stdlib.h>
#include "./crossprodmat.h"
#include "./Polynomial.h"
#include "./cstat.h"


/*****************************************************************************************************************
  typedefs
*****************************************************************************************************************/

typedef double(*pt2margFun)(int *, int *,
                            struct marginalPars *);  //pointer to function to compute marginal densities & prior prob (used for model selection)

typedef std::vector<int *> intptrvec; //vector where each element is a pointer to an integer





/*****************************************************************************************************************
  Define structures
*****************************************************************************************************************/


// Structure containing parameters needed to calculate the integrated likelihood
struct marginalPars {
//    int *family;
//    int *priorcode;
    int *sel;
    int *nsel;
    int *n;        //number of observations
//    int *nuncens;  //number of uncensored observations
    int *p;
    double *y;
//    int *uncens;
    double *sumy2;
//    double *sumy;
//    double *sumlogyfact; //sum(log(y!)), used in Poisson regression
//    double *x;
//    double *colsumsx;   //column sums of x
    crossprodmat *XtX;  //t(x) %*% x using all observations
//    crossprodmat *XtXuncens; //t(x) %*% x using uncensored observations
//    covariancemat *V0inv;  // covariance matrix for coef and groups priors
    double *ytX;             //t(x) %*% y using all observations
//    double *ytXuncens;       //t(x) %*% y using uncensored observations
//    double *m;  //Sinv * Xty   (needed by mom and emom)
//    double **S;  //XtX + I/tau  (needed by mom and emom)
    int *method; //method==0 for Laplace; method==1 for Monte Carlo; method==2 for ALA (method== -1 for automatic choice)
//    int *adjoverdisp; //Used only for ALA. 0 for no adjustment; 1 to estimate overdispersion from intercept-only model, 2 to estimate from residuals
//    int *hesstype; //for asymmetric Laplace residuals hess=1 means using asymptotic hessian, hess=2 means using diagonal adjustment to asymp hessian
//    int *optimMethod; //optimization method to find mode
//    int *optim_maxit; //maximum number of iterations
//    int *usethinit; //usethinit==1 tells optimization algorithms to store the optimal model parameters at thinit; usethinit==2 to initialize at thinit upon entry and store optimal value at thinit upon exit; usethinit==3 to initialize at thinit but not to change thinit; usethinit==0 to ignore thinit
//    double *thinit; //thinit[sel[j]] stores initial values for model parameters to be used by optimization algorithms
    int *B;      //number of Monte Carlo samples
//    double *alpha;    //prior for residual variance is IG(.5*alpha,.5*lambda)
//    double *lambda;
//    int *knownphi; //should dispersion parameter be considered known, e.g. error var in Gaussian regression, or phi=1 in logistic/poisson regression
    double *phi;      //residual variance
    double *tau;      //dispersion parameter in prior for regression coefficients. Also used to store the penalty parameter when using an info criteria rather than a prior
//    double *taugroup; //prior dispersion parameter on grouped coefficients, e.g. the block Zellner prior is prod_j N(delta_j; 0, (taugroup/ncol(X_j)) (X_j'X_j)^{-1})
//    double *taualpha; //dispersion parameter in prior for asymmetry parameter in two-piece Normal or two-piece Laplace residuals
//    double *fixatanhalpha; //fixed value for asymmetry parameter (usedful for quantile regression at fixed quantile levels)
    int *r;           //MOM power parameter for prior on coefficients
//    double *prDeltap; //For Binomial prior on model space, prDeltap is the prob of success. For complexity prior, the power parameter in the exponential
    double *parprDeltap; //For Beta-Binomial prior on model space, parprDeltap[0],parprDeltap[1] are the prior parameters
//    double *prConstrp; //idem for prior on number of included groups under hierarchical constraints
    double *parprConstrp;
    int *logscale;
//    double *offset;
    int *groups;  //group that each variable belongs to
//    int *isgroup; //isgroup[j]==1 indicates that variable j is in a group
    int *ngroups; //total number of groups
    int *ngroupsconstr; //number of groups that have a hierarchical constraint
    int *nvaringroup; //number of coefficients in group[0],...,group[ngroups-1]
    int *nconstraints; //number of constraints in group[0],...,group[ngroups-1]
//    int *ninvconstraints; //number of inverse constraints (number of groups depending on group[0],...,group[ngroups-1])
};





//*************************************************************************************
//Setting prior & marginals
//*************************************************************************************

int mspriorCode(int *prCoef, int *prGroup, struct marginalPars *pars);



//*************************************************************************************
//General Algebra
//*************************************************************************************

void Asym_xsel(double **A, int fi, double *x, int *sel,
               double *ans);  //multiply symmetric A[1..fi][1..fi] * x[sel[0]..sel[fi-1]]; Return in ans[1..fi]

void addct2XtX(double *ct, crossprodmat *XtX, int *sel, int *nsel, int *p,
               double **V); //add constant to diagonal elem of XtX



void set_f2opt_pars(double *m, double **S, double *sumy2, crossprodmat *XtX, double *ytX, double *alpha, double *lambda,
                    double *phi, double *tau, int *r, int *n, int *p, int *sel, int *nsel);


bool checkConstraints(int *addgroups, int *naddgroups, int *dropgroups, int *ndropgroups, intptrvec *constraints,
                      int *nconstraints, intptrvec *invconstraints, int *ninvconstraints, int *groups, int *nvaringroup,
                      int *sel, int *nsel);


void
nselConstraints(int *ngroups0, int *ngroups1, int *sel, int *nsel, int *group, int *nconstraints, int *nvaringroup);


double betabinPrior(int *sel, int *nsel, struct marginalPars *pars);


// pMOM on all coef
// TODO remove?
double pmomMarginalKC(int *sel, int *nsel, struct marginalPars *pars);


// piMOM on all coef
double pimomMarginalKC(int *sel, int *nsel, struct marginalPars *pars);

// TODO remove?
double pmomMarginalKC(int *sel, int *nsel, struct marginalPars *pars);


double f2opt_imom(double *th);

double fimomNegC_non0(double *th, crossprodmat *XtX, double *ytX, double *phi, double *tau, int *n, int *p, int *sel,
                      int *nsel);


void
fppimomNegC_non0(double **ans, double *th, crossprodmat *XtX, double *ytX, double *phi, double *tau, int *n, int *p,
                 int *sel, int *nsel);

void imomModeK(double *th, PolynomialRootFinder::RootStatus_T *status, crossprodmat *XtX, double *ytX, double *phi,
               double *tau, int *sel, int *nsel, int *p);


double IS_imom(double *thopt, double **Voptinv, int *sel, int *nsel, int *n, int *p, crossprodmat *XtX, double *ytX,
               double *phi, double *tau, int *B);


#endif /* MODELSEL_H */

