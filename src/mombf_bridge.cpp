//
// Created by jakob on 4/13/25.
//
#include "mombf_bridge.h"

namespace bisam {
    /* -------------------------------------------------------------------------- */
    /*                               modelSelection                               */
    /* -------------------------------------------------------------------------- */
    // This is a (very much simplified) version of the R function. Many parameters are hard coded since their values
    // are expected to remain the way they are for now in order for us not having to translate every single function
    // used within this one.

    arma::Col<int> modelSelection(const arma::vec &y,
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
                                  double prDeltap,
                                  arma::vec thinit,
                                  InitType initpar_type,
                                  int method,
                                  int hesstype,
                                  int optimMethod,
                                  int optim_maxit,
                                  int B,
                                  int knownphi,
                                  int r,
                                  double alpha,
                                  double lambda) {
        int p = x.n_cols;
        int n = y.n_elem;

        /* -------------------------------- postMode -------------------------------- */
        // Since familyint == 1 for now -> postMode is a vector of length p filled with zeros
        arma::vec postMode = arma::vec().zeros(p);

        /* ------------------------------ postModeProb ------------------------------ */
        double postModeProb = 0.0; // Initialized to 0.0

        /* -------------------------------- knownphi -------------------------------- */
        // We give phi as parameter -> phi is known
        // TODO potentially change to bool
        // Replaced with smiple Parameter
        // int knownphi = 1;

        /* ------------------------------- familyint -------------------------------- */
        // Set to 1 since we have family == "normal" and !issurvival (in formatFamily())
        int familyint = 1;

        /* --------------------------------- prior ---------------------------------- */
        // Set to 1. Defined in formatmsPriorsMarg. Equal to prior = piMOM
        int prior = 1;

        /* -------------------------------- priorgr --------------------------------- */
        // Set to 1. Defined in formatmsPriorsMarg. Equal to prior = piMOM
        int priorgr = 1;

        /* -------------------------- ndeltaini & deltaini -------------------------- */
        // Translation of the respective R code:
        // ndeltaini <- as.integer(sum(deltaini | includevars))
        // deltaini <- as.integer(which(deltaini | includevars) - 1)
        arma::Col<int> includevars(p); // = arma::vec().zeros(p);
        includevars.fill(arma::fill::zeros);

        int ndeltaini           = arma::sum(deltaini_input || includevars);
        arma::uvec indices      = arma::find(deltaini_input || includevars);
        arma::Col<int> deltaini = arma::conv_to<arma::Col<int> >::from(indices);

        /* ------------------- ystd, sumy2, sumy, xstd, colsumsx -------------------- */
        // The following calculations depend on one another and an actual translation is needed.
        // They have all been put in this one section since they are related.
        // TODO This part may include that is not of interest if scaling is never of interest?

        // Pre define variables so they are available after the if blocks aswell
        double my;
        double sy;

        // mx <- colMeans(x)
        arma::rowvec mx = arma::mean(x, 0);

        // sx <- sqrt(colMeans(x^2) - mx^2) * sqrt(n / (n - 1))
        // TODO Theoretically the case where n = 1 would need to be caught and handled to avoid div by zero
        arma::rowvec sx = (arma::sqrt(arma::mean(arma::pow(x, 2), 0) - arma::pow(mx, 2)))
                          * std::sqrt(static_cast<double>(n) / (n - 1));

        // ct <- (sx == 0)
        // TODO very much not sure if this is correct?
        arma::uvec ct = arma::conv_to<arma::uvec>::from(sx == 0);

        // if (!center) {
        //      my <- 0
        //      mx <- rep(0, p)
        // } else {
        //      my <- mean(y)
        // }
        if (!center) {
            my = 0.0;
            mx.fill(arma::fill::zeros);
        } else {
            my = arma::mean(y);
        }

        // if (!scale) {
        //      sy <- 1
        //      sx <- rep(1,p)
        // } else {
        //      sy <- sd(y)
        // }
        if (!scale) {
            sy = 1.0;
            sx.fill(arma::fill::ones);
        } else {
            sy = arma::stddev(y);
        }

        // Since typeofvar is always "numeric" in our case so far, this code block is skipped.
        // mx[typeofvar == "factor"] <- 0
        // sx[typeofvar == "factor"] <- 1

        // Since for now always false, this section of code is skipped
        // if (!(outcometype %in% c("Continuous", "Survival"))) {
        //   my <- 0
        //   sy <- 1
        // }

        // ystd <- (y - my) / sy
        arma::vec ystd = (y - my) / sy;

        // xstd <- x
        arma::mat xstd = x;

        // xstd[, !ct] <- t((t(x[, !ct]) - mx[!ct]) / sx[!ct])
        // In R, this standardizes columns with non-zero standard deviation
        // In Armadillo, we'll loop through columns and only standardize those with ct(j) == 0
        for (size_t j = 0; j < p; ++j) {
            if (ct(j) == 0) {
                // If sx(j) is not 0 (equivalent to !ct in R)
                xstd.col(j) = (x.col(j) - mx(j)) / sx(j);
            }
        }


        // sumy2 <- as.double(sum(ystd^2))
        // Using dot product which is more efficient for sum of squares
        double sumy2 = arma::dot(ystd, ystd);

        // sumy <- as.double(sum(ystd))
        double sumy = arma::sum(ystd);

        // ytX <- as.vector(matrix(ystd, nrow = 1) %*% xstd)
        // In Armadillo, we can multiply a column vector's transpose with a matrix to get row vector
        arma::rowvec ytX = ystd.t() * xstd;

        // colsumsx <- as.double(colSums(xstd))
        // 0 indicates summing along rows to get column sums
        arma::vec colsumsx = arma::sum(xstd, 0).t();

        // if (XtXprecomp) { XtX <- t(xstd) %*% xstd; hasXtX <- as.logical(TRUE) }
        // else { XtX <- double(0); hasXtX <- as.logical(FALSE) }
        arma::mat XtX;
        bool hasXtX;
        if (XtXprecomp) {
            // Compute the cross-product matrix
            XtX    = xstd.t() * xstd;
            hasXtX = true;
        } else {
            // Create an empty matrix instead of a single double
            XtX.reset();
            hasXtX = false;
        }

        /* ------------------------------ sumlogyfact ------------------------------- */
        // Hard coded since family is not 22 sumlogyfact <- as.double(0)
        double sumlogyfact = 0;

        /* --------------------------------- uncens --------------------------------- */
        // Hard coded to 0 since "formula" %in% class(y) = False and "Surv" %in% class(y) in
        // formatInputData
        int uncens = 0;

        /* --------------------------------- method --------------------------------- */
        // Fixed to int 0 in formats method. Because opti method = auto and outcome type = glm
        // Replaced with simple input parameter
        // int method = 0;

        /* ------------------------------ adj_overdisp ------------------------------ */
        // Fixed to 1 because of formats method. adj.overdisp <- as.integer(ifelse(adj.overdisp == "none", 0,
        // ifelse(adj.overdisp == "intercept", 1, 2))). We set adj.overdisp='intercept'.
        int adj_overdisp = 1;

        /* -------------------------------- hesstype -------------------------------- */
        // Fixed to 1 because of formats method. hesstype <- as.integer(ifelse(hess == "asympDiagAdj", 2, 1))
        // Replace with simple input parameter
        // int hesstype = 1;

        /* ------------------------------ optimMethod ------------------------------- */
        // Fixed to 2 because of formats method. We set optimMethod="auto".
        // Replaced with simple Paramter
        // int optimMethod = 2;

        /* ------------------------------ optim_maxit ------------------------------- */
        // Default is set to 10? Not really a default but we set it in our list of "default" parameters?
        // Review mit Lucas: setzen wir mal auf 0
        // Replaced with simple parameter
        // int optim_maxit = 0;

        /* -------------------------- thinit and usethinit -------------------------- */
        // Thinit uses a lasso for initialization - put into its own function for clearity
        int usethinit;

        std::pair<arma::vec, int> res = getthinit(y, x, familyint, thinit, false, initpar_type);
        thinit                        = res.first;
        usethinit                     = res.second;


        /* ----------------------------------- B ------------------------------------ */
        // Default parameter set to 100000
        // Replace with simple Parameter
        // int B = 100000;

        /* --------------------------------- alpha ---------------------------------- */
        // Fixed to 0.01 - parameter set by us. alpha <- as.double(priorVar@priorPars["alpha"])
        // Replace with simple paramater
        // double alpha = 0.01;

        /* --------------------------------- lambda --------------------------------- */
        // Fixed to 0.01 - parameter set by us. lambda <- as.double(priorVar@priorPars["lambda"])
        // Replaced with simple parameter
        // double lambda = 0.01;

        /* ---------------------------------- phi ----------------------------------- */
        // Set to be s2_i by us. Parameter of this function no further calculations.

        /* ---------------------------------- tau ----------------------------------- */
        // priorCoef = imomprior(tau = tau); originally. We just add a function paramter to get tau

        /* -------------------------------- taugroup -------------------------------- */
        // We have missing(priorgroup) == True ---> priorGroup <- priorCoef
        // Thus taugroup == tau? Hastaugroupstd == false ---> taugroup <- as.double(priorGroup@priorPars["tau"])
        double taugroup = tau;

        /* -------------------------------- taualpha -------------------------------- */
        // Since "msPriorSpec" %in% class(priorSkew) == true ---> taualpha <- as.double(priorSkew@priorPars["tau"])
        double taualpha = priorSkew;

        /* ----------------------------- fixatanhalpha ------------------------------ */
        // Since "msPriorSpec" %in% class(priorSkew) == true ---> fixatanhalpha <- as.double(-10000)
        int fixatanhalpha = -10000;

        /* ----------------------------------- r ------------------------------------ */
        // Fixed to 1. In formatmsPriorsMarg for piMom
        // Replaced with simple parameter
        // int r = 1;

        /* -------------------------------- prDelta --------------------------------- */
        int prDelta = 1;

        // TODO INVESTIGATION NEEDED
        // prDelta = 2;

        /* -------------------------------- prDeltap -------------------------------- */
        // Parameter priorDelta@priorPars[["p"]]

        /* ------------------------------ parprDeltap ------------------------------- */
        // parprDeltap <- as.double(length(prDeltap))
        int parprDeltap = 1;

        /* -------------------------------- prConstr -------------------------------- */
        // Fixed to one in formatmsPriorsMarg
        int prConstr = 1;

        // TODO INVESTIGATION NEEDED
        // prConstr = 2;

        /* ------------------------------- prConstrp -------------------------------- */
        // priorConstraints < -defaultpriorConstraints(priorDelta)
        // prConstrp <- as.double(priorConstraints@priorPars[["p"]])
        // Fixed to 0.5 in our case
        double prConstrp = 0.5;

        /* ------------------------------ parprConstrp ------------------------------ */
        // parprConstrp <- as.double(length(prConstrp))
        int parprConstrp = 1;

        /* --------------------------------- groups --------------------------------- */
        // Vector with numbers 0 to 20. Groups is modiefied (compared to default parameter) in
        // codeGroupsAndConstraints
        arma::Col<int> groups = arma::regspace<arma::Col<int> >(0, p - 1);

        /* -------------------------------- ngroups --------------------------------- */
        // = p
        int ngroups = p;

        /* ------------------------------ nvaringroup ------------------------------- */
        // nvaringroup <- as.integer(rep(1, p))
        // arma::uvec nvaringroup = arma::regspace<arma::uvec>(1, p);
        arma::Col<int> nvaringroup = arma::Col<int>().ones(p);

        /* ------------------------------ constraints ------------------------------- */
        // I am not 100% sure about constraints and invconstraints. Further checking may be needed
        arma::ivec constraints = arma::ivec().zeros(p);

        /* ----------------------------- invconstraints ----------------------------- */
        // I am not 100% sure about constraints and invconstraints. Further checking may be needed
        arma::ivec invconstraints = arma::ivec().zeros(p);

        /* -------------------------------- verbose --------------------------------- */
        // Has to be set to 0 otherwise some mombf code (calling Rprintf?) segfaults
        int verbose = 0;

        arma::Col<int> output = modelSelectionGibbsCI(postMode,
                                                      postModeProb,
                                                      knownphi,
                                                      familyint,
                                                      prior,
                                                      priorgr,
                                                      niter,
                                                      thinning,
                                                      burnin,
                                                      ndeltaini,
                                                      deltaini,
                                                      includevars,
                                                      n,
                                                      p,
                                                      ystd,
                                                      uncens,
                                                      sumy2,
                                                      sumy,
                                                      sumlogyfact,
                                                      xstd,
                                                      colsumsx,
                                                      hasXtX,
                                                      XtX,
                                                      ytX,
                                                      method,
                                                      adj_overdisp,
                                                      hesstype,
                                                      optimMethod,
                                                      optim_maxit,
                                                      thinit,
                                                      usethinit,
                                                      B,
                                                      alpha,
                                                      lambda,
                                                      phi,
                                                      tau,
                                                      taugroup,
                                                      taualpha,
                                                      fixatanhalpha,
                                                      r,
                                                      prDelta,
                                                      prDeltap,
                                                      parprDeltap,
                                                      prConstr,
                                                      prConstrp,
                                                      parprConstrp,
                                                      groups.memptr(),
                                                      ngroups,
                                                      nvaringroup,
                                                      constraints,
                                                      invconstraints,
                                                      verbose);


        return output;
    }

    /* -------------------------------------------------------------------------- */
    /*                           ModelSelectionGibbsCI                            */
    /* -------------------------------------------------------------------------- */
    // The following code is a modified version of the original modelSelectionGibbsCI. The original function
    // prepares the R (SEXP) objects for use with cpp. Since we already have the cpp objects ready not much
    // is needed here. We mainly make sure all our parameters are in a form to call modelSelectionGibbs.

    // This is an "output object" for modelSelectionGibbsCI mimicing what the original version outputs,
    // but since we are now only returning postSample this is kind of redundant now. I have to keep it for now doe.
    struct GibbsOutput {
        int *postSample;
        double *margpp;
        int *postMode;
        double *postModeProb;
        double *postProb;
        int mcmc2save;
        int mycols;
        int mycols2;

        GibbsOutput() : postSample(nullptr), margpp(nullptr), postMode(nullptr), postModeProb(nullptr),
                        postProb(nullptr), mcmc2save(0), mycols(0), mycols2(0) {
        }

        ~GibbsOutput() {
            delete[] postSample;
            delete[] margpp;
            delete[] postMode;
            delete[] postModeProb;
            delete[] postProb;
        }
    };


    arma::Col<int> modelSelectionGibbsCI(const arma::vec &SpostModeini,
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
                                         double SparprDeltap,
                                         int SpriorConstr,
                                         double SprConstrp,
                                         double SparprConstrp,
                                         int *Sgroups,
                                         int Sngroups,
                                         arma::Col<int> &Snvaringroup,
                                         arma::ivec &Sconstraints,
                                         arma::ivec &Sinvconstraints,
                                         int Sverbose) {
        GibbsOutput output;

        int i, j, idxj, logscale = 1, mcmc2save, *postSample, *postMode, mycols, mycols2, *nconstraints, *
                ninvconstraints, nuncens, ngroupsconstr = 0, *isgroup, usethinit = Susethinit, priorcode;

        double offset = 0, *margpp, *postModeProb, *postProb, *ytXuncens = NULL, *thinit;
        intptrvec constraints, invconstraints;
        crossprodmat *XtX, *XtXuncens = NULL;
        struct marginalPars pars;

        output.mcmc2save = floor((Sniter - Sburnin + .0) / (Sthinning + .0));
        if (Sfamily != 0) {
            output.mycols = output.mycols2 = Sp;
        } else {
            output.mycols  = 2 + Sp;
            output.mycols2 = output.mycols + 2;
        }
        mcmc2save = output.mcmc2save;
        mycols    = output.mycols;
        mycols2   = output.mycols2;

        thinit = dvector(0, mycols2 + 1);
        if (usethinit != 3) {
            for (j = 0; j <= mycols2 + 1; j++) { thinit[j] = 0; }
        } else {
            for (j = 0; j <= Sp; j++) { thinit[j] = Sthinit[j]; }
        }

        output.postSample = new int[mcmc2save * mycols];
        postSample        = output.postSample;
        for (j = 0; j < (mcmc2save * mycols); j++) postSample[j] = 0;

        output.margpp = new double[mycols2];
        margpp        = output.margpp;

        output.postMode = new int[mycols];
        postMode        = output.postMode;
        for (j = 0; j < mycols; j++) { postMode[j] = SpostModeini[j]; }

        output.postModeProb = new double[1];
        postModeProb        = output.postModeProb;
        postModeProb[0]     = SpostModeiniProb;

        output.postProb = new double[mcmc2save];
        postProb        = output.postProb;

        isgroup         = ivector(0, Sp);
        nconstraints    = ivector(0, Sngroups);
        ninvconstraints = ivector(0, Sngroups);

        // This is our modified version of countConstraints which takes cpp objects instead of SEXP
        countConstraints(nconstraints,
                         &constraints,
                         ninvconstraints,
                         &invconstraints,
                         &ngroupsconstr,
                         isgroup,
                         &Sngroups,
                         (int *) Snvaringroup.memptr(),
                         Sconstraints,
                         Sinvconstraints);

        if (ShasXtX) {
            XtX = new crossprodmat(SXtX.memptr(), Sn, Sp, true);
        } else {
            XtX = new crossprodmat(Sx.memptr(), Sn, Sp, false);
        }

        if (Suncens > 0) {
            //if there's censoring, also store t(x) %*% x and t(x) %*% y computed over uncensored observations
            int n       = Sn, *uncens       = &Suncens; // JAKOB added & in front of Suncens since its already int.
            double *pty = Sy.memptr(), *ptx = Sx.memptr();
            for (nuncens = 0; (nuncens < n) && (uncens[nuncens] == 1); nuncens++) {
            }
            XtXuncens = new crossprodmat(Sx.memptr(), Sn, Sp, false, nuncens, 0);
            ytXuncens = dvector(0, Sp);
            for (j = 0; j < Sp; j++) {
                for (i = 0, ytXuncens[j] = 0, idxj = j * n; i < nuncens; i++) {
                    ytXuncens[j] += pty[i] * ptx[i + idxj];
                }
            }
        } else { nuncens = Sn; }

        set_marginalPars(&pars,
                         &Sfamily,
                         &Sn,
                         &nuncens,
                         &Sp,
                         Sy.memptr(),
                         &Suncens,
                         &Ssumy2,
                         &Ssumy,
                         &Ssumlogyfact,
                         Sx.memptr(),
                         Scolsumsx.memptr(),
                         XtX,
                         SytX.memptr(),
                         &Smethod,
                         &Sadjoverdisp,
                         &Shesstype,
                         &SoptimMethod,
                         &Soptim_maxit,
                         &usethinit,
                         thinit,
                         &SB,
                         &Salpha,
                         &Slambda,
                         &Sknownphi,
                         &Sphi,
                         &Stau,
                         &Staugroup,
                         &Staualpha,
                         &Sfixatanhalpha,
                         &Sr,
                         &SprDeltap,
                         &SparprDeltap,
                         &SprConstrp,
                         &SparprConstrp,
                         &logscale,
                         &offset,
                         Sgroups,
                         isgroup,
                         &Sngroups,
                         &ngroupsconstr,
                         Snvaringroup.memptr(),
                         nconstraints,
                         ninvconstraints,
                         XtXuncens,
                         ytXuncens);

        priorcode      = mspriorCode(&SpriorCoef, &SpriorGroup, &pars);
        pars.priorcode = &priorcode;


        modelSelectionGibbs(postSample,
                            margpp,
                            postMode,
                            postModeProb,
                            postProb,
                            &SpriorDelta,
                            &SpriorConstr,
                            &Sniter,
                            &Sthinning,
                            &Sburnin,
                            &Sndeltaini,
                            Sdeltaini.memptr(),
                            (int *) Sincludevars.memptr(),
                            &constraints,
                            &invconstraints,
                            &Sverbose,
                            &pars);

        free_dvector(thinit, 0, mycols2 + 1);
        free_ivector(isgroup, 0, Sp);
        free_ivector(nconstraints, 0, Sngroups);
        free_ivector(ninvconstraints, 0, Sngroups);
        delete XtX;
        delete XtXuncens; // Added deallocation for XtXuncens


        // TODO this extra copy here is really not all that nice!
        arma::Col<int> post_sample = arma::Col<int>(mcmc2save * mycols);
        for (j = 0; j < (mcmc2save * mycols); j++) post_sample[j] = postSample[j];
        return post_sample;
    }

    /* -------------------------------------------------------------------------- */
    /*                                 getthinit                                  */
    /* -------------------------------------------------------------------------- */
    // The following code provides 2 functions and an enum to emulate the getthinit and initParameters
    // R functions.


    // The original function used the initpar as a string (always a string and never vec here!). Instead
    // we use our enum. For now we assume family to always be 1 (=normal!). Thus some code is not translated
    arma::vec initParameters(const arma::vec &y, const arma::mat &x, int family,
                             InitType initpar) {
        int n = y.n_elem;
        int p = x.n_cols;

        if (initpar == InitType::AUTO) {
            if (p <= n / 2) {
                initpar = InitType::MLE;
            } else {
                initpar = InitType::L1;
            }
        }

        // I think these lines should never be needed. Setting the family here should not do anything since we
        // are hard coding it later. And initpar should never be none here.
        // if (!(family % in % c('binomial', 'poisson'))) family = 'gaussian'
        // if (initpar == 'none') {
        //     ans = rep(0, p)
        // }

        // TODO if needed the lasso code can probably easily be adapted to accept arma types directly
        // Convert Armadillo data to the format expected by LassoRegression
        std::vector<std::vector<double> > samples(n, std::vector<double>(p));
        std::vector<double> target(n);

        for (int i = 0; i < n; ++i) {
            target[i] = y(i);
            for (int j = 0; j < p; ++j) {
                samples[i][j] = x(i, j);
            }
        }

        // Furthermore I think that we ALWAYS hit the L1 case? So we will simply use a lasso here not matter what
        // Create LassoRegression object
        LassoRegression lasso(samples, target);

        // Set tolerance and alpha for coordinate descent
        double tolerance = 0.001;
        double alpha     = 1; // This corresponds to the regularization strength


        // Run cyclical coordinate descent to get the weights
        double *weights = lasso.cyclicalCoordinateDescent(tolerance, alpha);

        // Convert the result (double*) to arma::vec
        arma::vec arma_weights(weights, p, true, true);

        // Clean up the dynamically allocated weights array
        delete[] weights;

        return arma_weights;
    }


    // Translation of the R function found in initParameters.R
    // Since mixing types in C++ is not really a thing and the original function uses the initpar parameter
    // as either a string or a vector, instead here initpar is only a vec. If it is NOT empty, the the value
    // is used directly. If it is empty, the type in initpar_type shall be used to determin what value
    // usethinit should be set to
    std::pair<arma::vec, int> getthinit(const arma::vec &y, const arma::mat &x, int family,
                                        const arma::vec &initpar, bool enumerate, InitType initpar_type) {
        arma::vec thinit;
        int usethinit;

        if (!initpar.empty()) {
            thinit    = initpar;
            usethinit = 1;
            if (thinit.n_elem != x.n_cols) {
                throw std::runtime_error(
                    "x has " + std::to_string(x.n_cols) + " columns but initpar has " + std::to_string(thinit.n_elem) +
                    " elements");
            }
        } else {
            if (initpar_type == InitType::NONE) {
                if (enumerate) {
                    usethinit = 0;
                } else {
                    usethinit = 1;
                }
            } else {
                usethinit = 3;
            }
            thinit = initParameters(y, x, family, initpar_type);


            // Passing an empty string as initpar to match the original logic for missing/auto
        }

        return {thinit, usethinit};
    }


    /* ---------------------------- countConstraints ---------------------------- */
    // Since count constraints in mombf accepts FOR WHATEVER REASON SEXP objects and we are not using them,
    // we need a version without them. Other than that this should be an almost 1 to 1 translation of the
    // original function

    void countConstraints(int *nconstraints,
                          intptrvec *constraints,
                          int *ninvconstraints,
                          intptrvec *invconstraints,
                          int *ngroupsconstr,
                          int *isgroup,
                          int *ngroups,
                          int *nvaringroup,
                          arma::ivec &Sconstraints,
                          arma::ivec &Sinvconstraints) {
        /*Count number of constraints, number of groups with constraints, determine which variables are in a group*/
        int i, j, jj;
        int offset_constraints    = 0;
        int offset_invconstraints = 0;
        *ngroupsconstr            = 0;

        constraints->clear();
        invconstraints->clear();

        // Create temporary storage for constraint data that intptrvec can point to
        int **temp_constraints    = new int *[*ngroups];
        int **temp_invconstraints = new int *[*ngroups];

        for (j = 0, jj = 0; j < *ngroups; j++) {
            // Get number of constraints for this group
            nconstraints[j] = Sconstraints(j);

            // Add pointer to constraints for this group
            if (nconstraints[j] > 0) {
                // Allocate memory for constraints
                temp_constraints[j] = new int[nconstraints[j]];

                // Copy constraint values
                for (i = 0; i < nconstraints[j]; i++) {
                    temp_constraints[j][i] = Sconstraints(*ngroups + offset_constraints + i);
                }

                constraints->push_back(temp_constraints[j]);
                offset_constraints += nconstraints[j];
                (*ngroupsconstr)++;
            } else {
                temp_constraints[j] = nullptr;
                constraints->push_back(nullptr);
            }

            // Get number of inverse constraints for this group
            ninvconstraints[j] = Sinvconstraints(j);

            // Add pointer to inverse constraints for this group
            if (ninvconstraints[j] > 0) {
                // Allocate memory for inverse constraints
                temp_invconstraints[j] = new int[ninvconstraints[j]];

                // Copy inverse constraint values
                for (i = 0; i < ninvconstraints[j]; i++) {
                    temp_invconstraints[j][i] = Sinvconstraints(*ngroups + offset_invconstraints + i);
                }

                invconstraints->push_back(temp_invconstraints[j]);
                offset_invconstraints += ninvconstraints[j];
            } else {
                temp_invconstraints[j] = nullptr;
                invconstraints->push_back(nullptr);
            }

            // Mark which variables belong to a group
            isgroup[jj] = ((int) (nvaringroup[j] + 0.1)) > 1;
            jj++;
            for (i = 1; i < nvaringroup[j]; i++, jj++) {
                isgroup[jj] = isgroup[jj - 1];
            }
        }

        // Note: We're not freeing temp_constraints and temp_invconstraints arrays here
        // as the pointers are now stored in constraints and invconstraints vectors.
        // They should be freed when no longer needed elsewhere in the code.
        // delete[] temp_constraints;
        // delete[] temp_invconstraints;
    }


    void rnlpPost_lm(double *ans, int niter, int burnin, int thinning, double *y, double *x, int n, int p, int r,
                     double tau, double a_phi, double b_phi, int prior, arma::dvec &thinit, bool use_thinit,
                     double phiinit, bool use_phiinit) {
        bool posdef;
        int i, j, k, isave, nsave;
        double *m, *mortho, *alpha, **S, **Sinv, **cholSinv, **inv_cholSinv, **K, **D, tauinv = 1.0 / tau, *Xty, *thcur,
                phicur, phinew, sqrtphi, th2sum, th2invsum, apost, bpost, *linpred, ssr;
        //Pre-compute stuff
        nsave        = (int) floor((niter - burnin + .0) / (thinning + .0));
        m            = dvector(1, p);
        mortho       = dvector(1, p);
        alpha        = dvector(1, p);
        thcur        = dvector(1, p);
        linpred      = dvector(0, n - 1);
        S            = dmatrix(1, p, 1, p);
        Sinv         = dmatrix(1, p, 1, p);
        cholSinv     = dmatrix(1, p, 1, p);
        inv_cholSinv = dmatrix(1, p, 1, p);
        K            = dmatrix(1, p, 1, p);
        D            = dmatrix(1, p, 1, p);

        AvectBvec(x, n, p, x, n, p, S); //S= t(x) %*% x + 1/tau
        for (i = 1; i <= p; i++) S[i][i] += tauinv;
        inv_posdef(S, p, Sinv, &posdef);
        choldc(Sinv, p, cholSinv, &posdef);
        choldc_inv(Sinv, p, inv_cholSinv, &posdef); //inverse of chol(Sinv)

        Xty = dvector(1, p);
        Atvecx(x, y, Xty + 1, 0, p - 1, 0, n - 1); //m= solve(S) %*% t(x) %*% y
        Ax(Sinv, Xty, m, 1, p, 1, p);
        Ax(inv_cholSinv, m, mortho, 1, p, 1, p);
        free_dvector(Xty, 1, p);

        if (prior == 0) apost = .5 * (a_phi + n + 3 * p);
        else if (prior == 1) apost = .5 * (a_phi + n - p);
        else apost                 = .5 * (a_phi + n + p);
        //Initialize
        th2sum = 0;
        phicur = sqrtphi = 1.0;

        if (!use_thinit) {
            for (j = 1; j <= p; j++) {
                thcur[j] = m[j];
                th2sum += thcur[j] * thcur[j];
            }
        } else {
            // TODO can we get rid of loop?
            // -1 because thinit is arma vec with 0 based index
            for (j = 1; j <= p; j++) {
                thcur[j] = thinit[j - 1];
                th2sum += thcur[j] * thcur[j];
            }
        }

        // std::cout << "thinit:\n" << thinit << std::endl;

        //Ax(cholSinv, thcur, Dthcur, 1, p, 1, p);
        isave = 0;
        for (i = 1; i <= niter; i++) {
            //for (j=1; j<=p; j++) Dthcur[j] = Dthcur[j] / sqrtphi;
            Avecx(x, thcur + 1, linpred, 0, n - 1, 0, p - 1);
            ssr = 0;
            for (j = 0; j < n; j++) ssr += pow(y[j] - linpred[j], 2.0);
            if (!use_phiinit) {
                if (prior == 0) {
                    bpost   = .5 * (b_phi + th2sum / tau + ssr);
                    phicur  = 1.0 / rgammaC(apost, bpost);
                    sqrtphi = sqrt(phicur);
                } else {
                    if (prior == 1) bpost = .5 * (b_phi + ssr);
                    else bpost            = .5 * (b_phi + th2sum / tau + ssr);
                    phinew    = 1.0 / rgammaC(apost, bpost);
                    th2invsum = 0;
                    for (j = 1; j <= p; j++) th2invsum += 1 / (thcur[j] * thcur[j]);
                    if (runif() < exp((phicur - phinew) * tau * th2invsum)) {
                        phicur  = phinew;
                        sqrtphi = sqrt(phicur);
                    }
                }
            } else {
                phicur  = phiinit;
                sqrtphi = sqrt(phicur);
            }
            for (j = 1; j <= p; j++) {
                alpha[j] = mortho[j] / sqrtphi;
                //Dthcur[j]= Dthcur[j] * sqrtphi;
                for (k = 1; k <= j; k++) {
                    D[j][k] = cholSinv[j][k] * sqrtphi;
                    K[j][k] = inv_cholSinv[j][k] / sqrtphi;
                }
            }
            rnlp_Gibbs(thcur, p, alpha, D, K, &tau, &phicur, r, prior);
            if (i > burnin && ((i - burnin) % thinning) == 0) {
                for (j = 1; j <= p; j++) ans[isave + (j - 1) * nsave] = thcur[j];
                ans[isave + p * nsave] = phicur;
                isave++;
            }
        }
        free_dvector(m, 1, p);
        free_dvector(mortho, 1, p);
        free_dvector(alpha, 1, p);
        free_dvector(thcur, 1, p);
        free_dvector(linpred, 0, n - 1);
        free_dmatrix(S, 1, p, 1, p);
        free_dmatrix(Sinv, 1, p, 1, p);
        free_dmatrix(cholSinv, 1, p, 1, p);
        free_dmatrix(inv_cholSinv, 1, p, 1, p);
        free_dmatrix(K, 1, p, 1, p);
        free_dmatrix(D, 1, p, 1, p);
    }
}
