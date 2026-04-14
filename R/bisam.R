#' Bayesian Indicator Saturation Model
#'
#' @param data A matrix with columns for unit indices, time indices, y variable, and X variables
#' @param i_index Column index for unit identifier (1-based R indexing)
#' @param t_index Column index for time identifier (1-based R indexing)
#' @param y_index Column index for dependent variable (1-based R indexing)
#' @param Ndraw Number of MCMC draws
#' @param Nburn Number of burn-in iterations
#' @param b_prior Prior type for coefficients ("g" or "f" or "hs")
#' @param lambda_b Parameter for g/f prior
#' @param c0 Parameter for inverse gamma prior on variance
#' @param C0 Parameter for inverse gamma prior on variance
#' @param va Parameter
#' @param vb Parameter
#' @param tau Parameter for model selection
#' @param priorDelta Prior object for model selection (modelbbprior or modelbinomprior)
#' @param geweke Boolean for Geweke test
#' @param use_phiinit Boolean for using initial phi values. TRUE -> Use our sampling, FALSE use rnlp sampling
#' @param const_val Boolean for including a constant
#' @param ife Boolean for individual fixed effects
#' @param tfe Boolean for time fixed effects
#' @param iis Boolean for indicator saturation
#' @param sis Boolean for stepshift saturation
#' @param computation_strategy Computation strategy (STANDARD, SPLIT_SEQUENTIAL, SPLIT_PARALLEL)
#' @param new_par_method New parameter method
#' @param new_par_hesstype New parameter hessian type
#' @param new_par_optim_method New parameter optimization method
#' @param new_par_optim_maxit New parameter optimization max iterations
#' @param new_par_B New parameter B
#' @param new_par_knownphi New parameter known phi
#' @param new_par_r New parameter r
#' @param new_par_alpha New parameter alpha
#' @param new_par_lambda New parameter lambda
#'
#' @return A list containing MCMC results and summary statistics
#' @export
estimate_model <- function(
    data,
    i_index = 1,
    t_index = 2,
    y_index = 3,
    Ndraw = 5000,
    Nburn = 500,
    b_prior = "g",
    lambda_b = 100.0,
    c0 = 0.001,
    C0 = 0.001,
    va = 1.0,
    vb = 1.0,
    tau = 1.0,
    priorDelta = modelbbprior(1, 1),
    geweke = FALSE,
    use_phiinit = TRUE,
    const_val = FALSE,
    ife = FALSE,
    tfe = FALSE,
    iis = TRUE,
    sis = TRUE,
#     new_par_include_vars = c(),
    new_par_method = 1,
    new_par_hesstype = 1,
    new_par_optim_method = 1,
    new_par_optim_maxit = 100,
    new_par_B = 10,
    new_par_knownphi = 0,
    new_par_r = 1,
    new_par_alpha = 0.05,
    new_par_lambda = 1.0,
    computation_strategy = SPLIT_SEQUENTIAL()  # Default to SPLIT_SEQUENTIAL
) {
    # Convert data to matrix if it's a data frame
    if (is.data.frame(data)) {
        data <- as.matrix(data)
    }

    # Adjust indices for R to C++ (0-based indexing)
    i_index_cpp <- i_index - 1
    t_index_cpp <- t_index - 1
    y_index_cpp <- y_index - 1

    a <- 1
    b <- 1
    p <- 0.5
    # prDelta says whether to use modelbinom or modelbb
    # modelbinom --> 1
    # modelbb --> 2
    #
    # 1 as default

    prDelta <- 1

    # Extract the prior variables depending on the type
    if (inherits(priorDelta, "modelbbprior")) {
          a <- priorDelta$a
          b <- priorDelta$b
          prDelta <- 2
      } else if (inherits(priorDelta, "modelbinomprior")) {
          p <- priorDelta$p
          prDelta <- 1
      } else {
        stop("priorDelta must be either a modelbbprior or modelbinomprior object")
      }




    # Call the C++ function
    result <- rcpp_estimate_model(
        data, i_index_cpp, t_index_cpp, y_index_cpp, Ndraw, Nburn, b_prior,
        lambda_b, c0, C0, va, vb, tau, use_phiinit,
        const_val, ife, tfe, iis, sis,
#          new_par_include_vars,
        new_par_method, new_par_hesstype, new_par_optim_method, new_par_optim_maxit,
        new_par_B, new_par_knownphi, new_par_r, new_par_alpha, new_par_lambda,
        prDelta, p, a, b,
        computation_strategy
    )

    return(result)
}

#' Constants for computation strategy
#'
#' @export
STANDARD <- function() {
  comp_strategy_standard()
}

#' @export
SPLIT_SEQUENTIAL <- function() {
  comp_strategy_split_sequential()
}

#' @export
SPLIT_PARALLEL <- function() {
  comp_strategy_split_parallel()
}



#' Create a Beta-Binomial prior object
#'
#' @param a First parameter of the Beta-Binomial prior (shape1)
#' @param b Second parameter of the Beta-Binomial prior (shape2)
#' @return A modelbbprior object
#' @export
modelbbprior <- function(a = 1, b = 1) {
  if (!is.numeric(a) || !is.numeric(b) || length(a) != 1 || length(b) != 1) {
    stop("Parameters 'a' and 'b' must be single numeric values")
  }
  if (a <= 0 || b <= 0) {
    stop("Parameters 'a' and 'b' must be positive")
  }

  structure(
    list(a = as.double(a), b = as.double(b)),
    class = "modelbbprior"
  )
}

#' Create a Binomial prior object
#'
#' @param p Probability parameter for the Binomial prior
#' @return A modelbinomprior object
#' @export
modelbinomprior <- function(p = 0.5) {
  if (!is.numeric(p) || length(p) != 1) {
    stop("Parameter 'p' must be a single numeric value")
  }
  if (p < 0 || p > 1) {
    stop("Parameter 'p' must be between 0 and 1")
  }

  structure(
    list(p = as.double(p)),
    class = "modelbinomprior"
  )
}

# Updated Wrapper

#' Fast Bayesian model selection (optimized backend)
#'
#' Exposes the package's optimized C++ model-selection step (a Gibbs-style update for
#' the model indicator \eqn{\delta}) so you can prepare inputs in R and avoid the extra
#' BISAM preparation done inside `estimate_model()`.
#'
#' This is intended as a near drop-in for the *matrix-interface* of `mombf::modelSelection()`
#' (i.e., `y` numeric vector and `x` numeric matrix), but note the return type differs:
#' this function returns a single sampled inclusion vector (0/1), not an `msfit` object.
#'
#' @param y Numeric vector of length n. Typically the current residual-like target, e.g.
#'   `y_hat = y - X %*% beta` in your BISAM Gibbs loop.
#' @param x Numeric matrix n x p. Candidate regressors to be selected (e.g., your Z matrix).
#'
#' @param center,scale Logical. Whether to center/scale internally in the model-selection code.
#' @param XtXprecomp Logical. Whether to precompute X'X internally.
#'
#' @param niter Integer. Number of internal Gibbs iterations.
#' @param thinning Integer. Thinning factor.
#' @param burnin Integer. Burn-in. Default `round(niter/10)`.
#'
#' @param deltaini Optional initial inclusion vector (logical/integer) length p.
#'
#' @param priorCoef Either NULL, a scalar `tau`, or an object/list with element `$tau`.
#'   Only `tau` is used by this wrapper.
#' @param tau Numeric. Prior scale parameter (used if `priorCoef` is NULL).
#'
#' @param priorDelta Model-space prior: `modelbbprior(a,b)` or `modelbinomprior(p)`.
#'
#' @param phi Numeric. Residual variance passed to the underlying routine.
#' @param knownphi Integer 0/1. Passed through.
#' @param priorSkew Numeric. Passed through to your translated mombf internals.
#'
#' @param thinit Optional numeric vector length p for initialization of continuous parameters
#'   inside the marginal-likelihood computations.
#' @param initpar_type Integer. Cast to your internal `InitType` enum (0 typically means AUTO).
#'
#' @param method,hesstype,optimMethod,optim_maxit,B,r,alpha,lambda Passed through to the C++
#'   internals (these are your integer-coded options, not mombf character strings).
#'
#' @param computation_strategy One of `STANDARD()`, `SPLIT_SEQUENTIAL()`, `SPLIT_PARALLEL()`.
#' @param n_units Integer. Number of partitions/units used by split strategies (often the number
#'   of individuals `n` in your panel setting).
#' @param max_threads Integer. Only used for `SPLIT_PARALLEL()`. 0 means "use all available".
#'
#' @return Integer vector (length p) with 0/1 inclusion indicators (one sampled draw).
#' @export
fast_model_selection <- function(
    y, x,
    center = TRUE,
    scale = TRUE,
    XtXprecomp = ifelse(ncol(x) < 1e4, TRUE, FALSE),
    niter = 5000L,
    thinning = 1L,
    burnin = round(niter / 10),
    deltaini = NULL,
    priorCoef = NULL,
    tau = 0.348,
    priorDelta = modelbbprior(1, 1),
    phi = 1.0,
    knownphi = 0L,
    priorSkew = 0.348,
    thinit = NULL,
    initpar_type = 0L,
    method = 1L,
    hesstype = 1L,
    optimMethod = 1L,
    optim_maxit = 10L,
    B = 100000L,
    r = 1L,
    alpha = 0.05,
    lambda = 1.0,
    computation_strategy = STANDARD(),
    n_units = 1L,
    max_threads = 0L
) {
  if (is.data.frame(x)) x <- as.matrix(x)
  y <- as.numeric(y)
  storage.mode(x) <- "double"

  p <- ncol(x)

  # allow priorCoef to mimic mombf-ish usage
  if (!is.null(priorCoef)) {
    if (is.numeric(priorCoef) && length(priorCoef) == 1) {
      tau <- as.double(priorCoef)
    } else if (is.list(priorCoef) && !is.null(priorCoef$tau)) {
      tau <- as.double(priorCoef$tau)
    } else {
      stop("priorCoef must be NULL, a scalar tau, or a list/object with element $tau.")
    }
  }

  # decode priorDelta exactly like your estimate_model wrapper
  prDelta <- 1L
  prDeltap <- 0.5
  a <- 1.0
  b <- 1.0
  if (inherits(priorDelta, "modelbbprior")) {
    a <- as.double(priorDelta$a)
    b <- as.double(priorDelta$b)
    prDelta <- 2L
  } else if (inherits(priorDelta, "modelbinomprior")) {
    prDeltap <- as.double(priorDelta$p)
    prDelta <- 1L
  } else {
    stop("priorDelta must be either a modelbbprior() or modelbinomprior() object.")
  }

  if (is.null(deltaini)) deltaini <- integer(p)
  if (length(deltaini) != p) stop("deltaini must have length ncol(x).")
  deltaini <- as.integer(deltaini)

  if (!is.null(thinit) && length(thinit) != p) {
    stop("thinit must be NULL or have length ncol(x).")
  }

  rcpp_fast_model_selection(
    y = y,
    x = x,
    niter = as.integer(niter),
    thinning = as.integer(thinning),
    burnin = as.integer(burnin),
    deltaini = deltaini,
    center = isTRUE(center),
    scale = isTRUE(scale),
    XtXprecomp = isTRUE(XtXprecomp),
    phi = as.double(phi),
    tau = as.double(tau),
    priorSkew = as.double(priorSkew),
    thinit = thinit,
    initpar_type = as.integer(initpar_type),
    method = as.integer(method),
    hesstype = as.integer(hesstype),
    optimMethod = as.integer(optimMethod),
    optim_maxit = as.integer(optim_maxit),
    B = as.integer(B),
    knownphi = as.integer(knownphi),
    r = as.integer(r),
    alpha = as.double(alpha),
    lambda = as.double(lambda),
    prDelta = as.integer(prDelta),
    prDeltap = as.double(prDeltap),
    prDelta_a = as.double(a),
    prDelta_b = as.double(b),
    computation_strategy = as.integer(computation_strategy),
    n_units = as.integer(n_units),
    max_threads = as.integer(max_threads)
  )
}

#' Minimal MOM prior spec (tau only)
#' @param tau numeric scalar
#' @export
momprior <- function(tau = 0.348) {
  structure(list(tau = as.double(tau)), class = "momprior")
}