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
#' @param geweke Boolean for Geweke test
#' @param use_phiinit Boolean for using initial phi values
#' @param const_val Boolean for including a constant
#' @param ife Boolean for individual fixed effects
#' @param tfe Boolean for time fixed effects
#' @param iis Boolean for indicator saturation
#' @param sis Boolean for stepshift saturation
#'
#' @return A list containing MCMC results and summary statistics
#' @export
b_ism <- function(
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
    geweke = FALSE,
    use_phiinit = TRUE,
    const_val = FALSE,
    ife = FALSE,
    tfe = FALSE,
    iis = TRUE,
    sis = TRUE
) {
    # Convert data to matrix if it's a data frame
    if (is.data.frame(data)) {
        data <- as.matrix(data)
    }

    # Adjust indices for R to C++ (0-based indexing)
    i_index_cpp <- i_index - 1
    t_index_cpp <- t_index - 1
    y_index_cpp <- y_index - 1

    # Call the C++ function
    result <- b_ism_rcpp(
        data, i_index_cpp, t_index_cpp, y_index_cpp, Ndraw, Nburn, b_prior,
        lambda_b, c0, C0, va, vb, tau, geweke, use_phiinit,
        const_val, ife, tfe, iis, sis
    )

    return(result)
}