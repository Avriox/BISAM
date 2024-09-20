b_ism <- function(data, include_constant = TRUE, tfe = FALSE, ife = FALSE, 
                  iis = FALSE, sis = FALSE, y_index = 1, i_index = 2, t_index = 3,
                  Ndraw = 1000, Nburn = 100, b_prior = "normal", lambda_b = 0,
                  c0 = 0.01, C0 = 100, geweke = FALSE, 
                  model_sel_optimization = "No Optimization") {
  
  # Check if data is a matrix or data frame
  if (!is.matrix(data) && !is.data.frame(data)) {
    stop("'data' must be a matrix or data frame")
  }
  
  # Convert data to matrix if it's a data frame
  if (is.data.frame(data)) {
    data <- as.matrix(data)
  }
  
  # Call the C++ wrapper function
  result <- b_ism_wrapper(data, include_constant, tfe, ife, iis, sis, 
                          y_index, i_index, t_index, Ndraw, Nburn, 
                          b_prior, lambda_b, c0, C0, geweke, 
                          model_sel_optimization)
  
  # Return the result
  return(result)
}