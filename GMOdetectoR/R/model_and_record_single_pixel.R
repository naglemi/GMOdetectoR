#' Run regression for a single pixel in hyperspectral image and record results
#'
#' @param i An integer representing the row of data to run regression over
#' @param Xmatrix The design matrix to be used for regression
#' @param full_spectrum_to_regress_over The Y matrix from which Yi=y, the vector for a specific pixel to regress over
#' @param fluorophore_ID_vector A vector of names for fluorophores in the design matrix
#' @param A Pre-calculated once for each design matrix and input to save time
#' @param Z Pre-calculated once for each design matrix and input to save time
#' @param backend Either \code{ByHand} or \code{FastLmPure}; the former is slighty faster and allows for further optimization in the future
#' @param record_residuals A boolean
#' @param intercept An integer set to 1 if the design matrix includes an intercept, or 0 otherwise
#'
#' @return Nothing is returned. Outputs are written to the global object \code{CLS_table}
#' @export
#'
#' @examples \dontrun{model_and_record_single_pixel(i, fluorophore_ID_vector = fluorophore_ID_vector,A = A, Z = Z, backend = backend,record_residuals = record_residuals)}
model_and_record_single_pixel <- function(i, Xmatrix = mm, fluorophore_ID_vector,
                                          A, Z, backend, record_residuals, intercept, CLS_table, spectrum_formatted){

  if(backend=="ByHand"){
    Beta <- Z %*% spectrum_formatted
    resid <- spectrum_formatted - (Xmatrix %*% Beta)
    p <- ncol(Xmatrix) - intercept # minus one if we have intercept because we don't count the intercept as a parameter
    n <- nrow(Xmatrix)
    MSE <- (t(resid) %*% resid)/(n-p-1)
    Beta.covar.matrix <- as.vector(MSE)*A
    Beta.se <- sqrt(diag(Beta.covar.matrix))
  }

  # Make sure colnames are set properly
  # df is n - p where n is # wavelengths and is parameters
  n <- nrow(Xmatrix)
  p <- ncol(Xmatrix)
  df <- n - p - 1

  for(j in 1:length(fluorophore_ID_vector)){

    if(intercept==1){
      k <- j+1
    }
    if(intercept==0){
      k <- j
    }

    if(backend=="FastLmPure"){
      this_t_stat <- model_out$coefficients[k]/model_out$se[k]
    }
    if(backend=="ByHand"){

      this_t_stat <- Beta[k]/Beta.se[k]
    }

    set(CLS_table, i, eval(fluorophore_ID_vector[j]), this_t_stat)
  }

  if(record_residuals == TRUE){
    #J <- colnames(CLS_table)[-(1:2+length(fluorophore_ID_vector))]
    set(CLS_table,
        i,
        ( (3+length(fluorophore_ID_vector) ):ncol(CLS_table)),
        as.list(resid))
  }

  return(CLS_table)
}
