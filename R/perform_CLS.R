#' Run regression over selected pixels of an image
#'
#' @param pixels_to_test A list of integer indices for the rows of a spectrum from a hyperSpec object for wish we wish to run
#' @param verbose
#' @param fluorophore_ID_vector
#' @param backend
#' @param record_residuals
#' @param intercept
#'
#' @return
#' @export
#'
#' @examples
perform_CLS <- function(pixels_to_test,
                        verbose=FALSE,
                        fluorophore_ID_vector,
                        backend="ByHand",
                        record_residuals,
                        intercept,
                        spectrum_in_to_CLS,
                        CLS_table,
                        Xmatrix,
                        t_spectrum_in_to_CLS,
                        vectorized = FALSE){

  #start_time <- Sys.time()
  if(backend=="ByHand"){
    #stop("Debug here")
    A <- solve(t(Xmatrix) %*% Xmatrix)
    Z <- A %*% t(Xmatrix)
  }else{
    A <- NA
    Z <- NA
  }

  if(vectorized == TRUE){
    BETA <- Z %*% t_spectrum_in_to_CLS
    RESID <- round(t_spectrum_in_to_CLS - (Xmatrix %*% BETA))
    P <- ncol(Xmatrix) - intercept # minus one if we have intercept because we don't count the intercept as a parameter
    N <- nrow(Xmatrix)
    ptm <- proc.time()
    profvis({

      use_condaenv("r-reticulate")

      np <- import("numpy", convert=FALSE)
      tf <- import("tensorflow", convert=FALSE)

      gpu <- import("gnumpy", convert=FALSE)

      (RESIDt <- np$transpose(RESID))
      #(RESIDtRESID <- np$matmul(RESIDt, RESID))

      (RESIDtRESID <- tf$mul(RESIDt, RESID))

      py_to_r(RESIDtRESID)

      MSE <- diag(RESIDtRESID/(N-P-1))

      #MSE <- diag((t(RESID) %*% RESID)/(N-P-1))
    })

    print((proc.time()-ptm))
    Beta.covar.matrix <- as.vector(MSE)*A
    Beta.se <- sqrt(diag(Beta.covar.matrix))
  }

  # Need to get rid of the row and col columns I used for indexing earlier
  #spectrum_in_to_CLS$rows <- NULL
  #spectrum_in_to_CLS$cols <- NULL
  if (shiny::isRunning() == "PROGRESSBARSLOWSDOWN"){
    #if (shiny::isRunning() == TRUE){
    withProgress(message = "Running CLS over all pixels passing denoising",
                 max=length(pixels_to_test),
                 value = 0, {
                   for(i in pixels_to_test){
                     #print(i)
                     spectrum_formatted <- t_spectrum_in_to_CLS[,i]

                     CLS_table <- model_and_record_single_pixel(i, fluorophore_ID_vector = fluorophore_ID_vector,
                                                                A = A,
                                                                Z = Z,
                                                                backend = backend,
                                                                record_residuals = record_residuals,
                                                                intercept = intercept,
                                                                CLS_table = CLS_table,
                                                                Xmatrix = Xmatrix,
                                                                spectrum_formatted = spectrum_formatted)
                     incProgress(1)
                   }
                 })
  }else{
    for(i in pixels_to_test){
      #browser()
      spectrum_formatted <- t_spectrum_in_to_CLS[,i]
      CLS_table <- model_and_record_single_pixel(i, fluorophore_ID_vector = fluorophore_ID_vector,
                                                 A = A,
                                                 Z = Z,
                                                 backend = backend,
                                                 record_residuals = record_residuals,
                                                 intercept = intercept,
                                                 CLS_table = CLS_table,
                                                 Xmatrix = Xmatrix,
                                                 spectrum_formatted = spectrum_formatted)
    }
  }

  options(warn=0)
  return(CLS_table)
}
