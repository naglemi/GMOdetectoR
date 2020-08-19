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
                        t_spectrum_in_to_CLS){

  #start_time <- Sys.time()
  if(backend=="ByHand"){
    #stop("Debug here")
    A <- solve(t(Xmatrix) %*% Xmatrix)
    Z <- A %*% t(Xmatrix)
  }else{
    A <- NA
    Z <- NA
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
