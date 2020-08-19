#' Prepare a table to record regression results to
#'
#' @param spectrum_in_to_CLS The spectrum from a hyperSpec object, which we need to know the # of rows (pixels) in
#' @param fluorophore_ID_vector A vector of strings representing each fluorophore regression is run over
#'
#' @return An empty CLS_table which results are recorded to during regression
#' @export
#'
#' @examples \dontrun{CLS_table <- prepare_CLS_table(spectrum_in_to_CLS, fluorophore_ID_vector)}
prepare_CLS_table <- function(spectrum_in_to_CLS, fluorophore_ID_vector, record_residuals){
  #stop("Debug making CLS table")
  #browser()
  if(record_residuals==TRUE){
    x <- ncol(spectrum_in_to_CLS)-2
  }
  if(record_residuals==FALSE){
    x <- 0
  }

  CLS_table <- data.table(cbind(as.numeric(spectrum_in_to_CLS$rows),
                                as.numeric(spectrum_in_to_CLS$cols),
                                matrix(as.numeric(NA),
                                       nrow=dim(spectrum_in_to_CLS)[1],
                                       ncol=length(fluorophore_ID_vector)+x)))

  # FIX THIS
  colnames(CLS_table)[1:2] <- c("rows", "cols")

  colnames(CLS_table)[3:(2+length(fluorophore_ID_vector))] <- fluorophore_ID_vector

  if(record_residuals==TRUE){
    colnames(CLS_table)[(3+length(fluorophore_ID_vector)):ncol(CLS_table)] <- colnames(spectrum_in_to_CLS)[3:ncol(spectrum_in_to_CLS)]
  }

  return(CLS_table)
}
