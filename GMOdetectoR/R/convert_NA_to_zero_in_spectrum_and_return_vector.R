#' Convert NA to 0 in a vector
#'
#' @param a vector representing the spectrum of a fluorophore
#' @param verbose
#'
#' @return a vector equal to the input but with all NA replaced with 0
#' @export
#'
#' @examples
convert_NA_to_zero_in_spectrum_and_return_vector <- function(spectrum_in, verbose=FALSE){
  #stop("break here to debug")
  vector <- as.vector(spectrum_in$'Normalized intensity (fitted)')
  if(verbose==TRUE) print(vector[1:5])
  vector[which(is.na(vector)==TRUE)] <- 0
  if(verbose==TRUE) print(vector[1:5])
  return(as.numeric(as.character(vector)))

}
