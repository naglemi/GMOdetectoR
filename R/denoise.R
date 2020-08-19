#' Denoise a reduced image spectrum table
#'
#' @param image_spectrum_table an object from \code{reduce_image} containing x and y positions of pixels and the r, g, b values for each pixel
#' @param threshold_FP the threshold at which fluorescent protein (labeled as green) is cut off
#' @param threshold_Chl the threshold at which chlorophyll (labeled as red) is cut off
#' @param verbose A Boolean; if true, the minimum chlorophyll prior to cutoff is printed.
#'
#' @return An \code{image_spectrum_table} with values below thresholds set to zero
#' @export
#'
#' @examples \dontrun{image_spectrum_table <- denoise(image_spectrum_table, threshold_FP=denoise_threshold_FP, threshold_Chl=denoise_threshold_Chl)}
denoise <- function(image_spectrum_table, threshold_FP, threshold_Chl, verbose=FALSE){

  if(verbose==TRUE){
    print("denoising....")
    print(min(image_spectrum_table$r))
  }
  image_spectrum_table$r[which(image_spectrum_table$r < threshold_Chl)] <- threshold_Chl
  image_spectrum_table$g[which(image_spectrum_table$g < threshold_FP)] <- threshold_FP
  return(image_spectrum_table)
}
