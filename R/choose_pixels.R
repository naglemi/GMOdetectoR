#' Title
#'
#' @param FP_lambda
#' @param Chl_lambda
#' @param full_spectrum
#' @param CLS_threshold_FP
#' @param CLS_threshold_Chl
#' @param Mode
#' @param Indices_submitted_to_subset_to
#' @param verbose
#'
#' @return
#' @export
#'
#' @examples
choose_pixels <- function(FP_lambda='582.6952', Chl_lambda='682.8465', full_spectrum=this_full_spectrum, CLS_threshold_FP, CLS_threshold_Chl,
                          Mode = mode, Indices_submitted_to_subset_to = indices_submitted_to_subset_to, verbose=FALSE){

  all_pixel_indices <- seq.int(nrow(full_spectrum))

  if(verbose==TRUE){
    print("Length of all_pixel_indices")
    print(length(all_pixel_indices))
  }

  if (Mode == "threshold"){

    # Need to chop off x and y columns up front in this mode
    # Actually, I don't think this is needed since all this function outputs is a
    # ... list of pixel indices
    full_spectrum=full_spectrum[,3:ncol(full_spectrum)]

    if(max(full_spectrum[,get(FP_lambda)])<CLS_threshold_FP){
      print(paste0("Maximum signal for FP is ", max(full_spectrum[,get(FP_lambda)])))
      warning("Need a lower threshold for FP to get these pixels")
    }
    if(max(full_spectrum[,get(Chl_lambda)])<CLS_threshold_Chl){
      print(paste0("Maximum signal for Chl is ", max(full_spectrum[,get(Chl_lambda)])))
      warning("Need a lower threshold for Chl to get these pixels")
    }

    pixels_list_first_pass <- which(arr.ind = TRUE, x = full_spectrum[,get(FP_lambda)] > CLS_threshold_FP)

    pixels_with_significant_Chl <- all_pixel_indices[which(full_spectrum[,get(Chl_lambda)] > CLS_threshold_Chl)]

    pixels_with_something <- unique(c(pixels_list_first_pass, pixels_with_significant_Chl))
  }

  if (Mode == "integrate"){
    pixels_with_something <- all_pixel_indices[which(Indices_submitted_to_subset_to[,4] > 0.05)]
  }

  return(pixels_with_something)
}
