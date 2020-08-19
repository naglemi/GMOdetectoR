#' Categorize an explant as transgenic or not based on whether X pixels meet Y intensity threshold of fluorescent protein signal
#'
#' @param CLS_table_path
#' @param npixel
#' @param intensity
#'
#' @return
#' @export
#'
#' @examples \dontrun{categorize_via_thresholds(files[i], 10, 5)}

categorize_via_thresholds <- function(CLS_table_path, npixel, intensity){
  CLS_table_in <- fread(CLS_table_path)
  n_pixels_passing_threshold <- length(which(CLS_table_in$DsRed > intensity, arr.ind = TRUE))
  if(n_pixels_passing_threshold >= npixel){
    transgenic_regeneration <- TRUE
  } else {
    transgenic_regeneration <- FALSE
  }
  return(c(CLS_table_path, transgenic_regeneration))
}

#for file in list.files("/scratch2/NSF_GWAS/")
