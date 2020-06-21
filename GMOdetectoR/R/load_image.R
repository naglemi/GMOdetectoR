#' Load a hyperspectral image
#'
#' @param image_path the path to a hyperspectral image to be loaded as a hyperSpec object
#'
#' @return A hyperSpec object with wavelengths added as column names for the spectra
#' @export
#'
#' @examples \dontrun{load_image(image_path)}
load_image <- function(image_path){
  image_in <- read.ENVI(image_path,
                        headerfile=paste0(tools::file_path_sans_ext(image_path), ".hdr"))
  colnames(image_in$spc) <- image_in@wavelength
  return(image_in)
}
