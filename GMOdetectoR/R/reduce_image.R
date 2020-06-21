#' Extract specific wavelengths from hyperSpec object to prepare for false color image
#'
#' @param image_in A hyperSpec object obtained from \code{load_image}
#' @param FP An optional parameter (set to "unspecified" by default) that allows the user to select specific FPs with known fluorescence peaks rather than inputting the peak wavelength itself
#' @param wavelength_make_green If option FP is not submitted, the wavelength to which green false color will be applied
#' @param wavelength_make_red If option FP is not submitted, the wavelength to which red false color will be applied
#'
#' @return An image spectrum table containing five columns: rows, cols, r, g, b
#' @export
#'
#' @examples \dontrun{reduce_image(image_in, FP="DsRed")}
reduce_image <- function(image_in = image_spectrum, FP="unspecified",
                         wavelength_make_green, wavelength_make_red){

  # THIS LINE EATS MEMORY
  image_spectrum_loaded <- image_in$spc
  wavelengths <- image_in@wavelength
  colnames(image_spectrum_loaded) <- image_in@wavelength

  # Find the closest wavelengths to those specified by user
  # http://adomingues.github.io/2015/09/24/finding-closest-element-to-a-number-in-a-list/

  if(FP=="unspecified"){
    selected_FP_wavelength <- which.min(abs(wavelengths - wavelength_make_green))
    selected_Chl_wavelength <- which.min(abs(wavelengths - wavelength_make_red))

    image_table <- as.data.table(cbind(image_in$y, image_in$x,
                                       image_spectrum_loaded[,which(colnames(image_spectrum_loaded)==selected_Chl_wavelength)],
                                       image_spectrum_loaded[,which(colnames(image_spectrum_loaded)==selected_FP_wavelength)],
                                       rep(0, length(image_in$x))))
  } else {
    if(FP=="DsRed"){
      FPlambda <- '582.6952'
    }
    if(FP=="GFP"){
      FPlambda <- '508.763'
    }

    image_table <- as.data.table(cbind(image_in$y, image_in$x,
                                       image_spectrum_loaded[,'682.8465'],
                                       image_spectrum_loaded[,which(colnames(image_spectrum_loaded)==FPlambda)],
                                       rep(0, length(image_in$x))))
  }




  colnames(image_table) <- c("rows", "cols", "r", "g", "b")
  return(image_table)
}
