# fitting a smooth curve https://stackoverflow.com/questions/3480388/how-to-fit-a-smooth-curve-to-my-data-in-r
#' Title
#'
#' @param spectra_path an absolute or relative data path to a csv file containing spectra with columns `emission wavelength (nm)` and `Normalized emission`
#' @param wavelengths a vector of wavelengths (for user-specific camera) we wish to fit a spectra to... can be taken from 'wavelengths' attribute of hyperSpec object (e.g. hyperimage@wavelengths)
#' @param plot TRUE or FALSE - whether to make plot of fitted spectra in addition to outputting it as a dataframe
#'
#' @return a data frame with spectra fitted to the specific wavelengths on the user's camera
#' @export
#'
#' @examples read_plot_spectra(path, wavelengths = wavelengths)
#'
#'
read_plot_spectra <- function(spectra_path, wavelengths, plot=FALSE){
  pub_spectra_in <- fread(spectra_path)
  pub_emission_spectrum <- data.frame(wavelength=pub_spectra_in$`emission wavelength (nm)`,
                                      intensity=pub_spectra_in$`Normalized emission`)
  fit <- loess(intensity~wavelength, data=pub_emission_spectrum, span=0.1)

  if(plot==TRUE){
    plot(pub_emission_spectrum,
         main=(paste0(basename(tools::file_path_sans_ext(spectra_path)), " (Published)")),
         xlab="Wavelength",
         ylab="Normalized emission")

    lines(pub_emission_spectrum$wavelength,
          predict(fit,pub_emission_spectrum$wavelength), col='red', lwd=2)
  }

  predict_from <- data.frame(wavelength=wavelengths)
  predictions <- predict(fit, predict_from)
  scaled_emission_spectra <- cbind(wavelengths, predictions)

  if(plot==TRUE){
    plot(scaled_emission_spectra,
         main=(paste0(basename(tools::file_path_sans_ext(spectra_path)), " (Fitted for macroPhor Array camera)")),
         xlab="Wavelength",
         ylab="Normalized emission")

    # Note here we are plotting the published lines over our data (should be the same as if we used our lines)
    lines(pub_emission_spectrum$wavelength,
          predict(fit,pub_emission_spectrum$wavelength), col='red', lwd=2)
  }

  #return(predictions)
  fitted_spectra_dataframe <- as.data.frame(cbind(wavelengths, predictions))
  colnames(fitted_spectra_dataframe) <- c("Wavelength", "Normalized intensity (fitted)")

  return(fitted_spectra_dataframe)
}
