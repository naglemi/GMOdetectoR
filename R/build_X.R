#' Prepare a design matrix (X) for regression
#'
#' @param fluorophore_ID_vector a vector of fluorophores in .csv with spectra contained in \code{spectra_library}
#' @param intercept Boolean indicating whether the model should have a y-intercept that is fitted, rather than constrained to 0
#'
#' @return X (mm), a design matrix
#' @export
#'
#' @examples
build_X <- function(fluorophore_ID_vector, intercept, wavelengths, desired_wavelength_range){

  for (i in 1:length(fluorophore_ID_vector)){
    path <- paste0("spectra_library/", fluorophore_ID_vector[i], ".csv")
    if(!file.exists(path)){
      stop(paste0("Error: Spectrum for ", fluorophore_ID_vector[i], " not found in spectra_library folder"))
    }
    #stop("Debug here")
    spectra <- convert_NA_to_zero_in_spectrum_and_return_vector(read_plot_spectra(path,
                                                                                  wavelengths = wavelengths))
    #browser()
    # assign(paste0(fluorophore_ID_vector[i], "_vector"), spectra)
    #
    #
    # fluorophore_lambda_vector_to_retrieve <- paste0(fluorophore_ID_vector[i], "_vector")
    # this_fluorophore_lambda_vector <- get(fluorophore_lambda_vector_to_retrieve)
    # If it is the first one, initialize
    if(i==1 & intercept==1){
      mm <- as.matrix(cbind(intercept, spectra))
    }
    if(i==1 & intercept==0){
      mm <- as.matrix(spectra)
    }
    if(i>1){
      mm <- cbind(mm, spectra)
    }
  }

  if(intercept==1){
    colnames(mm) <- c("Intercept", paste0(fluorophore_ID_vector, "_vector"))
  }
  if(intercept==0){
    colnames(mm) <- c(paste0(fluorophore_ID_vector, "_vector"))
  }

  mm <- mm[which(wavelengths>=min(desired_wavelength_range)&wavelengths<=max(desired_wavelength_range)),]



  # browser()
  # x <- colnames(img_in_backend$spc)
  # x <- x[x>=517]
  # x <- x[x<=722]
  # rownames(mm) <- format(as.numeric(x), digits = 1)
  # colnames(mm) <- c("Intercept", "DsRed", "ZsYellow", "Chlorophyll", "Diffraction")
  # y <- as.matrix((22*mm[,3])+(28*mm[,2])+(8*mm[,4]))
  # y <- as.matrix(y[nrow(y):1])
  # plot(mm[nrow(mm):1,2:4], border=NA, col=viridis(200), )
  # rownames(y) <- rownames(mm)
  # plot(as.matrix(y[nrow(y):1]), border=NA, col=viridis(200), )

  return(mm)
}
