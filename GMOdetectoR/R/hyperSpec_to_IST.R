hyperSpec_to_IST <- function(hyperSpec_object, desired_wavelength_range){
  # These lines should be moved to another function where more appropriate
  wavelengths <- hyperSpec_object@wavelength
  desired_wavelengths <- which( wavelengths>=min(desired_wavelength_range) & wavelengths<=max(desired_wavelength_range) ,
                                arr.ind = TRUE)

  hyperSpec_object <- data.table(cbind(hyperSpec_object$x, hyperSpec_object$y, hyperSpec_object[[]][,desired_wavelengths]))
  colnames(hyperSpec_object) <- c("cols", "rows", wavelengths[desired_wavelengths])
  return(hyperSpec_object)
}
