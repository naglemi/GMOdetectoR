cap_signal_in_table <- function(image_spectrum_table, max_intensity_FP, max_intensity_Chl){
  image_spectrum_table$g[which(image_spectrum_table$g > max_intensity_FP)] <- max_intensity_FP
  image_spectrum_table$r[which(image_spectrum_table$r > max_intensity_Chl)] <- max_intensity_Chl
  return(image_spectrum_table)
}
