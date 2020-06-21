scale_image_spectrum_table <- function(image_spectrum_table,
                                       denoise_threshold_Chl,
                                       denoise_threshold_FP,
                                       max_intensity_Chl,
                                       max_intensity_FP,
                                       standardize_rgb = TRUE){
  if(standardize_rgb==TRUE){
    image_spectrum_table$r[2] <- denoise_threshold_Chl
    image_spectrum_table$g[2] <- denoise_threshold_FP
    image_spectrum_table$r[1] <- max_intensity_Chl
    image_spectrum_table$g[1] <- max_intensity_FP
  }
  image_spectrum_table$r <- rescale(image_spectrum_table$r,
                                    from=range(c(denoise_threshold_Chl,max_intensity_Chl)),
                                    to=c(0,255))
  image_spectrum_table$g <- rescale(image_spectrum_table$g,
                                    from=range(c(denoise_threshold_FP,max_intensity_FP)),
                                    to=c(0,255))
  return(image_spectrum_table)
}
