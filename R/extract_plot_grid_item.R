extract_plot_grid_item <- function(image_spectrum_table,
                                   chroma_table,
                                   FC=TRUE,
                                   cap=TRUE,
                                   max_intensity_FP,
                                   max_intensity_Chl,
                                   dim_chlorophyll = numeric(1),
                                   scale=TRUE,
                                   normalize=TRUE,
                                   to_denoise=TRUE,
                                   denoise_threshold_FP,
                                   denoise_threshold_Chl,
                                   FP,
                                   min_for_rgb_scaling=0,
                                   standardize_rgb=TRUE,
                                   cropping_option,
                                   #filename=filename,
                                   filename,
                                   verbose=FALSE,
                                   grid_item){

  # This function takes an input of an image spectrum table, which already has false color applied,
  # and depending on options input by user, normalizes, denoises, caps, scales, dims chlorophyll channel
  # then outputs a new image spectrum table

  # Must normalize before scaling since normalization depends on cropping scaling depends on first pixel not being cropped out
  # just moved this to top in v0.19
  if(normalize==TRUE){
    image_spectrum_table <- normalize_image(image_spectrum_table, chroma_table)
  }

  ## Denoising
  if(to_denoise==TRUE){
    image_spectrum_table <- denoise(image_spectrum_table, threshold_FP=denoise_threshold_FP, threshold_Chl=denoise_threshold_Chl)
  }

  if(cap==TRUE){
    image_spectrum_table <- cap_signal_in_table(image_spectrum_table = image_spectrum_table,
                                                max_intensity_FP = max_intensity_FP,
                                                max_intensity_Chl = max_intensity_Chl)
  }



  # do this AFTER CROPPING and BEFORE SCALING so it doesn't screw up later steps or get nullified by prior steps

  if(scale==TRUE){
    image_spectrum_table <- scale_image_spectrum_table(image_spectrum_table = image_spectrum_table,
                                                       denoise_threshold_Chl = denoise_threshold_Chl,
                                                       denoise_threshold_FP = denoise_threshold_FP,
                                                       max_intensity_Chl = max_intensity_Chl,
                                                       max_intensity_FP = max_intensity_FP,
                                                       standardize_rgb = standardize_rgb)
  }

  # Will probably get rid of this later
  if(dim_chlorophyll<1){
    image_spectrum_table$r <- image_spectrum_table$r * dim_chlorophyll
  }

  # commenting out this line for GUI because we know img we're looking at and function was split up so
  # no cropping in this function any more
  # Need to uncomment in v0.19 for PhenotypeAssistant
  # Uncommenting didn't fix object'FullID' not found error.
  # I AM CURRENTLY DEFINING THIS BOTH IN FUNCTIONS extract_plot_grid_item and crop_and_plot
  # IF DEFINED IN THE FORMER, ONLY APPEARS FOR WHOLE PLATE
  FullID <<- paste0(parse_trayplateID(name_being_parsed = filename), "_exp", grid_item)
  # Debugging statement for v0.19

  if(verbose==TRUE){
    print(paste0("Full ID is ", FullID))
  }


  if(FC==TRUE){
    image_spectrum_table$color <- rgb(red= image_spectrum_table$r,
                                      green= image_spectrum_table$g,
                                      blue= image_spectrum_table$b,
                                      maxColorValue = max(image_spectrum_table$g))
  }
  return(image_spectrum_table)
}
