crop_and_plot <- function(image_spectrum_table, mode="whole_plate", grid_item, image_type="hyperspectral",
                          first_pass_FP_threshold, first_pass_Chl_threshold, image_to_crop = img_in_backend, input_toggle = TRUE,
                          name_to_parse, CLSmode = "threshold", indices_submitted_to_subset_to, verbose=FALSE,
                          fluorophore_ID_vector, tissue_type = "unsegmented", record_residuals, grid_type, job_id,
                          intercept, desired_wavelength_range, plotting = "TRUE", ...){
  # In addition to cropping and plotting, this function runs CLS between cropping and plotting if in CLS mode

  #browser()
  if(input_toggle==FALSE){
    return(NA)
  }

  if(image_type=="crop_segment_output"){
    warning('No need to use crop_and_plot for crop_segment_output. Returning input image_spectrum_table')
    return(image_spectrum_table)
  }
  #browser()
  image_to_crop <- crop(image_being_cropped = image_to_crop,
                               mode = mode,
                               grid_file = NA,
                               image_type = image_type,
                               grid_type = grid_type,
                               grid_item = grid_item,
                               desired_wavelength_range = desired_wavelength_range,
                               ...)
  gc()

  if (image_type=="CLS"){
    #stop("Debug in crop_and_plot")
    image_to_crop <- CLS_workflow(spectrum_in_to_CLS = image_to_crop,
                                         pass_FP_threshold_from_input = first_pass_FP_threshold,
                                         pass_Chl_threshold_from_input = first_pass_Chl_threshold,
                                         mode = CLSmode,
                                         fluorophore_ID_vector = fluorophore_ID_vector,
                                         filename = name_to_parse,
                                         grid_item = grid_item,
                                         tissue_type = tissue_type,
                                         record_residuals = record_residuals,
                                         intercept = intercept,
                                         desired_wavelength_range = desired_wavelength_range,
                                         job_id = job_id)

    if(plotting==FALSE){
      return("Not plotting")
    }
    cumulative_t_stats <- colSums(image_to_crop[,3:(3+length(fluorophore_ID_vector)-1)], na.rm = TRUE)
  }

  #return("Not plotting to save time")

  #stop("Debug crop_and_plot")

  if(mode=="whole_plate" | mode=="both"){

    p <- generate_plots(image_spectrum_table = image_to_crop,
                        FullID = FullID,
                        fluorophore_ID_vector = fluorophore_ID_vector,
                        name_to_parse = name_to_parse,
                        mode_information = paste0("FullPlate"),
                        grid_type = grid_type,
                        image_type = image_type,
                        job_id = job_id,
                        grid_item = grid_item,
                        cumulative_t_stats = cumulative_t_stats)

  }
  if(mode=="single_explant" | mode=="both"){

    q <- generate_plots(image_spectrum_table = image_to_crop,
                        FullID = FullID,
                        fluorophore_ID_vector = fluorophore_ID_vector,
                        name_to_parse = name_to_parse,
                        mode_information = paste0("SingleExplantNo",grid_item),
                        grid_type = grid_type,
                        image_type = image_type,
                        job_id = job_id,
                        grid_item = grid_item,
                        cumulative_t_stats = cumulative_t_stats)
  }
  if(mode=="both"){
    return(grid.arrange(p,q))
    #return(list(p,q))
  }
  if(mode=="whole_plate"){
    return(p)
  }
  if(mode=="single_explant"){
    return(q)
  }
}
