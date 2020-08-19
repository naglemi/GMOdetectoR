CLS_workflow <- function(spectrum_in_to_CLS = image_full_spectrum_cropped,
                         pass_FP_threshold_from_input,
                         pass_Chl_threshold_from_input,
                         mode = "threshold",
                         fluorophore_ID_vector,
                         filename,
                         grid_item,
                         tissue_type,
                         record_residuals,
                         intercept=0,
                         desired_wavelength_range = c(0,1000),
                         job_id){

  print(paste0("Record residuals: ", record_residuals))
  #browse()
  #stop("Stop here to debug")

  # CLS_table <- prepare_CLS_table(spectrum_in_to_CLS = spectrum_in_to_CLS,
  #                                fluorophore_ID_vector = fluorophore_ID_vector)

  #browser()
  mm <- build_X(fluorophore_ID_vector = fluorophore_ID_vector,
                intercept = intercept,
                wavelengths = colnames(spectrum_in_to_CLS)[3:ncol(spectrum_in_to_CLS)],
                desired_wavelength_range = desired_wavelength_range)

  print("Done building mm")
  #browser()

  chosen_pixels <- choose_pixels(full_spectrum=spectrum_in_to_CLS,
                                 CLS_threshold_FP=pass_FP_threshold_from_input,
                                 CLS_threshold_Chl=pass_Chl_threshold_from_input,
                                 Mode = mode)

  print("Done choosing pixels")
  CLS_table <- prepare_CLS_table(spectrum_in_to_CLS = spectrum_in_to_CLS,
                                 fluorophore_ID_vector = fluorophore_ID_vector,
                                 record_residuals = record_residuals)


  print("CLS table ready")
  # NEED TO ADD THIS TO GET RID OF BOTTLENECK OF CONVERTING ONE COLUMN AT A TIME TO MATRIX
  # v0.27
  # ALSO NEED TO TRANSPOSE SINCE SELECTION OF COLS FROM MATRIX IS SUPERFAST
  t_spectrum_in_to_CLS <- t(as.matrix(spectrum_in_to_CLS[,3:ncol(spectrum_in_to_CLS)]))

  #browser()
  #full_spectrum <- t(as.matrix(spectrum_in_to_CLS))
  print("CLS_spectrum transposed")
  # For now, make this conditional because residual tables still depend on non-transposed version of this.
  if(record_residuals==FALSE){
    spectrum_in_to_CLS <- NULL
  }
  gc()
  print("Run CLS now")
  CLS_table <-perform_CLS(pixels_to_test = chosen_pixels,
                          fluorophore_ID_vector = fluorophore_ID_vector,
                          record_residuals = record_residuals,
                          intercept = intercept,
                          CLS_table = CLS_table,
                          Xmatrix = mm,
                          t_spectrum_in_to_CLS = t_spectrum_in_to_CLS)
  print("CLS done.")
  # CLS_table <- finish_processing_CLS_table(CLS_table,
  #                                          fluorophore_ID_vector = fluorophore_ID_vector,
  #                                          filename = filename,
  #                                          grid_item = grid_item,
  #                                          tissue_type = tissue_type,
  #                                          record_residuals = record_residuals,
  #                                          job_id = job_id)

  #browser()
  #options(warn=1)
  print("About to go into record residuals plot?")
  if(record_residuals==TRUE){
    print("Yes.")
    make_residual_plots(CLS_table = CLS_table,
                        t_spectrum_in_to_CLS = t_spectrum_in_to_CLS,
                        spectrum_in_to_CLS = spectrum_in_to_CLS,
                        bins = 200,
                        job_id = job_id,
                        filename = filename,
                        grid_item = grid_item)
  }

  # Currently this is one twice, once here and again for plots
  cumulative_t_stats <- colSums(CLS_table[,3:(2+length(fluorophore_ID_vector))], na.rm = TRUE)

  write_tables(table_subfolder = "CLS_tables",
               table = CLS_table,
               filename = filename,
               grid_item = grid_item,
               fluorophore_ID_vector = fluorophore_ID_vector,
               job_id = job_id)

  # if(record_residuals==TRUE){
  #   write_tables(table_subfolder = "residual_tables",
  #                table = list_CLS_residual_tables[[2]],
  #                filename = filename,
  #                grid_item = grid_item)
  # }

  write_sum_stats(job_id = job_id,
                  fluorophore_ID_vector = fluorophore_ID_vector,
                  filename = filename,
                  grid_item = grid_item,
                  tissue_type = tissue_type,
                  cumulative_t_stats = cumulative_t_stats)

  return(CLS_table)
}
