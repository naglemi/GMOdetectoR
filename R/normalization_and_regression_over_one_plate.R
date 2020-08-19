normalization_and_regression_over_one_plate <- function(files_to_loop, i, chroma_in_backend, job_id, by_explant=TRUE, grid_type=12,
                                                        intercept, FP_threshold, Chl_threshold, record_residuals, plotting,
                                                        desired_wavelength_range, fluorophore_ID_vector){
  gc() # Check if this line helps with memory issues, crashing
  filename <- files_to_loop[i]

  print(paste0("Running for plate ", files_to_loop[i]))

  img_in_backend <- load_image(filename)

  # pixels_chroma_in <- nrow(chroma_in_backend)
  #
  # number_rows_to_include_in_chroma_colmeans <- 100
  # middle_of_chroma <- round(nrow(chroma_in_backend$spc)/2, digits = 0)
  # chroma_colmeans <- colMeans(chroma_in_backend$spc[(middle_of_chroma-(number_rows_to_include_in_chroma_colmeans/2)):(middle_of_chroma+(number_rows_to_include_in_chroma_colmeans/2)),])
  #
  # scaling_factor <- mean(chroma_colmeans)
  #
  # img_in_backend$spc <- (img_in_backend$spc/chroma_colmeans)*scaling_factor

  input <- list()
  input$hys_CLS_PCA <- 2
  input$denoise_threshold_FP <- FP_threshold
  input$denoise_threshold_Chl <- Chl_threshold

  if(by_explant==TRUE){
    what_to_plot <- decide_what_to_plot(mode = c(2))
    for(j in 1:grid_type){
      input$grid_position <- j
      crop_and_plot(mode=what_to_plot,
                    image_to_crop = img_in_backend,
                    input_toggle=(2 %in% input$hys_CLS_PCA),
                    grid_item = input$grid_position,
                    image_type = "CLS",
                    first_pass_FP_threshold = FP_threshold,
                    first_pass_Chl_threshold = Chl_threshold,
                    #first_pass_FP_threshold = input$denoise_threshold_FP,
                    #first_pass_Chl_threshold = input$denoise_threshold_Chl,
                    name_to_parse = filename,
                    fluorophore_ID_vector = fluorophore_ID_vector,
                    grid_type = grid_type,
                    job_id = job_id,
                    intercept = intercept,
                    desired_wavelength_range = desired_wavelength_range,
                    record_residuals = record_residuals,
                    plotting = plotting)
    }
  }else{
    what_to_plot <- decide_what_to_plot(mode = c(1))

    crop_and_plot(mode=what_to_plot,
                  image_to_crop = img_in_backend,
                  input_toggle=(2 %in% input$hys_CLS_PCA),
                  # The below line is currently the same for any type of plot
                  #image_spectrum_table = image_spectrum_table_colored(),
                  #image_spectrum_table = image_spectrum_table_colored,
                  grid_item = "WholePlate",
                  image_type = "CLS",
                  first_pass_FP_threshold = input$denoise_threshold_FP,
                  first_pass_Chl_threshold = input$denoise_threshold_Chl,
                  name_to_parse = filename,
                  fluorophore_ID_vector = fluorophore_ID_vector,
                  grid_type=grid_type,
                  job_id = job_id,
                  intercept = intercept,
                  desired_wavelength_range = desired_wavelength_range,
                  record_residuals = record_residuals,
                  plotting = plotting,
                  left_edge = 0,
                  right_edge = 2000,
                  bottom_edge = 2000,
                  top_edge = 0,)
  }
}
