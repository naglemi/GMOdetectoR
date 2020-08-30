normalization_and_regression_over_all_data_one_timepoint <- function(files_to_loop, maximum_CPU_number, job_id, by_explant,
                                                                     intercept, FP_threshold, Chl_threshold, grid_type,
                                                                     record_residuals, plotting, desired_wavelength_range,
                                                                     fluorophore_ID_vector, parallel_mode = TRUE){

  #try(dev.off())
  gc()

  # chroma_standard_path <- files_to_loop[1]
  # chroma_in_backend <<- load_image(image_path = chroma_standard_path)
  #browser()
  print(paste0('To loop over ', length(files_to_loop), ' files'))

  if(parallel_mode == TRUE){
    cl <- makeCluster(min(maximum_CPU_number, length(files_to_loop)), type = "FORK", outfile="output/cluster.outfile")
    registerDoParallel(cl)
    on.exit(stopCluster(cl))
  }

  files_to_loop <- files_to_loop[!grepl("chroma", files_to_loop)]


  # Create a log so we can print within dopar -------------------------------
  # https://stackoverflow.com/questions/10903787/how-can-i-print-when-using-dopar
  # log.socket <- make.socket(port=4000)

  # Run parallel ------------------------------------------------------------

  foreach(i=1:length(files_to_loop),
          .packages = c("shiny", #"RcppEigen",
                        "hyperSpec", "data.table",
                        "scales", "gridExtra",
                        "tools", "stringr"),
          .export = c("normalization_and_regression_over_one_plate"),
          .errorhandling = 'stop') %dopar% {
            #for(i in 1:length(files_to_loop)){
            print(paste0("Running for file ", i, " of ", length(files_to_loop)))
            print(files_to_loop[i])
            normalization_and_regression_over_one_plate(files_to_loop = files_to_loop,
                                                        i = i,
                                                        chroma_in_backend = chroma_in_backend,
                                                        job_id = job_id,
                                                        by_explant = by_explant,
                                                        intercept = intercept,
                                                        FP_threshold = FP_threshold,
                                                        Chl_threshold = Chl_threshold,
                                                        grid_type = grid_type,
                                                        record_residuals = record_residuals,
                                                        plotting = plotting,
                                                        desired_wavelength_range = desired_wavelength_range,
                                                        fluorophore_ID_vector = fluorophore_ID_vector)
            #cat("\n\n")
          }
  # Redundant. Remove this line if including on.exit(stopCluster(cl))
  #stopCluster(cl)
}
