run_parallel <- function(timepoint,
                         maximum_CPU_number,
                         intercept,
                         FP_threshold,
                         Chl_threshold,
                         by_explant,
                         grid_type,
                         record_residuals,
                         plotting,
                         desired_wavelength_range = c(545, 722),
                         fluorophore_ID_vector,
                         parallel_mode = TRUE
                         ){

  print(paste0('Running for timepoint ', timepoint))

  files_to_loop <- list.files(timepoint,
                              pattern=".raw",
                              full.names = TRUE)

  job_id <- paste0(str_split_fixed(timepoint, "/", 5)[5], "/intercept", intercept,
                   "/",
                   Sys.Date())

  print(paste0("Starting workflow for ", job_id))

  normalization_and_regression_over_all_data_one_timepoint(files_to_loop = files_to_loop,
                                                           maximum_CPU_number = maximum_CPU_number,
                                                           job_id = job_id,
                                                           by_explant = by_explant,
                                                           FP_threshold = FP_threshold,
                                                           Chl_threshold = Chl_threshold,
                                                           intercept = intercept,
                                                           grid_type = grid_type,
                                                           record_residuals = record_residuals,
                                                           plotting,
                                                           desired_wavelength_range = desired_wavelength_range,
                                                           fluorophore_ID_vector = fluorophore_ID_vector,
                                                           parallel_mode = parallel_mode)

  output <- check_success_high_throughput(files_to_loop = files_to_loop,
                                          sum_stats_file = paste0("/scratch2/NSF_GWAS/GMOdetectoR/output/",
                                                                  job_id,
                                                                  "/sum_stats.csv"),
                                          grid_type = grid_type,
                                          by_explant = by_explant)

  files_failed_first_time <- unique(c(files_to_loop[1], output[is.na(output$V3),]$file))
  length(files_failed_first_time)
  # 18

  if(length(files_failed_first_time)>1){
    stop("Not all jobs completed.")
  }
  print("Complete. On to the next one...")
}
