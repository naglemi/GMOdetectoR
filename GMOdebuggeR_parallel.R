library(shiny)
library(RcppEigen)
library(hyperSpec)
library(data.table)
library(scales)
library(gridExtra)
library(tools)
library(stringr)
library(foreach)
library(doParallel)
tempdir("/scratch2/NSF_GWAS/Rtmp/")
setwd("/scratch2/NSF_GWAS/GMOdetectoR/")
source("/scratch2/NSF_GWAS/GMOdetectoR/GMOdetectoRv0.30cool.R")

closeAllConnections()
gc()

normalization_and_regression_over_one_plate <- function(files_to_loop, i, chroma_in_backend){
  gc() # Check if this line helps with memory issues, crashing
  filename <- files_to_loop[i]
  
  img_in_backend <- load_image(filename)
  
  pixels_chroma_in <- nrow(chroma_in_backend)
  
  number_rows_to_include_in_chroma_colmeans <- 100
  middle_of_chroma <- round(nrow(chroma_in_backend$spc)/2, digits = 0)
  chroma_colmeans <- colMeans(chroma_in_backend$spc[(middle_of_chroma-(number_rows_to_include_in_chroma_colmeans/2)):(middle_of_chroma+(number_rows_to_include_in_chroma_colmeans/2)),])
  
  scaling_factor <- mean(chroma_colmeans)
  
  img_in_backend$spc <- (img_in_backend$spc/chroma_colmeans)*scaling_factor
  
  # The use of sweep on hyperspec objects is described in the hyperspec manual under 11.9: Normalization
  #img_in_backend$spc <- sweep (img_in_backend$spc, 2, chroma_colmeans, "/")*scaling_factor
  
  #img_in_backend$spc[1:pixels_chroma_in,] <- img_in_backend$spc[1:pixels_chroma_in,] / chroma_in_backend$spc
  
  for(i in 1:12){
    input$grid_position <- i
    crop_and_plot(mode=what_to_plot,
                  image_to_crop = img_in_backend,
                  input_toggle=(2 %in% input$hys_CLS_PCA),
                  # The below line is currently the same for any type of plot
                  #image_spectrum_table = image_spectrum_table_colored(),
                  #image_spectrum_table = image_spectrum_table_colored,
                  grid_item = input$grid_position,
                  image_type = "CLS",
                  first_pass_FP_threshold = input$denoise_threshold_FP,
                  first_pass_Chl_threshold = input$denoise_threshold_Chl,
                  name_to_parse = filename,
                  fluorophore_ID_vector = c("DsRed", "ZsYellow", "ChlA", "ChlB"),
                  record_residuals = TRUE,
                  grid_type=12)
  }
}

normalization_and_regression_over_all_data_one_timepoint <- function(files_to_loop){
  
  #try(dev.off())
  gc()
  
  chroma_standard_path <- files_to_loop[1]
  chroma_in_backend <<- load_image(image_path = chroma_standard_path)
  
  cl <<- makeCluster(min(12, length(files_to_loop)), type = "FORK", outfile="output/cluster.outfile")
  registerDoParallel(cl)
  
  files_to_loop[!grepl("chroma", files_to_loop)]
  #files_to_loop <- files_to_loop[41]
  
  foreach(i=1:length(files_to_loop),
          .packages = c("shiny", "RcppEigen", "hyperSpec", "data.table",
                        "scales", "gridExtra", "tools", "stringr"),
          .export = c("normalization_and_regression_over_one_plate")) %dopar% {
            normalization_and_regression_over_one_plate(files_to_loop = files_to_loop,
                                                        i = i,
                                                        chroma_in_backend = chroma_in_backend)
          }
  stopCluster(cl)
}

what_to_plot <- decide_what_to_plot(mode = c(2))
tempdir("/scratch2/NSF_GWAS/Rtmp/")

input <- list()
input$hys_CLS_PCA <- 2
#input$grid_position <- 18
input$denoise_threshold_FP <- 0
input$denoise_threshold_Chl <- 0

files_to_loop <- list.files("/scratch2/NSF_GWAS/macroPhor_Array/T16_DEV_genes/EA/wk1/", pattern=".raw",
                            full.names = TRUE)

Start_time_parallel <- Sys.time()
start_time_parallel <- proc.time()

normalization_and_regression_over_all_data_one_timepoint(files_to_loop = files_to_loop)

finish_time_parallel <- proc.time() - start_time_parallel
Finish_time_parallel <- Sys.time()

check_success_high_throughput <- function(files_to_loop, sum_stats_file, grid_type, pattern_to_exclude = "chroma"){
  
  # This line assumes that all chroma standard files contain "chroma" unless pattern_to_exclude option is 
  # modified by user. These will be excluded.
  
  files_to_loop <- files_to_loop[!grepl(pattern_to_exclude, files_to_loop)]
  
  # We will create a dt of all possible grid items for the specified grid type, then compare this to the output
  # so we can see if anything is missing.
  sum_stats_file_in <- fread(sum_stats_file)
  colnames(sum_stats_file_in)[1:2] <- c("file", "grid_item")
  comprehensive_dt <- data.table(file = c(rep(files_to_loop, times = 1, each = grid_type)),
                                 grid_item = c(rep(1:12, times = length(files_to_loop))))
  comprehensive_dt_with_results <- merge(sum_stats_file_in,
                                         comprehensive_dt, by=c("file", "grid_item"),
                                         all = TRUE)
  return(comprehensive_dt_with_results)
}



output <- check_success_high_throughput(files_to_loop = files_to_loop,
                              sum_stats_file = "/scratch2/NSF_GWAS/GMOdetectoR/output/sum_stats/2020-01-09.csv",
                              grid_type = 12)
