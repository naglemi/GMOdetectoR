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
#source("/scratch2/NSF_GWAS/GMOdetectoR/GMOdetectoRv0.34cool.R")

closeAllConnections()
gc()

normalization_and_regression_over_one_plate <- function(files_to_loop, i, chroma_in_backend, job_id, by_explant=TRUE, grid_type=12,
                                                        intercept){
  gc() # Check if this line helps with memory issues, crashing
  filename <- files_to_loop[i]
  
  img_in_backend <- load_image(filename)
  
  pixels_chroma_in <- nrow(chroma_in_backend)
  
  number_rows_to_include_in_chroma_colmeans <- 100
  middle_of_chroma <- round(nrow(chroma_in_backend$spc)/2, digits = 0)
  chroma_colmeans <- colMeans(chroma_in_backend$spc[(middle_of_chroma-(number_rows_to_include_in_chroma_colmeans/2)):(middle_of_chroma+(number_rows_to_include_in_chroma_colmeans/2)),])
  
  scaling_factor <- mean(chroma_colmeans)
  
  img_in_backend$spc <- (img_in_backend$spc/chroma_colmeans)*scaling_factor
  
  if(by_explant==TRUE){
    what_to_plot <- decide_what_to_plot(mode = c(2))
    for(j in 1:grid_type){
      input$grid_position <- j
      crop_and_plot(mode=what_to_plot,
                    image_to_crop = img_in_backend,
                    input_toggle=(2 %in% input$hys_CLS_PCA),
                    grid_item = input$grid_position,
                    image_type = "CLS",
                    first_pass_FP_threshold = 0,
                    first_pass_Chl_threshold = 0,
                    #first_pass_FP_threshold = input$denoise_threshold_FP,
                    #first_pass_Chl_threshold = input$denoise_threshold_Chl,
                    name_to_parse = filename,
                    fluorophore_ID_vector = c("DsRed", "ZsYellow", "Chl", "Diffraction"),
                    record_residuals = TRUE,
                    grid_type = 12,
                    job_id = job_id,
                    intercept = intercept,
                    desired_wavelength_range = c(545, 722))
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
                  fluorophore_ID_vector = c("DsRed", "ZsYellow", "Chl", "Diffraction"),
                  record_residuals = TRUE,
                  grid_type=12,
                  job_id = job_id,
                  intercept = intercept)
  }
}

normalization_and_regression_over_all_data_one_timepoint <- function(files_to_loop, maximum_CPU_number, job_id, by_explant,
                                                                     intercept){
  
  #try(dev.off())
  gc()
  
  chroma_standard_path <- files_to_loop[1]
  chroma_in_backend <<- load_image(image_path = chroma_standard_path)
  
  cl <<- makeCluster(min(maximum_CPU_number, length(files_to_loop)), type = "FORK", outfile="output/cluster.outfile")
  registerDoParallel(cl)
  
  files_to_loop[!grepl("chroma", files_to_loop)]
  #files_to_loop <- files_to_loop[41]
  
  foreach(i=1:length(files_to_loop),
          .packages = c("shiny", "RcppEigen", "hyperSpec", "data.table",
                        "scales", "gridExtra", "tools", "stringr"),
          .export = c("normalization_and_regression_over_one_plate")) %dopar% {
            normalization_and_regression_over_one_plate(files_to_loop = files_to_loop,
                                                        i = i,
                                                        chroma_in_backend = chroma_in_backend,
                                                        job_id,
                                                        by_explant = by_explant,
                                                        intercept = intercept)
          }
  stopCluster(cl)
}

check_success_high_throughput <- function(files_to_loop, sum_stats_file, grid_type, pattern_to_exclude = "chroma", by_explant=TRUE){
  
  # This line assumes that all chroma standard files contain "chroma" unless pattern_to_exclude option is 
  # modified by user. These will be excluded.
  
  files_to_loop <- files_to_loop[!grepl(pattern_to_exclude, files_to_loop)]
  
  # We will create a dt of all possible grid items for the specified grid type, then compare this to the output
  # so we can see if anything is missing.
  sum_stats_file_in <- fread(sum_stats_file)
  
  if(by_explant==TRUE){
    colnames(sum_stats_file_in)[1:2] <- c("file", "grid_item")
    comprehensive_dt <- data.table(file = c(rep(files_to_loop, times = 1, each = grid_type)),
                                   grid_item = c(rep(1:12, times = length(files_to_loop))))
    comprehensive_dt_with_results <- merge(sum_stats_file_in,
                                           comprehensive_dt, by=c("file", "grid_item"),
                                           all.x = FALSE,
                                           all.y = TRUE)
  }else{
    colnames(sum_stats_file_in)[1:2] <- c("file")
    comprehensive_dt <- data.table(file = c(files_to_loop))
    comprehensive_dt_with_results <- merge(sum_stats_file_in,
                                           unique(comprehensive_dt$file), by=c("file"),
                                           all.x = FALSE,
                                           all.y = TRUE)
  }

  return(comprehensive_dt_with_results)
}

run_parallel <- function(timepoint, maximum_CPU_number, intercept){
  files_to_loop <- list.files(timepoint,
                              pattern=".raw",
                              full.names = TRUE)
  
  job_id <- paste0(str_split_fixed(timepoint, "/", 5)[5], "_intercept", intercept)
  print(paste0("Starting workflow for ", job_id))
  
  normalization_and_regression_over_all_data_one_timepoint(files_to_loop = files_to_loop,
                                                           maximum_CPU_number = 10,
                                                           job_id = job_id,
                                                           by_explant = FALSE)
  
  output <- check_success_high_throughput(files_to_loop = files_to_loop,
                                          sum_stats_file = paste0("/scratch2/NSF_GWAS/GMOdetectoR/output/",
                                                                  Sys.Date(),
                                                                  "/",
                                                                  job_id,
                                                                  "/sum_stats.csv"),
                                          grid_type = 12,
                                          by_explant = FALSE)
  
  files_failed_first_time <- unique(c(files_to_loop[1], output[is.na(output$V3),]$file))
  length(files_failed_first_time)
  # 18
  
  if(length(files_failed_first_time)>1){
    stop("Not all jobs completed.")
  }
  print("Complete. On to the next one...")
}

what_to_plot <- decide_what_to_plot(mode = c(2))
tempdir("/scratch2/NSF_GWAS/Rtmp/")

input <- list()
input$hys_CLS_PCA <- 2
input$grid_position <- 2
input$denoise_threshold_FP <- 0
input$denoise_threshold_Chl <- 0

## If we are to run data for a single timepoint...

data_path <- "/scratch2/NSF_GWAS/macroPhor_Array/T16_DEV_genes/EA/wk7/"

files_to_loop <- list.files(data_path,
                            pattern=".raw",
                            full.names = TRUE)
i <- 2
job_id <- "test"
intercept <- 0
# job_id <- str_split_fixed(data_path, "/", 5)[5]
# 
# Start_time_parallel <- Sys.time()
# start_time_parallel <- proc.time()
# 
# normalization_and_regression_over_all_data_one_timepoint(files_to_loop = files_to_loop,
#                                                          maximum_CPU_number = 10,
#                                                          job_id = job_id)
# 
# finish_time_parallel <- proc.time() - start_time_parallel
# Finish_time_parallel <- Sys.time()

# > finish_time_parallel/60
#      user    system   elapsed 
# 32.640650  1.228967 61.551750 
# 
# 
# output <- check_success_high_throughput(files_to_loop = files_to_loop,
#                               sum_stats_file = "/scratch2/NSF_GWAS/GMOdetectoR/output/2020-01-17/T16_DEV_genes/EA/wk7/sum_stats.csv",
#                               grid_type = 12)
# 
# files_failed_first_time <- unique(c(files_to_loop[1], output[is.na(output$V3),]$file))
# length(files_failed_first_time)
# # 18
# 
# if(length(files_failed_first_time)>1){
#   stop("Not all jobs completed.")
# }

## What if we want to run over multiple timepoints back-to-back
# Note: Here we are parallelizing within timepoints, not across. Could squeeze out
# a little more performance in cases where we're left with a number of plates
# not divisible by number of cores if we change this... but not much.

all_timepoints <- list.files("/scratch2/NSF_GWAS/macroPhor_Array/T16_DEV_genes/EA",
                           full.names = TRUE)

for(timepoint in all_timepoints){
  
run_parallel(timepoint = timepoint,
             maximum_CPU_number = 10,
             intercept = 1)
  

}

Sys.time() # Started at 3:47pm

## Back to looking at individual timepoints... Now to make some plots.
## Parse from filename all the info for the plate, by splitting filename, referencing randomization data sheet

for(i in 1:nrow(output)){
  output$ID[i] <- parse_trayplateID(name_being_parsed = output$file[i])
}

library(readxl)
randomization_datasheet <- read_excel("/scratch2/NSF_GWAS/macroPhor_Array/T16_DEV_genes/EA_randomized.xlsx")

randomization_datasheet$ID <- paste0(randomization_datasheet$`Tray ID`,
                                     "_",
                                     randomization_datasheet$`Image #`)

head(randomization_datasheet$ID)
colnames(output)[4:7] <- c("DsRed", "ZsYellow", "Chl", "Diffraction")

combined_data <- merge(output, randomization_datasheet, by="ID")

combined_data$Genotype_ID <- as.factor(combined_data$Genotype_ID)

bp_X_factors(data = combined_data,
             plot_together = TRUE,
             factor1col = "Treatment name",
             factor2col = "Genotype_ID",
             yname = "DsRed",
             factor2colname = "DEV gene",
             y = "DsRed")

model <- lm(combined_data$DsRed ~ combined_data$`Treatment name` + combined_data$Genotype_ID)
summary(model)


normalization_and_regression_over_all_data_one_timepoint(files_to_loop = files_failed_first_time)

# REPEAT

output2 <- check_success_high_throughput(files_to_loop = files_failed_first_time,
                              sum_stats_file = "/scratch2/NSF_GWAS/GMOdetectoR/output/sum_stats/2020-01-09.csv",
                              grid_type = 12)

files_failed_second_time <- unique(c(files_failed_first_time[1], output[is.na(output$V3),]$file))
length(files_failed_second_time)
# 18

## WHY FAILING FOR THESE SAMPLES?
# Try without parallelization
closeAllConnections()
gc()
normalization_and_regression_over_all_data_one_timepoint_single_core <- function(files_to_loop){
  
  #try(dev.off())
  gc()
  
  chroma_standard_path <- files_to_loop[1]
  chroma_in_backend <<- load_image(image_path = chroma_standard_path)
  
  files_to_loop[!grepl("chroma", files_to_loop)]
  #files_to_loop <- files_to_loop[41]
  
  foreach(i=1:length(files_to_loop),
          .packages = c("shiny", "RcppEigen", "hyperSpec", "data.table",
                        "scales", "gridExtra", "tools", "stringr"),
          .export = c("normalization_and_regression_over_one_plate")) %do% {
            normalization_and_regression_over_one_plate(files_to_loop = files_to_loop,
                                                        i = i,
                                                        chroma_in_backend = chroma_in_backend)
            dev.off()
          }
}
normalization_and_regression_over_all_data_one_timepoint_single_core(files_to_loop = files_failed_first_time)

## Did it work without parallelization?
output3 <- check_success_high_throughput(files_to_loop = files_failed_first_time,
                                         sum_stats_file = "/scratch2/NSF_GWAS/GMOdetectoR/output/sum_stats/2020-01-12.csv",
                                         grid_type = 12)
# Seems to have worked but I stopped it early. Finish the rest with parallelization.

files_not_complete_second_time <- unique(c(files_failed_first_time[1], output3[is.na(output3$V3),]$file))
length(files_not_complete_second_time)

normalization_and_regression_over_all_data_one_timepoint(files_to_loop = files_not_complete_second_time)

# Finally all complete this time?
output4 <- check_success_high_throughput(files_to_loop = files_failed_first_time,
                                         sum_stats_file = "/scratch2/NSF_GWAS/GMOdetectoR/output/sum_stats/2020-01-12.csv",
                                         grid_type = 12)

files_not_complete_third_time <- unique(c(files_failed_first_time[1], output4[is.na(output4$V3),]$file))
length(files_not_complete_third_time)

# Working now!
# Run over EB data at wk7
files_to_loop <- list.files("/scratch2/NSF_GWAS/macroPhor_Array/T16_DEV_genes/EB/wk7/", pattern=".raw",
                            full.names = TRUE)

Start_time_parallel <- Sys.time()
start_time_parallel <- proc.time()

normalization_and_regression_over_all_data_one_timepoint(files_to_loop = files_to_loop)

finish_time_parallel <- proc.time() - start_time_parallel
Finish_time_parallel <- Sys.time()
# Finished in 68min

