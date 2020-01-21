filename <- "CV3_F1.9_I5.0_L100_cyan_001552_9_1_4.raw"

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
library(Brobdingnag)
tempdir("/scratch2/NSF_GWAS/Rtmp/")
setwd("/scratch2/NSF_GWAS/GMOdetectoR/")

source("/scratch2/NSF_GWAS/GMOdetectoR/GMOdetectoRv0.32cool.R")

# Define server logic required to draw a histogram
# server <- function(input, output, session) {
img_in_backend <<- load_image(image_path = paste0("/scratch2/NSF_GWAS/macroPhor_Array/CT_CU_CV_raw/wk6/",
                                                  filename))

reduced_img_in_backend <- reduce_image(img_in_backend)

plotting_mode <- "hyperspectral"


image_spectrum_table_colored <- extract_plot_grid_item(image_spectrum_table = reduced_img_in_backend,
                                                       chroma_table = chroma_table,
                                                       FC = TRUE,
                                                       cap = TRUE,
                                                       max_intensity_FP = 800,
                                                       max_intensity_Chl = 800,
                                                       scale = TRUE,
                                                       normalize = FALSE,
                                                       to_denoise = TRUE,
                                                       standardize_rgb = TRUE,
                                                       denoise_threshold_FP = 400,
                                                       denoise_threshold_Chl = 400,
                                                       dim_chlorophyll = as.numeric(1),
                                                       FP = input$FP,
                                                       filename = filename,
                                                       grid_item = 18)


#what_to_plot <- decide_what_to_plot(mode = c(1, 2))
what_to_plot <- decide_what_to_plot(mode = c(1))
tempdir("/scratch2/NSF_GWAS/Rtmp/")

input <- list()
input$hys_CLS_PCA <- 2

crop_and_plot(mode=what_to_plot,
              input_toggle=TRUE,
              # The below line is currently the same for any type of plot
              image_spectrum_table = image_spectrum_table_colored,
              grid_item = 18,
              image_type = "hyperspectral",
              first_pass_FP_threshold = 200,
              first_pass_Chl_threshold = 200,
              name_to_parse = filename,
              grid_type = 20)

input$grid_position <- 18
input$denoise_threshold_FP <- 0
input$denoise_threshold_Chl <- 0


#source("/scratch2/NSF_GWAS/GMOdetectoR/GMOdetectoRv0.29cool.R")

#crop_and_plot(mode=what_to_plot(),
#Rprof(filename = "Rprof_a3.out")
crop_and_plot(image_type = "CLS",
              mode=what_to_plot,
              input_toggle=(2 %in% input$hys_CLS_PCA),
              # The below line is currently the same for any type of plot
              #image_spectrum_table = image_spectrum_table_colored(),
              image_spectrum_table = img_in_backend,
              grid_item = input$grid_position,
              first_pass_FP_threshold = input$denoise_threshold_FP,
              first_pass_Chl_threshold = input$denoise_threshold_Chl,
              name_to_parse = filename,
              fluorophore_ID_vector = c("DsRed", "ZsYellow", "Chl", "Diffraction"),
              record_residuals = TRUE,
              grid_type = 20)
#Rprof(filename = NULL)

files_to_loop <- list.files("/scratch2/NSF_GWAS/macroPhor_Array/CT_CU_CV_raw/wk6/", pattern=".raw")[2:11]

cl <- makeCluster(2)
registerDoParallel(cl)

start_time_single_core <- proc.time()
foreach(i=1:2) %do% {
  
  img_in_backend <<- load_image(image_path = paste0("/scratch2/NSF_GWAS/macroPhor_Array/CT_CU_CV_raw/wk6/",
                                                    files_to_loop[i]))
  
  reduced_img_in_backend <- reduce_image(img_in_backend)
  
  image_spectrum_table_colored <- extract_plot_grid_item(image_spectrum_table = reduced_img_in_backend,
                                                         chroma_table = chroma_table,
                                                         FC = TRUE,
                                                         cap = TRUE,
                                                         max_intensity_FP = 800,
                                                         max_intensity_Chl = 800,
                                                         scale = TRUE,
                                                         normalize = FALSE,
                                                         to_denoise = TRUE,
                                                         standardize_rgb = TRUE,
                                                         denoise_threshold_FP = 400,
                                                         denoise_threshold_Chl = 400,
                                                         dim_chlorophyll = as.numeric(1),
                                                         FP = input$FP,
                                                         filename = filename,
                                                         grid_item = 18)
  
  crop_and_plot(mode=what_to_plot,
                input_toggle=(2 %in% input$hys_CLS_PCA),
                # The below line is currently the same for any type of plot
                #image_spectrum_table = image_spectrum_table_colored(),
                image_spectrum_table = image_spectrum_table_colored,
                grid_item = input$grid_position,
                image_type = "CLS",
                first_pass_FP_threshold = input$denoise_threshold_FP,
                first_pass_Chl_threshold = input$denoise_threshold_Chl,
                name_to_parse = filename,
                fluorophore_ID_vector = c("DsRed", "ZsYellow", "ChlA", "ChlB"),
                record_residuals = TRUE)
  
}
finish_time_single_core <- proc.time() - start_time_single_core

start_time_parallel <- proc.time()
foreach(i=1:2) %dopar% {
  
  img_in_backend <<- load_image(image_path = paste0("/scratch2/NSF_GWAS/macroPhor_Array/CT_CU_CV_raw/wk6/",
                                                    files_to_loop[i]))
  
  reduced_img_in_backend <- reduce_image(img_in_backend)
  
  image_spectrum_table_colored <- extract_plot_grid_item(image_spectrum_table = reduced_img_in_backend,
                                                         chroma_table = chroma_table,
                                                         FC = TRUE,
                                                         cap = TRUE,
                                                         max_intensity_FP = 800,
                                                         max_intensity_Chl = 800,
                                                         scale = TRUE,
                                                         normalize = FALSE,
                                                         to_denoise = TRUE,
                                                         standardize_rgb = TRUE,
                                                         denoise_threshold_FP = 400,
                                                         denoise_threshold_Chl = 400,
                                                         dim_chlorophyll = as.numeric(1),
                                                         FP = input$FP,
                                                         filename = filename,
                                                         grid_item = 18)
  
  crop_and_plot(mode=what_to_plot,
                input_toggle=(2 %in% input$hys_CLS_PCA),
                # The below line is currently the same for any type of plot
                #image_spectrum_table = image_spectrum_table_colored(),
                image_spectrum_table = image_spectrum_table_colored,
                grid_item = input$grid_position,
                image_type = "CLS",
                first_pass_FP_threshold = input$denoise_threshold_FP,
                first_pass_Chl_threshold = input$denoise_threshold_Chl,
                name_to_parse = filename,
                fluorophore_ID_vector = c("DsRed", "ZsYellow", "ChlA", "ChlB"),
                record_residuals = TRUE)
  
}
finish_time_parallel <- proc.time() - start_time_single_core
