# PhenotypeAssistant.R
# Requires GMOdetectoR backend, originally version 0.18

# Paths set meant to be run from macroPhor_Array Windows PC, not server

source("C:/Users/OSU/code/GMOdetectoR/GMOdetectoRv0.19.R")

folder_path_all_images <- "G:/Transformation/T12_Different_CIMs_TAG_TAH_TAI/TAI/wk3/"

chroma_standard_path <- "G:/Transformation/T12_Different_CIMs_TAG_TAH_TAI/TAI/wk3/chroma_I0.85_F2.2_L100_092210_0_0_0.raw"

chroma_in <- load_image(image_path = chroma_standard_path)

destination_path <- paste0(getwd(), "/", str_split_fixed(folder_path_all_images, "/", 2)[,2])
if(!dir.exists(destination_path)){
  dir.create(destination_path,
             recursive = TRUE)
}

setwd(destination_path)

phenotype_assistant <- function(image_path, chroma_standard_preloaded, grid_type=20){
  
  img_in_backend <<- load_image(image_path = image_path)
  
  reduced_img_in_backend <- reduce_image(img_in_backend)
  
  plotting_mode <- "hyperspectral"
  
  grid_item <<- " Full Plate"
  
  image_spectrum_table_colored <- extract_plot_grid_item(image_spectrum_table = reduced_img_in_backend,
                                                         chroma_table = chroma_table,
                                                         FC = TRUE,
                                                         cap = TRUE,
                                                         max_intensity_FP = 300,
                                                         max_intensity_Chl = 300,
                                                         scale = TRUE,
                                                         normalize = FALSE,
                                                         to_denoise = TRUE,
                                                         standardize_rgb = TRUE,
                                                         denoise_threshold_FP = 100,
                                                         denoise_threshold_Chl = 100,
                                                         dim_chlorophyll = as.numeric(1),
                                                         FP = input$FP)
  
  
  what_to_plot <- decide_what_to_plot(mode = c(1))
  
  # First let's print the whole plate
  whole_plate_plot <- crop_and_plot(mode=what_to_plot,
                                    input_toggle=TRUE,
                                    # The below line is currently the same for any type of plot
                                    image_spectrum_table = image_spectrum_table_colored,
                                    grid_item = 18,
                                    image_type = "hyperspectral",
                                    first_pass_FP_threshold = 100,
                                    first_pass_Chl_threshold = 100)
  
  ggsave(paste0(tools::file_path_sans_ext(basename(image_path)),".png"), plot = last_plot())
  
  # Now each individual explant
  if(grid_type!=20){
    stop("Support not added yet for grid types other than 20 explant")
  }
  
  what_to_plot <- decide_what_to_plot(mode = c(2))
  
  if(grid_type==20){
    for(i in 1:20){
      this_explant_plot <- crop_and_plot(mode=what_to_plot,
                                         input_toggle=TRUE,
                                         # The below line is currently the same for any type of plot
                                         image_spectrum_table = image_spectrum_table_colored,
                                         grid_item = i,
                                         image_type = "hyperspectral",
                                         first_pass_FP_threshold = 100,
                                         first_pass_Chl_threshold = 100)
      print(this_explant_plot)
      
      ggsave(paste0(FullID, ".png") , plot = last_plot())
    }
  }
  
}

all_image_paths <- list.files("G:/Transformation/T12_Different_CIMs_TAG_TAH_TAI/TAI/wk3/", pattern=".raw", full.names = TRUE)

filename <- all_image_paths[10]
phenotype_assistant(image_path = filename,
                    chroma_standard_preloaded = chroma_in)
