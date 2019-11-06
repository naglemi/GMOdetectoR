# PhenotypeAssistant.R
# Requires GMOdetectoR backend, originally version 0.18

# Paths set meant to be run from macroPhor_Array Windows PC, not server

source("C:/Users/OSU/code/GMOdetectoR/GMOdetectoRv0.19.R")

phenotype_assistant_loop <- function(folder_path_all_images,
                                     chroma_standard_path,
                                     grid_type,
                                     Chl_cap,
                                     FP_cap,
                                     denoise_thr_Chl,
                                     denoise_thr_FP){
  
  chroma_in <- reduce_image(load_image(image_path = chroma_standard_path))
  
  destination_path <- paste0("G:/PhenotypeAssistant/", str_split_fixed(folder_path_all_images, "/", 2)[,2])
  if(!dir.exists(destination_path)){
    dir.create(destination_path,
               recursive = TRUE)
  }
  
  setwd(destination_path)
  
  phenotype_assistant <- function(image_path, chroma_standard_preloaded, grid_type=20, Chl_cap, FP_cap, denoise_thr_Chl, denoise_thr_FP){
    
    img_in_backend <<- load_image(image_path = image_path)
    
    reduced_img_in_backend <- reduce_image(img_in_backend)
    
    plotting_mode <- "hyperspectral"
    
    grid_item <<- " Full Plate"
    
    image_spectrum_table_colored <- extract_plot_grid_item(image_spectrum_table = reduced_img_in_backend,
                                                           chroma_table = chroma_standard_preloaded,
                                                           FC = TRUE,
                                                           cap = TRUE,
                                                           max_intensity_FP = FP_cap,
                                                           max_intensity_Chl = Chl_cap,
                                                           scale = TRUE,
                                                           normalize = TRUE,
                                                           to_denoise = TRUE,
                                                           standardize_rgb = TRUE,
                                                           denoise_threshold_FP = denoise_thr_FP,
                                                           denoise_threshold_Chl = denoise_thr_Chl,
                                                           dim_chlorophyll = as.numeric(1),
                                                           FP = input$FP,
                                                           filename = image_path)
    
    
    what_to_plot <- decide_what_to_plot(mode = c(1))
    
    # First let's print the whole plate
    whole_plate_plot <- crop_and_plot(mode=what_to_plot,
                                      input_toggle=TRUE,
                                      # The below line is currently the same for any type of plot
                                      image_spectrum_table = image_spectrum_table_colored,
                                      grid_item = 18,
                                      image_type = "hyperspectral",
                                      first_pass_FP_threshold = 100,
                                      first_pass_Chl_threshold = 100,
                                      name_to_parse = image_path)
    
    # The variable dirname will refer both to the actual directory name (containing images
    # for individual explants) and the name of the whole-plate file
    
    dirname <- paste0(str_split_fixed(FullID, "_", 3)[,1], "_plate", str_split_fixed(FullID, "_", 3)[,2])
    if(!dir.exists(dirname)){
      dir.create(dirname)
    }
    
    ggsave(paste0(tools::file_path_sans_ext(basename(dirname)),".png"), plot = last_plot())
    
    # Now each individual explant
    if(grid_type!=20){
      stop("Support not added yet for grid types other than 20 explant")
    }
    
    what_to_plot <- decide_what_to_plot(mode = c(2))
    
    if(grid_type==20){
      for(j in 1:20){
        print(paste0("On explant ", j))
        
        this_explant_plot <- crop_and_plot(mode=what_to_plot,
                                           input_toggle=TRUE,
                                           # The below line is currently the same for any type of plot
                                           image_spectrum_table = image_spectrum_table_colored,
                                           grid_item = j,
                                           image_type = "hyperspectral",
                                           first_pass_FP_threshold = denoise_thr_FP,
                                           first_pass_Chl_threshold = denoise_thr_Chl,
                                           name_to_parse = image_path)
        print(this_explant_plot)
        
        ggsave(paste0(dirname, "/", FullID, ".png") , plot = last_plot())
      }
    }
    
  }
  
  all_image_paths <- list.files(folder_path_all_images, pattern=".raw", full.names = TRUE)
  
  for(i in 1:length(all_image_paths)){
    print(paste0("Running workflow for image #", i, " of ", length(all_image_paths), ": ", all_image_paths[i]))
    phenotype_assistant(image_path = all_image_paths[i],
                        chroma_standard_preloaded = chroma_in,
                        FP_cap = FP_cap,
                        Chl_cap = Chl_cap,
                        denoise_thr_Chl = denoise_thr_Chl,
                        denoise_thr_FP = denoise_thr_FP)
  }
  
}

phenotype_assistant_loop(folder_path_all_images = "G:/Transformation/T12_Different_CIMs_TAG_TAH_TAI/TAI/wk3/",
                         chroma_standard_path = "G:/Transformation/T12_Different_CIMs_TAG_TAH_TAI/TAI/wk3/chroma_I0.85_F2.2_L100_092210_0_0_0.raw",
                         FP_cap = 0.15,
                         Chl_cap = 2.5,
                         denoise_thr_Chl = 0.5,
                         denoise_thr_FP = 0.04)

# Started at 4:15pm
Sys.time()

#print("Image processing complete! This window can now be closed.")

#Sys.sleep(60*24*3)