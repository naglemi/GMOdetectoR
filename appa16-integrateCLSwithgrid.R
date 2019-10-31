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

get_row_col_for_grid_item <- function(grid_item, verbose=FALSE){
  if(verbose==TRUE) print(paste0("Grid item: ", grid_item))
  if (grid_item == 20){
    row <- 4
    col <- 2
  }
  if (grid_item == 19){
    row <- 4
    col <- 3
  }
  if (grid_item == 18){
    row <- 4
    col <- 4
  }
  if (grid_item == 17){
    row <- 4
    col <- 5
  }
  if (grid_item == 16){
    row <- 3
    col <- 1
  }
  if (grid_item == 15){
    row <- 3
    col <- 2
  }
  if (grid_item == 14){
    row <- 3
    col <- 3
  }
  if (grid_item == 13){
    row <- 3
    col <- 4
  }
  if (grid_item == 12){
    row <- 3
    col <- 5
  }
  if (grid_item == 11){
    row <- 3
    col <- 6
  }
  if (grid_item == 10){
    row <- 2
    col <- 1
  }
  if (grid_item == 9){
    row <- 2
    col <- 2
  }
  if (grid_item == 8){
    row <- 2
    col <- 3
  }
  if (grid_item == 7){
    row <- 2
    col <- 4
  }
  if (grid_item == 6){
    row <- 2
    col <- 5
  }
  if (grid_item == 5){
    row <- 2
    col <- 6
  }
  if (grid_item == 4){
    row <- 1
    col <- 2
  }
  if (grid_item == 3){
    row <- 1
    col <- 3
  }
  if (grid_item == 2){
    row <- 1
    col <- 4
  }
  if (grid_item == 1){
    row <- 1
    col <- 5
  }
  return(c(row, col))
}

normalize_image <- function(image_spectrum_table, chroma_table, verbose=FALSE){
  # Crop chroma and pre-normalized image to same size
  image_spectrum_table <- image_spectrum_table[image_spectrum_table$rows <= max(chroma_table$rows),]
  #image_spectrum_table[cols < max(chroma_table$cols)]
  if(verbose==TRUE){
    print(dim(image_spectrum_table))
    print(max(image_spectrum_table$rows))
    print(dim(chroma_table))
    print(max(chroma_table$rows))
  }
  chroma_table <- chroma_table[chroma_table$rows <= max(image_spectrum_table$rows),]
  #image_spectrum_table[cols < max(chroma_table$cols)]
  
  image_spectrum_table$r <- image_spectrum_table$r / chroma_table$r
  image_spectrum_table$g <- image_spectrum_table$g / chroma_table$g
  
  return(image_spectrum_table)
}

load_image <- function(image_path, FP){
  image_in <- (read.ENVI(image_path))
  # CAN I PULL THIS OFF
  colnames(image_in$spc) <- image_in@wavelength
  return(image_in)
}

reduce_image <- function(image_in = image_spectrum, FP="DsRed"){
  image_spectrum_loaded <- image_in$spc
  ## Applying false color
  if(FP=="DsRed"){
    FPlambda <- '582.6952'
  }
  if(FP=="GFP"){
    FPlambda <- '508.763'
  }
  
  #colnames(image_spectrum_loaded) <-image_in@wavelength
  image_table <- as.data.table(cbind(image_in$y, image_in$x,
                                     image_spectrum_loaded[,'672.6443'],
                                     image_spectrum_loaded[,which(colnames(image_spectrum_loaded)==FPlambda)],
                                     rep(0, length(image_in$x))))
  
  colnames(image_table) <- c("rows", "cols", "r", "g", "b")
  return(image_table)
}

denoise <- function(image_spectrum_table, threshold_FP, threshold_Chl, verbose=FALSE){
  if(verbose==TRUE){
    print("denoising....")
    print(min(image_spectrum_table$r))
  }
  image_spectrum_table$r[which(image_spectrum_table$r < threshold_Chl)] <- threshold_Chl
  image_spectrum_table$g[which(image_spectrum_table$g < threshold_FP)] <- threshold_FP
  return(image_spectrum_table)
}

assign_ID_index_from_row_column_on_tray <- function(data_to_parse, components_list, mode="table"){
  dictionary <- fread("/scratch2/NSF_GWAS/macroPhor_Array/row_column_key.csv")
  # Get the ID of position in tray in according to row and column
  dictionary$row_column <- paste0(dictionary$row, "_", dictionary$column)
  dictionary[,1:2] <- NULL
  if(mode=="table"){
    # Set colnames for spectral components if multiple are same
    colnames(data_to_parse)[1:length(components_list)] <- components_list
    #colnames(data_in)[1:2] <- c("ChlA", "ChlB")
    data_merged <- merge(data_to_parse, dictionary, by="row_column", all.x = TRUE, all.y = TRUE)
    return(data_merged)
  }
  if(mode=="filename"){
    row <- str_split_fixed(basename(file_path_sans_ext(filename)), "_", 9)[8]
    col <- str_split_fixed(basename(file_path_sans_ext(filename)), "_", 9)[9]
    row_col <- paste0(row, "_", col)
    ID <- dictionary[which(dictionary$row_column == row_col),]$ID
    return(ID)
  }
}

parse_trayplateID <- function(image_path){
  filename <- file_path_sans_ext(basename(image_path))
  trayID <- str_split_fixed(filename, "_", 2)[1]
  plateID <- assign_ID_index_from_row_column_on_tray(data_to_parse = filename, mode="filename")
  trayplateID <- paste0(trayID, "_", plateID)
  return(trayplateID)
}

CLS_and_plot <- function(image_spectrum_table,
                                   chroma_table,
                                   normalize=TRUE, 
                                   to_denoise=TRUE,
                                   denoise_threshold_FP,
                                   denoise_threshold_Chl,
                                   FP,
                                   cropping_option){
  
  ## Denoising
  if(to_denoise==TRUE){
    image_spectrum_table <- denoise(image_spectrum_table, threshold_FP=denoise_threshold_FP, threshold_Chl=denoise_threshold_Chl)
  }
  
  if(cap==TRUE){
    image_spectrum_table$g[which(image_spectrum_table$g > max_intensity_FP)] <- max_intensity_FP
    image_spectrum_table$r[which(image_spectrum_table$r > max_intensity_Chl)] <- max_intensity_Chl
  }
  
  # Must normalize before scaling since normalization depends on cropping scaling depends on first pixel not being cropped out
  if(normalize==TRUE){
    image_spectrum_table <- normalize_image(image_spectrum_table, chroma_table)
  }
  
}

extract_plot_grid_item <- function(image_spectrum_table,
                                   chroma_table,
                                   FC=TRUE,
                                   cap=TRUE,
                                   max_intensity_FP,
                                   max_intensity_Chl,
                                   dim_chlorophyll,
                                   scale=TRUE,
                                   normalize=TRUE, 
                                   to_denoise=TRUE,
                                   denoise_threshold_FP,
                                   denoise_threshold_Chl,
                                   FP,
                                   min_for_rgb_scaling=0,
                                   standardize_rgb=TRUE,
                                   cropping_option){
  
  ## Denoising
  if(to_denoise==TRUE){
    image_spectrum_table <- denoise(image_spectrum_table, threshold_FP=denoise_threshold_FP, threshold_Chl=denoise_threshold_Chl)
  }
  
  if(cap==TRUE){
    image_spectrum_table$g[which(image_spectrum_table$g > max_intensity_FP)] <- max_intensity_FP
    image_spectrum_table$r[which(image_spectrum_table$r > max_intensity_Chl)] <- max_intensity_Chl
  }
  
  # Must normalize before scaling since normalization depends on cropping scaling depends on first pixel not being cropped out
  if(normalize==TRUE){
    image_spectrum_table <- normalize_image(image_spectrum_table, chroma_table)
  }
  
  # do this AFTER CROPPING and BEFORE SCALING so it doesn't screw up later steps or get nullified by prior steps
  if(standardize_rgb==TRUE){
    image_spectrum_table$r[2] <- denoise_threshold_Chl
    image_spectrum_table$g[2] <- denoise_threshold_FP
    image_spectrum_table$r[1] <- max_intensity_Chl
    image_spectrum_table$g[1] <- max_intensity_FP
  }
  
  if(scale==TRUE){
    image_spectrum_table$r <- rescale(image_spectrum_table$r,
                                      from=range(c(denoise_threshold_Chl,max_intensity_Chl)),
                                      to=c(0,255))
    image_spectrum_table$g <- rescale(image_spectrum_table$g,
                                      from=range(c(denoise_threshold_FP,max_intensity_FP)),
                                      to=c(0,255))
  }
  
  if(dim_chlorophyll<1){
    image_spectrum_table$r <- image_spectrum_table$r * dim_chlorophyll
  }
  
  # commenting out this line for GUI because we know img we're looking at and function was split up so
  # no cropping in this function any more
  #FullID <- paste0(parse_trayplateID(filename), "_exp", grid_item)
  
  
  if(FC==TRUE){
    image_spectrum_table$color <- rgb(red= image_spectrum_table$r,
                                      green= image_spectrum_table$g,
                                      blue= image_spectrum_table$b,
                                      maxColorValue = max(image_spectrum_table$g))
  }
  return(image_spectrum_table)
}

crop_and_plot <- function(mode="whole_plate", grid_item, image_spectrum_table = image_spectrum_table_colored, image_type="hyperspectral",
                          first_pass_FP_threshold, first_pass_Chl_threshold, image_to_crop = img_in_backend){
  
    ## Now for cropping
    # First, set vars for just cropping down to the whole grid (crop out edges outside of all grid spaces)
    y_crop_leftside <- 14
    y_crop_rightside <- 1406
    x_crop_leftside <- 1256
    x_crop_rightside <- 226
    
    if(is.na(mode)==TRUE){
      return("Pick a mode")
    }
    
    if(mode=="whole_plate" | mode=="both"){
      crop_bottom <- (x_crop_leftside)
      crop_top <- (x_crop_rightside)
      crop_left <- (y_crop_leftside)
      crop_right <- (y_crop_rightside)
      
      if(image_type=="hyperspectral"){
        image_spectrum_table_all_cropped <- image_spectrum_table[which(image_spectrum_table$cols < crop_bottom & image_spectrum_table$cols > crop_top & image_spectrum_table$rows>crop_left & image_spectrum_table$rows<crop_right),]
      }
      if (image_type=="CLS"){
        image_full_spectrum_cropped <- data.table(cbind(image_to_crop$x, image_to_crop$y, image_to_crop[[]]))
        # WHOA MAYBE I SHOULDNT CHANGE THIS BUT I THINK I SHOULD
        colnames(image_full_spectrum_cropped) <- c("cols", "rows", image_to_crop@wavelength)
        # Need to change from colnames to indices because of data structure.. NOT. Change to data.table earlier instead.
        image_full_spectrum_cropped <- image_full_spectrum_cropped[which(image_full_spectrum_cropped$cols < crop_bottom & image_full_spectrum_cropped$cols > crop_top & image_full_spectrum_cropped$rows>crop_left & image_full_spectrum_cropped$rows<crop_right),]
        
        image_spectrum_table_all_cropped <- CLS_workflow(spectrum_in_to_CLS = image_full_spectrum_cropped,
                                             pass_FP_threshold_from_input = first_pass_FP_threshold,
                                             pass_Chl_threshold_from_input = first_pass_Chl_threshold)
      }
      
      
      p <- ggplot(data=image_spectrum_table_all_cropped, aes(x=cols, y=rows, fill=color)) +
        coord_equal() + geom_tile() + scale_fill_identity() +
        theme(axis.ticks = element_blank(), axis.title = element_blank(),
              panel.background = element_blank(), axis.text = element_blank(),
              panel.border = element_blank(), panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(), axis.line = element_blank())
      #+ggtitle(FullID)
      
    }
    if(mode=="single_explant" | mode=="both"){
      y_newsize <- y_crop_rightside - y_crop_leftside
      x_newsize <- x_crop_leftside - x_crop_rightside
      
      # Calculate size of each grid item
      y_grid_rows <- 4
      x_grid_columns <- 6
      x_grid_item_size <- x_newsize/y_grid_rows
      y_grid_item_size <- y_newsize/x_grid_columns
      
      # Find row and column of desired grid item
      row_col <- get_row_col_for_grid_item(grid_item = grid_item)
      row <- row_col[1]
      col <- row_col[2]
      
      # Crop to desired grid item
      crop_bottom <- (x_crop_leftside - ((row-1)*x_grid_item_size))
      crop_top <- (x_crop_leftside - (row*x_grid_item_size))
      crop_left <- (y_crop_leftside + ((col-1)*y_grid_item_size))
      crop_right <- (y_crop_leftside + (col*y_grid_item_size))
      
      if(image_type=="hyperspectral"){
        image_spectrum_table_single_cropped <- image_spectrum_table[which(image_spectrum_table$cols < crop_bottom & image_spectrum_table$cols > crop_top & image_spectrum_table$rows>crop_left & image_spectrum_table$rows<crop_right),]
      }
      
      if (image_type=="CLS"){
        image_full_spectrum_cropped <- data.table(cbind(image_to_crop$x, image_to_crop$y, image_to_crop[[]]))
        colnames(image_full_spectrum_cropped) <- c("cols", "rows", image_to_crop@wavelength)
        image_full_spectrum_cropped <- image_full_spectrum_cropped[which(image_full_spectrum_cropped$cols < crop_bottom & image_full_spectrum_cropped$cols > crop_top & image_full_spectrum_cropped$rows>crop_left & image_full_spectrum_cropped$rows<crop_right),]
        
        image_spectrum_table_single_cropped <- CLS_workflow(spectrum_in_to_CLS = image_full_spectrum_cropped,
                                                         pass_FP_threshold_from_input = first_pass_FP_threshold,
                                                         pass_Chl_threshold_from_input = first_pass_Chl_threshold)
      }
      
      q <- ggplot(data=image_spectrum_table_single_cropped, aes(x=cols, y=rows, fill=color)) +
        coord_equal() + geom_tile() + scale_fill_identity() +
        theme(axis.ticks = element_blank(), axis.title = element_blank(),
              panel.background = element_blank(), axis.text = element_blank(),
              panel.border = element_blank(), panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(), axis.line = element_blank())
      #+ggtitle(FullID)
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

decide_what_to_plot <- function(mode=NA){
  print(mode)
  if(is.na(mode)==TRUE){
    return("Pick a mode")
  }
  if(1 %in% mode){
    to_do <- "whole_plate"
  }
  if(2 %in% mode){
    to_do <- "single_explant"
  }
  if(length(mode) == 2){
    to_do <- "both"
  }
  print(to_do)
  return(to_do)
}

CLS_workflow <- function(spectrum_in_to_CLS = image_full_spectrum_cropped, pass_FP_threshold_from_input, pass_Chl_threshold_from_input, plot=FALSE){
  
  # Get rid of the rows and column columns (x and y position)
  
  choose_pixels <- function(FP_lambda='582.6952', Chl_lambda='508.763', full_spectrum=this_full_spectrum, CLS_threshold_FP, CLS_threshold_Chl){
    all_pixel_indices <- seq.int(nrow(full_spectrum))
    
    print("Length of all_pixel_indices")
    print(length(all_pixel_indices))
    print(FP_lambda)
    print(Chl_lambda)
    #print(colnames(full_spectrum))
    #print("Head of full spectrum:")
    #print(head(full_spectrum))
    print("Dimensions of full spectrum:")
    print(dim(full_spectrum))
    #print("Head of full spectrum evaluated at FP wavelength:")
    # HAD TO GET RID OF EVAL AND CALL BY PASSING VARIABLE DIRECTLY INTO COLUMN NAME!!!
    # https://stackoverflow.com/questions/19730806/r-access-data-frame-column-using-variable
    # No actually this is a data frame and I had to use get instead of eval
    #print(head(full_spectrum[,get(FP_lambda)]))
    #print("length of full spectrum evaluated at FP wavelength:")
    #print(length(full_spectrum[,get(FP_lambda)]))
    
    print("Are there pixels in the full spectrum that have intensity above the reporter protein threshold?")
    print(head(which(full_spectrum[,get(FP_lambda)] > CLS_threshold_FP))) # This should give a vector of TRUE/FALSE
    print("Max:")
    print(max(full_spectrum[,get(FP_lambda)]))
    if(max(full_spectrum[,get(FP_lambda)])<CLS_threshold_FP){
      print(stop("Need a lower max for FP"))
    }
    if(max(full_spectrum[,get(Chl_lambda)])<CLS_threshold_Chl){
      print(stop("Need a lower max for Chl"))
    }
    # Maybe I can replace the whole use of the all_pixel_indices vector with arr.ind=TRUE option in which
    # Will that fix bug?
    
    #pixels_with_significant_FP <- all_pixel_indices[which(full_spectrum[,get(FP_lambda)] > CLS_threshold_FP)]
    pixels_with_significant_FP <- which(arr.ind = TRUE, x = full_spectrum[,get(FP_lambda)] > CLS_threshold_FP)
    print(paste0("Significant FP is at least ", CLS_threshold_FP))
    print("Number of pixels with significant FP:")
    
    print(length(pixels_with_significant_FP))
    print("Head of pixels with significant FP:")
    print(head(pixels_with_significant_FP))
    pixels_with_significant_Chl <- all_pixel_indices[which(full_spectrum[,get(Chl_lambda)] > CLS_threshold_Chl)]
    pixels_with_something <- unique(c(pixels_with_significant_FP, pixels_with_significant_Chl))
    return(pixels_with_something)
  }
  
  model_and_record_single_pixel <- function(i, Xmatrix = mm, full_spectrum_to_regress_over){
    # This should not take in the image from the backend, but something which is passed through
    # and can be either cropped or uncropped image
    # print(spectrum_in_to_CLS[10:10])
    spectrum_formatted <- as.numeric(as.matrix(unlist(spectrum_in_to_CLS[i,]), nrow=1))

    model_out <- fastLmPure(Xmatrix, spectrum_formatted)

    # Make sure colnames are set properly
    CLS_table[i, DsRed := model_out$coefficients[2]]
    CLS_table[i, ZsYellow := model_out$coefficients[3]]
    CLS_table[i, ChlA := model_out$coefficients[4]]
    CLS_table[i, ChlB := model_out$coefficients[5]]
  }
  
  perform_CLS <- function(pixels_to_test){
    options(warn=-1)
    pb <- txtProgressBar(min = 0, max = length(pixels_to_test), initial = 0)
    print("About to start CLS over all pixels passing threshold. Here is table:")
    print(head(CLS_table))
    #CLS_table_out <- foreach(i=pixels_to_test) %do% {
    print(paste0("Number of pixels to test: ", length(pixels_to_test)))
    print("First 5 pixels to test: ")
    print(head(pixels_to_test))
    # Need to get rid of the row and col columns I used for indexing earlier
    spectrum_in_to_CLS$rows <<- NULL
    spectrum_in_to_CLS$cols <<- NULL
    for(i in pixels_to_test){
      model_and_record_single_pixel(i)
      setTxtProgressBar(pb,i)
    }
    options(warn=0)
    print("Done. Here is table again:")
    print(head(CLS_table))
    
  }
  
  finish_processing_CLS_table <- function(CLS_table, chosen_pixels_in = chosen_pixels){
    # THIS IS SLOW
    # Faster way: https://stackoverflow.com/questions/38226323/replace-all-values-in-a-data-table-given-a-condition?rq=1
    print(Sys.time())
    CLS_table[is.na(CLS_table)] <- 0
    print(Sys.time())
    #CLS_table[!chosen_pixels_in, DsRed := 0]
    #CLS_table[!chosen_pixels_in, ZsYellow := 0]
    #CLS_table[!chosen_pixels_in, ChlA := 0]
    #CLS_table[!chosen_pixels_in, ChlB := 0]
    
    CLS_table$ChlB <- rescale(CLS_table$ChlB,
                              #from=range(c(denoise_threshold_Chl,max_intensity_Chl)),
                              to=c(0,255))
    CLS_table$DsRed <- rescale(CLS_table$DsRed,
                               #from=range(c(denoise_threshold_FP,max_intensity_FP)),
                               to=c(0,255))
    CLS_table$color <- rgb(red= CLS_table$ChlB,
                           green= CLS_table$DsRed,
                           blue= rep(0, nrow(CLS_table)),
                           maxColorValue = 255)
    print(dim(CLS_table))
    print(head(CLS_table))
    return(CLS_table)
  }
  
  # fitting a smooth curve https://stackoverflow.com/questions/3480388/how-to-fit-a-smooth-curve-to-my-data-in-r
  read_plot_spectra <- function(spectra_path, plot=TRUE){
    pub_spectra_in <- fread(spectra_path)
    pub_emission_spectrum <- data.frame(wavelength=pub_spectra_in$`emission wavelength (nm)`,
                                        intensity=pub_spectra_in$`Normalized emission`)
    fit <- loess(intensity~wavelength, data=pub_emission_spectrum, span=0.1)
    
    if(plot==TRUE){
      plot(pub_emission_spectrum,
           main=(paste0(basename(tools::file_path_sans_ext(spectra_path)), " (Published)")),
           xlab="Wavelength",
           ylab="Normalized emission")
      
      lines(pub_emission_spectrum$wavelength,
            predict(fit,pub_emission_spectrum$wavelength), col='red', lwd=2)
    }
    
    wavelengths <- colnames(spectrum_in_to_CLS)[3:ncol(spectrum_in_to_CLS)]
    predict_from <- data.frame(wavelength=wavelengths)
    predictions <- predict(fit, predict_from)
    scaled_emission_spectra <- cbind(wavelengths, predictions)
    
    if(plot==TRUE){
      plot(scaled_emission_spectra,
           main=(paste0(basename(tools::file_path_sans_ext(spectra_path)), " (Fitted for macroPhor Array camera)")),
           xlab="Wavelength",
           ylab="Normalized emission")
      
      lines(pub_emission_spectrum$wavelength,
            predict(fit,pub_emission_spectrum$wavelength), col='red', lwd=2)
    }
    
    #return(predictions)
    fitted_spectra_dataframe <- as.data.frame(cbind(wavelengths, predictions))
    colnames(fitted_spectra_dataframe) <- c("Wavelength", "Normalized intensity (fitted)")
    
    return(fitted_spectra_dataframe)
  }
  ChlB <- read_plot_spectra("/scratch2/NSF_GWAS/macroPhor_Array/Fluorophore_spectra/ChlB.csv")
  ChlA <- read_plot_spectra("/scratch2/NSF_GWAS/macroPhor_Array/Fluorophore_spectra/ChlA.csv")
  DsRed <- read_plot_spectra("/scratch2/NSF_GWAS/macroPhor_Array/Fluorophore_spectra/DsRed_Clontech.csv")
  ZsYellow <- read_plot_spectra("/scratch2/NSF_GWAS/macroPhor_Array/Fluorophore_spectra/ZsYellow_cutoff.csv")
  
  convert_NA_to_zero_in_spectrum_and_return_vector <- function(spectrum_in){
    print("KOONICHIWA")
    vector <- as.vector(spectrum_in$'Normalized intensity (fitted)')
    print(vector[1:5])
    vector[which(is.na(vector)==TRUE)] <- 0
    print(vector[1:5])
    return(as.numeric(as.character(vector)))
    
  }
  ChlA_vector <- convert_NA_to_zero_in_spectrum_and_return_vector(ChlA)
  ChlB_vector <- convert_NA_to_zero_in_spectrum_and_return_vector(ChlB)
  DsRed_vector <- convert_NA_to_zero_in_spectrum_and_return_vector(DsRed)
  ZsYellow_vector <- convert_NA_to_zero_in_spectrum_and_return_vector(ZsYellow)
  print(ChlA_vector[1:5])
  
  mm <<- as.matrix(cbind(1, DsRed_vector, ZsYellow_vector, ChlA_vector, ChlB_vector))
  print("Show me all of mm:")
  print(mm)
  
  
  
  print("Making CLS table")
  CLS_table_server_side <- data.table(cbind(spectrum_in_to_CLS$rows, spectrum_in_to_CLS$cols, matrix(nrow=dim(spectrum_in_to_CLS)[1], ncol=4)))
  assign("CLS_table", CLS_table_server_side, envir=globalenv())
  
  colnames(CLS_table) <- c("rows", "cols", "DsRed", "ZsYellow", "ChlA", "ChlB")
  print("About to perform CLS")
  print("Pixels we will perform CLS over will first be determined.")
  print(paste0("Threshold for FP: ", pass_FP_threshold_from_input))
  print(paste0("Treshold for Chl: ", pass_Chl_threshold_from_input))
  chosen_pixels <<- choose_pixels(full_spectrum=spectrum_in_to_CLS[,3:ncol(spectrum_in_to_CLS)],
                                 CLS_threshold_FP=pass_FP_threshold_from_input,
                                 CLS_threshold_Chl=pass_Chl_threshold_from_input)
  print(paste0("Number of chosen pixels: ", length(chosen_pixels)))
  print("First 5: ")
  print(head(chosen_pixels))
  perform_CLS(pixels_to_test = chosen_pixels)
  print("CLS done for whole image. Finish data processing.")
  CLS_table <- finish_processing_CLS_table(CLS_table)
  print("CLS data processing finished.")
  return(CLS_table)
}

# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("GMOdetectoR"),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      position=c("left"),
      sidebarPanel(
        # Copy the line below to make a file upload manager
        fileInput("file", label = h3("Chroma standard")),
        
        fileInput("file", label = h3("Sample image")),
        
        #hr(),
        fluidRow(column(4, verbatimTextOutput("value"))),
        
        # Copy the line below to make a set of radio buttons
        radioButtons("radio", label = h3("Reporter protein"),
                     choices = list("DsRed" = 1, "ZsYellow" = 2, "GFP" = 3), 
                     selected = 1),
        
        #hr(),
        
        checkboxGroupInput("checkGroup", label = h3("View"), 
                           choices = list("Whole plate" = 1, "Single explant" = 2),
                           selected = c(1, 2)),
        
        numericInput("grid_position", label = h3("Grid position"), value = 1),
        
        
         #sliderInput("bins",
         #            "Number of bins:",
         #            min = 1,
         #            max = 50,
        #             value = 30),
         sliderInput("denoise_threshold_Chl",
                     "Denoising threshold for Chlorophyll:",
                     min = 0,
                     max = 200,
                     value = 100),
         sliderInput("denoise_threshold_FP",
                     "Denoising threshold for reporter protein:",
                     min = 0,
                     max = 200,
                     value = 120),
         sliderInput("max_intensity_Chl",
                     "Maximum intensity for Chlorophyll:",
                     min = 1,
                     max = 1000,
                     value = 200),
         sliderInput("max_intensity_FP",
                     "Maximum intensity for reporter protein:",
                     min = 1,
                     max = 1000,
                     value = 300),
         sliderInput("dim_chlorophyll",
                     "Chlorophyll signal",
                     min = 0,
                     max = 1,
                     value = 1),
        actionButton("runCLS", label = "CLS"),
        actionButton("runPCA", label = "PCA"),
        textInput("text",
                  #label = h3("Text input"),
                  label = NULL,
                  value = "Number of components")
      ),
      
      
      # Show a plot of the generated distribution
      mainPanel(
         #plotOutput("distPlot"),
         plotOutput("image")
         #plotOutput("single_explant_image")
      )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
   
   output$distPlot <- renderPlot({
      # generate bins based on input$bins from ui.R
      x    <- faithful[, 2] 
      bins <- seq(min(x), max(x), length.out = input$bins + 1)
      
      # draw the histogram with the specified number of bins
      hist(x, breaks = bins, col = 'darkgray', border = 'white')
   })
   ###################################################################################
   
   img_in_backend <<- load_image(image_path = paste0("/scratch2/NSF_GWAS/macroPhor_Array/CT_CU_CV_raw/wk6/",
                                                    filename))
   
   reduced_img_in_backend <- reduce_image(img_in_backend)
   
   plotting_mode <- "hyperspectral"
   
   observeEvent(input$runCLS, {
     session$sendCustomMessage(type = 'testmessage',
                               message = 'Thank you for clicking')
     plotting_mode <<- "CLS"
     print("Event observed!!!")
   })
   
   observeEvent(input$runPCA, {
     session$sendCustomMessage(type = 'testmessage',
                               message = 'Thank you for clicking')
   })
   
   filename <- "CV3_F1.9_I5.0_L100_cyan_001552_9_1_4.raw"
   
   image_spectrum_table_colored <- reactive({
     image_spectrum_table_colored <- extract_plot_grid_item(image_spectrum_table = reduced_img_in_backend,
                                                            chroma_table = chroma_table,
                                                            FC = TRUE,
                                                            cap = TRUE,
                                                            max_intensity_FP = input$max_intensity_FP,
                                                            max_intensity_Chl = input$max_intensity_Chl,
                                                            scale = TRUE,
                                                            normalize = FALSE,
                                                            to_denoise = TRUE,
                                                            standardize_rgb = TRUE,
                                                            denoise_threshold_FP = input$denoise_threshold_FP,
                                                            denoise_threshold_Chl = input$denoise_threshold_Chl,
                                                            dim_chlorophyll = as.numeric(input$dim_chlorophyll),
                                                            FP = input$FP)
   })
   
   what_to_plot <- reactive({
     what_to_plot <- decide_what_to_plot(mode = input$checkGroup)
   })
   
   output$image <- renderPlot({
     # parenthesis are needed to avoid error because our table object is really kind of a function
     # https://stackoverflow.com/questions/40623749/what-is-object-of-type-closure-is-not-subsettable-error-in-shiny/40623750
     crop_and_plot(mode=what_to_plot(),
                   # The below line is currently the same for any type of plot
                   image_spectrum_table = image_spectrum_table_colored(),
                   grid_item = input$grid_position,
                   image_type = plotting_mode,
                   first_pass_FP_threshold = input$denoise_threshold_FP,
                   first_pass_Chl_threshold = input$denoise_threshold_Chl)
   }, height = 1200, width=1200)
   
   output$CLSplot <- renderPlot({
     crop_and_plot(mode=what_to_plot(),
                   image_spectrum_table = )
   })

   ###################################################################################
   
   output$FP <- renderPrint({ input$radio })
   
   output$normalize_now <- renderPrint({ input$action })
   
   output$griditem <- renderPrint({ input$num })
   
}

# Run the application 
shinyApp(ui = ui, server = server)

