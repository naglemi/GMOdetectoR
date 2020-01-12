# GMOdetectoR complete backend v0.18 (a18)

library(RcppEigen)
library(hyperSpec)
library(data.table)
library(scales)
library(gridExtra)
library(tools)
library(stringr)
library(foreach)
library(doParallel)
library(wesanderson)
library(Brobdingnag)

get_row_col_for_grid_item <- function(grid_item, verbose=FALSE, grid_type){
  if(verbose==TRUE) print(paste0("Grid item: ", grid_item))
  if(grid_type==20){
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
  }
  if(grid_type==12){
    if (grid_item == 12){
      row <- 1
      col <- 4
    }
    if (grid_item == 11){
      row <- 2
      col <- 4
    }
    if (grid_item == 10){
      row <- 3
      col <- 4
    }
    if (grid_item == 9){
      row <- 4
      col <- 4
    }
    if (grid_item == 8){
      row <- 2
      col <- 1
    }
    if (grid_item == 7){
      row <- 2
      col <- 2
    }
    if (grid_item == 6){
      row <- 2
      col <- 3
    }
    if (grid_item == 5){
      row <- 2
      col <- 4
    }
    if (grid_item == 4){
      row <- 1
      col <- 1
    }
    if (grid_item == 3){
      row <- 1
      col <- 2
    }
    if (grid_item == 2){
      row <- 1
      col <- 3
    }
    if (grid_item == 1){
      row <- 1
      col <- 4
    }
  }
  return(c(row, col))
}

normalize_image <- function(image_spectrum_table, chroma_table, verbose=FALSE){
  # Crop chroma and pre-normalized image to same size
  
  image_spectrum_table <- image_spectrum_table[image_spectrum_table$rows <= max(chroma_table$rows),]
  #image_spectrum_table[cols < max(chroma_table$cols)]
  if(verbose==TRUE){
    print("Head of r and g channels for image and chroma:")
    print(head(image_spectrum_table$r))
    print(head(image_spectrum_table$g))
    print(head(chroma_table$r))
    print(head(chroma_table$g))
    print("Dimensions and maximum row index for image and chroma:")
    print(dim(image_spectrum_table))
    print(max(image_spectrum_table$rows))
    print(dim(chroma_table))
    print(max(chroma_table$rows))
  }
  chroma_table <- chroma_table[chroma_table$rows <= max(image_spectrum_table$rows),]
  #image_spectrum_table[cols < max(chroma_table$cols)]
  
  # Added re-scaling to mean of original to keep denoise thresholds working as they do without normalization in 
  # v0.19
  
  image_spectrum_table$r <- as.numeric(as.character(image_spectrum_table$r))
  image_spectrum_table$g <- as.numeric(as.character(image_spectrum_table$g))
  
  min_image_spectrum_table_r_before_normalizing <- min(image_spectrum_table$r)
  min_image_spectrum_table_g_before_normalizing <- min(image_spectrum_table$g)
  max_image_spectrum_table_r_before_normalizing <- max(image_spectrum_table$r)
  max_image_spectrum_table_g_before_normalizing <- max(image_spectrum_table$g)
  
  
  
  # Think I might be losing detail for green channel because chroma standard not on same scale
  # Due to limits of precision?
  #chroma_table$r <- rescale(chroma_table$r,
  #                                  to=c(min_image_spectrum_table_r_before_normalizing,
  #                                       max_image_spectrum_table_r_before_normalizing))
  #chroma_table$g <- rescale(chroma_table$g,
  #                                  to=c(min_image_spectrum_table_g_before_normalizing,
  #                                      max_image_spectrum_table_g_before_normalizing))
  

  
  image_spectrum_table$r <- image_spectrum_table$r / chroma_table$r
  image_spectrum_table$g <- image_spectrum_table$g / chroma_table$g
  
  if(verbose==TRUE){
    print(paste0("Max for R channel: ", max(image_spectrum_table$r)))
    print(paste0("Max for G channel: ", max(image_spectrum_table$g)))
    print(paste0("Min for R channel: ", min(image_spectrum_table$r)))
    print(paste0("Min for G channel: ", min(image_spectrum_table$g)))
    
    print(paste0("Max for R channel after normalizing: ", max(image_spectrum_table$r)))
    print(paste0("Max for G channel after normalizing: ", max(image_spectrum_table$g)))
    
    #image_spectrum_table$r <- rescale(image_spectrum_table$r,
    #                                  to=c(min_image_spectrum_table_r_before_normalizing,
    #                                       max_image_spectrum_table_r_before_normalizing))
    #image_spectrum_table$g <- rescale(image_spectrum_table$g,
    #                                  to=c(min_image_spectrum_table_g_before_normalizing,
    #                                       max_image_spectrum_table_g_before_normalizing))
    
    print(paste0("Max for R channel after rescaling: ", max(image_spectrum_table$r)))
    print(paste0("Max for G channel after rescaling: ", max(image_spectrum_table$g)))
    print(paste0("Min for R channel after rescaling: ", min(image_spectrum_table$r)))
    print(paste0("Min for G channel after rescaling: ", min(image_spectrum_table$g)))
    
    print(head(image_spectrum_table$r))
    print(head(image_spectrum_table$g))
  }
  

  
  return(image_spectrum_table)
}

load_image <- function(image_path, FP){
  image_in <- read.ENVI(image_path,
              headerfile=paste0(tools::file_path_sans_ext(image_path), ".hdr"))
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

# v0.19 fixed redundancy in how both filename and data_to_parse refer to same thing. Make both "data_to_parse" within this scope
assign_ID_index_from_row_column_on_tray <- function(data_to_parse = filename, components_list, mode="table", verbose=FALSE){
  #dictionary <- fread("/scratch2/NSF_GWAS/macroPhor_Array/row_column_key.csv")
  # Hardcoding dictionary in v0.19
  dictionary <- cbind(c(0,0,0,0,0,0,0,
                        1,1,1,1,1,1,1,
                        2,2,2,2,2,2,2),
                      c(0,1,2,3,4,5,6,
                        0,1,2,3,4,5,6,
                        0,1,2,3,4,5,6),
                      c(1:21))
  dictionary <- as.data.table(dictionary)
  colnames(dictionary) <- c("row", "column", "ID")
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
    
    # Patch added in v0.19 for compatibility regardless of whether "_cyan" is at end of filename
    if(grepl("cyan", data_to_parse)==1){
      ndelimiters=9
    }else{
      ndelimiters=8
    }
    
    
    row <- str_split_fixed(basename(file_path_sans_ext(data_to_parse)), "_", 9)[ndelimiters-1]
    # Changed in v0.19 along with patch above
    col <- str_split_fixed(basename(file_path_sans_ext(data_to_parse)), "_", ndelimiters)[ndelimiters]
    row_col <- paste0(row, "_", col)
    ID <- dictionary[which(dictionary$row_column == row_col),]$ID
    # Debugging lines added in v0.19
    if(verbose==TRUE){
      print(paste0("This row is ", row))
      print(paste0("This col is ", col))
      print(paste0("This row_col is ", row_col))
      print(paste0("This filename (stripped) is ", basename(file_path_sans_ext(data_to_parse))))
      print(paste0("This ID about to be returned from assign_ID_index_from_roW_column_on_tray is ", ID))
    }

    return(ID)
  }
}

parse_trayplateID <- function(name_being_parsed){
  pass_to_dodge_error <- name_being_parsed
  imgpath_stripped <- file_path_sans_ext(basename(pass_to_dodge_error))
  trayID <- str_split_fixed(imgpath_stripped, "_", 2)[1]
  plateID <- assign_ID_index_from_row_column_on_tray(data_to_parse = imgpath_stripped, mode="filename")
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
    stop("Still need to modify function for FULL normalization of all wavelengths, not just those used in false color")
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
                                   cropping_option,
                                   #filename=filename,
                                   filename,
                                   verbose=FALSE,
                                   grid_item){
  if(verbose==TRUE){
    print(paste0("Filename is ", filename))
  }
  
  # Must normalize before scaling since normalization depends on cropping scaling depends on first pixel not being cropped out
  # just moved this to top in v0.19
  if(normalize==TRUE){
    image_spectrum_table <- normalize_image(image_spectrum_table, chroma_table)
  }
  
  ## Denoising
  if(to_denoise==TRUE){
    image_spectrum_table <- denoise(image_spectrum_table, threshold_FP=denoise_threshold_FP, threshold_Chl=denoise_threshold_Chl)
  }
  
  if(cap==TRUE){
    image_spectrum_table$g[which(image_spectrum_table$g > max_intensity_FP)] <- max_intensity_FP
    image_spectrum_table$r[which(image_spectrum_table$r > max_intensity_Chl)] <- max_intensity_Chl
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
  # Need to uncomment in v0.19 for PhenotypeAssistant
  # Uncommenting didn't fix object'FullID' not found error.
  # I AM CURRENTLY DEFINING THIS BOTH IN FUNCTIONS extract_plot_grid_item and crop_and_plot
  # IF DEFINED IN THE FORMER, ONLY APPEARS FOR WHOLE PLATE
  FullID <<- paste0(parse_trayplateID(name_being_parsed = filename), "_exp", grid_item)
  # Debugging statement for v0.19
  
  if(verbose==TRUE){
    print(paste0("Full ID is ", FullID))
  }

  
  if(FC==TRUE){
    image_spectrum_table$color <- rgb(red= image_spectrum_table$r,
                                      green= image_spectrum_table$g,
                                      blue= image_spectrum_table$b,
                                      maxColorValue = max(image_spectrum_table$g))
  }
  return(image_spectrum_table)
}

crop_and_plot <- function(mode="whole_plate", grid_item, image_spectrum_table = image_spectrum_table_colored, image_type="hyperspectral",
                          first_pass_FP_threshold, first_pass_Chl_threshold, image_to_crop = img_in_backend, input_toggle,
                          sum_DsRed_imported_for_some_reason = sum_DsRed_global, name_to_parse,
                          CLSmode = "threshold", indices_submitted_to_subset_to, verbose=FALSE, fluorophore_ID_vector,
                          tissue_type = "unsegmented", record_residuals, grid_type){
  
  # Credit for below function: https://stackoverflow.com/questions/31993704/storing-ggplot-objects-in-a-list-from-within-loop-in-r
  plot_data_column = function (data, column, filename, mode_information, fluorophore_ID_vector, grid_type) {
    
    name_for_global_variable_with_cumulative_t_statistic <- paste0("sum_", column, "_global")
    
    colfunc <- colorRampPalette(c("white", "blue"))
    
    plot_out <- ggplot(data=data)+
      # Fill is different depending on if we are in p-val or z-score plotting mode
      # Note redundancy with explicitly calling data
      geom_tile(aes(y=rows, x=cols, fill=data[,get(column)]))+ #https://stackoverflow.com/questions/32184252/how-to-select-columns-in-data-table-using-a-character-vector-of-certain-column-n
      scale_fill_gradientn(colours = colfunc(10))+
      theme(axis.ticks = element_blank(), axis.title = element_blank(),
            panel.background = element_blank(), axis.text = element_blank(),
            panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_blank())+
      labs(title = paste0(column, ": ", format(get(name_for_global_variable_with_cumulative_t_statistic)/1000, digits=1)),
           subtitle = expression('x' ~10^3 ~ '=' ~ Sigma^'pixels' ~ Tau['FP signal']),
           fill = paste0("T(", column,")"))
    
    if(column==fluorophore_ID_vector[1]){
      # Save on the first iteration of the loop (bad practice)
      # Now need to bring the vector into this scope
      output_dir <- paste0("cooloutput/", gsub("-","_",Sys.Date()), "/")
      if(dir.exists(output_dir)==FALSE) dir.create(output_dir, recursive = TRUE)
      
      output_prefix <- paste0(output_dir,
                              file_path_sans_ext(basename(filename)),
                              "_", mode_information)
      
      if(file.exists(paste0(output_prefix, ".CLS"))){
        output_prefix <- paste0(output_prefix,
                                "_",
                                Sys.time())
      }
      
      print(paste0("Saving plot and data for for ", output_prefix))
      ggsave(paste0(output_prefix, ".png"), plot = plot_out)
      fwrite(data, paste0(output_prefix, ".CLS"))
    }
    

    
    return(plot_out)}
  
  if(input_toggle==FALSE){
    return(NA)
  }
  
  ## Now for cropping
  # First, set vars for just cropping down to the whole grid (crop out edges outside of all grid spaces)
  # v0.21 do this only if cropping down to single explants. Leave edges otherwise so we can overlay with RGB images later.
  y_crop_leftside <- 14
  y_crop_rightside <- 1406
  x_crop_leftside <- 1256
  x_crop_rightside <- 226
  
  if(mode=="whole_plate" | mode=="both"){
    crop_bottom <- (x_crop_leftside)
    crop_top <- (x_crop_rightside)
    crop_left <- (y_crop_leftside)
    crop_right <- (y_crop_rightside)
    
    # All I need to do is comment this out? (also see 13 lines down)
    # What about the variable image_spectrum_all_cropped?
    #if(image_type=="hyperspectral"){
    #  image_spectrum_table_processed <- image_spectrum_table[which(image_spectrum_table$cols < crop_bottom & image_spectrum_table$cols > crop_top & image_spectrum_table$rows>crop_left & image_spectrum_table$rows<crop_right),]
    #}
    
    if (image_type=="CLS"){
      image_spectrum_table <- data.table(cbind(image_to_crop$x, image_to_crop$y, image_to_crop[[]]))
      # WHOA MAYBE I SHOULDNT CHANGE THIS BUT I THINK I SHOULD (pre-v0.18)
      colnames(image_spectrum_table) <- c("rows", "cols", image_to_crop@wavelength)
      # Patch added in v.22
      if(CLSmode=="integrate"){
        colnames(image_spectrum_table)[1:2] <- c("rows", "cols")
      }
      # Try changing back to see if it fixed mismatch pixels (stretched plants) in plot (v0.19)
      # The problem of stretched plants (lines) came back when I changed fixed the plot cropping
      # (Was cropping x and y in reverse for segment CLS images only)
      # So now change back and see what happens (v0.20)
      #colnames(image_full_spectrum_cropped) <- c("rows", "cols", image_to_crop@wavelength)
      # Need to change from colnames to indices because of data structure.. NOT. Change to data.table earlier instead.
      
      # And comment this out?
      #image_full_spectrum_cropped <- image_full_spectrum_cropped[which(image_full_spectrum_cropped$cols < crop_bottom & image_full_spectrum_cropped$cols > crop_top & image_full_spectrum_cropped$rows>crop_left & image_full_spectrum_cropped$rows<crop_right),]
      
      image_spectrum_table <- CLS_workflow(spectrum_in_to_CLS = image_spectrum_table,
                                           pass_FP_threshold_from_input = first_pass_FP_threshold,
                                           pass_Chl_threshold_from_input = first_pass_Chl_threshold,
                                           mode = CLSmode, indices_submitted_to_subset_to = indices_submitted_to_subset_to,
                                           fluorophore_ID_vector = fluorophore_ID_vector,
                                           filename = name_to_parse, grid_item = grid_item, tissue_type = tissue_type,
                                           record_residuals = record_residuals)
    }
    if(image_type=="crop_segment_output"){
      # This mode is different, outputs the cropped (post-transformed by alignment) segmentation output
      # ...instead of a plot like the other two modes.
      # I realized I had to flip x and y here 11.22.19
      # Or not. This made the lines come back.
      # I need to revert the WHOA... but only for segment...?
      # Comment this out? v0.22
      #image_spectrum_table <- image_spectrum_table[which(image_spectrum_table$x < crop_bottom & image_spectrum_table$x > crop_top & image_spectrum_table$y>crop_left & image_spectrum_table$y<crop_right),]
      return(image_spectrum_table)
    }
    
    #pal <- wes_palette("Zissou1", max(image_spectrum_table[,3:5]), type = "continuous") #https://www.datanovia.com/en/blog/ggplot-colors-best-tricks-you-will-love/
    
    # New debugging statement in v0.19 because of error Error in +ggtitle(FullID) : invalid argument to unary operator
    # Is this variable defined within the proper scope?
    # Note: I think that error only happened because I had the + before ggtitle on the same line
    # instead of the line before
    if(verbose==TRUE){
      print(paste0("Full ID right before plotting is ", FullID))
    }
    
    
    # This was commented out for some reason in v0.18 but I am uncommenting in v0.19 because of
    # error when using backend for PhenotypeAssistant
    # Error in crop_and_plot(mode = what_to_plot, input_toggle = TRUE, image_spectrum_table = image_spectrum_table_colored,  : 
    # object 'p' not found
    p <- ggplot(data=image_spectrum_table, aes(x=cols, y=rows, fill=color)) +
      coord_equal() + geom_tile() + scale_fill_identity() +
      theme(axis.ticks = element_blank(), axis.title = element_blank(),
            panel.background = element_blank(), axis.text = element_blank(),
            panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_blank()) +
      ggtitle(FullID)
    
    if (image_type=="CLS"){
      if(verbose==TRUE){
        print("Total DsRed in final image about to be plotted: ")
        print(sum_DsRed_imported_for_some_reason)
      }

      p <- p + ggtitle(paste0("Total DsRed coefficient over area: ", sum_DsRed_imported_for_some_reason))

      myplots <- lapply(fluorophore_ID_vector, plot_data_column, data=image_spectrum_table, filename=name_to_parse,
                        mode_information="Uncropped", fluorophore_ID_vector = fluorophore_ID_vector, grid_type = grid_type)
      p <- grid.arrange(grobs=myplots)
      
    }
    
  }
  if(mode=="single_explant" | mode=="both"){
    # Adding this in v0.19, not sure why it needs to be here suddenly
    FullID <<- paste0(parse_trayplateID(name_to_parse), "_exp", grid_item)
    y_newsize <- y_crop_rightside - y_crop_leftside
    x_newsize <- x_crop_leftside - x_crop_rightside
    
    # Calculate size of each grid item
    if(grid_type==20){
      y_grid_rows <- 4
      x_grid_columns <- 6
    }
    if(grid_type==12){
      y_grid_rows <- 3
      x_grid_columns <- 4
    }

    x_grid_item_size <- x_newsize/y_grid_rows
    y_grid_item_size <- y_newsize/x_grid_columns
    
    # Find row and column of desired grid item
    row_col <- get_row_col_for_grid_item(grid_item = grid_item, grid_type = grid_type)
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
                                                          pass_Chl_threshold_from_input = first_pass_Chl_threshold,
                                                          fluorophore_ID_vector = fluorophore_ID_vector,
                                                          filename = name_to_parse, grid_item = grid_item,
                                                          tissue_type = tissue_type,
                                                          record_residuals = record_residuals)
    }
 
    if(verbose==TRUE){
      print(paste0("Full ID right before plotting is ", FullID))
    }
    
    q <- ggplot(data=image_spectrum_table_single_cropped, aes(x=cols, y=rows, fill=color)) +
      coord_equal() + geom_tile() + scale_fill_identity() +
      theme(axis.ticks = element_blank(), axis.title = element_blank(),
            panel.background = element_blank(), axis.text = element_blank(),
            panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_blank())+
      # Uncommented in v0.19
    ggtitle(FullID)
    
    if (image_type=="CLS"){
      if(verbose==TRUE){
        print("Total DsRed in final image about to be plotted: ")
        print(sum_DsRed_imported_for_some_reason)
      }

      q <- q + ggtitle(paste0("Total DsRed coefficient over area: ", sum_DsRed_imported_for_some_reason))
      
      image_spectrum_table_single_cropped[image_spectrum_table_single_cropped < 0] <- 0
      
      #myplots <- vector('list', ncol(imagr_spectrum_table_single_cropped))
      
      #fluorophore_pval_ID_vector <- paste0(fluorophore_ID_vector, "_pval")
      
      myplots <- lapply(fluorophore_ID_vector, plot_data_column,
                        data=image_spectrum_table_single_cropped,
                        filename=name_to_parse,
                        mode_information=paste0("SingleExplantNo",grid_item),
                        fluorophore_ID_vector = fluorophore_ID_vector,
                        grid_type = grid_type)
      q <- grid.arrange(grobs=myplots)
      #print(q)
      #stop("Break here")
      
      # b <- ggplot(data=image_spectrum_table_single_cropped)+
      #   geom_tile(aes(y=rows, x=cols, fill=ZsYellow))+
      #   scale_fill_gradientn(colours = heat.colors(10))+
      #   theme(axis.ticks = element_blank(), axis.title = element_blank(),
      #         panel.background = element_blank(), axis.text = element_blank(),
      #         panel.border = element_blank(), panel.grid.major = element_blank(),
      #         panel.grid.minor = element_blank(), axis.line = element_blank())+
      #   ggtitle(paste0("Total ZsYellow in image:\n", 
      #                  format(sum_ZsYellow_global/1000, digits=1), "k"))
      # 
      # c <- ggplot(data=image_spectrum_table_single_cropped)+
      #   geom_tile(aes(y=rows, x=cols, fill=ChlA))+
      #   scale_fill_gradientn(colours = heat.colors(10))+
      #   theme(axis.ticks = element_blank(), axis.title = element_blank(),
      #         panel.background = element_blank(), axis.text = element_blank(),
      #         panel.border = element_blank(), panel.grid.major = element_blank(),
      #         panel.grid.minor = element_blank(), axis.line = element_blank())+
      #   ggtitle(paste0("Total Chl-A in image:\n", 
      #                  format(sum_ChlA_global/1000, digits=1), "k"))
      # 
      # d <- ggplot(data=image_spectrum_table_single_cropped)+
      #   geom_tile(aes(y=rows, x=cols, fill=ChlB))+
      #   scale_fill_gradientn(colours = heat.colors(10))+
      #   theme(axis.ticks = element_blank(), axis.title = element_blank(),
      #         panel.background = element_blank(), axis.text = element_blank(),
      #         panel.border = element_blank(), panel.grid.major = element_blank(),
      #         panel.grid.minor = element_blank(), axis.line = element_blank())+
      #   ggtitle(paste0("Total Chl-B in image:\n", 
      #                  format(sum_ChlB_global/1000, digits=1), "k"))
      # 
      # q <- grid.arrange(a, b, c, d)
    }
    
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

decide_what_to_plot <- function(mode=NA, verbose=FALSE){
  if(verbose==TRUE) print(mode)
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
  if(verbose==TRUE) print(to_do)
  return(to_do)
}

CLS_workflow <- function(spectrum_in_to_CLS = image_full_spectrum_cropped, pass_FP_threshold_from_input, pass_Chl_threshold_from_input, plot=FALSE,
                         mode = "threshold", indices_submitted_to_subset_to, verbose=FALSE, fluorophore_ID_vector,
                         filename, grid_item, tissue_type, record_residuals){
  
  
  choose_pixels <- function(FP_lambda='582.6952', Chl_lambda='508.763', full_spectrum=this_full_spectrum, CLS_threshold_FP, CLS_threshold_Chl,
                            Mode = mode, Indices_submitted_to_subset_to = indices_submitted_to_subset_to, verbose=FALSE){
    
    
    
    all_pixel_indices <- seq.int(nrow(full_spectrum))
    
    if(verbose==TRUE){
      print("Length of all_pixel_indices")
      print(length(all_pixel_indices))
    }
    
    if (Mode == "threshold"){
      
      # Need to chop off x and y columns up front in this mode
      # Actually, I don't think this is needed since all this function outputs is a 
      # ... list of pixel indices
      full_spectrum=full_spectrum[,3:ncol(full_spectrum)]
      
      if(verbose==TRUE){
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
      }

      # Maybe I can replace the whole use of the all_pixel_indices vector with arr.ind=TRUE option in which
      # Will that fix bug?
      
      #pixels_list_first_pass <- all_pixel_indices[which(full_spectrum[,get(FP_lambda)] > CLS_threshold_FP)]
      # First subset by FP, then Chl intensity
      pixels_list_first_pass <- which(arr.ind = TRUE, x = full_spectrum[,get(FP_lambda)] > CLS_threshold_FP)
      
      if(verbose==TRUE){
        print(paste0("Significant FP is at least ", CLS_threshold_FP))
        print("Number of pixels with significant FP:")
        
        print(length(pixels_list_first_pass))
        print("Head of pixels with significant FP:")
        
        print(head(pixels_list_first_pass))
      }

      pixels_with_significant_Chl <- all_pixel_indices[which(full_spectrum[,get(Chl_lambda)] > CLS_threshold_Chl)]
      
      pixels_with_something <- unique(c(pixels_list_first_pass, pixels_with_significant_Chl))
    }
    
    if (Mode == "integrate"){
      if(verbose==TRUE){
        print("Start subsetting to selected pixels of this tissue at...")
        print(Sys.time())
        
        # THIS (0.3) SHOULD NOT BE HARDCODED
        print("What is the format of the object Indices_submitted_to_subset_to ?")
        print(dim(Indices_submitted_to_subset_to))
        print(head(Indices_submitted_to_subset_to))
        #pixels_with_something <- which(arr.ind = TRUE, x = Indices_submitted_to_subset_to[,4] > 0.3)
        print("What is the length of the object all_pixel_indices?")
        print(dim(all_pixel_indices))
        print("What does all_pixel_indices it look like?")
        print(head(all_pixel_indices))
        print(tail(all_pixel_indices))
        print("Head, dim, head factors, min and max of Indices_submitted_to_subset_to:")
        print(head(Indices_submitted_to_subset_to))
        print(dim(Indices_submitted_to_subset_to))
        print(head(levels(factor(Indices_submitted_to_subset_to[,4]))))
        print(min(Indices_submitted_to_subset_to[,4]))
        print(max(Indices_submitted_to_subset_to[,4]))
      }


      pixels_with_something <- all_pixel_indices[which(Indices_submitted_to_subset_to[,4] > 0.05)]
      
      if(verbose==TRUE){
        print("To help troubleshoot cropping, maximum value (for pixel index) in pixels_with_something is")
        print(max(pixels_with_something))
        print("... and the dimensions of pixels_with_something are: ")
        print(dim(pixels_with_something))
        print("Finish at")
        print(Sys.time())
        #stop("Stopping here now")
      }
    }
    
    if(verbose==TRUE) {
      print("Done with function choose_pixels, now returning pixels_with_something")
      cat("\n\n")
    }

    return(pixels_with_something)
  }
  
  model_and_record_single_pixel <- function(i, Xmatrix = mm, full_spectrum_to_regress_over, fluorophore_ID_vector,
                                            A, Z, backend, profiling = TRUE, record_residuals){
    #stop("Break here to profile")
    #if(profiling==TRUE){
    #  profile_outpath <- paste0("profiling_outputs/single_pixel_", Sys.time(), ".Rprof")
    #  Rprof(filename = profile_outpath)
    #}
    # This should not take in the image from the backend, but something which is passed through
    # and can be either cropped or uncropped image
    # print(spectrum_in_to_CLS[10:10])
    
    # Putting statements in here a bad idea if working w/ many pixels
    
    #print("Unformatted spectrum: ")
    #print(spectrum_in_to_CLS[200:210,1:10])
    #print("Dimensions of spectrum_in_to_CLS within model_and_record_single pixel:")
    #print(dim(spectrum_in_to_CLS))
    
    
    
    #spectrum_formatted <- as.numeric(as.matrix(unlist(spectrum_in_to_CLS[i,]), nrow=1))
    #spectrum_formatted <- spectrum_in_to_CLS[.I[i]]
    
    spectrum_formatted <- t_spectrum_in_to_CLS[,i]
    
    #cat("\n\n")
    
    #print("Xmatrix: ")
    #print(Xmatrix[200:210,1:5])
    #print("Is Xmatrix numeric?")
    #print(is.numeric(Xmatrix))
    #cat("\n\n")
    
    if(backend=="ByHand"){
      Beta <- Z %*% spectrum_formatted
      resid <- spectrum_formatted - (mm %*% Beta)
      p <- ncol(mm) - 1 # minus one because we don't count the intercept as a parameter
      n <- nrow(mm)
      MSE <- (t(resid) %*% resid)/(n-p-1)
      Beta.covar.matrix <- as.vector(MSE)*A
      Beta.se <- sqrt(diag(Beta.covar.matrix))
    }
    
    #model_out <- fastLmPure(Xmatrix, spectrum_formatted)
    #print(model_out$coefficients)
    
    #print("Dimensions of CLS_table within model_and_record_single_pixel")
    #print(dim(CLS_table))
    # Make sure colnames are set properly
    # df is n - p where n is # wavelengths and is parameters
    n <- nrow(Xmatrix)
    p <- ncol(Xmatrix)
    df <- n - p - 1
    #stop("Check n and p")
    for(j in 1:length(fluorophore_ID_vector)){
      
      # WHY are t-stats so massive? Over 1500 for DsRed...
      if(backend=="FastLmPure"){
        this_t_stat <- model_out$coefficients[j+1]/model_out$se[j+1]
      }
      if(backend=="ByHand"){
        this_t_stat <- Beta[j+1]/Beta.se[j+1]
      }
      #this_p_val <- as.brob(pnorm(-abs(this_t_stat)))
      #this_p_val <- as.numeric(-log(this_p_val, base=10))
      #this_p_val <- pt(-abs(this_t_stat), df=df)
      #this_p_val <- -log(this_p_val, base=10)
      #stop("Break here")
      #CLS_table[i, eval(fluorophore_ID_vector[j]) := this_t_stat]
      set(CLS_table, i, eval(fluorophore_ID_vector[j]), this_t_stat)
      #CLS_table[i, eval(paste0(fluorophore_ID_vector[j], "_pval")) := this_p_val]
      #colMeans(CLS_table[,3:10], na.rm = TRUE)
      #stop("Stop inside model")
      
      if(record_residuals == TRUE){
        set(CLS_table, i, length(fluorophore_ID_vector)+j, resid[j])
      }
      
    }
    
    #stop("Break inside model and record single pixel")
    
    # Our test is one-sided and we are only interested in enrichment of protein:
    # https://www.biostars.org/p/17227/
    #CLS_table[i, DsRed_pval := model_out$coefficients[2]/model_out$se[2]]
    #CLS_table[i, ZsYellow_pval := model_out$coefficients[3]/model_out$se[3]]
    #CLS_table[i, ChlA_pval := model_out$coefficients[4]/model_out$se[4]]
    #CLS_table[i, ChlB_pval := model_out$coefficients[5]/model_out$se[5]]
    
    #print("This model:")
    #print(model_out)
    #stop("Stopping here")
    
    #if(profiling == TRUE){
    #  Rprof(filename = NULL)
    #}
  }
  
  perform_CLS <- function(pixels_to_test,
                          verbose=FALSE,
                          fluorophore_ID_vector,
                          backend="ByHand",
                          record_residuals){
    

    
    start_time <- Sys.time()
    if(backend=="ByHand"){
      A <- solve(t(mm) %*% mm)
      Z <- A %*% t(mm)
    }else{
      A <- NA
      Z <- NA
    }
    
    #options(warn=-1)
    #pb <- txtProgressBar(min = 0, max = length(pixels_to_test), initial = 0)
    if(verbose==TRUE){
      print("About to start CLS over all pixels passing threshold. Here is table:")
      print(head(CLS_table))
      #CLS_table_out <- foreach(i=pixels_to_test) %do% {
      print(paste0("Number of pixels to test: ", length(pixels_to_test)))
      print("First 5 pixels to test: ")
      print(head(pixels_to_test))
    }

    # Need to get rid of the row and col columns I used for indexing earlier
    spectrum_in_to_CLS$rows <<- NULL
    spectrum_in_to_CLS$cols <<- NULL
    if(verbose==TRUE) print(head(CLS_table))
    if (shiny::isRunning() == "PROGRESSBARSLOWSDOWN"){
    #if (shiny::isRunning() == TRUE){
      withProgress(message = "Running CLS over all pixels passing denoising",
                   max=length(pixels_to_test),
                   value = 0, {
                     for(i in pixels_to_test){
                       #print(i)
                       model_and_record_single_pixel(i, fluorophore_ID_vector = fluorophore_ID_vector,
                                                     A = A, Z = Z, backend = backend,
                                                     record_residuals = record_residuals)
                       incProgress(1)
                     }      
                   })
    }else{
      for(i in pixels_to_test) model_and_record_single_pixel(i, fluorophore_ID_vector = fluorophore_ID_vector,
                                                             A = A, Z = Z, backend = backend,
                                                             record_residuals = record_residuals)
    }
    
    options(warn=0)
    
    if(verbose==TRUE){
      print("Done. Here is table again:")
      print(head(CLS_table))
    }
    run_time <- Sys.time() - start_time
    print(paste0("Runtime with backend ", backend, " was ", run_time))
  }
  
  finish_processing_CLS_table <- function(CLS_table, chosen_pixels_in = chosen_pixels, fluorophore_ID_vector,
                                          filename, grid_item, tissue_type, record_residuals){
    # THIS IS SLOW
    # Faster way: https://stackoverflow.com/questions/38226323/replace-all-values-in-a-data-table-given-a-condition?rq=1
    
    print("Printing fluorophore ID vector inside finishing_processing_CLS_table:")
    print(fluorophore_ID_vector)
    
    if(verbose==TRUE) print(Sys.time())
    
    #CLS_table[is.na(CLS_table)] <- 0
    if(verbose==TRUE) print(Sys.time())
    #CLS_table[!chosen_pixels_in, DsRed := 0]
    #CLS_table[!chosen_pixels_in, ZsYellow := 0]
    #CLS_table[!chosen_pixels_in, ChlA := 0]
    #CLS_table[!chosen_pixels_in, ChlB := 0]
    
    if(verbose==TRUE){
      print(paste0("Maximum DsRed: ", max(CLS_table$DsRed)))
      print(paste0("Maximum ChlB: ", max(CLS_table$ChlB)))
    }
    
    #Somehow lost colnames and need to get them back. Possibly related to a similar issue noted in coment in
    #FastLM.R on Github (RcppCore/RcppArmadillo)???
    colnames(CLS_table)[1:2] <- c("rows", "cols")
    colname_vector <- c(fluorophore_ID_vector, paste0(fluorophore_ID_vector, "_pval"))
    colnames(CLS_table)[3:ncol(CLS_table)] <- colname_vector
    print(paste0("Colnames for CLS_table in finishing_processing_CLS_table"))
    print(colnames(CLS_table))
    #colnames(CLS_table)[1:6] <- c("rows", "cols", "DsRed", "ZsYellow", "ChlA", "ChlB")
    
    cumulative_t_stats <- colSums(CLS_table[,3:(3+length(fluorophore_ID_vector))], na.rm = TRUE)
    print("Cumulative t-stats are:")
    print(cumulative_t_stats)
    
    for(i in 1:length(fluorophore_ID_vector)){
      name_for_this_global_variable <- paste0("sum_", fluorophore_ID_vector[i], "_global")
      print(paste0("Making global variable ", name_for_this_global_variable,
                   " equal to ", cumulative_t_stats[i]))
      assign(name_for_this_global_variable, cumulative_t_stats[i], envir=globalenv())
    }
      
    CLS_table_output_folder_path <- paste0("output/CLS_tables/", Sys.Date(), "/")
    
    if(!dir.exists(CLS_table_output_folder_path)){
      dir.create(CLS_table_output_folder_path, recursive = TRUE)
    }
    
    write.table(CLS_table,
                file = paste0(paste0(CLS_table_output_folder_path,
                                     basename(file_path_sans_ext(filename)),
                                     ".csv")),
                sep = ",",
                quote = FALSE,
                col.names = TRUE,
                row.names = FALSE)
    
    
    if(!dir.exists(paste0("output/sum_stats/"))){
      dir.create("output/sum_stats/")
    }
    
    
    
    if(!file.exists(paste0("output/sum_stats/", Sys.Date(), ".csv"))){
      file.create(paste0("output/sum_stats/", Sys.Date(), ".csv"))
      # Check this on another day to make sure it works properly
      #header_to_write <- matrix(c("Filename", "Grid_Item", "Tissue_Type", fluorophore_ID_vector))
      #write.table(header_to_write, file = paste0("output/sum_stats/",
      #                                           Sys.Date(),
      #                                           ".csv"),
      #            sep = ",",
      #            append = TRUE,
      #            quote = FALSE,
      #            col.names = FALSE,
      #            row.names = FALSE)
    }
    
    
    line_to_write <- matrix(c(filename, grid_item, tissue_type, cumulative_t_stats), nrow=1)
    
    write.table(line_to_write,
                file = paste0("output/sum_stats/",
                              Sys.Date(),
                              ".csv"),
                sep = ",",
                append = TRUE,
                quote = FALSE,
                col.names = FALSE,
                row.names = FALSE)
    
    #writeLines(line_to_write, sum_stat_file_connection)
    
    CLS_table$ChlB[which(CLS_table$ChlB<0)] <- 0
    CLS_table$DsRed[which(CLS_table$DsRed<0)] <- 0
    CLS_table$ChlA[which(CLS_table$ChlA<0)] <- 0
    CLS_table$ZsYellow[which(CLS_table$ZsYellow<0)] <- 0
    
    #for(i in 1:length(fluorophore_ID_vector)){
    #  cumulative_t_stat <- sum(CLS_table[,eval(fluorophore_ID_vector[i])])
    #  name_for_this_global_variable <- paste0("sum_", fluorophore_ID_vector[i], "_global")
    #  print(paste0("Making global variable ", name_for_this_global_variable,
    #               " equal to ", cumulative_t_stat))
    #  assign(name_for_this_global_variable, cumulative_t_stat, envir=globalenv())
    #}
    
    #sum_DsRed_local <- sum(CLS_table$DsRed)
    #sum_ZsYellow_local <- sum(CLS_table$ZsYellow)
    #sum_ChlA_local <- sum(CLS_table$ChlA)
    #sum_ChlB_local <- sum(CLS_table$ChlB)
    #assign("sum_DsRed_global", sum_DsRed_local, envir=globalenv())
    #assign("sum_ZsYellow_global", sum_ZsYellow_local, envir=globalenv())
    #assign("sum_ChlA_global", sum_ChlA_local, envir=globalenv())
    #assign("sum_ChlB_global", sum_ChlB_local, envir=globalenv())
    
    if(verbose==TRUE){
      print("!!!!!!!!!!!!!!SUMS!!!!!!!!!!!!!!!!!")
      print(sum_DsRed_local)
      print(sum_ZsYellow_local)
      print(sum_ChlA_local)
      print(sum_ChlB_local)
      print(sum_DsRed_global)
      print(sum_ZsYellow_global)
      print(sum_ChlA_global)
      print(sum_ChlB_global)
    }

    
    #CLS_table$ChlBfalsecolor <- rescale(CLS_table$ChlB,
                                        #from=range(c(denoise_threshold_Chl,max_intensity_Chl)),
    #                                    to=c(0,255))
    #CLS_table$DsRedfalsecolor <- rescale(CLS_table$DsRed,
                                         #from=range(c(denoise_threshold_FP,max_intensity_FP)),
    #                                     to=c(0,255))
    #CLS_table$color <- rgb(red= CLS_table$ChlBfalsecolor,
    #                       green= CLS_table$DsRedfalsecolor,
    #                       blue= rep(0, nrow(CLS_table)),
    #                       maxColorValue = 255)
    if(verbose==TRUE){
      print(dim(CLS_table))
      print(head(CLS_table))
    }

    return(CLS_table)
  }
  
  # fitting a smooth curve https://stackoverflow.com/questions/3480388/how-to-fit-a-smooth-curve-to-my-data-in-r
  read_plot_spectra <- function(spectra_path, plot=FALSE){
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
      
      # Note here we are plotting the published lines over our data (should be the same as if we used our lines)
      lines(pub_emission_spectrum$wavelength,
            predict(fit,pub_emission_spectrum$wavelength), col='red', lwd=2)
    }
    
    #return(predictions)
    fitted_spectra_dataframe <- as.data.frame(cbind(wavelengths, predictions))
    colnames(fitted_spectra_dataframe) <- c("Wavelength", "Normalized intensity (fitted)")
    
    # These two lines may be superfluous
    # Confirmed they are superfluous in integration_a2 jupyter notebook
    #fitted_spectra_dataframe[,2] <- as.numeric(as.character(fitted_spectra_dataframe[,2]))
    #fitted_spectra_dataframe[is.na(fitted_spectra_dataframe)] <- 0
    
    
    if(verbose==TRUE) print(fitted_spectra_dataframe)
    
    return(fitted_spectra_dataframe)
  }
  ChlB <- read_plot_spectra("/scratch2/NSF_GWAS/macroPhor_Array/Fluorophore_spectra/ChlB.csv")
  ChlA <- read_plot_spectra("/scratch2/NSF_GWAS/macroPhor_Array/Fluorophore_spectra/ChlA.csv")
  DsRed <- read_plot_spectra("/scratch2/NSF_GWAS/macroPhor_Array/Fluorophore_spectra/DsRed_Clontech.csv")
  ZsYellow <- read_plot_spectra("/scratch2/NSF_GWAS/macroPhor_Array/Fluorophore_spectra/ZsYellow_cutoff.csv")
  
  convert_NA_to_zero_in_spectrum_and_return_vector <- function(spectrum_in, verbose=FALSE){
    vector <- as.vector(spectrum_in$'Normalized intensity (fitted)')
    if(verbose==TRUE) print(vector[1:5])
    vector[which(is.na(vector)==TRUE)] <- 0
    if(verbose==TRUE) print(vector[1:5])
    return(as.numeric(as.character(vector)))
    
  }
  
  ChlA_vector <- convert_NA_to_zero_in_spectrum_and_return_vector(ChlA)
  ChlB_vector <- convert_NA_to_zero_in_spectrum_and_return_vector(ChlB)
  DsRed_vector <- convert_NA_to_zero_in_spectrum_and_return_vector(DsRed)
  ZsYellow_vector <- convert_NA_to_zero_in_spectrum_and_return_vector(ZsYellow)
  
  if(verbose==TRUE) print(ChlA_vector[1:5])
  
  #_.-._.-._.-._.-._.-._.-._.-._.-._.-._.-._.-._.-._.-._.-._
  ##,'_.-._.-._.-._.-._.-._.-._.-._.-._.-._.-._.-._.-._.-._.-._`.
  #( (                                                         ) )
  ##) )                                                       ( (
  #( (                                                         ) )
  ##) )                                                       ( (
  #( (             STOP HARDWIRING                            ) )
  ##) )                                                       ( (
  #( (                                                         ) )
  ##) )                                                       ( (
  #( (_.-._.-._.-._.-._.-._.-._.-._.-._.-._.-._.-._.-._.-._.-._) )
  #`._.-._.-._.-._.-._.-._.-._.-._.-._.-._.-._.-._.-._.-._.-._,'
  
  fluorophore_ID_vector <- c("DsRed", "ZsYellow", "ChlA", "ChlB")
  
  n_fluorophores <- length(fluorophore_ID_vector)
  
  for (i in 1:length(fluorophore_ID_vector)){
    fluorophore_lambda_vector_to_retrieve <- paste0(fluorophore_ID_vector[i], "_vector")
    this_fluorophore_lambda_vector <- get(fluorophore_lambda_vector_to_retrieve)
    # If it is the first one, initialize
    if(i==1){
      mm <<- as.matrix(cbind(1, this_fluorophore_lambda_vector))
    }else{
      mm <<- cbind(mm, this_fluorophore_lambda_vector)
    }
    if(i==length(fluorophore_ID_vector)){
      colnames(mm) <<- c("Intercept", paste0(fluorophore_ID_vector, "_vector"))
    }

    
  }
  
  print(head(mm))
  
  #_.-._.-._.-._.-._.-._.-._.-._.-._.-._.-._.-._.-._.-._.-._
  ##,'_.-._.-._.-._.-._.-._.-._.-._.-._.-._.-._.-._.-._.-._.-._`.
  #( (                                                         ) )
  ##) )                                                       ( (
  #( (                                                         ) )
  ##) )                                                       ( (
  #( (             STOP HARDWIRING                            ) )
  ##) )                                                       ( (
  #( (                                                         ) )
  ##) )                                                       ( (
  #( (_.-._.-._.-._.-._.-._.-._.-._.-._.-._.-._.-._.-._.-._.-._) )
  #`._.-._.-._.-._.-._.-._.-._.-._.-._.-._.-._.-._.-._.-._.-._,'
  
  # HERE WE CAN SET THE INTERCEPT AS ZERO OR ONE. SHOULD BE HARDWIRED?
  #mm <<- as.matrix(cbind(1, DsRed_vector, ZsYellow_vector, ChlA_vector, ChlB_vector))
  
  if(verbose==TRUE){
    print("Making CLS table")
    print("Dimensions of spectrum_in_to_CLS being used to make CLS table: ")
    print(dim(spectrum_in_to_CLS))
  }

  CLS_table_server_side <- data.table(cbind(as.numeric(spectrum_in_to_CLS$rows),
                                            as.numeric(spectrum_in_to_CLS$cols),
                                            matrix(nrow=dim(spectrum_in_to_CLS)[1],
                                                   ncol=length(fluorophore_ID_vector)*2)))
  
  CLS_table_server_side[] <- lapply(CLS_table_server_side, as.numeric)
  
  if(verbose==TRUE){
    print(CLS_table_server_side[1:10,1:6])
    print("type of row col")
    print(typeof(spectrum_in_to_CLS$rows))
    print("Dimensions of CLS table after being made:")
    print(dim(CLS_table_server_side))
    
    # CLS_table_server_side[,1] <- as.numeric(as.character(CLS_table_server_side[,1]))
    # CLS_table_server_side[,2] <- as.numeric(as.character(CLS_table_server_side[,2]))
    
    print("Converting columns to numeric")
  }

  #CLS_table_server_side[,3] <- as.numeric(as.character(CLS_table_server_side[,3]))
  #CLS_table_server_side[,4] <- as.numeric(as.character(CLS_table_server_side[,4]))
  #CLS_table_server_side[,5] <- as.numeric(as.character(CLS_table_server_side[,5]))
  #CLS_table_server_side[,6] <- as.numeric(as.character(CLS_table_server_side[,6]))
  
  if(verbose==TRUE){
    print("Done converting")
    
    print("Assigning CLS table to a global object")
  }

  assign("CLS_table", CLS_table_server_side, envir=globalenv())
  
  # FIX THIS
  colnames(CLS_table)[1:2] <<- c("rows", "cols")
  
  colname_vector <- c(fluorophore_ID_vector, paste0(fluorophore_ID_vector, "_pval"))
  
  colnames(CLS_table)[3:ncol(CLS_table)] <<- colname_vector
  
  #colnames(CLS_table)[1:6] <- c("rows", "cols", "DsRed", "ZsYellow", "ChlA", "ChlB")
  
  if(verbose==TRUE){
    print("About to perform CLS")
    print("Pixels we will perform CLS over will first be determined.")
    print(paste0("Threshold for FP: ", pass_FP_threshold_from_input))
    print(paste0("Treshold for Chl: ", pass_Chl_threshold_from_input))
  }

  chosen_pixels <<- choose_pixels(full_spectrum=spectrum_in_to_CLS,
                                  CLS_threshold_FP=pass_FP_threshold_from_input,
                                  CLS_threshold_Chl=pass_Chl_threshold_from_input)
  
  if(verbose==TRUE){
    print("Type of output from choose_pixels for chosen pixels is: ")
    print(typeof(chosen_pixels))
    print("...and it has dimensions: ")
    print(dim(chosen_pixels))
    print(paste0("Number of chosen pixels: ", length(chosen_pixels)))
    print("First 5: ")
    print(head(chosen_pixels))
  }

  # NEED TO ADD THIS TO GET RID OF BOTTLENECK OF CONVERTING ONE COLUMN AT A TIME TO MATRIX
  # v0.27
  # ALSO NEED TO TRANSPOSE SINCE SELECTION OF COLS FROM MATRIX IS SUPERFAST
  t_spectrum_in_to_CLS <<- t(as.matrix(spectrum_in_to_CLS[,3:ncol(spectrum_in_to_CLS)]))
  #full_spectrum <- t(as.matrix(spectrum_in_to_CLS))
  
  perform_CLS(pixels_to_test = chosen_pixels, fluorophore_ID_vector = fluorophore_ID_vector,
              record_residuals = record_residuals)
  
  if(verbose==TRUE) print("CLS done for whole image. Finish data processing.")
  
  CLS_table <- finish_processing_CLS_table(CLS_table, fluorophore_ID_vector = fluorophore_ID_vector,
                                           filename = filename, grid_item = grid_item, tissue_type = tissue_type,
                                           record_residuals = record_residuals)

  
  if(verbose==TRUE){
    print("Head of CLS_table right after it comes out of finish_processing_CLS_table workflow")
    print(head(CLS_table))
    print("total dsred is.......")
    print(sum(CLS_table$DsRed))
    print("CLS data processing finished.")
  }
  
  return(CLS_table)
}































extract_pixels_of_tissue <- function(segmented_image, tissue_type, verbose=FALSE){
  segmented_image[,,,-(tissue_type)] = 0
  
  if(tissue_type==2){
    segmented_image[,,,2][which(segmented_image[,,,2] < 0.3)] = 0
  }
  if(tissue_type==3){
    segmented_image[,,,3][which(segmented_image[,,,3] < 0.3)] = 0
  }
  plot(segmented_image)
  if(verbose==TRUE){
    print("Head of levels of signal for tissue type on color channel ", tissue_type )
    print(head(levels(factor(segmented_image[,,,tissue_type]))))
  }
  
  return(as.data.table(as.data.frame(as.cimg(segmented_image))))
}

stack_segments_and_hyperspectral <- function(segmented_image_preprocessed, hyperspectral_data, verbose=FALSE){
  segmented_image_preprocessed$x = segmented_image_preprocessed$x - 1
  segmented_image_preprocessed$y = segmented_image_preprocessed$y - 1
  
  # Need to crop to match sizes since we chop off a few pixels due to varying hyperspectral camera?
  x_mutual <- min(max(hyperspectral_data$x), max(segmented_image_preprocessed$x))
  y_mutual <- min(max(hyperspectral_data$y), max(segmented_image_preprocessed$y))
  
  hyperspectral_data <- hyperspectral_data[which(hyperspectral_data$x <= x_mutual & which(hyperspectral_data$y <= y_mutual)),]
  
  segmented_image_preprocessed <- segmented_image_preprocessed[which(segmented_image_preprocessed$x <= x_mutual & segmented_image_preprocessed$y <= y_mutual),]
  #segmented_image_preprocessed <- segmented_image_preprocessed[which(segmented_image_preprocessed$cc == 2)]
  
  # Need to sort x then y for both hyperspectral and deep segmentation output
  # This is default on hyperspectral
  keycol <- c("x", "y")
  setorderv(segmented_image_preprocessed, keycol)
  
  order_column_desired <- as.data.table(cbind(paste0(hyperspectral_data$x, "_", hyperspectral_data$y), as.numeric(seq(1, nrow(hyperspectral_data)))))
  colnames(order_column_desired) = c("order_column", "index")
  
  # We need to take a roundabout approach to reordering the rows of the deep segmentation output to the same order as in the hyperspectral output
  # Not so straightforward since it seems to be a strange order, need to make a dummy column
  segmented_image_preprocessed$order_column <- paste0(segmented_image_preprocessed$x, "_", segmented_image_preprocessed$y)
  
  merged = merge(segmented_image_preprocessed, order_column_desired, by="order_column")
  
  merged$index <- as.numeric(merged$index)
  
  setorder(merged, index)
  
  segmented_image_preprocessed <- merged[,2:5]
  
  if(verbose==TRUE){
    print(typeof(segmented_image_preprocessed))
    print(typeof(hyperspectral_data))
    cat("\n\n")
    #print(segmented_image_preprocessed[,1][1000:2000] == hyperspectral_data[,1][1000:2000])
    print("Head of segmented image and hyperspectral data, respectively:")
    print(head(segmented_image_preprocessed))
    print(head(hyperspectral_data))
  }
  
  
  output <- list()
  output[['segmented_image_preprocessed']]<- segmented_image_preprocessed
  output[['hyperspectral_data']]<- hyperspectral_data
  return(output)
}

post_deep_segmentation_CLS <- function(hyperspectral_file_path,
                                       deep_segmentation_output_file_path,
                                       desired_tissue_type, verbose=FALSE, fluorophore_ID_vector,
                                       grid_type){
  
  segmented_image_subset <- extract_pixels_of_tissue(segmented_image=load.image(deep_segmentation_output_file_path),
                                                     tissue_type=desired_tissue_type)
  
  img_in_backend <<- load_image(image_path = hyperspectral_file_path)
  
  stacked_segments_and_hyperspectral <- stack_segments_and_hyperspectral(segmented_image_preprocessed = segmented_image_subset,
                                                                         hyperspectral_data = img_in_backend)
  
  reduced_stacked_hyperspectral <- reduce_image(stacked_segments_and_hyperspectral[2])
  
  image_spectrum_table_colored <- extract_plot_grid_item(image_spectrum_table = reduced_stacked_hyperspectral,
                                                         FC = FALSE,
                                                         cap = FALSE,
                                                         max_intensity_FP = NA,
                                                         max_intensity_Chl = NA,
                                                         scale = FALSE,
                                                         normalize = FALSE,
                                                         to_denoise = FALSE,
                                                         standardize_rgb = FALSE,
                                                         denoise_threshold_FP = NA,
                                                         denoise_threshold_Chl = NA,
                                                         dim_chlorophyll = FALSE,
                                                         FP = "DsRed",
                                                         filename = basename(hyperspectral_file_path),
                                                         verbose = FALSE,
                                                         grid_item = "dummy")
  
  cropped_segment <- crop_and_plot(mode="whole_plate",
                                   grid_item = NA,
                                   image_spectrum_table = stacked_segments_and_hyperspectral[[1]],
                                   image_type="crop_segment_output",
                                   first_pass_FP_threshold = NA,
                                   first_pass_Chl_threshold = NA,
                                   image_to_crop = stacked_segments_and_hyperspectral[[2]],
                                   input_toggle = TRUE,
                                   sum_DsRed_imported_for_some_reason = sum_DsRed_global,
                                   name_to_parse,
                                   CLSmode = NA,
                                   indices_submitted_to_subset_to = NA,
                                   fluorophore_ID_vector = fluorophore_ID_vector,
                                   name_to_parse = paste0(basename(file_path_sans_ext(hyperspectral_file_path)),
                                                          "_tissue_",desired_tissue_type),
                                   grid_type = grid_type)
  if(verbose==TRUE){
    print("LEVELS OF COLOR CHANNEL:")
    print(levels(factor(cropped_segment$cc)))
    
    print("Is the value column for cropped_segment numeric?")
    print(is.numeric(cropped_segment[,4]))
    print("What type is it?")
    print(typeof(cropped_segment[,4]))
    
    print("Why is it (not) numeric? Show the head and head levels of it, respectively.")
    
    print(head(cropped_segment))
    print(head(levels(factor(cropped_segment[,4]))))
    
    print("Now make it numberic")
    cropped_segment[,4] <- as.numeric(as.character(unlist(cropped_segment[,4])))
    
    print("Now look at the head and head levels again...")
    print(head(cropped_segment))
    print(head(levels(factor(cropped_segment[,4]))))
  }
  
  crop_and_plot(mode="whole_plate",
                grid_item = NA,
                image_spectrum_table = image_spectrum_table_colored,
                image_type="CLS",
                first_pass_FP_threshold = NA,
                first_pass_Chl_threshold = NA,
                image_to_crop = img_in_backend,
                input_toggle = TRUE,
                sum_DsRed_imported_for_some_reason = sum_DsRed_global,
                name_to_parse = paste0(basename(file_path_sans_ext(hyperspectral_file_path)),
                                       "_tissue_",desired_tissue_type),
                CLSmode = "integrate",
                indices_submitted_to_subset_to = cropped_segment[which(cropped_segment$cc == desired_tissue_type)],
                fluorophore_ID_vector = fluorophore_ID_vector,
                tissue_type = desired_tissue_type,
                grid_type = grid_type)
  cat("\n\n")
}