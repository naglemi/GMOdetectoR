#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
filename <- "CV3_F1.9_I5.0_L100_cyan_001552_9_1_4.raw"
library(shiny)

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
  ## Applying false color
  if(FP=="DsRed"){
    FPlambda <- '582.8245'
  }
  if(FP=="GFP"){
    FPlambda <- '509.6405'
  }
  
  image_in <- (read.ENVI(image_path))
  image_spectrum <- image_in$spc
  colnames(image_spectrum) <- format(seq(from=399.8645,to=799.8529, by=((799.8529-399.8645)/317)), digits=7)
  image_table <- as.data.table(cbind(image_in$y, image_in$x,
                                     image_spectrum[,'672.4118'],
                                     image_spectrum[,which(colnames(image_spectrum)==FPlambda)],
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

extract_plot_grid_item <- function(grid_item,
                                   image_spectrum_table,
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
                                   standardize_rgb=TRUE){
  
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
  
  ## Now for cropping
  # First, set vars for just cropping down to the whole grid (crop out edges outside of all grid spaces)
  y_crop_leftside <- 14
  y_crop_rightside <- 1406
  x_crop_leftside <- 1256
  x_crop_rightside <- 226
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
  
  image_spectrum_table <- image_spectrum_table[which(image_spectrum_table$cols < crop_bottom & image_spectrum_table$cols > crop_top & image_spectrum_table$rows>crop_left & image_spectrum_table$rows<crop_right),]
  
  # do this AFTER CROPPING and BEFORE SCALING so it doesn't screw up later steps or get nullified by prior steps
  if(standardize_rgb==TRUE){
    image_spectrum_table$r[2] <- denoise_threshold_Chl
    image_spectrum_table$g[2] <- denoise_threshold_FP
    #image_spectrum_table$r[ncol(image_spectrum_table)] <- denoise_threshold_Chl
    #image_spectrum_table$g[ncol(image_spectrum_table)] <- denoise_threshold_FP
    
    image_spectrum_table$r[1] <- max_intensity_Chl
    image_spectrum_table$g[1] <- max_intensity_FP
    #image_spectrum_table$r[ncol(image_spectrum_table)-1] <- max_intensity_Chl
    #image_spectrum_table$g[ncol(image_spectrum_table)-1] <- max_intensity_FP
    #print("Done standardizing")
  }
  
  if(scale==TRUE){
    image_spectrum_table$r <- rescale(image_spectrum_table$r,
                                      from=range(c(denoise_threshold_Chl,max_intensity_Chl)),
                                      to=c(0,255))
    image_spectrum_table$g <- rescale(image_spectrum_table$g,
                                      from=range(c(denoise_threshold_FP,max_intensity_FP)),
                                      to=c(0,255))
  }
  
  if(dim_chlorophyll!=FALSE){
    image_spectrum_table$r <- image_spectrum_table$r * dim_chlorophyll
  }
  
  FullID <- paste0(parse_trayplateID(filename), "_exp", grid_item)
  
  
  if(FC==TRUE){
    image_spectrum_table$color <- rgb(red= image_spectrum_table$r,
                                      green= image_spectrum_table$g,
                                      blue= image_spectrum_table$b,
                                      maxColorValue = max(image_spectrum_table$g))
    
    p <- ggplot(data=image_spectrum_table, aes(x=cols, y=rows, fill=color)) +
      coord_equal() + geom_tile() + scale_fill_identity() +
      theme(axis.ticks = element_blank(), axis.title = element_blank(),
            panel.background = element_blank(), axis.text = element_blank(),
            panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_blank()) +
      ggtitle(FullID)
    # return(list(FullID, p))
    # this is the only line in the function that is being changed
    return(p)
  } else {
    image(matrix(spectrum_cropped[,"g"], nrow=y_grid_item_size-1))
  }
  
  
}

manage_grid_item_database <- function(database_in=NA, new_grid_item=NA){
  grid_item_append <- c()
}

# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("GMOdetectoR"),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
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
        numericInput("num", label = h3("Grid position"), value = 1),
        
        
         sliderInput("bins",
                     "Number of bins:",
                     min = 1,
                     max = 50,
                     value = 30),
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
                     value = 1)
      ),
      
      
      # Show a plot of the generated distribution
      mainPanel(
         plotOutput("distPlot"),
         plotOutput("hyperspectral_false_color")
      )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
   
   output$distPlot <- renderPlot({
      # generate bins based on input$bins from ui.R
      x    <- faithful[, 2] 
      bins <- seq(min(x), max(x), length.out = input$bins + 1)
      
      # draw the histogram with the specified number of bins
      hist(x, breaks = bins, col = 'darkgray', border = 'white')
   })
   ###################################################################################
   
   img_in_backend <- load_image(image_path = paste0("/scratch2/NSF_GWAS/macroPhor_Array/CT_CU_CV_raw/wk6/",
                                                    filename),
                                FP = "DsRed")
   
   filename <- "CV3_F1.9_I5.0_L100_cyan_001552_9_1_4.raw"
   output$hyperspectral_false_color <- renderPlot({
     extract_plot_grid_item(18,
                          image_spectrum_table = img_in_backend,
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
                          dim_chlorophyll = input$dim_chlorophyll,
                          FP = "DsRed")
   })
   ###################################################################################
   
   output$FP <- renderPrint({ input$radio })
   
   output$normalize_now <- renderPrint({ input$action })
   
   output$griditem <- renderPrint({ input$num })
   
}

# Run the application 
shinyApp(ui = ui, server = server)

