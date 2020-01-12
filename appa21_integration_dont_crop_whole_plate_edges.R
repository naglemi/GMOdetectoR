#filename <- "CV1_F1.9_I5.0_L100_cyan_231158_4_0_4.raw"
filename <- "CV3_F1.9_I5.0_L100_cyan_001552_9_1_4.raw" # the one with the good explant on grid item 19
library(shiny)

source("/scratch2/NSF_GWAS/GMOdetectoR/GMOdetectoRv0.25cool.R")

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
        
        numericInput("grid_position", label = h3("Grid position"), value = 1),
        
        # Copy the line below to make a set of radio buttons
        radioButtons("radio", label = h3("Reporter protein"),
                     choices = list("DsRed" = 1, "ZsYellow" = 2, "GFP" = 3), 
                     selected = 1, inline = TRUE),
        
        #hr(),
        
        checkboxGroupInput("checkGroup", label = h3("Plot cropping"), 
                           choices = list("Whole plate" = 1, "Single explant" = 2),
                           selected = c(1, 2), inline = TRUE),
        
        
        checkboxGroupInput("hys_CLS_PCA", label = h3("Plots to build"),
                           choices = list("Hyperspectral"=1,
                                          "CLS"=2,
                                          "PCA"=3),
                           selected = NULL,
                           inline = TRUE),
        
        #h3("Plots to build:"),
        #checkboxInput("hyperspectral_on", "Hyperspectral", FALSE),
        #checkboxInput("CLS_on", "CLS", FALSE),
        #checkboxInput("PCA_on", "PCA", FALSE),
        
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
                     value = 1)
      ),
      
      
      # Show a plot of the generated distribution
      mainPanel(
        tabsetPanel(tabPanel("Hyperspectral Plots", plotOutput("image")),
                    tabPanel("CLS plots", plotOutput("image2")),
                    tabPanel("PCA plots", plotOutput("image3")))
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
                                                            FP = input$FP,
                                                            filename = filename,
                                                            grid_item = input$grid_position)
   })
   
   what_to_plot <- reactive({
     what_to_plot <- decide_what_to_plot(mode = input$checkGroup)
   })
   
   output$image <- renderPlot({
         # parenthesis are needed to avoid error because our table object is really kind of a function
         # https://stackoverflow.com/questions/40623749/what-is-object-of-type-closure-is-not-subsettable-error-in-shiny/40623750
         crop_and_plot(mode=what_to_plot(),
                       input_toggle=(1 %in% input$hys_CLS_PCA),
                       # The below line is currently the same for any type of plot
                       image_spectrum_table = image_spectrum_table_colored(),
                       grid_item = input$grid_position,
                       image_type = "hyperspectral",
                       first_pass_FP_threshold = input$denoise_threshold_FP,
                       first_pass_Chl_threshold = input$denoise_threshold_Chl,
                       name_to_parse = filename,
                       fluorophore_ID_vector = c("DsRed", "ZsYellow", "ChlA", "ChlB"))
       }, height = 800, width=800)
   
   output$image2 <- renderPlot({
     # parenthesis are needed to avoid error because our table object is really kind of a function
     # https://stackoverflow.com/questions/40623749/what-is-object-of-type-closure-is-not-subsettable-error-in-shiny/40623750
     crop_and_plot(mode=what_to_plot(),
                   input_toggle=(2 %in% input$hys_CLS_PCA),
                   # The below line is currently the same for any type of plot
                   image_spectrum_table = image_spectrum_table_colored(),
                   grid_item = input$grid_position,
                   image_type = "CLS",
                   first_pass_FP_threshold = input$denoise_threshold_FP,
                   first_pass_Chl_threshold = input$denoise_threshold_Chl,
                   name_to_parse = filename,
                   fluorophore_ID_vector = c("DsRed", "ZsYellow", "ChlA", "ChlB"))
   }, height = 600, width=600)

   PCA_workflow <- function(image_in=img_in_backend, input_toggle=FALSE){
     print("Starting PCA workflow")
     if(input_toggle==FALSE){
       return(NA)
     }
     pca <- prcomp(~ spc,
                   data = image_in$.,
                   center = FALSE)
     loadings <- decomposition(image_in,
                               t(pca$rotation),
                               scores = FALSE,
                               label.spc = "loading I / a.u.")
     
     scores <- decomposition(image_in,
                             pca$x,
                             label.wavelength = "PC",
                             label.spc = "score / a.u.")
     
     layout.matrix <- matrix(c(1,2, 3, 4), nrow = 2, ncol = 2)
     
     #layout(mat = layout.matrix,
     #       heights = c(2, 2), # Heights of the two rows
     #       widths = c(2, 2)) # Widths of the two columns
     
     #layout.show(4)
     print("Make plots")
     a <- plot(loadings[1:3])
     b <- plotmap (scores [, , 3])
     c <- plotmap (scores [, , 2])
     d <- plotmap (scores [, , 1])
     print("Done PCA workflow")
     return(list(a,b,c,d))
   }
   
  output$image3 <- renderPlot({
    PCA_workflow(image_in=img_in_backend,
                 input_toggle=(3 %in% input$hys_CLS_PCA))
  })
     

   ###################################################################################
   
   output$FP <- renderPrint({ input$radio })
   
   output$normalize_now <- renderPrint({ input$action })
   
   output$griditem <- renderPrint({ input$num })
   
}

# Run the application 
shinyApp(ui = ui, server = server)

