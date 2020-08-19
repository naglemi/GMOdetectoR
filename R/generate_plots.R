generate_plots <- function(image_spectrum_table, FullID, fluorophore_ID_vector, name_to_parse, mode_information, grid_type,
                           image_type, job_id, grid_item, cumulative_t_stats){

  #stop("debug generate_plots")

  plot_directory <- paste0("output/",
                           job_id,
                           "/plots/")

  if(!dir.exists(plot_directory)){
    dir.create(plot_directory, recursive = TRUE)
  }

  plot_output <- paste0(plot_directory,
                        file_path_sans_ext(basename(name_to_parse)),
                        "_GridItem",
                        grid_item)

  if (image_type=="hyperspectral"){
    p <- ggplot(data=image_spectrum_table, aes(x=cols, y=rows, fill=color)) +
      coord_equal() + geom_tile() + scale_fill_identity() +
      theme(axis.ticks = element_blank(), axis.title = element_blank(),
            panel.background = element_blank(), axis.text = element_blank(),
            panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_blank()) +
      ggtitle(FullID)
    ggsave(paste0(plot_output, ".png"),
           plot = p)
  }

  if (image_type=="CLS"){
    fluorophore_t_stats_table <- data.table(cbind(fluorophore_ID_vector, cumulative_t_stats))
    myplots <- lapply(X = fluorophore_ID_vector,
                      FUN = plot_data_column,
                      table = image_spectrum_table,
                      filename = name_to_parse,
                      mode_information = "Uncropped",
                      grid_type = grid_type,
                      fluorophore_t_stats_table = fluorophore_t_stats_table)
    #p <- grid.arrange(grobs=myplots)

    # Cannot ggsave results from grid.arrange, but can save to pdf
    # https://stackoverflow.com/questions/17059099/saving-grid-arrange-plot-to-file
    # pdf(paste0(plot_output, ".pdf"))
    # grid.arrange(grobs=myplots)
    # dev.off()
    g <-arrangeGrob(grobs=myplots)
    ggsave(file=paste0(plot_output, ".png"),
           g)
  }

  print(paste0("Saved plot ",
               plot_output))

  return(last_plot())
}
