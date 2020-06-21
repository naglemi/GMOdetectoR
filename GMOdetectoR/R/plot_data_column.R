# Credit for below function: https://stackoverflow.com/questions/31993704/storing-ggplot-objects-in-a-list-from-within-loop-in-r
plot_data_column <- function (fluorophore_ID,
                              table,
                              filename,
                              mode_information,
                              grid_type,
                              fluorophore_t_stats_table) {

  #stop("debug plot_table_column")

  cumulative_t_stat <- fluorophore_t_stats_table[which(fluorophore_t_stats_table[,1] == fluorophore_ID), 2]

  colfunc <- colorRampPalette(c("white", "blue"))

  plot_out <- ggplot(data=table)+
    # Fill is different depending on if we are in p-val or z-score plotting mode
    # Note redundancy with explicitly calling table
    geom_tile(aes(y=rows, x=cols, fill=table[,get(fluorophore_ID)]))+ #https://stackoverflow.com/questions/32184252/how-to-select-fluorophore_IDs-in-table-table-using-a-character-vector-of-certain-fluorophore_ID-n
    scale_fill_gradientn(colours = colfunc(10))+
    theme(axis.ticks = element_blank(),
          axis.title = element_blank(),
          panel.background = element_blank(),
          axis.text = element_blank(),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_blank())+
    labs(title = paste0(fluorophore_ID, ": ",
                        format(as.numeric(cumulative_t_stat)/1000,
                                          digits=1)),
         subtitle = expression('x' ~10^3 ~ '=' ~ Sigma^'pixels' ~ Tau['FP signal']),
         fill = paste0("T(", fluorophore_ID,")"))

  return(plot_out)}
