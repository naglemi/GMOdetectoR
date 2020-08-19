write_sum_stats <- function(job_id, fluorophore_ID_vector, filename, grid_item, tissue_type, cumulative_t_stats){
  sum_stats_path <- paste0("output/",
                           job_id,
                           "/sum_stats.csv")

  if(!file.exists(sum_stats_path)){
    file.create(sum_stats_path)
    header_to_write <- matrix(c("Filename", "Grid_Item", "Tissue_Type", fluorophore_ID_vector))
    write.table(t(header_to_write), file = sum_stats_path,
                sep = ",",
                append = TRUE,
                quote = FALSE,
                col.names = FALSE,
                row.names = FALSE)
  }


  line_to_write <- matrix(c(filename, grid_item, tissue_type, cumulative_t_stats), nrow=1)

  #stop("Break here to debug failure of output to write for certain images")

  write.table(line_to_write,
              file = sum_stats_path,
              sep = ",",
              append = TRUE,
              quote = FALSE,
              col.names = FALSE,
              row.names = FALSE)
}
