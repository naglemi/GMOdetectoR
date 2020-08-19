write_tables <- function(table_subfolder, table, filename, grid_item, job_id, fluorophore_ID_vector){
  table_output_folder_path <- paste0("output/",
                                     job_id,
                                     "/",
                                     table_subfolder,
                                     "/")

  if(!dir.exists(table_output_folder_path)){
    dir.create(table_output_folder_path, recursive = TRUE)
  }

  fwrite(na.omit(table[,1:(2+length(fluorophore_ID_vector))]),
              file = paste0(paste0(table_output_folder_path,
                                   basename(file_path_sans_ext(filename)),
                                   "_GridItem",
                                   grid_item,
                                   ".csv")),
              sep = ",",
              quote = FALSE,
              col.names = TRUE,
              row.names = FALSE)
}
