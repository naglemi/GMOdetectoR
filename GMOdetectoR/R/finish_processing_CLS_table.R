#' Title
#'
#' @param CLS_table
#' @param chosen_pixels_in
#' @param fluorophore_ID_vector
#' @param filename
#' @param grid_item
#' @param tissue_type
#' @param job_id
#'
#' @return
#' @export
#'
#' @examples
finish_processing_CLS_table <- function(CLS_table, chosen_pixels_in = chosen_pixels, fluorophore_ID_vector,
                                        filename, grid_item, tissue_type, job_id){

  write_tables(table_subfolder = "CLS_tables",
               table = CLS_table,
               filename = filename)

  write_sum_stats(job_id = job_id,
                  fluorophore_ID_vector = fluorophore_ID_vector,
                  filename = filename,
                  grid_item = grid_item,
                  tissue_type = tissue_type,
                  cumulative_t_stats = cumulative_t_stats)

  CLS_table[CLS_table < 0] <- 0

  return(CLS_table)
}
