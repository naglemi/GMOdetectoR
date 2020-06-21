check_success_high_throughput <- function(files_to_loop, sum_stats_file, grid_type, pattern_to_exclude = "chroma", by_explant=TRUE){

  # This line assumes that all chroma standard files contain "chroma" unless pattern_to_exclude option is
  # modified by user. These will be excluded.

  files_to_loop <- files_to_loop[!grepl(pattern_to_exclude, files_to_loop)]

  # We will create a dt of all possible grid items for the specified grid type, then compare this to the output
  # so we can see if anything is missing.
  sum_stats_file_in <- fread(sum_stats_file)

  if(by_explant==TRUE){
    colnames(sum_stats_file_in)[1:2] <- c("file", "grid_item")
    comprehensive_dt <- data.table(file = c(rep(files_to_loop, times = 1, each = grid_type)),
                                   grid_item = c(rep(1:grid_type, times = length(files_to_loop))))
    comprehensive_dt_with_results <- merge(sum_stats_file_in,
                                           comprehensive_dt, by=c("file", "grid_item"),
                                           all.x = FALSE,
                                           all.y = TRUE)
  }else{
    colnames(sum_stats_file_in)[1] <- c("file")
    comprehensive_dt <- data.table(file = c(files_to_loop))
    comprehensive_dt_with_results <- merge(sum_stats_file_in,
                                           comprehensive_dt[!duplicated(comprehensive_dt[,c('file')]),], by=c("file"),
                                           all.x = FALSE,
                                           all.y = TRUE)
  }

  return(comprehensive_dt_with_results)
}
