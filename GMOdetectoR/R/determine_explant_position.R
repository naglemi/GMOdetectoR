determine_explant_position <- function(grid_type,
                                       grid_item,
                                       left_edge = 14,
                                       right_edge = 1406,
                                       bottom_edge = 1256,
                                       top_edge = 226){
  y_newsize <- right_edge - left_edge
  x_newsize <- bottom_edge - top_edge

  # Calculate size of each grid item
  if(grid_type==20){
    y_grid_rows <- 4
    x_grid_columns <- 6
  }
  if(grid_type==12){
    y_grid_rows <- 3
    x_grid_columns <- 4
  }
  if(grid_type!= 12 & grid_type!= 20){
    if(is.na(grid_file_path)){
      stop("Error: Need to provide grid file if not using standard 12 or 20 section grid")
    }
    # Finish up this functionality later
    grid_file <- fread(grid_file_path)
    y_grid_rows <- max(grid_file$row)
    x_grid_columns <- max(grid_file$col)
  }

  x_grid_item_size <- x_newsize/y_grid_rows
  y_grid_item_size <- y_newsize/x_grid_columns

  # Find row and column of desired grid item
  row_col <- get_row_col_for_grid_item(grid_item = grid_item, grid_type = grid_type)
  row <- row_col[1]
  col <- row_col[2]

  # Crop to desired grid item
  #stop("debug here")
  cropped_bottom_edge <- (bottom_edge - ((row-1)*x_grid_item_size))
  cropped_top_edge <- (bottom_edge - (row*x_grid_item_size))
  cropped_left_edge <- (left_edge + ((col-1)*y_grid_item_size))
  cropped_right_edge <- (left_edge + (col*y_grid_item_size))

  return(list(cropped_bottom_edge,
              cropped_top_edge,
              cropped_left_edge,
              cropped_right_edge))
}
