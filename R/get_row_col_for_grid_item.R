#' Find the row and column for a specific grid item in sample
#'
#' @param grid_item_in An integer representing the specific item for which the row and column in grid are desired
#' @param grid_type An integer representing the type of grid, for grid types already in \code{GMOdetectoR/grids/}
#' @param grid_file A path to a grid file if using a grid type not already in \code{GMOdetectoR/grids/}
#'
#' @return A list with values for the row and column of the grid item
#' @export
#'
#' @examples \dontrun{get_row_col_for_grid_item(grid_item = grid_item, grid_type = grid_type)}
get_row_col_for_grid_item <- function(grid_item_in, grid_type, grid_file){

  if(grid_type==20){
    grid_file <- fread("sections20.txt")
  }
  if(grid_type==12){
    grid_file <- fread("sections12.txt")
  }
  # Column name can't be same as variable name for this to work. Hence naming of variable as grid_item_in
  this_grid_item_position <- grid_file[which(grid_file$grid_item == grid_item_in),]
  row <- this_grid_item_position$row
  col <- this_grid_item_position$col
  return(c(row, col))
}

