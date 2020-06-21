crop <- function(image_being_cropped,
                 mode,
                 grid_file = NA,
                 grid_type,
                 grid_item,
                 image_type,
                 left_edge = 14,
                 right_edge = 1406,
                 bottom_edge = 1256,
                 top_edge = 226,
                 desired_wavelength_range){

  #browser()

  if(image_type=="CLS"){
    image_being_cropped <- hyperSpec_to_IST(hyperSpec_object = image_being_cropped,
                                            desired_wavelength_range = desired_wavelength_range)
  }

  if(mode=="single_explant"){

    if(grid_item>grid_type){
      stop(paste0("Grid item #", grid_item, " not found in grid type #", grid_type))
    }

    explant_position <- determine_explant_position(grid_type = grid_type,
                                                   grid_item = grid_item,
                                                   left_edge = left_edge,
                                                   right_edge = right_edge,
                                                   bottom_edge = bottom_edge,
                                                   top_edge = top_edge)

    bottom_edge <- explant_position[[1]]
    top_edge <- explant_position[[2]]
    left_edge <- explant_position[[3]]
    right_edge <- explant_position[[4]]
  }

  image_being_cropped <- image_being_cropped[which(image_being_cropped$cols < bottom_edge & image_being_cropped$cols > top_edge & image_being_cropped$rows>left_edge & image_being_cropped$rows<right_edge),]
  return(image_being_cropped)
}
