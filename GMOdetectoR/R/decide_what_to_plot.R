decide_what_to_plot <- function(mode=NA, verbose=FALSE){
  if(verbose==TRUE) print(mode)
  if(is.na(mode)==TRUE){
    stop("Pick a mode")
  }
  if(1 %in% mode){
    to_do <- "whole_plate"
  }
  if(2 %in% mode){
    to_do <- "single_explant"
  }
  if(length(mode) == 2){
    to_do <- "both"
  }
  if(verbose==TRUE) print(to_do)
  return(to_do)
}
