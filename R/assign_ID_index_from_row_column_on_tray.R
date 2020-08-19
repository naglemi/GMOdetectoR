#' Title
#'
#' @param data_to_parse
#' @param components_list
#' @param mode
#' @param verbose
#'
#' @return
#' @export
#'
#' @examples
assign_ID_index_from_row_column_on_tray <- function(data_to_parse = filename, components_list, mode="table", verbose=FALSE){
  dictionary <- cbind(c(0,0,0,0,0,0,0,
                        1,1,1,1,1,1,1,
                        2,2,2,2,2,2,2),
                      c(0,1,2,3,4,5,6,
                        0,1,2,3,4,5,6,
                        0,1,2,3,4,5,6),
                      c(1:21))
  dictionary <- as.data.table(dictionary)
  colnames(dictionary) <- c("row", "column", "ID")
  # Get the ID of position in tray in according to row and column
  dictionary$row_column <- paste0(dictionary$row, "_", dictionary$column)
  dictionary[,1:2] <- NULL
  if(mode=="table"){
    # Set colnames for spectral components if multiple are same
    #colnames(data_to_parse)[1:length(components_list)] <- components_list
    data_merged <- merge(data_to_parse, dictionary, by="row_column", all.x = TRUE, all.y = TRUE)
    return(data_merged)
  }
  if(mode=="filename"){

    # Patch added in v0.19 for compatibility regardless of whether "_cyan" is at end of filename
    if(grepl("cyan", data_to_parse)==1){
      ndelimiters=9
    }else{
      ndelimiters=8
    }

    row <- str_split_fixed(basename(file_path_sans_ext(data_to_parse)), "_", 9)[ndelimiters-1]
    # Changed in v0.19 along with patch above
    col <- str_split_fixed(basename(file_path_sans_ext(data_to_parse)), "_", ndelimiters)[ndelimiters]
    row_col <- paste0(row, "_", col)
    ID <- dictionary[which(dictionary$row_column == row_col),]$ID
    # Debugging lines added in v0.19
    if(verbose==TRUE){
      print(paste0("This row is ", row))
      print(paste0("This col is ", col))
      print(paste0("This row_col is ", row_col))
      print(paste0("This filename (stripped) is ", basename(file_path_sans_ext(data_to_parse))))
      print(paste0("This ID about to be returned from assign_ID_index_from_roW_column_on_tray is ", ID))
    }

    return(ID)
  }
}
