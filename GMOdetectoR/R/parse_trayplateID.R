#' Given macroPhor Array output filename, parse out tray and plate IDs
#'
#' @param filename for macroPhor Array output, with file naming as used in Strauss Lab
#'
#' @return A character string with tray ID and plate ID delimited by "_"
#' @export
#'
#' @examples
parse_trayplateID <- function(name_being_parsed){
  pass_to_dodge_error <- name_being_parsed
  imgpath_stripped <- file_path_sans_ext(basename(pass_to_dodge_error))
  trayID <- str_split_fixed(imgpath_stripped, "_", 2)[1]
  plateID <- assign_ID_index_from_row_column_on_tray(data_to_parse = imgpath_stripped, mode="filename")
  trayplateID <- paste0(trayID, "_", plateID)
  return(trayplateID)
}
