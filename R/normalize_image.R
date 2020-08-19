#' Normalize a reduced image
#'
#' @param image_spectrum_table an image spectrum table processed by \code{reduce_image}
#' @param chroma_table an image spectrum table for a standard, also processed by \code{reduce_image}
#'
#' @return an image spectrum table which has been normalized according to the standard
#' @export
#'
#' @examples \dontrun{image_spectrum_table <- normalize_image(image_spectrum_table, chroma_table)}
normalize_image <- function(image_spectrum_table, chroma_table){

  image_spectrum_table <- image_spectrum_table[image_spectrum_table$rows <= max(chroma_table$rows),]

  chroma_table <- chroma_table[chroma_table$rows <= max(image_spectrum_table$rows),]

  image_spectrum_table$r <- as.numeric(as.character(image_spectrum_table$r))
  image_spectrum_table$g <- as.numeric(as.character(image_spectrum_table$g))

  image_spectrum_table$r <- image_spectrum_table$r / chroma_table$r
  image_spectrum_table$g <- image_spectrum_table$g / chroma_table$g

  return(image_spectrum_table)
}
