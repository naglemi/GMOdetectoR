# This function allows us to print within dopar and is taken entirely from stackoverflow.
# https://stackoverflow.com/questions/10903787/how-can-i-print-when-using-dopar
Log <- function(text, ...) {
  msg <- sprintf(paste0(as.character(Sys.time()), ": ", text, "\n"), ...)
  cat(msg)
  write.socket(log.socket, msg)
}
