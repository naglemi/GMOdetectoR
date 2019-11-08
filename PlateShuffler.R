# Load GMOdetectoR backend, needed for: 1) assign_tray_plate_ID, 2)get_row_col_for_grid_item, ...

library(stringr)

source("/scratch2/NSF_GWAS/GMOdetectoR/GMOdetectoRv0.19.R")
#source("C:/Users/OSU/code/GMOdetectoR/GMOdetectoRv0.19.R")

# Need exifr in order to read exif data from images and keep only the most recent one of each
#install.packages("exifr")
library(exifr)

source_directory <- "/scratch2/NSF_GWAS/macroPhor_Array/Randomization_Datasheets"
source_files <- list.files(source_directory, full.names = TRUE)

randomization_database <- data.table(matrix(ncol = ncol(fread(source_files[1]))))
colnames(randomization_database) <- colnames(fread(source_files[1]))
colnames(randomization_database)[1] <- "PlateNumber"

for(file in source_files){
  print(file)
  this_file <- fread(file)[,1:ncol(randomization_database)]
  colnames(this_file) <- colnames(randomization_database)
  randomization_database <- rbind(randomization_database, this_file)
}

all_images_directory <- "/scratch2/NSF_GWAS/macroPhor_Array/all_jpgs_as_of_11.6.19/"
image_files <- list.files(all_images_directory, full.names=TRUE, pattern="", recursive = TRUE)

image_metadata <- read_exif(image_files)
image_creation_dates <- image_metadata$CreateDate
image_db <- data.table(cbind(image_files, image_creation_dates))

image_db$timepoint <- str_split_fixed(image_files, "/", 9)[,8]

image_basename <- basename(file_path_sans_ext(image_files))

image_db$TrayID <- str_split_fixed(image_basename, "_", 2)[,1]

image_db$row_col <- paste0(str_split_fixed(image_basename, "_", 9)[,7],
                           "_",
                           str_split_fixed(image_basename, "_", 9)[,8])

image_basename <- gsub("_rgb", "", image_basename)

image_db$PlateID <- rep(NA, nrow(image_db))
for(i in length(image_basename)){
  this_ID <- assign_ID_index_from_row_column_on_tray(data_to_parse = image_basename[i], mode="filename")
  image_db$PlateID[i] <- this_ID
}



image_files <- basename(tools::file_path_sans_ext(image_files))
