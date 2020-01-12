# Load GMOdetectoR backend, needed for: 1) assign_tray_plate_ID, 2)get_row_col_for_grid_item, ...

library(stringr)

# Need lubridate for getting rid of earlier photos which have problems
library(lubridate)

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
image_files <- list.files(all_images_directory, full.names=TRUE, pattern="GW", recursive = TRUE)

# Pull out metadata for image creation date/time and put this together with
# Path in a data table
image_metadata <- read_exif(image_files)
image_creation_dates <- image_metadata$CreateDate
image_db <- data.table(cbind(image_files, image_creation_dates))

# Parse out timepoint and row_col from file path
image_db$timepoint <- str_split_fixed(image_files, "/", 9)[,8]
image_basename <- basename(file_path_sans_ext(image_files))
image_db$TrayID <- str_split_fixed(image_basename, "_", 2)[,1]
image_db$row_col <- paste0(str_split_fixed(image_basename, "_", 9)[,7],
                           "_",
                           str_split_fixed(image_basename, "_", 9)[,8])
image_basename <- gsub("_rgb", "", image_basename)

# Assign IDs based on row_col
image_db$PlateID <- rep(NA, nrow(image_db))
for(i in 1:length(image_basename)){
  this_ID <- assign_ID_index_from_row_column_on_tray(data_to_parse = image_basename[i], mode="filename")
  image_db$PlateID[i] <- this_ID
}

# Any time there are two images for a specific plate at a specific timepoint, throw out the older one
image_db$tray_plate_timepoint <- paste0(image_db$TrayID, "_",
                                        image_db$PlateID, "_",
                                        image_db$timepoint)

# This should be rewritten into a while loop since there are many tray_plate_timepoint_combo with >2 images
for(tray_plate_timepoint_combo in levels(factor(image_db$tray_plate_timepoint))){
  indices_for_this_tray_plate_timepoint_combo <- which(image_db$tray_plate_timepoint == tray_plate_timepoint_combo)
  if(nrow(image_db[indices_for_this_tray_plate_timepoint_combo,]) > 1){
    print(paste0("Multiple images exist for ", tray_plate_timepoint_combo))
    print("Taken at:")
    print(image_db$image_creation_dates[indices_for_this_tray_plate_timepoint_combo])
    print("The earlier one is: ")
    earlier_time <- image_db$image_creation_dates[indices_for_this_tray_plate_timepoint_combo][1]
    print(earlier_time)
    image_db$tray_plate_timepoint[which(image_db$image_creation_dates == earlier_time & image_db$tray_plate_timepoint == tray_plate_timepoint_combo)] <- NA
    image_db$image_files[which(image_db$image_creation_dates == earlier_time & image_db$tray_plate_timepoint == tray_plate_timepoint_combo)] <- NA
  }
}

# Now to match up randomization datasheets with image file data so we can pick X images
# Actually, not really needed... Don't want the genotype information, just tray and plate ID
# which is all we really need to match image labels up to filenames
# while manually checking labels during scoring

image_db <- na.omit(image_db)
nrow(image_db)

# Shuffle https://stackoverflow.com/questions/6422273/how-to-randomize-or-permute-a-dataframe-rowwise-and-columnwise
image_db_shuffled <- image_db[sample(nrow(image_db)),]

head(image_db)

head(image_db_shuffled)

# Subset to desired number of images
image_db_shuffled_subset <- image_db_shuffled[1:200,]

# Copy these images to a new folder
getwd()
dir.create("/scratch2/NSF_GWAS/macroPhor_Array/all_jpgs_as_of_11.6.19/subset/")
setwd("/scratch2/NSF_GWAS/macroPhor_Array/all_jpgs_as_of_11.6.19/subset")

for(file in image_db_shuffled_subset$image_files){
  file.copy(file, basename(file))
}

fwrite(image_db_shuffled_subset, "Subset_images_to_score.csv", sep=",", quote = FALSE)
