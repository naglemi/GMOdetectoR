# Load GMOdetectoR backend, needed for: 1) assign_tray_plate_ID, 2)get_row_col_for_grid_item, ...

source("C:/Users/OSU/code/GMOdetectoR/GMOdetectoRv0.19.R")

source_directory <- "Z:/Groups/Tgerc/RESEARCH PROJECTS/NSF project/In Vitro GWAS study/Randomization_Datasheets"
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

all_images_directory <- "J:/Downloads/copy_all_jpgs_on_Box_as_of_11.5.19_sync/"
image_files <- list.files(all_images_directory, full.names=TRUE)

image_files <- basename(tools::file_path_sans_ext(image_files))
