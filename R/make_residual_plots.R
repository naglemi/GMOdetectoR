make_residual_plots <- function(CLS_table, t_spectrum_in_to_CLS, spectrum_in_to_CLS, bins, job_id, filename, grid_item){
  #browser()
  print("Now inside make_residual_plots")
  print("Properties of CLS table (colnames and dimensions):")
  print(colnames(CLS_table))
  print(dim(CLS_table))
  print("Wavelengths being included (# wavelengths, head, and tail)")
  the_wavelengths <- rownames(t_spectrum_in_to_CLS)
  print(length(the_wavelengths))
  print(head(the_wavelengths))
  print(tail(the_wavelengths))
  #browser()
  molten_residual_table <- melt(data = CLS_table,
                                id.vars = c("rows", "cols"),
                                measure.vars = rownames(t_spectrum_in_to_CLS))

  print("Melted residual table, now to melt y_obs")
  molten_y_obs <- melt(data = spectrum_in_to_CLS,
                       id.vars = c("rows", "cols"),
                       measure.vars = rownames(t_spectrum_in_to_CLS))

  print("Done melting")
  colnames(molten_residual_table)[3:4] <- c("wavelength", "residual")
  colnames(molten_y_obs)[3:4] <- c("wavelength", "observed")

  molten_combined <- cbind(molten_residual_table, molten_y_obs[,4])

  molten_residual_table <- molten_y_obs <- NULL

  molten_combined$predicted <- molten_combined$observed - molten_combined$residual

  molten_combined$residual_normalized <- (molten_combined$residual-mean(molten_combined$residual))/sd(molten_combined$residual)

  print("Done calculating normalized residuals")

  molten_combined <- molten_combined[!which(molten_combined$predicted==0)]

  molten_combined$wavelength <- as.numeric(as.character(molten_combined$wavelength))

  table_output_folder_path <- paste0("output/",
                                     job_id,
                                     "/plots/")

  if(!dir.exists(table_output_folder_path)){
    dir.create(table_output_folder_path, recursive = TRUE)
  }

  print("Making plot")
  p <- ggplot(molten_combined, aes(x=predicted, y=residual_normalized)) +
    geom_bin2d(bins=200)

  print(paste0("Saving residual plots ", basename(file_path_sans_ext(filename)),
                                                "GridItem",
                                                grid_item))

  ggsave(paste0(table_output_folder_path,
                paste0(basename(file_path_sans_ext(filename)),
                       "GridItem",
                       grid_item,
                       "_res1.png")),
         p,
         width = 10,
         height = 3)

  p <- NULL
  gc()

  p <- ggplot(molten_combined, aes(x=wavelength, y=residual_normalized)) +
    geom_bin2d(bins=200) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

  ggsave(paste0(table_output_folder_path,
                paste0(basename(file_path_sans_ext(filename)),
                       "GridItem",
                       grid_item,
                       "_res2.png")),
         p,
         width = 10,
         height = 3)
}
