#tempdir("/scratch2/NSF_GWAS/Rtmp/")

library(GMOdetectoR)
message(paste0("Running GMOdetectoR version ", packageVersion('GMOdetectoR')))

#!/usr/bin/env Rscript
library("optparse")

option_list = list(
  make_option(c("-F", "--folder"), type="character", default=NULL,
              help="dataset folder containing .raw and .hdr files, for one timepoint", metavar="character"),
  make_option(c("-g", "--grid_type"), type="numeric", default=12,
              help="output file name [default= %default]", metavar="numeric"),
  make_option(c("-m", "--max_cores"), type="numeric", default=20,
              help="maximum number of CPU cores for parallelization", metavar="numeric"),
  make_option(c("-r", "--record_residuals"), type="numeric", default=0,
              help="whether to record and plot residuals for each sample (plate or grid item, depending on by_explant)", metavar="numeric"),
  make_option(c("-p", "--plotting"), type="numeric", default=0,
              help="whether to make plots of regression outputs", metavar="numeric"),
  make_option(c("-e", "--explant_level"), type="numeric", default=0,
              help="whether to produce outputs for grid items (explants) individually", metavar="numeric"),
  make_option(c("-f", "--fluorophore_list"), type="character", default="fluorophore_list.txt",
              help="Path to a text file with name of one fluorophore per line (must be in spectra_library folder)", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#what_to_plot <- decide_what_to_plot(mode = c(2))
print(opt)
print(opt$max_cores)
run_parallel(timepoint = opt$folder,
             maximum_CPU_number = opt$max_cores[1],
             intercept = 1,
             FP_threshold = 0,
             Chl_threshold = 0,
             by_explant = opt$explant_level,
             grid_type = opt$grid_type,
             record_residuals = opt$record_residuals,
             plotting = opt$plotting,
             fluorophore_ID_vector = as.vector(unlist(read.csv(opt$fluorophore_list, header = FALSE))))
