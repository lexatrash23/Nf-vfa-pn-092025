#!/usr/bin/env Rscript
library(dplyr)
args <- commandArgs(trailingOnly = TRUE)
#list files 

args <- commandArgs(trailingOnly = TRUE)
Transdf <- args[1]
species_name  <- args[2]
Sample_name <- args[3]

transdf <- read.csv(Transdf, header = TRUE )

#filter for those that are complete + SP(Sec/SPI) + TMHMM is false 
transdf_filtered <- transdf[
  (transdf$ORF_type == "complete (+)" | transdf$ORF_type == "complete (-)") &
    transdf$SP_Prediction == "SP(Sec/SPI)" &
    transdf$TMHMM == "FALSE",
]
#add two columns to the start. One species and one sample name 
transdf_filtered$Species <- species_name
transdf_filtered$Sample_name <- Sample_name

transdf_filtered <- transdf_filtered[c("Species", "Sample_name", setdiff(names(transdf_filtered), c("Species", "Sample_name")))]


write.csv(transdf_filtered, paste0(Sample_name, "_filtered_nomasspec_.csv"), row.names = FALSE)

#the dataframe for the invidual ones 

write.csv(transdf, file=gzfile(paste0(Sample_name, "_distinct_nomasspec.csv.gz")), row.names = FALSE)
