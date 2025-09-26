#!/usr/bin/env Rscript
library(dplyr)
#list files




args <- commandArgs(trailingOnly = TRUE)
Transdf <- args[1]
mass_spec <- args[2]
species_name  <- args[3]
Sample_name <- args[4]

transdf <- read.csv(Transdf, header = TRUE )
mass_spec <- read.csv(mass_spec, header = TRUE)
colnames(mass_spec)[which(names(mass_spec) == "Accession")] <- "Transdecoder_ID"

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

#merge with the mass spec data 
filtered_mass_spec <- left_join(transdf_filtered,mass_spec, by = "Transdecoder_ID")
filtered_select <- filtered_mass_spec[, c(1:32, 36)]
write.csv(filtered_select, paste0(Sample_name, "_filtered_masspec_select.csv"), row.names = FALSE)
#the dataframe for the invidual ones
distinct_mass_spec <- left_join(transdf,mass_spec, by = "Transdecoder_ID")

write.csv(distinct_mass_spec, file=gzfile(paste0(Sample_name, "_distinct_masspec.csv.gz")), row.names = FALSE)
