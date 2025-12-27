#!/usr/bin/env Rscript
library(dplyr)



#Allow for command line arguments
args <- commandArgs(trailingOnly = TRUE)
Transdf <- args[1] #transdf_distinct
mass_spec <- args[2]
species_name  <- args[3]
Sample_name <- args[4]


#read in transdf
transdf <- read.csv(Transdf, header = TRUE )
mass_spec <- read.csv(mass_spec, header = TRUE)
colnames(mass_spec)[which(names(mass_spec) == "Accession")] <- "Transdecoder_ID" # rename Accession column for later left join

#filter for those that are complete + SP(Sec/SPI) + TMHMM is false 
transdf_filtered <- transdf[
  (transdf$ORF_type == "complete (+)" | transdf$ORF_type == "complete (-)") &
    transdf$SP_Prediction == "SP(Sec/SPI)" &
    transdf$TMHMM == "FALSE",
]



#add two columns to the start. One species and one sample name
transdf_filtered$Species <- species_name
transdf_filtered$Sample_name <- Sample_name

#move those to the end
transdf_filtered <- transdf_filtered[c("Species", "Sample_name", setdiff(names(transdf_filtered), c("Species", "Sample_name")))]

#merge with the mass spec data
transdf_massspec <- left_join(transdf,mass_spec, by = "Transdecoder_ID")


#for the final R app for species with no genome but with mass spec data.gz to reduce file size
write.csv(transdf_massspec, file=gzfile(paste0(Sample_name, "_transdf_distinct_masspec.csv.gz")), row.names = FALSE)

write.csv(transdf_massspec, paste0(Sample_name, "_transdf_distinct_masspec.csv"), row.names = FALSE)
filtered_mass_spec <- left_join(transdf_filtered,mass_spec, by = "Transdecoder_ID")
filtered_select <- filtered_mass_spec[, c(1:32, 36)]
write.csv(filtered_select, paste0(Sample_name, "_filtered_masspec_select.csv"), row.names = FALSE)

