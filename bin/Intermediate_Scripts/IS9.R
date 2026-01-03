#!/usr/bin/env Rscript
library(dplyr)


#Allow for command line arguments
args <- commandArgs(trailingOnly = TRUE)
Transdf <- args[1] #transdf_distinct
species_name  <- args[2] #species name
Sample_name <- args[3] #sample name 


transdf <- read.csv(Transdf, header = TRUE )


#add two columns to the start. One species and one sample name
transdf$Species <- species_name
transdf$Sample_name <- Sample_name

#move those to the end
transdf <- transdf[c("Species", "Sample_name", setdiff(names(transdf), c("Species", "Sample_name")))]



#Fulldf # all transdf distinct no mass spec 
write.csv(transdf, paste0(Sample_name, "_transdf_distinct_nomasspec_full.csv"), row.names = FALSE)

#Simplified df #only complete and signalp transdf distinct and select mass spec columns 

#filter for those that are complete + SP(Sec/SPI) + TMHMM is false 
transdf_filtered <- transdf[
  (transdf$ORF_type == "complete (+)" | transdf$ORF_type == "complete (-)") &
    transdf$SP_Prediction == "SP(Sec/SPI)" &
    transdf$TMHMM == "FALSE",
]


write.csv(transdf_filtered, paste0(Sample_name, "_transdf_distinct_nomasspec_simplified.csv"), row.names = FALSE)

