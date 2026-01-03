#!/usr/bin/env Rscript
library(dplyr)



#Allow for command line arguments
args <- commandArgs(trailingOnly = TRUE)
Transdf <- args[1] #transdf_distinct
mass_spec <- args[2] #massspec csv 
species_name  <- args[3] #species name
Sample_name <- args[4] #sample name 


#read in files 
transdf <- read.csv(Transdf, header = TRUE )
mass_spec <- read.csv(mass_spec, header = TRUE)
# rename Accession column for later left join
colnames(mass_spec)[which(names(mass_spec) == "Accession")] <- "Transdecoder_ID" 


#add two columns to the start. One species and one sample name
transdf$Species <- species_name
transdf$Sample_name <- Sample_name

#move those to the end
transdf <- transdf[c("Species", "Sample_name", setdiff(names(transdf), c("Species", "Sample_name")))]


#Fulldf # all transdf distinct all mass spec columns 
#merge with the mass spec data
transdf_massspec <- left_join(transdf_filtered,mass_spec, by = "Transdecoder_ID")

write.csv(transdf_massspec, paste0(Sample_name, "_transdf_distinct_masspec_full.csv"), row.names = FALSE)


#Simplified df #only complete and signalp transdf distinct and select mass spec columns 
#filter for those that are complete + SP(Sec/SPI) + TMHMM is false 
transdf_filtered <- transdf[
  (transdf$ORF_type == "complete (+)" | transdf$ORF_type == "complete (-)") &
    transdf$SP_Prediction == "SP(Sec/SPI)" &
    transdf$TMHMM == "FALSE",
]



massspec_select <- mass_spec %>% dplyr::select(Top, Transdecoder_ID, X.10LgP, Coverage..., X.Peptides, X.Unique, PTM )
transdf_massspec_select <- left_join(transdf_filtered,massspec_select, by = "Transdecoder_ID")
write.csv(transdf_massspec_select, paste0(Sample_name, "_transdf_distinct_masspec_simplified.csv"), row.names = FALSE)



