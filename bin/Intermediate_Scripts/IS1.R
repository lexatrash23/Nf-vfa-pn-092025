#!/usr/bin/env Rscript
#Intermediate_Script_1
#IS1

library(Biostrings)
library(tidyr)
library(dplyr)


# variables defined by command line arguments
args <- commandArgs(trailingOnly = TRUE)
Trinity_file <- args[1] #Trinity fasta file
blastxunitox6_file <- args[2] #output version 6 with following fields "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qframe qcovs"
kallisto_file <- args[3] #Kallisto csv file with following columns "rownames(unnamed), target_id, length, eff_length, est_counts, tpm, percent, cumulativepercent"
sample_name <- args[4] #basename


#loading in the trinity file as a DNAstringset object
Trinity_Fasta_set_sequences <- readDNAStringSet(Trinity_file)

#making a dataframe out of this object by getting the name, sequence and length
Trinity_Fasta_df <- data.frame (
  Trinity_ID = names(Trinity_Fasta_set_sequences),
  Sequence = as.character(Trinity_Fasta_set_sequences),
  Length = width(Trinity_Fasta_set_sequences),
  stringsAsFactors = FALSE)
#Edit name of the Trinity IDs to remove everything after .
Trinity_Fasta_df$Trinity_ID <- sub(" .*", "", Trinity_Fasta_df$Trinity_ID)



#load in the unitox txt file
Unitox_blastx_6 <- read.table(blastxunitox6_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
#name the columns of the file
colnames(Unitox_blastx_6) <- c("Trinity_ID", "Hit", "Percentage_Identity", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "ssend" "E_value", "BitScore", "Frame", "query_coverage") #name columns

#Default hit column is sp|Code|Hit_Species, this separates that to make three separate columns one with just the Code, one with the Hit and one that has Hit_Species
#this removes the sp| that precedes each Hit
Unitox_blastx_6$Hit <- gsub("^sp\\|", "", Unitox_blastx_6$Hit) #remove sp|
#This sep the hit and the code column by the | separating them. so now the Code column has the code and the Hit column has the Hit_Species.
Unitox_blastx_6 <- separate(Unitox_blastx_6,
                            col = "Hit", 
                            into = c("Code", "Hit"), 
                            sep = "\\|", 
                            extra = "merge", 
                            fill = "right")

#Creating a new column called Hit_species and this will have the Hit_species values
Unitox_blastx_6$Hit_species <- Unitox_blastx_6$Hit
#For the original Hit column we want to remove the species so we are removing anything after the "_.*)
Unitox_blastx_6$Hit <- gsub("_.*", "", Unitox_blastx_6$Hit)
#for the merged data, we are going to drop some of the blastx6 columns to reduce file size
Unitox_blastx_6 <- select(Unitox_blastx_6, -length, -mismatch, -gapopen, -qstart, -qend, -sstart, -ssend)


##Kallisto data
#Reading in kallisto datasets
Trinity_Kallisto_data <- read.csv(kallisto_file, stringsAsFactors = FALSE)
#rename the target_id column in kallisto so that we can merge by column name in the next step
colnames(Trinity_Kallisto_data)[colnames(Trinity_Kallisto_data) == 'target_id'] <- 'Trinity_ID'
#left_join the Trinity kallisto and trinity one. left_join retains all the rows of the trinity_fasta and just adds in the column data from the kallisto where it is available
Trinity_Kallisto_merge_all <- left_join(Trinity_Fasta_df, Trinity_Kallisto_data, by = "Trinity_ID")
#left join the trinity_kallisto merge. merge with all.x=TRUE is basically a left join, i just wanted to try another way lol.
combined_kallisto_unitox_trinity <- merge(Trinity_Kallisto_merge_all, Unitox_blastx_6, by = "Trinity_ID", all.x = TRUE)


#now we save this dataframe as a csv file
write.csv(combined_kallisto_unitox_trinity, paste0(basename + "_TBK.csv"), row.names = FALSE)


#Distinct list where each trinity ID is only present once, retaining the result with the highest blast bitscore if any
combined_kallisto_unitox_trinity_ordered <- combined_kallisto_unitox_trinity[order(combined_kallisto_unitox_trinity$BitScore, decreasing = TRUE), ]
combined_kallisto_unitox_trinity_distinct <- distinct(combined_kallisto_unitox_trinity_ordered,Trinity_ID, .keep_all = TRUE )


#this is saved as a .csv.gz, this file is used in the R shiny app to parse the dataframe. zipped and filtered for distinct to reduce file size

write.csv(combined_kallisto_unitox_trinity_distinct, file=gzfile(paste0(sample_name + "_TBK_distinct.csv.gz")), row.names = FALSE)
