#!/usr/bin/env Rscript
#Intermediate_Script_5
#IS_5
#userinputversion 
#20250121

# Function to check if a package is installed and install if not
install_if_missing <- function(package) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package, dependencies = TRUE, repos = "https://cloud.r-project.org")
    library(package, character.only = TRUE)
  }
}
BiocManager::install("Biostrings")
library(Biostrings)
packages <- c("tidyr", "ggplot2", "gridExtra", "data.table", "dplyr")
lapply(packages, install_if_missing)

#print the current working directory as a check #The original command should be run from the directory containing all input files 
args <- commandArgs(trailingOnly = TRUE)


transdecoder_pep_file <- args[1]
transdecoder_cds_file <- args[2]
blastpunitox6_file <- args[3]
signalpmature_file <- args[4]
signalplist_file <- args[5]
signalpseq_file <- args[6]
Interproscan_file <- args[7]
kallistotrans <- args[8]
sample_name <- args[9]



Interproscan <- read.csv(file = Interproscan_file, header = TRUE, sep = ",")

#kallisto 

#read in transdecoder files
transdecoder_pep <- readAAStringSet(transdecoder_pep_file)
transdecoder_pep_df <- data.frame(
  Transdecoder_ID = names(transdecoder_pep),  # Extracts the sequence names
  Sequence = as.character(transdecoder_pep),   # Converts AAStringSet to a character vector
  Length = width(transdecoder_pep),            # Extracts the sequence lengths (works for AAStringSet)
  stringsAsFactors = FALSE,
  row.names = NULL  # Ensure no row names
)

transdecoder_pep_df$ORF_type <- sub(".*ORF type:([^,]+).*", "\\1", transdecoder_pep_df$Transdecoder_ID)
transdecoder_pep_df$Transdecoder_ID <- sub(" .*", "", transdecoder_pep_df$Transdecoder_ID)
colnames(transdecoder_pep_df)[colnames(transdecoder_pep_df) == "Sequence"] <- "PEP_Sequence"
colnames(transdecoder_pep_df)[colnames(transdecoder_pep_df) == "Length"] <- "PEP_Length"

transdecoder_cds <- readDNAStringSet(transdecoder_cds_file)
transdecoder_cds_df <- data.frame (
  Transdecoder_ID = names(transdecoder_cds),
  Sequence = as.character(transdecoder_cds),
  Length = width(transdecoder_cds),
  stringsAsFactors = FALSE,
  row.names = NULL 
)
transdecoder_cds_df$Transdecoder_ID <- sub(" .*", "", transdecoder_cds_df$Transdecoder_ID)
colnames(transdecoder_cds_df)[colnames(transdecoder_cds_df) == "Sequence"] <- "CDS_Sequence"
colnames(transdecoder_cds_df)[colnames(transdecoder_cds_df) == "Length"] <- "CDS_Length"

#merge transdecoder files 
Transdecoder_pep_cds_merge <- full_join(transdecoder_pep_df, transdecoder_cds_df, by = "Transdecoder_ID")

# readin blasatpunitox6 file 
blastpunitox6 <- read.table(blastpunitox6_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
blastpunitox6 <- blastpunitox6[, c(1, 2, 3, ncol(blastpunitox6)-1, ncol(blastpunitox6))]

colnames(blastpunitox6) <- c("Transdecoder_ID", "Hit", "Percentage_Identity", "E_value", "BitScore") #name columns
blastpunitox6$Hit <- gsub("^sp\\|", "", blastpunitox6$Hit) #remove sp|
blastpunitox6 <- separate(blastpunitox6,  #seperate hit and code column 
                          col = "Hit", 
                          into = c("Code", "Hit"), 
                          sep = "\\|", 
                          extra = "merge", 
                          fill = "right")
blastpunitox6$Hit_species <- blastpunitox6$Hit
blastpunitox6$Hit <- gsub("_.*", "", blastpunitox6$Hit)

#merge transdecoder with blastpunitox6hits 
Transdecoder_blastp <- full_join(Transdecoder_pep_cds_merge, blastpunitox6, by = "Transdecoder_ID")

Transdecoder_blastp <- Transdecoder_blastp[, c(setdiff(names(Transdecoder_blastp), c("PEP_Sequence", "CDS_Sequence")), "PEP_Sequence", "CDS_Sequence")]

#read in signalp files 
signalplist <- read.table(signalplist_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)

signalplist <- signalplist[, c(1, 2)]
colnames(signalplist) <- c("Transdecoder_ID", "SP_Prediction") #name columns

signalpmature <- readAAStringSet(signalpmature_file)
signalpmature_df <- data.frame(
  Transdecoder_ID = names(signalpmature),  # Extracts the sequence names
  Sequence = as.character(signalpmature),   # Converts AAStringSet to a character vector
  Length = width(signalpmature),            # Extracts the sequence lengths (works for AAStringSet)
  stringsAsFactors = FALSE,
  row.names = NULL  # Ensure no row names
)
signalpmature_df$Transdecoder_ID <- sub(" .*", "", signalpmature_df$Transdecoder_ID)
colnames(signalpmature_df)[colnames(signalpmature_df) == "Sequence"] <- "mature_sequence"
colnames(signalpmature_df)[colnames(signalpmature_df) == "Length"] <- "mature_length"

#merge signalplist and signalmature 

signalplist_mature_df <- full_join(signalplist, signalpmature_df, by = "Transdecoder_ID")
#merge
transdecoder_blastp_signalplist_mature <-  full_join(Transdecoder_blastp, signalplist_mature_df, by = "Transdecoder_ID")
transdecoder_blastp_signalplist_mature <- transdecoder_blastp_signalplist_mature[, c(setdiff(names(transdecoder_blastp_signalplist_mature), c("PEP_Sequence", "CDS_Sequence")), "PEP_Sequence", "CDS_Sequence")]
# need to add signal sequence 
signalpsignalseq <- readAAStringSet(signalpseq_file)
signalpmature_df <- data.frame(
  Transdecoder_ID = names(signalpsignalseq),  # Extracts the sequence names
  Signal_Sequence = as.character(signalpsignalseq),   # Converts AAStringSet to a character vector
  Signal_Length = width(signalpsignalseq),            # Extracts the sequence lengths (works for AAStringSet)
  stringsAsFactors = FALSE,
  row.names = NULL  # Ensure no row names
)

signalplist_mature_signal_df <- full_join(signalpmature_df, transdecoder_blastp_signalplist_mature, by = "Transdecoder_ID")

#transdecoder kallisto 
kallistotrans_file <-read.table(file = kallistotrans, sep = ',', header = TRUE)
colnames(kallistotrans_file)[colnames(kallistotrans_file) == "target_id"] <- "Transdecoder_ID"
fin_signalplist_mature_signal_kallisto_df <- full_join(kallistotrans_file, signalplist_mature_signal_df, by = "Transdecoder_ID")

#Cleaning up csv before saving 
keeps <- c("X", "Transdecoder_ID", "ORF_type", "PEP_Length", "CDS_Length", "SP_Prediction", "Signal_Length", "mature_length", "percent", "cumulativepercent", "Code", "Hit", "Percentage_Identity", "E_value", "BitScore", "Hit_species", "Signal_Sequence", "mature_sequence", "PEP_Sequence", "CDS_Sequence")
FINAL_CSV <- fin_signalplist_mature_signal_kallisto_df[keeps]



# View result

FINAL_CSV <- left_join(FINAL_CSV,Interproscan, by="Transdecoder_ID")
#SAVE 
output_path <- file.path("_transdf.csv")
write.csv(FINAL_CSV, output_path, row.names = FALSE)
#this file is to generate the tables for the htmls 

FINAL_CSV_distinct <- FINAL_CSV[order(FINAL_CSV$BitScore, decreasing = TRUE), ]
FINAL_CSV_distinct <- distinct(FINAL_CSV_distinct, Transdecoder_ID, .keep_all = TRUE)

output_path <- file.path("_transdf_distinct.csv")
write.csv(FINAL_CSV_distinct, output_path, row.names = FALSE)
#this file is for the search and download function on R shiny 

#export my annotated putative venom sequences to fasta sequence 

#filter for ORF_type column including the word complete, The SP_prediction column including the word SP and the TMHMM column being FALSE
 #35658
FINAL_CSV_distinct_filtered <- FINAL_CSV_distinct %>%
  filter(
    grepl("complete", ORF_type, ignore.case = TRUE),
    grepl("SP", SP_Prediction, ignore.case = TRUE),
    TMHMM == FALSE
  ) #1127 transcripts vv
FINAL_CSV_distinct_filtered_min <- subset(FINAL_CSV_distinct_filtered, select = c(Transdecoder_ID, CDS_Length,PEP_Length, Hit_species,Percentage_Identity,E_value,BitScore,InterPro_accession_Names,GO_name,Phobius_Name,Panther_ID_Name,PEP_Sequence))
colnames(FINAL_CSV_distinct_filtered)
FINAL_CSV_distinct_filtered_min[is.na(FINAL_CSV_distinct_filtered_min)] <- ""
# Define the fixed sample name (replace with your desired sample name)

# Open the file for writing
fasta_file <- file.path("FINAL_CSV_distinct_filtered_putative_toxins.fasta")
file_conn <- file(fasta_file, open = "w")

# Loop through each row to write custom FASTA format
for (i in 1:nrow(FINAL_CSV_distinct_filtered_min)) {
  row <- FINAL_CSV_distinct_filtered_min[i, ]
  
  # Format lengths and percentage
  cds_len <- paste0(row$CDS_Length, "n")
  pep_len <- paste0(row$PEP_Length, "aa")
  perc_id <- paste0(row$Percentage_Identity, "%")
  
  # Build the description parts, including the fixed sample name
  desc_parts <- c(
    sample_name,  # Add the fixed sample name
    cds_len,
    pep_len,
    row$Hit_species,
    perc_id,
    row$E_value,
    row$BitScore,
    row$InterPro_accession_Names,
    row$GO_name,
    row$Phobius_Name,
    row$Panther_ID_Name
  )
  
  # Combine description parts with " | "
  description <- paste(desc_parts, collapse = " | ")
  
  # Create FASTA header line
  header <- paste0(">", row$Transdecoder_ID, " | ", description)
  
  # Write header and sequence
  writeLines(c(header, row$PEP_Sequence), file_conn)
}

# Close file connection
close(file_conn)
