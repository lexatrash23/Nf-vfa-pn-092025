#!/usr/bin/env Rscript
#Intermediate_Script_5
#IS_5

#library
library(Biostrings)
library(dplyr)
library(tidyr)

#allow for command line arguments
args <- commandArgs(trailingOnly = TRUE)


transdecoder_pep_file <- args[1] #Transdecoder pep output raw
transdecoder_cds_file <- args[2] #Transdecoder cds output raw
blastpunitox6_file <- args[3]    #Blastp6 output raw
signalpmature_file <- args[4]    #SignalP output raw mature fasta raw
signalplist_file <- args[5]      #SignalP output raw summary file raw
signalpseq_file <- args[6]       #SignalSequence output fasta file from IS2 (Proces 2)
Interproscan_file <- args[7]     #Interproscan file output from IS4 (Process 5)
kallistotrans <- args[8]         #Kallisto transdecoder csv output from kallistoanalysistrans.py (Process 2)
sample_name <- args[9]          #sample name



#read in transdecoder pep files
transdecoder_pep <- readAAStringSet(transdecoder_pep_file)
transdecoder_pep_df <- data.frame(
  Transdecoder_ID = names(transdecoder_pep),  # save sequence_id
  Sequence = as.character(transdecoder_pep),  # save sequence to as character
  Length = width(transdecoder_pep),            # save sequence length
  stringsAsFactors = FALSE,
  row.names = NULL  # ensures no row names
)

#add ORF type column to df
transdecoder_pep_df$ORF_type <- sub(".*ORF type:([^,]+).*", "\\1", transdecoder_pep_df$Transdecoder_ID) #searches the transdecoderID captures [^<,]+ any characters are the 'ORF_type" and before the first comma after this comma as the first capturing groups, saving it in the ORF_type columm.
#rename the TransdecoderID  removing everything after the first space
transdecoder_pep_df$Transdecoder_ID <- sub(" .*", "", transdecoder_pep_df$Transdecoder_ID)
#rename sequence and length columns to be pep specific
colnames(transdecoder_pep_df)[colnames(transdecoder_pep_df) == "Sequence"] <- "PEP_Sequence"
colnames(transdecoder_pep_df)[colnames(transdecoder_pep_df) == "Length"] <- "PEP_Length"

#read in transdecoder cds file
transdecoder_cds <- readDNAStringSet(transdecoder_cds_file)
transdecoder_cds_df <- data.frame (
  Transdecoder_ID = names(transdecoder_cds),
  Sequence = as.character(transdecoder_cds),
  Length = width(transdecoder_cds),
  stringsAsFactors = FALSE,
  row.names = NULL 
)

#rename TransdecoderID to remove everything after the first space
transdecoder_cds_df$Transdecoder_ID <- sub(" .*", "", transdecoder_cds_df$Transdecoder_ID)
#rename sequence and length columns to be cds specific
colnames(transdecoder_cds_df)[colnames(transdecoder_cds_df) == "Sequence"] <- "CDS_Sequence"
colnames(transdecoder_cds_df)[colnames(transdecoder_cds_df) == "Length"] <- "CDS_Length"
qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qframe qcovs
#merge transdecoder files full join by ID column
Transdecoder_pep_cds_merge <- full_join(transdecoder_pep_df, transdecoder_cds_df, by = "Transdecoder_ID")

# readin blasatpunitox6 file 
blastpunitox6 <- read.table(blastpunitox6_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(blastpunitox6) <- c("Transdecoder_ID", "Hit", "Percentage_Identity", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "ssend" "E_value", "BitScore", "Frame", "query_coverage") #name columns
keeps <- c("Transdecoder_ID", "Hit", "Percentage_Identity","E_value", "BitScore")
blastpunitox6 <- blastpunitox6[keeps]
#similar to IS1, separating hit to hit, code, hit_species
blastpunitox6$Hit <- gsub("^sp\\|", "", blastpunitox6$Hit) #remove sp|
blastpunitox6 <- separate(blastpunitox6,  #seperate hit and code column 
                          col = "Hit", 
                          into = c("Code", "Hit"), 
                          sep = "\\|", 
                          extra = "merge", 
                          fill = "right")
blastpunitox6$Hit_species <- blastpunitox6$Hit
blastpunitox6$Hit <- gsub("_.*", "", blastpunitox6$Hit)

#merge transdecoder with blastpunitox6hits full join by ID
Transdecoder_blastp <- full_join(Transdecoder_pep_cds_merge, blastpunitox6, by = "Transdecoder_ID")

#moving the sequences columns right to the end
Transdecoder_blastp <- Transdecoder_blastp[, c(setdiff(names(Transdecoder_blastp), c("PEP_Sequence", "CDS_Sequence")), "PEP_Sequence", "CDS_Sequence")]

#read in signalp files 
signalplist <- read.table(signalplist_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
#just interested in the first two columns
signalplist <- signalplist[, c(1, 2)]
#name the two columns
colnames(signalplist) <- c("Transdecoder_ID", "SP_Prediction")


#read in signal p mature fasta file
signalpmature <- readAAStringSet(signalpmature_file)
#save data as a df
signalpmature_df <- data.frame(
  Transdecoder_ID = names(signalpmature),  # sequence names
  Sequence = as.character(signalpmature),   # sequence
  Length = width(signalpmature),            # sequence length
  stringsAsFactors = FALSE,
  row.names = NULL
)

#there should not be anything after the space, but just in case
signalpmature_df$Transdecoder_ID <- sub(" .*", "", signalpmature_df$Transdecoder_ID)
#rename sequence and length to be mature sequence related
colnames(signalpmature_df)[colnames(signalpmature_df) == "Sequence"] <- "mature_sequence"
colnames(signalpmature_df)[colnames(signalpmature_df) == "Length"] <- "mature_length"

#merge signalplist and signalmature full join
signalplist_mature_df <- full_join(signalplist, signalpmature_df, by = "Transdecoder_ID")

#merge signalp mature data with previous df
transdecoder_blastp_signalplist_mature <-  full_join(Transdecoder_blastp, signalplist_mature_df, by = "Transdecoder_ID")

# Add signal sequences
signalpsignalseq <- readAAStringSet(signalpseq_file)
signalpsignalseq_dfs <- data.frame(
  Transdecoder_ID = names(signalpsignalseq),  # ID
  Signal_Sequence = as.character(signalpsignalseq),   # Sequence
  Signal_Length = width(signalpsignalseq),            # Length
  stringsAsFactors = FALSE,
  row.names = NULL
)

#full join merge to add Signal Sequence Data
signalplist_mature_signal_df <- full_join(signalpsignalseq_dfs, transdecoder_blastp_signalplist_mature, by = "Transdecoder_ID")


#Add in kallisto data
kallistotrans_file <-read.table(file = kallistotrans, sep = ',', header = TRUE)
#rename column
colnames(kallistotrans_file)[colnames(kallistotrans_file) == "target_id"] <- "Transdecoder_ID"
#full join
fin_signalplist_mature_signal_kallisto_df <- full_join(kallistotrans_file, signalplist_mature_signal_df, by = "Transdecoder_ID")


#read in interproscan csv
Interproscan <- read.csv(file = Interproscan_file, header = TRUE, sep = ",")
#merge exisiting dataframe with interproscan
FINAL_CSV <- left_join(fin_signalplist_mature_signal_kallisto_df,Interproscan, by="Transdecoder_ID")
FINAL_CSV <- FINAL_CSV %>%
  mutate(SP = if_else(is.na(Signal_Sequence), "OTHER", "SP"))


#keep only columns of interest
keeps <- c("X", "Transdecoder_ID", "ORF_type", "PEP_Length", "CDS_Length", "SP", "Signal_Length", "mature_length","tpm", "percent", "cumulativepercent", "Code", "Hit", "Percentage_Identity", "E_value", "BitScore", "Hit_species", "InterPro_accession_Names","GO_name","Panther_ID_Name", "Phobius_Name", "TMHMM"," Signal_Sequence", "mature_sequence", "PEP_Sequence", "CDS_Sequence")
FINAL_CSV <- FINAL_CSV[keeps]


#write full df to csv
write.csv(FINAL_CSV, paste0(sample_name, "_transdf.csv"), row.names = FALSE)

#creating df of distinct transcripts only keeping the highest scoring hit for each transcript
FINAL_CSV_distinct <- FINAL_CSV[order(FINAL_CSV$BitScore, decreasing = TRUE), ]
FINAL_CSV_distinct <- distinct(FINAL_CSV_distinct, Transdecoder_ID, .keep_all = TRUE)

#saving this data frame, will be used to add to mass spec and genomeblastn data where available
write.csv(FINAL_CSV_distinct, paste0(sample_name, "_transdf_distinct.csv"), row.names = FALSE)

#Create fasta file of all secreted proteins with complete ORFs
FINAL_CSV_distinct_filtered <- FINAL_CSV_distinct %>%
  filter(
    grepl("complete", ORF_type, ignore.case = TRUE), #ORF has to be completed
    grepl("SP", SP_Prediction, ignore.case = TRUE), #only select those with positive signalp result
    TMHMM == FALSE #only select those where tmhmm equals to false as well
  )

#column data we want to include in the fasta sequence D
FINAL_CSV_distinct_filtered_min <- subset(FINAL_CSV_distinct_filtered, select = c(Transdecoder_ID, CDS_Length,PEP_Length, Hit_species,Percentage_Identity,E_value,BitScore,InterPro_accession_Names,GO_name,Phobius_Name,Panther_ID_Name,PEP_Sequence))
#replace na values with blank
FINAL_CSV_distinct_filtered_min[is.na(FINAL_CSV_distinct_filtered_min)] <- ""

# define output fasta file
pep_secreted_fasta_file <- file.path(paste0(sample_name, "_secreted_proteins.pep"))
cds_secreted_fasta_file <- file.path(paste0(sample_name, "_secreted_proteins.cds"))
#open
file_pep <- file(paste0(sample_name, "_secreted_proteins.pep"), open = "w")
file_cds <- file(paste0(sample_name, "_secreted_proteins.cds"), open = "w")

# Loop through each row to write fasta file
for (i in 1:nrow(FINAL_CSV_distinct_filtered_min)) {
  row <- FINAL_CSV_distinct_filtered_min[i, ]
  
  # add units for  lengths and percentage
  cds_len <- paste0(row$CDS_Length, "n")
  pep_len <- paste0(row$PEP_Length, "aa")
  perc_id <- paste0(row$Percentage_Identity, "%")
  
  # build the description parts, including the fixed sample name
  desc_parts <- c(
    sample_name,  # Add the fixed sample name
    cds_len,
    pep_len,
    row$Hit_species,
    row$tpm,
    perc_id,
    row$E_value,
    row$BitScore,
    row$InterPro_accession_Names,
    row$GO_name,
    row$Phobius_Name,
    row$Panther_ID_Name
  )
  
  # sep descriptors with pipe
  description <- paste(desc_parts, collapse = " | ")
  
  # (sequence ID defined)
  header <- paste0(">", row$Transdecoder_ID, " | ", description)
  
  # write
  writeLines(c(header, row$PEP_Sequence), file_pep)
  writeLines(c(header, row$CDS_Sequence), file_cds)
}

# close file
close(file_pep)
close(file_cds)
