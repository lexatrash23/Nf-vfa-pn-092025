#!/usr/bin/env Rscript
#Intermediate_Script_5
#IS_5

#library
library(Biostrings)
library(dplyr)
library(tidyr)
library(stringr)

#allow for command line arguments
args <- commandArgs(trailingOnly = TRUE)


transdecoder_pep_file <- args[1] #complete_pep
transdecoder_cds_file <- args[2] #complete_cds
blastpunitox6_file <- args[3]    #Blastp6 output raw
signalpmature_file <- args[4]    #SignalP output raw mature fasta raw or combined mature if DeepTMHMM is used 
signalpseq_file <- args[5]       #SignalSequence output fasta file from IS2 (Process 2)
Interproscan_file <- args[6]     #Interproscan file output from IS4 (Process 5)
kallistotrans <- args[7]         #Kallisto transdecoder csv output from kallistoanalysistrans.py (Process 2)
sample_name <- args[8]          #sample name



#read in transdecoder pep files
transdecoder_pep <- readAAStringSet(transdecoder_pep_file)
transdecoder_pep_df <- data.frame(
  Transdecoder_ID = names(transdecoder_pep),  # save sequence_id
PEP_Sequence = as.character(transdecoder_pep),  # save sequence to as character
PEP_Length = width(transdecoder_pep),            # save sequence length
  stringsAsFactors = FALSE,
  row.names = NULL  # ensures no row names
)

#add ORF type column to df
transdecoder_pep_df$ORF_type <- sub(".*ORF type:([^,]+).*", "\\1", transdecoder_pep_df$Transdecoder_ID) #searches the transdecoderID captures [^<,]+ any characters are the 'ORF_type" and before the first comma after this comma as the first capturing groups, saving it in the ORF_type columm.
#rename the TransdecoderID  removing everything after the first space
transdecoder_pep_df$Transdecoder_ID <- sub(" .*", "", transdecoder_pep_df$Transdecoder_ID)
#rename sequence and length columns to be pep specific

#read in transdecoder cds file
transdecoder_cds <- readDNAStringSet(transdecoder_cds_file)
transdecoder_cds_df <- data.frame (
  Transdecoder_ID = names(transdecoder_cds),
CDS_Sequence = as.character(transdecoder_cds),
CDS_Length = width(transdecoder_cds),
  stringsAsFactors = FALSE,
  row.names = NULL 
)

#rename TransdecoderID to remove everything after the first space
transdecoder_cds_df$Transdecoder_ID <- sub(" .*", "", transdecoder_cds_df$Transdecoder_ID)

#merge transdecoder files full join by ID column
Transdecoder_pep_cds_merge <- left_join(transdecoder_pep_df, transdecoder_cds_df, by = "Transdecoder_ID")

# readin blasatpunitox6 file 
blastpunitox6 <- read.table(blastpunitox6_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(blastpunitox6) <- c("Transdecoder_ID", "Hit", "Percentage_Identity", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "ssend", "E_value", "BitScore", "Frame", "query_coverage")
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
Transdecoder_blastp <- left_join(Transdecoder_pep_cds_merge, blastpunitox6, by = "Transdecoder_ID")

#moving the sequences columns right to the end
Transdecoder_blastp <- Transdecoder_blastp[, c(setdiff(names(Transdecoder_blastp), c("PEP_Sequence", "CDS_Sequence")), "PEP_Sequence", "CDS_Sequence")]


#read in signal p mature fasta file
signalpmature <- readAAStringSet(signalpmature_file)
#save data as a df
signalpmature_df <- data.frame(
  Transdecoder_ID = names(signalpmature),  # sequence names
Mature_Sequence = as.character(signalpmature),   # sequence
    Mature_Length = width(signalpmature),            # sequence length
  stringsAsFactors = FALSE,
  row.names = NULL
)

#there should not be anything after the space, but just in case
signalpmature_df$Transdecoder_ID <- sub(" .*", "", signalpmature_df$Transdecoder_ID)


#with the addition of DeepTMHMM the signalplist file is not needed, but kept on for now

#merge mature with previous df full join
mature_df <- left_join(Transdecoder_blastp, signalpmature_df, by = "Transdecoder_ID")


# Add signal sequences
signalpsignalseq <- readAAStringSet(signalpseq_file)
signalpsignalseq_dfs <- data.frame(
  Transdecoder_ID = names(signalpsignalseq),  # ID
  Signal_Sequence = as.character(signalpsignalseq),   # Sequence
  Signal_Length = width(signalpsignalseq),            # Length
  stringsAsFactors = FALSE,
  row.names = NULL
)
head(signalpsignalseq_dfs)
#full join merge to add Signal Sequence Data
signalplist_mature_signal_df <- left_join(mature_df, signalpsignalseq_dfs, by = "Transdecoder_ID")


#Add in kallisto data
kallistotrans_file <-read.table(file = kallistotrans, sep = ',', header = TRUE)
#rename column
colnames(kallistotrans_file)[colnames(kallistotrans_file) == "target_id"] <- "Transdecoder_ID"
#full join
fin_signalplist_mature_signal_kallisto_df <- left_join(signalplist_mature_signal_df,kallistotrans_file, by = "Transdecoder_ID")


#read in interproscan csv
Interproscan <- read.csv(file = Interproscan_file, header = TRUE, sep = ",")
#merge exisiting dataframe with interproscan
FINAL_CSV <- left_join(fin_signalplist_mature_signal_kallisto_df,Interproscan, by="Transdecoder_ID")
FINAL_CSV <- FINAL_CSV %>%
  mutate(SP = if_else(is.na(Signal_Sequence), "OTHER", "SP"))

FINAL_CSV <- FINAL_CSV %>%
  mutate(CysPer = ((str_count(PEP_Sequence, fixed("C"))) / (str_length(PEP_Sequence)-1)) * 100 )
#keep only columns of interest
keeps <- c("Transdecoder_ID", "ORF_type","SP", "PEP_Length", "CDS_Length", "Signal_Length", "Mature_Length","CysPer","est_counts","tpm", "percent", "cumulativepercent", "Code", "Hit", "Percentage_Identity", "E_value", "BitScore", "Hit_species", "InterPro_accession_Names","GO_name","Panther_ID_Name", "Phobius_Name", "TMHMM","Signal_Sequence", "Mature_Sequence", "PEP_Sequence", "CDS_Sequence")
FINAL_CSV <- FINAL_CSV[keeps]


#Create fasta file of all secreted proteins with complete ORFs
FINAL_CSV_filtered <- FINAL_CSV %>%
  filter(
    grepl("complete", ORF_type, ignore.case = TRUE), #ORF has to be completed
    grepl("SP", SP, ignore.case = TRUE), #only select those with positive signalp result
    TMHMM != TRUE #only select those where tmhmm doesnt equals to true as well
  )

#for the html reports complete orfs sectoin
CompleteORFsDf <- FINAL_CSV %>%
  filter(
    grepl("complete", ORF_type, ignore.case = TRUE) #ORF has to be completed
  )

#write full df to csv
write.csv(CompleteORFsDf, paste0(sample_name, "_transdf.csv"), row.names = FALSE)

#creating df of distinct transcripts only keeping the highest scoring hit for each transcript # this is only for secreted ORFs used used for interproscan figure generation.
FINAL_CSV_distinct <- FINAL_CSV_filtered[order(FINAL_CSV_filtered$BitScore, decreasing = TRUE), ]
FINAL_CSV_distinct <- distinct(FINAL_CSV_distinct, Transdecoder_ID, .keep_all = TRUE)

#saving this data frame, will be used to add to mass spec and genomeblastn data where available
write.csv(FINAL_CSV_distinct, paste0(sample_name, "_transdf_distinct.csv"), row.names = FALSE)


#column data we want to include in the fasta sequence D
FINAL_CSV_distinct_filtered_min <- subset(FINAL_CSV_distinct, select = c(Transdecoder_ID, CDS_Length,PEP_Length, Hit_species,Percentage_Identity,CysPer,E_value,BitScore,InterPro_accession_Names,GO_name,Phobius_Name,Panther_ID_Name,PEP_Sequence))
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
  CysPers <- paste0(row$CysPer, "%")
  
  
  # build the description parts, including the fixed sample name
  desc_parts <- c(
    paste0("sample=", sample_name),
    paste0("cds_len=", cds_len),
    paste0("pep_len=", pep_len),
    paste0("CysPercentage=", CysPers),
    paste0("hit_species=", row$Hit_species),
    paste0("tpm=", row$tpm),
    paste0("perc_id=", perc_id),
    paste0("evalue=", row$E_value),
    paste0("bitscore=", row$BitScore),
    paste0("interpro=", row$InterPro_accession_Names),
    paste0("go=", row$GO_name),
    paste0("phobius=", row$Phobius_Name),
    paste0("panther=", row$Panther_ID_Name)
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
