#!/usr/bin/env Rscript

#a rough overview of interproscan domains html graph + filtered sequence list of those that have at least one domain common with a known toxin
#installing and loading packages
library(Biostrings)
library(dplyr)
library(tidyr)
library(htmlwidgets)
library(plotly)


#command line arguments
args <- commandArgs(trailingOnly = TRUE)
Trans <- args[1] #transdf 
Toxin_domains_file <- args[2] #toxin domains file
fasta_file_pep <- args[3] # the previous secreted protein files that were made
fasta_file_cds <- args[4] # the previous secreted protein files that were made
sample <- args[4] #sample name

Trans <- "/Users/praveena/Desktop/PhD_all/2025/github/oldergithubscripts/_transdf.csv" #transdf_distinct
toxin_domains_file <- "/Users/praveena/Downloads/Toxin_domains_April26.tsv"
sample <- "DP3"

#Generating Interproscan plotly
#read in transdf_distinct
Trans <- read.csv(file = Trans, header = TRUE)  %>%
  dplyr::rename(SP =SP_Prediction) %>%
  mutate(CysPer = ((str_count(PEP_Sequence, fixed("C"))) / str_length(PEP_Sequence)) * 100 )
  

# Reordering columns and keeping only those of interest
keeps <- c("Transdecoder_ID", "ORF_type", "PEP_Length", "CDS_Length", "SP","Signal_Length", "mature_length", "CysPer", "percent", "cumulativepercent", "Code", "Hit", "Percentage_Identity", "E_value", "BitScore", "Hit_species","InterPro_accession_Names","GO_name","TMHMM","Signal_Sequence","mature_sequence", "PEP_Sequence","CDS_Sequence")
Trans <- Trans[keeps]


# sort by bitscore
Trans_sorted <- Trans[order(Trans$BitScore, decreasing = TRUE),]
#sort only for those complete with false tmhmm
FINAL_CSV_distinct_filtered <- Trans_sorted %>%
  filter(
    grepl("complete", ORF_type, ignore.case = TRUE),
    grepl("SP", SP, ignore.case = TRUE),
    TMHMM == FALSE
  )
#keep only columns of interest
keeps <- c("Transdecoder_ID", "InterPro_accession_Names","GO_name")
df_figures <-FINAL_CSV_distinct_filtered[keeps]



#separate the interproscan IPs by PIPs
annotation_list <- strsplit(df_figures$InterPro_accession_Names, split = "\\|")
#creates one long list of IPs and summarizes the counts for each
all_annotations <- unlist(annotation_list)
annotation_counts <- table(all_annotations)
annotation_counts_df <- as.data.frame(annotation_counts)
#sort by frequencing of IP
annotation_counts_df <- annotation_counts_df[order(-annotation_counts_df$Freq), ]

colnames(annotation_counts_df) <- c("Transdf_annotations", "Freq")

write.csv(annotation_counts_df, paste0(sample,"_annotation_counts_trandsf.csv"), row.names = FALSE)


#read in toxin domains tsv
Toxin_domains <- read.delim(file = Toxin_domains_file, header = TRUE, sep = "\t")
#same as above but different separator for this tsv
interpro_list <- strsplit(Toxin_domains$InterPro, split = ";")
all_interpro_ids <- unlist(interpro_list)
all_interpro_ids <- trimws(all_interpro_ids) #remove ws
unique_interpro_ids <- unique(all_interpro_ids) #only get unique IDs
unique_interpro_df <- data.frame(InterPro = unique_interpro_ids)
annotation_counts_toxin <- table(all_interpro_ids)
annotation_counts_toxin_df <- as.data.frame(annotation_counts_toxin)
annotation_counts_toxin_df <- annotation_counts_toxin[order(-annotation_counts_toxin$Freq), ]
colnames(annotation_counts_toxin_df) <- c("toxin_annotations", "Freq")

write.csv(annotation_counts_toxin_df, paste0(sample,"_annotation_counts_toxins.csv"), row.names = FALSE)

#iterate through the unique_interpro_df$InterPro, search df_figres_interpro_names df with this value as query, save the rows that have a match at least somewhere in  the row
matched_rows <- data.frame()
for (interpro_id in unique_interpro_df$InterPro) {
  # searches
  matches <- df_figures[apply(df_figures, 1, function(row) any(grepl(interpro_id, row, fixed = TRUE))), ]
  
  # add the whole df_Figures row row to a new dataframe
  if (nrow(matches) > 0) {
    matched_rows <- rbind(matched_rows, matches)
  }
}
#unique Transdecoder IDs only
matched_rows <- distinct(matched_rows, Transdecoder_ID, .keep_all = TRUE)

#remove any na rows
matched_rows <- matched_rows[!is.na(matched_rows$InterPro_accession_Names), ]

#Frequency list of toxin-related domains
annotation_list_matched <- strsplit(matched_rows$InterPro_accession_Names, split = "\\|")
all_annotations_matched  <- unlist(annotation_list_matched)
annotation_counts_matched <- table(all_annotations_matched)
annotation_counts_df_matched <- as.data.frame(annotation_counts_matched)
annotation_counts_df_matched <- annotation_counts_df_matched[order(-annotation_counts_df_matched$Freq), ]

colnames(annotation_counts_df_matched) <- c("matched_toxin_annotations", "Freq")

write.csv(annotation_counts_df_matched, paste0(sample,"_annotation_counts_trandf_matched_to_toxins.csv"), row.names = FALSE)




fasta <- readAAStringSet(file = fasta_file_pep)
seq_ids <- sapply(strsplit(names(fasta), " "), `[`, 1)
keep <- seq_ids %in% matched_rows$Transdecoder_ID
filtered_fasta <- fasta[keep]
writeXStringSet(filtered_fasta, file.path("filtered_sequences.fasta"))


group_totals <- annotation_counts_df_matched %>%
  group_by(Group) %>%
  summarise(TotalFreq = sum(Freq)) %>%
  mutate(GroupLabel = paste0(Group, " (", TotalFreq, ")"))
annotation_counts_df_matched <- annotation_counts_df_matched %>%
  left_join(group_totals, by = "Group")




