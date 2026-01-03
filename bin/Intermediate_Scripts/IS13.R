#!/usr/bin/env Rscript
#Intermediate_Script_1
#IS13
#userinputversion

#for samples with or without massspec and with  genome available
#same script for process 16 and process 17 except for change of file output name 
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

transdf_distinct_full <- args[1]
transdf_distinct_sim <- args[2]
Blastn_result <- args[3]
sample <- args[4]

transdf_distinct_full <- read.csv(transdf_distinct_full, header = TRUE)
transdf_distinct_sim <- read.csv(transdf_distinct_sim, header = TRUE)
Blastn_result_read <- read.table(Blastn_result)
colnames(Blastn_result_read) <- c("Transdecoder_ID", "sseqid", "genome_pident", "genome_length", "genome_mismatch","genome_gapopen","genome_qstart","genome_qend","genome_sstart","genome_send","genome_evalue", "genome_bitscore", "genome_qframe", "genome_qcovs")
Blastn_result_read_sortedbyBS <- Blastn_result_read[order(Blastn_result_read$bitscore, decreasing = TRUE), ]
#keep only strongest hit by bitscore
Blastn_result_read_sortedbyBS_distinct <- Blastn_result_read_sortedbyBS %>% distinct(Transdecoder_ID, .keep_all = TRUE) 


#full df 
#join transdf_distinct with Blastn_distinct 
transdf_distinct_blastn <- left_join(transdf_distinct_full,Blastn_result_read_sortedbyBS_distinct, by = "Transdecoder_ID")
write.csv(transdf_distinct_blastn, paste0(sample, "_transdf_distinct_masspec_blastn_full.csv"), row.names = FALSE)

#simplified df 
Blastn_sim <- Blastn_result_read_sortedbyBS_distinct %>% 
  select(Transdecoder_ID,sseqid,genome_pident,genome_length,genome_mismatch,genome_sstart,genome_send,genome_evalue,genome_bitscore,genome_qcovs)
transdf_distinct_blastn_sim <- left_join(transdf_distinct_sim,Blastn_sim, by = "Transdecoder_ID")
write.csv(transdf_distinct_blastn, paste0(sample, "_transdf_distinct_masspec_blastn_simplified.csv"), row.names = FALSE)

#onlythosewithmassspecevidenceandblastn 
transdf_distinct_blastn_filtered <- transdf_distinct_blastn %>%
  filter(!is.na(sseqid))  %>% filter(!is.na(Coverage...))
write.csv(transdf_distinct_blastn_filtered, paste0(sample, "_transdf_distinct_masspec_blastn_filtered_simplified.csv"), row.names = FALSE)
