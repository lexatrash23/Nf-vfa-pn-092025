#!/usr/bin/env Rscript
#Intermediate_Script_1
#IS12
#userinputversion

library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

Venn_diagram_union <- args[1]
transdf_distinct <- args[2]
Blastn_result <- args[3]




Venn_diagram_union_read <- read.csv(Venn_diagram_union, header = TRUE)


transdf_distinct_read <- read.csv(transdf_distinct, header = TRUE)
Blastn_result_read <- read.table(Blastn_result)
           
colnames(Blastn_result_read) <- c("Transdecoder_ID", "sseqid", "pident", "length", "mismatch","gapopen","qstart","qend","sstart","send","evalue", "bitscore")

Blastn_result_read_sortedbyBS <- Blastn_result_read[order(Blastn_result_read$bitscore, decreasing = TRUE), ]
Blastn_result_read_sortedbyBS_distinct <- Blastn_result_read_sortedbyBS %>% distinct(Transdecoder_ID, .keep_all = TRUE) 


#join transdf_distinct with Blastn_distinct 
transdf_distinct_blastn <- left_join(transdf_distinct_read,Blastn_result_read_sortedbyBS_distinct, by = "Transdecoder_ID")

write.csv(transdf_distinct_blastn, "transdf_distinct_blastn.csv", row.names = FALSE)
#join venn_overview with Blastn_distinct
venn_overview_blastn <- left_join(Venn_diagram_union_read,Blastn_result_read_sortedbyBS_distinct, by = "Transdecoder_ID")


write.csv(venn_overview_blastn, "venn_overview_blastn.csv", row.names = FALSE)

venn_overview_blastn_filtered <- venn_overview_blastn %>%
  filter(!is.na(sseqid))

write.csv(venn_overview_blastn_filtered, "venn_overview_blastn_filtered.csv", row.names = FALSE)
