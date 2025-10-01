#!/usr/bin/env Rscript
library(dpylr)
args <- commandArgs(trailingOnly = TRUE)
dataframe <- args[1]
toxin_data_csv <- args[2]


toxin_data <- read.csv(toxin_data_csv, sep = "\t", header = TRUE)

#read in the dataframe 
df <- read.csv(dataframe)
# order by bitscore, if bitscore is equal sort by evalue 
df <- df[order(-df$BitScore, df$E_value ), ]
df <- distinct(df, Transdecoder_ID, .keep_all = TRUE)
#filter for signalP present, no THMM and complete 

# mass spec group if present 
if (!"Coverage..." %in% colnames(df)) {
  set_B <- data.frame(Transdecoder_ID = character(0))
} else {
  set_B <- df[df[["Coverage..."]] > 50 & !is.na(df[["Coverage..."]]), ]
}

set_C <- df[df[["BitScore"]] > 250 & !is.na(df[["BitScore"]]), ]

set_D <-df[df[["percent"]] > 1 & !is.na(df[["percent"]]), ]

toxin_Domain_long <- toxin_data %>%
  separate_rows(InterPro, sep = ";") %>%    
  mutate(InterPro = str_trim(InterPro)) %>%  
  filter(InterPro != "") 

Toxin_IPs <- data.frame(unique(toxin_Domain_long$InterPro)) 

colnames(Toxin_IPs)[colnames(Toxin_IPs) == "unique.toxin_Domain_long.InterPro."] <- "IP"
pattern <- paste(Toxin_IPs$IP, collapse = "|") 
matching_rows <- df[grepl(pattern, df$InterPro_accession_Names), ]
set_A<-matching_rows

if (!"Coverage..." %in% colnames(df)) {
  venn_list <- list(
    TD = set_A$Transdecoder_ID,
    MS = set_B$Transdecoder_ID,
    TP = set_C$Transdecoder_ID,
    KE = set_D$Transdecoder_ID
  )
  p <- ggvenn(venn_list,c("TD", "MS","TP","KE"), fill_color = c("#E41A1C", "#377EB8", "#4DAF4A", "#EBAC4D"))
  Union_ABC <- Reduce(union, list(set_A$Transdecoder_ID, set_B$Transdecoder_ID, set_C$Transdecoder_ID, set_D$Transdecoder_ID))
  
} else {
  venn_list <- list(
    TD = set_A$Transdecoder_ID,
    TP = set_C$Transdecoder_ID,
    KE = set_D$Transdecoder_ID
    
  )
  p <- ggvenn(venn_list,c("TD","TP","KE"), fill_color = c("#E41A1C", "#4DAF4A", "#EBAC4D"))
  Union_ABC <- Reduce(union, list(set_A$Transdecoder_ID, set_C$Transdecoder_ID, set_D$Transdecoder_ID))
  
}

ggsave("Venn.png", plot = p, width = 6, height = 4, dpi = 300)


df_union <- df[df$Transdecoder_ID %in% Union_ABC, ]
write.csv(df_union, "Venn_Diagram_union.csv", row.names = FALSE)

