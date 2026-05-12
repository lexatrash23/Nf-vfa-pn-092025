#!/usr/bin/env Rscript

library(dplyr)
library(GO.db)
library(tidyr)
library(stringr)
library(AnnotationDbi)
library(archive)
library(readr)

args <- commandArgs(trailingOnly = TRUE)
ToxinDataTSV <- args[1] 
NonToxinDataTSV  <- args[2] 
NonToxinDataTSV <- "/Users/praveena/Desktop/PhD_all/2025/github/Nf-vfa-pn-092025/MetadataFiles/NonToxin_Domains.7z"
IPmetadata <- args[3] 
#Read in IPmetadatafile 
Interproscan_metadata_csv <- read.delim(file = IPmetadata, header = TRUE, sep = "\t") %>%
 dplyr::select(ENTRY_AC, ENTRY_TYPE,ENTRY_NAME ) %>%
  dplyr::rename(InterPro = ENTRY_AC )
#Read in ToxinDataTSV, only keep those with InterPro Values
toxin_data <- read.delim(ToxinDataTSV, header = TRUE) %>%
  separate_rows(InterPro, sep = ";") %>%    
  mutate(InterPro = str_trim(InterPro)) %>%  
  filter(!is.na(InterPro) & str_trim(InterPro) != "")

# Calculate percentage of toxin proteins that have each domain 
toxin_data <- toxin_data %>%
  count(InterPro) %>%
  mutate(ToxPercent = (n/(as.numeric(nrow(toxin_data))))*100) %>%
  dplyr::select(-n)

#Read in NonToxinDataTSV, only keep those with InterPro Values
nontoxin_data <- read_tsv(archive_read(NonToxinDataTSV), col_names = TRUE,show_col_types = FALSE) %>%
  separate_rows(InterPro, sep = ";") %>%    
  mutate(InterPro = str_trim(InterPro)) %>%  
  filter(!is.na(InterPro) & str_trim(InterPro) != "")


names(nontoxin_data) <- gsub(" ", ".", names(nontoxin_data))
names(toxin_data) <- gsub(" ", ".", names(toxin_data))

# Calculate percentage of nontoxin proteins that have each domain , full join with toxindata values, replace NA values with 0
#Calculate Relative expression, arrange by high to low, assign a rank number, add metadata info anddplyr::select columns of interest
summary_IP <- nontoxin_data %>%
  count(InterPro) %>%
  mutate(NonToxPercent = (n/(as.numeric(nrow(nontoxin_data))))*100) %>%
  dplyr::select(-n) %>%
  full_join(toxin_data, by = "InterPro") %>%
  mutate(NonToxPercent = ifelse(is.na(NonToxPercent), 0, NonToxPercent)) %>% 
  mutate(ToxPercent = ifelse(is.na(ToxPercent), 0, ToxPercent)) %>% # Replace NA with 0 in NonToxPercent
  mutate(RelativeExpression = ToxPercent / NonToxPercent) %>%
  arrange(desc(RelativeExpression)) %>%
  mutate(Rank = row_number()) %>%
  dplyr::select(Rank,InterPro,ToxPercent,NonToxPercent,RelativeExpression) %>%
  left_join(Interproscan_metadata_csv, by = "InterPro") %>%
  mutate(Name_short = str_extract(ENTRY_NAME, "^[^,]+")) %>%
  filter(!is.na(ENTRY_NAME)) %>%
  dplyr::select(Rank,InterPro,,RelativeExpression,ENTRY_TYPE, Name_short) 

#save this file for future processess
write.csv(summary_IP, "ToxNonTox_IP.csv", row.names = FALSE)

#Read in ToxinDataTSV, only keep those with GO Values
toxin_data <- read.delim(ToxinDataTSV, header = TRUE) %>%
  rename_with(~ make.names(.), everything()) %>%
  separate_rows(Gene.Ontology.IDs, sep = ";") %>%
  mutate(Gene.Ontology.IDs = str_trim(Gene.Ontology.IDs)) %>%  
  filter(!is.na(Gene.Ontology.IDs) & str_trim(Gene.Ontology.IDs) != "")

# Calculate percentage of toxin proteins that have each GO 
toxin_data <- toxin_data %>%
  count(Gene.Ontology.IDs) %>%
  mutate(ToxPercent = (n/(as.numeric(nrow(toxin_data))))*100) %>%
  dplyr::select(-n)

#Read in NonToxinDataTSV, only keep those with GO Values
nontoxin_data <- read.delim(NonToxinDataTSV, header = TRUE)  %>%
  rename_with(~ make.names(.), everything()) %>%
  separate_rows(Gene.Ontology.IDs, sep = ";") %>%
  mutate(Gene.Ontology.IDs = str_trim(Gene.Ontology.IDs)) %>%  
  filter(!is.na(Gene.Ontology.IDs) & str_trim(Gene.Ontology.IDs) != "")

# Calculate percentage of nontoxin proteins that have each GO 
#join with toxin data, make NA values 0, calculate relative expression and sort 
summary_IP <- nontoxin_data %>%
  count(Gene.Ontology.IDs) %>%
  mutate(NonToxPercent = (n/(as.numeric(nrow(nontoxin_data))))*100) %>%
  dplyr::select(-n) %>%
  full_join(toxin_data, by = "Gene.Ontology.IDs") %>%
  mutate(NonToxPercent = ifelse(is.na(NonToxPercent), 0, NonToxPercent)) %>% 
  mutate(ToxPercent = ifelse(is.na(ToxPercent), 0, ToxPercent)) %>% # Replace NA with 0 in NonToxPercent
  mutate(RelativeExpression = ToxPercent / NonToxPercent) %>%
  arrange(desc(RelativeExpression)) %>%
  dplyr::select(Gene.Ontology.IDs,ToxPercent,NonToxPercent,RelativeExpression) 


#pull GO metadata
#defining the columns of the dataframe
go_ids <- summary_IP$Gene.Ontology.IDs
term_names <- character(length(go_ids))
definitions <- character(length(go_ids))
ontologies <- character(length(go_ids))

# Loop through our unique df and to metadata
for (i in seq_along(go_ids)) {
  go_id <- go_ids[i]
  go_obj <- GOTERM[[go_id]]
  
  if (!is.null(go_obj)) {
    term_names[i] <- Term(go_obj)
    definitions[i] <- Definition(go_obj)
    ontologies[i] <- Ontology(go_obj)
  } else {
    term_names[i] <- NA
    definitions[i] <- NA
    ontologies[i] <- NA
  }
}

# Combine into a new data frame
go_metadata <- data.frame(
  GO_ID = go_ids,
  Name = term_names,
  Definition = definitions,
  Ontology = ontologies,
  stringsAsFactors = FALSE
)
#add Go_metadata remove rows that refer to cellular compartment
summary_IP <- summary_IP %>%
  dplyr::rename(GO_ID = Gene.Ontology.IDs) %>%
  left_join(go_metadata, by = "GO_ID") %>%
  filter(Ontology != "CC") 

#create MF summary and rank 
Summary_MF <- summary_IP %>%
  filter(Ontology == "MF") %>%
  mutate(Rank = row_number()) %>%
  dplyr::select(Rank,Ontology,GO_ID,Name,RelativeExpression) 

write.csv(Summary_MF, "ToxNonTox_MF.csv", row.names = FALSE)

#create BP summary and rank 
Summary_BP <- summary_IP %>%
  filter(Ontology == "BP") %>%
  mutate(Rank = row_number()) %>%
  dplyr::select(Rank,Ontology,GO_ID,Name,RelativeExpression) 

write.csv(Summary_BP, "ToxNonTox_BP.csv", row.names = FALSE)


