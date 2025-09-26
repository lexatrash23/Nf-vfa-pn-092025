#!/usr/bin/env Rscript
#load pacakages

install_if_missing <- function(package) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package, dependencies = TRUE, repos = "https://cloud.r-project.org/")
    library(package, character.only = TRUE)
  }
}

packages <- c("httr", "jsonlite", "dplyr", "tidyr", "data.table", "dplyr")
lapply(packages, install_if_missing)
BiocManager::install("GO.db")
BiocManager::install("biomaRt")

library(GO.db)
library(biomaRt)
#saving current directory as a variable

args <- commandArgs(trailingOnly = TRUE)


Interproscan <- args[1]
List <- args[2]
PANTHER_file <- args[3]
#Read in our interproscan output tsv file

#load file 
Interproscan<- read.delim(file = Interproscan, header = FALSE, sep = "\t")


#Naming columns based on the Interproscan manual
colnames(Interproscan) <- c ("Transdecoder_ID", "Sequence_MD5_digest","Sequence_length ", "Analysis","Signature_accession", "Signature_description","Start_location", "Stop_location","E-value","Status","Date","InterPro_annotations_accession", "InterPro_annotations_description", "GO_annotations", "Pathways_annotations")
#select only columns of most interest
data_subset <- subset(Interproscan, select = c(Transdecoder_ID, Analysis,Signature_accession, InterPro_annotations_accession,GO_annotations,Pathways_annotations))
# removing all the - characters marking the empty columns
data_subset$InterPro_annotations_accession<-gsub("-","",as.character(data_subset$InterPro_annotations_accession))
data_subset$GO_annotations<-gsub("-","",as.character(data_subset$GO_annotations))
data_subset$Pathways_annotations <- gsub("^-$", "", data_subset$Pathways_annotations)
data_subset <- data.frame(data_subset) #idk if this actually necessary but i always do this just be sure its a dataframe lol

#readingi in the IP tsv metadata file from Interproscan
IP_name <- read.delim(file = List, header = TRUE, sep = "\t")
# The entry name column is so long, so dividing it only take name as what is before the first comma
Lookup_table <- IP_name %>%
  separate(ENTRY_NAME, into = c("Name", "Type"), sep = ",") #separating that entry_name column cos i only want that name
#now we are only keeping the column with the IP id and the actual name
keep <- c("ENTRY_AC", "Name") #keeping just the ones we want to merge
Lookup_table<-Lookup_table[keep] #keeping just the ones we want to merge

#Renaming the columns for the full join
colnames(Lookup_table) = c("InterPro_annotations_accession","InterPro_Names")
#full joining with dataset based on the InterPro_annotations_accession
data_subset_with_IP_names <- full_join (data_subset, Lookup_table, by = "InterPro_annotations_accession")
#now just keeping columns of interest cause we want to later collapse the same trinity IDs into one row
keep <- c("Transdecoder_ID","InterPro_annotations_accession","InterPro_Names")
IP_set <- data_subset_with_IP_names[keep]
#creating a new column that combines the ID with the name  and just keeping that column
IP_set$InterPro_accession_Names <- paste0(IP_set$InterPro_annotations_accession, "(",IP_set$InterPro_Names ,")")
keep <- c("Transdecoder_ID","InterPro_accession_Names")
IP_set <- IP_set[keep]
#removing the NA values
IP_set$InterPro_accession_Names <- gsub("\\(NA\\)", "", IP_set$InterPro_accession_Names)
#collapsing back so we only have one row per Transdecoder_ID
IP_set_collapsed <- IP_set %>%
  group_by(Transdecoder_ID) %>%
  summarise(
    InterPro_accession_Names = paste(
      unique(InterPro_accession_Names[InterPro_accession_Names != ""]),
      collapse = "|"
    ),
    .groups = "drop"
  )




#Now basically the ame thing but now for the GO annotations
#obtaining our list of unique GO terms prsent in teh data file
go_list <- strsplit(data_subset$GO_annotations, split = "\\|") #split it
go_terms <- unlist(lapply(go_list, function(x) sub("\\(.*\\)", "", x))) #remove extra characters
unique_go_terms <- unique(go_terms) #find the unique go terms
unique_go_df <- data.frame(GO = unique_go_terms) #make into a dataframe


#from 81-112: extracting GO_ID names from numbers using
#defining the columns of the dataframe
go_ids <- unique_go_df$GO #
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

#for now we are only interested in the ID and the names
keep = c("GO_ID","Name")
go_metadata_subset <- go_metadata[keep]
head(go_metadata_subset)

#make a long version of teh the data subset so that theres only one GO per row  and remove teh characters in teh brackets
df_long <- data_subset %>%
  separate_rows(GO_annotations, sep = "\\|")
df_long$GO_annotations <- gsub("\\(.*\\)", "", df_long$GO_annotations)


#left join to keep all the rows of the subset but just add the ID names based on teh go IDs
df_joined <- df_long %>%
  left_join(go_metadata_subset, by = c("GO_annotations" = "GO_ID"))
#same as IP taking out columns of interest so we can join and then collapse later
keep <- c("Transdecoder_ID","GO_annotations","Name")
df_joined <-df_joined[keep]
#creating the column combining the go and the name  and just keeping that association
df_joined$GO_name <- paste0(df_joined$GO_annotations, "(",df_joined$Name, ")")
keep <- c("Transdecoder_ID", "GO_name")
df_joined <-df_joined[keep]
#get rid of NAs
df_joined$GO_name <- gsub("\\(NA\\)", "", df_joined$GO_name)

#collapse this so it's one Transdecoder_ID per row
df_collapsed_GO <- df_joined %>%
  group_by(Transdecoder_ID) %>%
  summarise(
    GO_name = paste(
      unique(GO_name[GO_name != ""]),
      collapse = "|"
    ),
    .groups = "drop"
  )

# Now same same but yo annotate the sig_accessions

#Panther
PANTHER <- read.csv(file = PANTHER_file, header = TRUE, row.names = NULL)
PANTHER <-PANTHER[2:3]

#rows with PANTHER analysis
#filter the rows in df1 where Analysis is "PANTHER" and merge with df2
Sig_accession_set <- subset(Interproscan, select = c(Transdecoder_ID, Analysis,Signature_accession))
Sig_accession_set_panther <- Sig_accession_set[Sig_accession_set$Analysis == "PANTHER",]
#we are gonna do a left join of the oanther_Unique and this set
Sig_accession_set_panther_merged <- left_join(Sig_accession_set_panther, PANTHER, by = c("Signature_accession" = "Panther_ID"))
head(Sig_accession_set_panther_merged)
Sig_accession_set_panther_merged$Panther_ID_Name <- paste0 (Sig_accession_set_panther_merged$Panther_Name,"(", Sig_accession_set_panther_merged$Analysis, ")")
                                          # View final dataframe

Sig_accession_set_panther_merged_subset <- subset(
  Sig_accession_set_panther_merged,
  select = c(Transdecoder_ID, Panther_ID_Name)
)
collapsed_panther <- Sig_accession_set_panther_merged_subset %>%
  group_by(Transdecoder_ID) %>%
  summarise(Panther_ID_Name = paste(unique(Panther_ID_Name), collapse = " | "), .groups = "drop")

#NexT Phobius
Sig_accession_set_phobius <- Sig_accession_set[Sig_accession_set$Analysis == "Phobius",]
Sig_accession_set_phobius$Phobius_Name <- paste0 (Sig_accession_set_phobius$Signature_accession,"(", Sig_accession_set_phobius$Analysis, ")")
Sig_accession_set_phobius_subset <- subset(
  Sig_accession_set_phobius,
  select = c(Transdecoder_ID, Phobius_Name)
)
collapsed_phobius <- Sig_accession_set_phobius_subset %>%
  group_by(Transdecoder_ID) %>%
  summarise(Phobius_Name = paste(unique(Phobius_Name), collapse = " | "), .groups = "drop")
head(collapsed_phobius)
#merged collapsed_phobius and collapsed_panther
phobius_panther <- full_join(collapsed_panther, collapsed_phobius, by = "Transdecoder_ID")
head(phobius_panther)



#subset with only Transdecoder ID, oathways and tmhmm
data_subset_2 <- subset(data_subset, select = c(Transdecoder_ID,Pathways_annotations))
data_subset_2 <- distinct(data_subset_2, Transdecoder_ID, .keep_all = TRUE)
#merge them all together
Final_interproscan_dataframe <- left_join(data_subset_2,IP_set_collapsed, by ="Transdecoder_ID")
Final_interproscan_dataframe <- left_join(Final_interproscan_dataframe,df_collapsed_GO, by ="Transdecoder_ID")
Final_interproscan_dataframe<-left_join(Final_interproscan_dataframe,phobius_panther, by ="Transdecoder_ID")
colnames(Final_interproscan_dataframe)
Final_interproscan_dataframe$TMHMM <- grepl("TRANSMEMBRANE\\(Phobius\\)", Final_interproscan_dataframe$Phobius_Name)
head(Final_interproscan_dataframe)
head(Final_interproscan_dataframe)
keep <- c("Transdecoder_ID","InterPro_accession_Names","GO_name","Panther_ID_Name","Phobius_Name","TMHMM")
Final_interproscan_dataframe <-Final_interproscan_dataframe[keep]

write.csv(Final_interproscan_dataframe, file = "Final_interproscan_dataframe.csv", row.names = FALSE)
