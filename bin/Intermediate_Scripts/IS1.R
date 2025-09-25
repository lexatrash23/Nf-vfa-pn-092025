#Intermediate_Script_4
#IS4
#userinputversion 
#20250202 

# Function to check if a package is installed and install if not
BiocManager::install("GenomeInfoDbData")

install_if_missing <- function(package) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package, dependencies = TRUE, repos = "https://cloud.r-project.org/")
    library(package, character.only = TRUE)
  }
}

packages <- c("tidyr", "ggplot2", "gridExtra", "data.table", "dplyr")
lapply(packages, install_if_missing)
BiocManager::install("Biostrings")
library(Biostrings)


#saves the current directory into the current_dir variable
args <- commandArgs(trailingOnly = TRUE)

Trinity_file <- args[1]
blastxunitox6_file <- args[2]
kallisto_file <- args[3]




#loading in the trinity file as a DNAstringset object
Trinity_Fasta_set_sequences <- readDNAStringSet(Trinity_file)

#making a dataframe out of this object by getting teh name, sequence and length
Trinity_Fasta_df <- data.frame (
  Trinity_ID = names(Trinity_Fasta_set_sequences),
  Sequence = as.character(Trinity_Fasta_set_sequences),
  Length = width(Trinity_Fasta_set_sequences),
  stringsAsFactors = FALSE)
#Edits the name of the Trinity IDs
Trinity_Fasta_df$Trinity_ID <- sub(" .*", "", Trinity_Fasta_df$Trinity_ID)


#load in the unitox txt file
Unitox_blastx_6 <- read.table(blastxunitox6_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
#name the columns of the file
colnames(Unitox_blastx_6) <- c("Trinity_ID", "Hit", "pident", "alignment_length", "mismatch", "E_value", "BitScore", "Frame") #name columns
#Default hit column is sp|Code|Hit_Species, this separates that to make three separate columns one with just the Code, one with the Hit and one that has Hit_Species

#this removes the sp| that precedes each Hit
Unitox_blastx_6$Hit <- gsub("^sp\\|", "", Unitox_blastx_6$Hit) #remove sp|
#This seprates the hit and the code column by the | separating them. so now the Code column has the code and the Hit column has the Hit_Species.
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


##Kallisto data
#Nnow we are reading in the kallisto data
Trinity_Kallisto_data <- read.csv(kallisto_file, stringsAsFactors = FALSE)
#rename the target_id column in kallisto so that we can merge by column name in the next step
colnames(Trinity_Kallisto_data)[colnames(Trinity_Kallisto_data) == 'target_id'] <- 'Trinity_ID'
#left_join the Trinity kallisto and trinity one. left_join retains all the rows of the trinity_fasta and just adds in the column data from the kallisto where it is available
Trinity_Kallisto_merge_all <- left_join(Trinity_Fasta_df, Trinity_Kallisto_data, by = "Trinity_ID")
#left join the trinity_kallisto merge. merge with all.x=TRUE is basically a left join, i just wanted to try another way lol.
combined_kallisto_unitox_trinity <- merge(Trinity_Kallisto_merge_all, Unitox_blastx_6, by = "Trinity_ID", all.x = TRUE)


#Distinct list where each trinity ID is only present once, retaining the result with the highest blast bitscore if any
combined_kallisto_unitox_trinity_ordered <- combined_kallisto_unitox_trinity[order(combined_kallisto_unitox_trinity$BitScore, decreasing = TRUE), ]
combined_kallisto_unitox_trinity_distinct <- distinct(combined_kallisto_unitox_trinity_ordered,Trinity_ID, .keep_all = TRUE )

#now we save this dataframe as a csv file in the intermediate_outputs.
write.csv(combined_kallisto_unitox_trinity, "_TBK.csv", row.names = FALSE)

#this is saved as a .csv.zip because that's what the R app is set up to read
write.csv(combined_kallisto_unitox_trinity_distinct, file=gzfile("_TBK_distinct.csv.gz"), row.names = FALSE)
