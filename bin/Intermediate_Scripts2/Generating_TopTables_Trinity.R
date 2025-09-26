#!/usr/bin/env Rscript
#installing and laoding packages

install_if_missing <- function(package) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package, dependencies = TRUE, repos = "https://cloud.r-project.org/",type = "binary")
    library(package, character.only = TRUE)
  }
}
packages <- c("dplyr", "DT", "knitr", "kableExtra","shiny")
lapply(packages, install_if_missing)
install.packages("seqinr", dependencies = TRUE, type = "source", repos = "https://cloud.r-project.org/")


#Takes three arguments, the first being the TBK csv , the second being the NCBI genome code and the third being the full species name
args <- commandArgs(trailingOnly = TRUE)
TBK <- args[1]
genome_code  <- args[2]
species_name <- args[3]

#check of current working directory. The script assumes the working directory is the Venomflowanalysis folder
#reading in our TBK file

TBK <- read.csv(TBK, header = TRUE)
# Reordering columns and keeping only those of interest
keeps <- c("Trinity_ID", "Length", "percent", "cumulativepercent", "Code", "Hit", "E_value", "BitScore", "Frame", "pident", "alignment_length", "mismatch", "Hit_species", "Sequence")
Kallisto_Blastx_Trinity <- TBK[keeps]
#Sort this table by decreasing BitScore
KBT_sortedbybitscore <- Kallisto_Blastx_Trinity[order(Kallisto_Blastx_Trinity$BitScore, decreasing = TRUE), ]
#include new column that will have external links to to Blastn against genome and a Tblastn against NCBI and create a hyperlink for uniprot ID 

KBT_sortedbybitscore_with_Hyperlinks <- KBT_sortedbybitscore %>%
  mutate(
    Code = paste0(
      "<a href=\"https://www.uniprot.org/uniprotkb/",
      Code,
      "/entry\" target=\"_blank\">",
      Code,
      "</a>"
    ),
    
    Sequence = paste0(
      "<a href=\"https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastx&PAGE_TYPE=BlastSearch&QUERY=",
      Sequence,
      "\" target=\"_blank\">BLASTx (NCBI Nucleotide DB) Search</a>",
      if (!is.null(genome_code)) {
        paste0(
          " | <a href=\"https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastx&BLAST_SPEC=GDH_",
          genome_code,
          "&PAGE_TYPE=BlastSearch&QUERY=",
          Sequence,
          "\" target=\"_blank\">BLASTn (", species_name, " Genome) Search</a>"
        )
      } else {
        ""
      }
    ),
    
    Full_Sequence = KBT_sortedbybitscore$Sequence
  )

#Table 1 is top 100 by bitscore regardless of hit/transcript repitition 
Table1 <- head(KBT_sortedbybitscore_with_Hyperlinks, n = 100)
#Table 2 is top 100 by bitscore each distinct transcript 
KBT_sortedbybitscore_with_Hyperlinks_Distinct_Transcripts <- KBT_sortedbybitscore_with_Hyperlinks %>%
  distinct(Trinity_ID, .keep_all = TRUE)
Table2 <- head(KBT_sortedbybitscore_with_Hyperlinks_Distinct_Transcripts, n = 100)
#Table 3 is top 100 by bitscore each distinct hit  
KBT_sortedbybitscore_with_Hyperlinks_Distinct_Hits <- KBT_sortedbybitscore_with_Hyperlinks %>%
  distinct(Hit, .keep_all = TRUE)
Table3 <- head(KBT_sortedbybitscore_with_Hyperlinks_Distinct_Hits, n = 100)
#Table 4 is top 25% of kallisto expression  

KBT_sortedbybitscore_with_Hyperlinks_Distinct_Transcripts_increasingcumulativepercent <- KBT_sortedbybitscore_with_Hyperlinks_Distinct_Transcripts[order(KBT_sortedbybitscore_with_Hyperlinks_Distinct_Transcripts$cumulativepercent, decreasing = FALSE), ]



Table4 <- KBT_sortedbybitscore_with_Hyperlinks_Distinct_Transcripts_increasingcumulativepercent[1:10, ]

#save all four tables to intermediate_outputs2
write.csv(Table1,file = file.path("Table1.csv"), row.names = FALSE)
write.csv(Table2,file = file.path("Table2.csv"), row.names = FALSE)
write.csv(Table3,file = file.path("Table3.csv"), row.names = FALSE)
write.csv(Table4,file = file.path("Table4.csv"), row.names = FALSE)

