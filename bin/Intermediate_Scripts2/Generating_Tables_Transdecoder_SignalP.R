#!/usr/bin/env Rscript
#installing and laoding packages
library(dplyr)
library(DT)


#this take two arguments, the first being the genome ID and the second being species name
args <- commandArgs(trailingOnly = TRUE)
Trans <- args[1]
genome_code   <- args[2]
species_name  <- args[3]


#check of current working directory. The script assumes the working directory is the

#read in our trandf file

Trans <- read.csv(Trans, header = TRUE)
# Reordering columns and keeping only those of interest
keeps <- c("Transdecoder_ID", "ORF_type", "PEP_Length", "CDS_Length", "SP_Prediction", "Signal_Length", "mature_length", "percent", "cumulativepercent", "Code", "Hit", "Percentage_Identity", "E_value", "BitScore", "Hit_species","InterPro_accession_Names","GO_name","TMHMM","Phobius_Name","Panther_ID_Name", "PEP_Sequence","Signal_Sequence","mature_sequence","CDS_Sequence")
Trans <- Trans[keeps]
#tables to save
#order by bitscore
Trans_sorted <- Trans[order(Trans$BitScore, decreasing = TRUE),]
#add hyperliniks

Trans_hyperlinks <- Trans_sorted %>%
  mutate(
    Code = paste0("<a href=\"https://www.uniprot.org/uniprotkb/", Code, "/entry\" target=\"_blank\">", Code, "</a>"), #mutates the Code column to make a hyperlink leading to the uniprot page for that code
    
    PEP_Sequence = paste0(
      "<a href=\"https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastp&PAGE_TYPE=BlastSearch&QUERY=",
      PEP_Sequence,
      "\" target=\"_blank\">BLASTp (NCBI Protein DB) Search</a> | ", #text that is hyperlinked
      "<a href=\"https://dtu.biolib.com/DeepTMHMM/?seq=",
      PEP_Sequence,
      "\" target=\"_blank\">DeepTMHMM Search</a> | ",
      "<a href=\"https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb#scrollTo=kOblAo-xetgx\" target=\"_blank\">AlphaFold2 Colab</a> | ",
      "<a href=\"https://www.ncbi.nlm.nih.gov/Structure/cdd/wrpsb.cgi?seqinput=",
      PEP_Sequence,
      "\" target=\"_blank\">NCBI CDD Search</a>",
        if (!is.null(genome_code)) {
          paste0(
            " | <a href=\"https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=tblastn&BLAST_SPEC=GDH_G",
            genome_code,
            "&PAGE_TYPE=BlastSearch&QUERY=",
            PEP_Sequence,
            "\" target=\"_blank\">TBLASTn (", species_name, " Genome) Search</a>"
          )
        } else {
          ""
        }
        
    ),
    
    PEP_Sequence_full = Trans_sorted$PEP_Sequence
  )

library(dplyr)
library(stringr)
#add the interproacession hyperlinks as well
Trans_hyperlinks <- Trans_hyperlinks %>%
  mutate(
    InterPro_accession_Names = str_split(InterPro_accession_Names, "\\|") %>%
      lapply(function(x) {
        sapply(x, function(entry) {
          ipr <- str_extract(entry, "IPR\\d+")
          name <- str_extract(entry, "\\(.*?\\)")  # gets the (name) part
          paste0('<a href="https://www.ebi.ac.uk/interpro/entry/InterPro/',
                 ipr, '" target="_blank">', ipr, '</a>', name)
        }) %>% paste(collapse = "|")
      }) %>% unlist()
  )


#add the GO name hyperlinks too
Trans_hyperlinks <- Trans_hyperlinks %>%
  mutate(
    GO_name = str_split(GO_name, "\\|") %>%
      lapply(function(x) {
        sapply(x, function(entry) {
          ipr <- str_extract(entry, "GO:\\d+")
          name <- str_extract(entry, "\\(.*?\\)")  # gets the (name) part
          paste0('<a href="https://www.ebi.ac.uk/QuickGO/term/',
                 ipr, '" target="_blank">', ipr, '</a>', name)
        }) %>% paste(collapse = "|")
      }) %>% unlist()
  )


#Table5 Top 100 by bitscore -ALL
Table5 <- head(Trans_hyperlinks, n = 100)
#Table6 Top 100 by bitscore -distinct_transcripts
Trans_hyperlinks_Distinct_Transcripts <- Trans_hyperlinks %>%
  distinct(Transdecoder_ID, .keep_all = TRUE)
Table6 <- head(Trans_hyperlinks_Distinct_Transcripts, n = 100)
#Table7 Top 100 by bitscore -distinct_hits
Trans_hyperlinks_Distinct_Hits <- Trans_hyperlinks %>%
  distinct(Hit, .keep_all = TRUE)
Table7 <- head(Trans_hyperlinks_Distinct_Hits, n = 100)
#Table8 Top 25% Kallisto Expression

Trans_hyperlinks_Distinct_Transcripts_increasingcumulativepercent <- Trans_hyperlinks_Distinct_Transcripts[order(Trans_hyperlinks_Distinct_Transcripts$cumulativepercent, decreasing = FALSE), ]


Table8 <- Trans_hyperlinks_Distinct_Transcripts_increasingcumulativepercent[1:10, ]

#NEXT set is signalp ones
#filter only for those with signalp
Trans_hyperlinks_signalp <- Trans_hyperlinks[
  Trans_hyperlinks$SP_Prediction == "SP(Sec/SPI)" &
    !is.na(Trans_hyperlinks$SP_Prediction) &
    grepl("complete", Trans_hyperlinks$ORF_type, ignore.case = TRUE),
]
#Table9
Table9 <-head(Trans_hyperlinks_signalp, n = 100)
#Table10
Trans_hyperlinks_signalp_Distinct_Transcripts <- Trans_hyperlinks_signalp %>%
  distinct(Transdecoder_ID, .keep_all = TRUE)
Table10 <- head(Trans_hyperlinks_signalp_Distinct_Transcripts, n = 100)
#Table11
Trans_hyperlinks_signalp_Distinct_Hits <- Trans_hyperlinks_signalp %>%
  distinct(Hit, .keep_all = TRUE)
Table11 <- head(Trans_hyperlinks_signalp_Distinct_Hits, n = 100)
#Table12 Top 10 Kallisto_Expression
Trans_hyperlinks_signalp_Distinct_Transcripts_sortK <- Trans_hyperlinks_signalp_Distinct_Transcripts[order(Trans_hyperlinks_signalp_Distinct_Transcripts$percent, decreasing = TRUE), ]
Table12 <- head(Trans_hyperlinks_signalp_Distinct_Transcripts_sortK, n = 10)



#Save tables to Intermediate_outputs2
write.csv(Table5,file = file.path("Table5.csv"), row.names = FALSE)
write.csv(Table6,file = file.path("Table6.csv"), row.names = FALSE)
write.csv(Table7,file = file.path("Table7.csv"), row.names = FALSE)
write.csv(Table8,file = file.path("Table8.csv"), row.names = FALSE)
write.csv(Table9,file = file.path("Table9.csv"), row.names = FALSE)
write.csv(Table10,file = file.path("Table10.csv"), row.names = FALSE)
write.csv(Table11,file = file.path("Table11.csv"), row.names = FALSE)
write.csv(Table12,file = file.path("Table12.csv"), row.names = FALSE)

