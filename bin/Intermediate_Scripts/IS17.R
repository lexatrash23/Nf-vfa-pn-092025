#!/usr/bin/env Rscript

# Combine our sample ToxinData with toxprot ToxinData for the annotation csv provided to protspace
library(dplyr)
library(tidyr)
library(readr)
library(stringr)


args <- commandArgs(trailingOnly = TRUE)
ProtSpaceAnnotatedCSV <- args[1]
ToxinToxinData <- args[2]
ToxinMetaddata <- "/Users/praveena/Desktop/PhD_all/2025/github/Nf-vfa-pn-092025/MetadataFiles/Toxin_Domains.tsv.gz"
ProtSpaceAnnotatedCSV <- "/Users/praveena/Desktop/PhD_all/2026/Todd/Pipelines/Analysis/results/Overview/Dataframes/Annotated/CV_TL_all_Annotated_df.csv"


#Enyme class labels 
ec_class_map <- c(
  "1" = "Oxidoreductase",
  "2" = "Transferase",
  "3" = "Hydrolase",
  "4" = "Lyase",
  "5" = "Isomerase",
  "6" = "Ligase",
  "7" = "Translocase"
)

SampleData <- read.csv(ProtSpaceAnnotatedCSV) %>%
  dplyr::select(Transdecoder_ID,UniqueSequenceName,PEP_Length,Signal_Length,Mature_Length,CysPer,tpm,Hit,Sample_name,Species,Protein_Name,Gene_Name,Enzyme_Class,Annotation_Source,Domain_Label, Molecular_Function, Biological_Process,AnnotationScore,ToxinDomainScore,ProteomicCoverageScore,KallistoExpressionScore,CysPerScore,BlastScore,GenomeSupportScore,Category,Filter, Minimap_mapq) %>%
  dplyr::rename(Sequence_ID =Transdecoder_ID,) %>%
  mutate(Source = "SampleData") %>%
  dplyr::select(Sequence_ID,Protein_Name,Gene_Name,Species,Source,Enzyme_Class,Hit,Category, Annotation_Source)


ToxinData <-read.csv(ToxinMetaddata, sep = "\t") 
names(ToxinData) <- gsub(" ", ".", names(ToxinData))
ToxinData <- ToxinData %>%
  dplyr::rename(AccessionNo = Entry.Name) %>%
  mutate(Protein.names = sub("\\s*\\(.*", "", Protein.names)) %>%
  dplyr::rename(Name = Protein.names) %>%
  tidyr::separate(Organism, into = c("Species", "CommonName"), sep = " \\(", remove = FALSE, extra = "merge", fill = "right") %>%
  tidyr::separate(AccessionNo, into = c("Hit", "extra"), sep = "\\_", remove = FALSE, extra = "merge", fill = "right") %>%
  mutate(Category = "ToxProtDatabaseProtein", Source = "ToxProtData", Annotation_Source = "UniProt"  ) %>%
  mutate(
    EC_class_number = str_extract(EC.number, "^[0-9]"),
    EC_class_name = ec_class_map[EC_class_number]
  ) %>%
  mutate(Gene_Name = Gene.Names) %>%
  mutate(Gene_Name = word(Gene_Name,1)) %>%
  select(AccessionNo, Name, Gene_Name, Species, Source, EC_class_name,Hit,Category,Source,Annotation_Source)%>%
  dplyr::rename(Sequence_ID = AccessionNo, Protein_Name = Name, Enzyme_Class = EC_class_name,) 

CombinedCsv <- rbind(SampleData, ToxinData)


write.csv(CombinedCsv, "ProtSpaceToxin.csv", row.names = FALSE)



