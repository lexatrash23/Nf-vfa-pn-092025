#!/usr/bin/env Rscript
#!/usr/bin/env Rscript

# Annotation Command line
library(rentrez)
library(xml2)
library(Biostrings)
library(stringr)
library(dplyr)
library(tidyr)
library(archive)
library(readr)
library(purrr)



args <- commandArgs(trailingOnly = TRUE)
transdf_final_filtered_lax_file <- args[1] 
toxprotblast6 <- args[2] 
toxprotblastmetadata <- args[3] 
nontoxprotblast6 <- args[4] 
nontoxprotmetadata <- args[5]
toxvsnontoxIP <- args[6]
toxvsnontoxMF <- args[7]
toxvsnontoxBP <- args[8]
sample <- args[9]
Diamondblast6 <- args[10]


#toxprotblast6, nontoxprotblast6, Diamondblast6, transdf_final_filtered_strict , transdf_final_filtered_lax 
#toxprot metadata #nontoxprotmetadata
#GO entry metadata 

#read in transdf_filtered file 
transdf_final_filtered_lax <-read.csv(transdf_final_filtered_lax_file, header = TRUE)

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
#Function to Parse UniProt blast database 
ParseBlast6UniProt <- function(blast6file, Metadata, prefix) {
  blast6 <- read.delim(blast6file, header = FALSE, sep = "\t")
  colnames(blast6) <- c("Transdecoder_ID", "AccessionNo", "pident", "length", "mismatch", "gapopen", "qstart","qend","sstart","send","evalue","bitscore","qframe", "qcovs")
  blast6[[2]] <- stringr::str_replace(blast6[[2]], "^[^|]*\\|[^|]*\\|", "")
  blast6 <- blast6 %>%
    arrange(desc(bitscore), desc(qcovs),desc(pident)) %>%
    distinct(Transdecoder_ID, .keep_all = TRUE) %>%
    dplyr::select(Transdecoder_ID, AccessionNo,pident,evalue,bitscore,qcovs) 
  metadata <-read.csv(Metadata, sep = "\t") 
    names(metadata) <- gsub(" ", ".", names(metadata))
    metadata <- metadata %>%
    dplyr::rename(AccessionNo = Entry.Name) %>%
    mutate(Protein.names = sub("\\s*\\(.*", "", Protein.names)) %>%
    dplyr::rename(Name = Protein.names) %>%
    tidyr::separate(Organism, into = c("Species", "CommonName"), sep = " \\(", remove = FALSE, extra = "merge", fill = "right") %>%
    mutate(Source = prefix) %>%
    mutate(
      EC_class_number = str_extract(EC.number, "^[0-9]"),
      EC_class_name = ec_class_map[EC_class_number]
    ) %>%
    select(-EC_class_number) %>%
    mutate(Gene_Name = Gene.Names) %>%
    mutate(Gene_Name = word(Gene_Name,1)) %>%
    select(AccessionNo, Name, Gene_Name, Species, Source, EC_class_name, EC.number)
  blast6_metadata <- left_join(blast6,metadata, by = "AccessionNo")
  return(blast6_metadata)
  
}

ParseBlast6UniProt7z <- function(blast6file, Metadata, prefix) {
  blast6 <- read.delim(blast6file, header = FALSE, sep = "\t")
  colnames(blast6) <- c("Transdecoder_ID", "AccessionNo", "pident", "length", "mismatch", "gapopen", "qstart","qend","sstart","send","evalue","bitscore","qframe", "qcovs")
  blast6[[2]] <- stringr::str_replace(blast6[[2]], "^[^|]*\\|[^|]*\\|", "")
  blast6 <- blast6 %>%
    arrange(desc(bitscore), desc(qcovs),desc(pident)) %>%
    distinct(Transdecoder_ID, .keep_all = TRUE) %>%
    dplyr::select(Transdecoder_ID, AccessionNo,pident,evalue,bitscore,qcovs) 
  metadata <-read_tsv(archive_read(Metadata), col_names = TRUE,show_col_types = FALSE) 
  names(metadata) <- gsub(" ", ".", names(metadata))
  metadata <- metadata %>%
    dplyr::rename(AccessionNo = Entry.Name) %>%
    mutate(Protein.names = sub("\\s*\\(.*", "", Protein.names)) %>%
    dplyr::rename(Name = Protein.names) %>%
    tidyr::separate(Organism, into = c("Species", "CommonName"), sep = " \\(", remove = FALSE, extra = "merge", fill = "right") %>%
    mutate(Source = prefix) %>%
    mutate(
      EC_class_number = str_extract(EC.number, "^[0-9]"),
      EC_class_name = ec_class_map[EC_class_number]
    ) %>%
    select(-EC_class_number) %>%
    mutate(Gene_Name = Gene.Names) %>%
    mutate(Gene_Name = word(Gene_Name,1)) %>%
    select(AccessionNo, Name, Gene_Name, Species, Source, EC_class_name, EC.number)
  blast6_metadata <- left_join(blast6,metadata, by = "AccessionNo")
  return(blast6_metadata)
  
}

#dataframe to combine with nontoxprot to select most significant match 
Toxprotblast <- ParseBlast6UniProt(toxprotblast6, toxprotblastmetadata, "ToxProt")
#dataframe to add to final annotation excel 
Toxprotblast_maindf <- Toxprotblast
prefix= "ToxProt"
Toxprotblast_maindf <- Toxprotblast_maindf %>%
  relocate(Name, .after = AccessionNo) %>%
  relocate(Species, .after = Name) %>%
  select(-Gene_Name,-Source,-EC_class_name,-EC.number) 
names(Toxprotblast_maindf)[names(Toxprotblast_maindf) == "AccessionNo"] <- paste0(prefix, "_AccessionNo")
names(Toxprotblast_maindf)[names(Toxprotblast_maindf) == "pident"] <- paste0(prefix, "_pident")
names(Toxprotblast_maindf)[names(Toxprotblast_maindf) == "evalue"] <- paste0(prefix, "_evalue")
names(Toxprotblast_maindf)[names(Toxprotblast_maindf) == "bitscore"] <- paste0(prefix, "_bitscore")
names(Toxprotblast_maindf)[names(Toxprotblast_maindf) == "qcovs"] <- paste0(prefix, "_qcovs")
names(Toxprotblast_maindf)[names(Toxprotblast_maindf) == "Name"] <- paste0(prefix, "_Name")
names(Toxprotblast_maindf)[names(Toxprotblast_maindf) == "Species"] <- paste0(prefix, "_Species")


#dataframe to combine with toxprot to select most significant match 
NonToxprotblast <- ParseBlast6UniProt7z(nontoxprotblast6, nontoxprotmetadata, "NonToxUniProt")
#dataframe to add to the final annotation excel file 
NonToxprotblast_maindf <- NonToxprotblast
prefix= "NonToxUniProt"
NonToxprotblast_maindf <- NonToxprotblast_maindf %>%
  relocate(Name, .after = AccessionNo) %>%
  relocate(Species, .after = Name) %>%
  select(-Gene_Name,-Source,-EC_class_name,-EC.number) 
names(NonToxprotblast_maindf)[names(NonToxprotblast_maindf) == "AccessionNo"] <- paste0(prefix, "_AccessionNo")
names(NonToxprotblast_maindf)[names(NonToxprotblast_maindf) == "pident"] <- paste0(prefix, "_pident")
names(NonToxprotblast_maindf)[names(NonToxprotblast_maindf) == "evalue"] <- paste0(prefix, "_evalue")
names(NonToxprotblast_maindf)[names(Toxprotblast_maindf) == "bitscore"] <- paste0(prefix, "_bitscore")
names(NonToxprotblast_maindf)[names(NonToxprotblast_maindf) == "qcovs"] <- paste0(prefix, "_qcovs")
names(NonToxprotblast_maindf)[names(NonToxprotblast_maindf) == "Name"] <- paste0(prefix, "_Name")
names(NonToxprotblast_maindf)[names(NonToxprotblast_maindf) == "Species"] <- paste0(prefix, "_Species")

#combine dataframes to add to final annotation excel 
blastresults2 <- full_join(Toxprotblast_maindf,NonToxprotblast_maindf, by ="Transdecoder_ID")

#df to find best reviewed hit per transript 
df_to_compare <- rbind(Toxprotblast,NonToxprotblast) %>%
  arrange(desc(bitscore), desc(qcovs),desc(pident)) %>%
  distinct(Transdecoder_ID, .keep_all = TRUE) %>%
  dplyr::select(Transdecoder_ID, Name,Gene_Name,Source,EC_class_name) %>%
  dplyr::rename(Protein_Name = Name, Enzyme_Class = EC_class_name, Annotation_Source= Source)


#Optional argument if available 
#add optional if diamond ends with addign all three blast results to the end of file, if no diamond end with two blast results added per file 

fetch_protein_metadata <- function(accessions, batch_size = 50, sleep_time = 0.4) {
  results <- list()
  n <- length(accessions)
  batches <- split(accessions, ceiling(seq_along(accessions)/batch_size))
  message("Starting retrieval of ", n, " accession numbers...")
  for (i in seq_along(batches)) {
    batch <- batches[[i]]
    message(sprintf("Processing batch %d/%d (%d records)...", 
                    i, length(batches), length(batch)))
    # Fetch XML from NCBI protein database
    xml_data <- tryCatch({
      entrez_fetch(db = "protein",
                   id = batch,
                   rettype = "gb",
                   retmode = "xml")
    }, error = function(e) {
      message("Error in batch ", i, ": ", e$message)
      return(NULL)
    })
    if (is.null(xml_data)) next
    # Parse XML
    xml_parsed <- read_xml(xml_data)
    records <- xml_find_all(xml_parsed, ".//GBSeq")
    batch_df <- map_dfr(records, parse_protein_xml)
    results[[i]] <- batch_df
    message(sprintf("Finished batch %d. Sleeping %.2f sec...", i, sleep_time))
    Sys.sleep(sleep_time)
  }
  message("Retrieval complete.")
  bind_rows(results)
}

parse_protein_xml <- function(xml_doc) {
  get_text <- function(node, xpath) {
    res <- xml_find_first(node, xpath)
    if (inherits(res, "xml_missing")) return(NA_character_)
    xml_text(res)
  }
  AN <- get_text(xml_doc, ".//GBSeq_accession-version")
  protein_name <- get_text(xml_doc, ".//GBSeq_definition")
  organism <- get_text(xml_doc, ".//GBSeq_organism")
  gene_name <- get_text(xml_doc, ".//GBQualifier[GBQualifier_name='gene']/GBQualifier_value")
  ec_number <- get_text(xml_doc, ".//GBQualifier[GBQualifier_name='EC_number']/GBQualifier_value")
  # Return as dataframe row
  data.frame(
  AccessionNo = AN,
    Name = protein_name,
    Species = organism,
    EC.number = ec_number,
    Gene_Name = gene_name,
    stringsAsFactors = FALSE
  )
}

if (!is.na(Diamondblast6) && Diamondblast6 != "NULL" && Diamondblast6 != "") {
  Diamond <- read.delim(Diamondblast6, sep = "\t", header = FALSE)
  colnames(Diamond) <- c("Transdecoder_ID", "AccessionNo", "pident", "length", "mismatch", "gapopen", "qstart","qend","sstart","send","evalue","bitscore","qcovs")
  Diamond <- Diamond %>%
    arrange(desc(bitscore), desc(qcovs),desc(pident)) %>%
    distinct(Transdecoder_ID, .keep_all = TRUE) %>%
    dplyr::select(Transdecoder_ID, AccessionNo,pident,evalue,bitscore,qcovs) %>%
    mutate(Soure = "NCBInr") 
  
  Diamond_AccessionNo_distinct <- Diamond %>%
    distinct(AccessionNo)
  
  AccessionNo_metadata <- fetch_protein_metadata(Diamond_AccessionNo_distinct$AccessionNo)
  Diamond_with_metadata <- left_join(Diamond,AccessionNo_metadata, by = "AccessionNo" )
  Diamond_with_metadata <- Diamond_with_metadata%>%
  mutate(Source = "NCBInr")  %>%
    mutate(
      EC_class_number = str_extract(EC.number, "^[0-9]"),
      EC_class_name = ec_class_map[EC_class_number]
    ) %>%
    select(-EC_class_number) 
  
  Diamond_maindf <- Diamond_with_metadata
  prefix= "NCBInr"
  Diamond_maindf <- Diamond_maindf %>%
    relocate(Name, .after = AccessionNo) %>%
    relocate(Species, .after = Name) %>%
    select(-EC.number)
  names(Diamond_maindf)[names(Diamond_maindf) == "AccessionNo"] <- paste0(prefix, "_AccessionNo")
  names(Diamond_maindf)[names(Diamond_maindf) == "pident"] <- paste0(prefix, "_pident")
  names(Diamond_maindf)[names(Diamond_maindf) == "evalue"] <- paste0(prefix, "_evalue")
  names(Diamond_maindf)[names(Diamond_maindf) == "bitscore"] <- paste0(prefix, "_bitscore")
  names(Diamond_maindf)[names(Diamond_maindf) == "qcovs"] <- paste0(prefix, "_qcovs")
  names(Diamond_maindf)[names(Diamond_maindf) == "Name"] <- paste0(prefix, "_Name")
  names(Diamond_maindf)[names(Diamond_maindf) == "Gene_Name"] <- paste0(prefix, "_Gene_Name")
  names(Diamond_maindf)[names(Diamond_maindf) == "EC_class_name"] <- paste0(prefix, "_EC_class_name")

  names(Diamond_maindf)[names(Diamond_maindf) == "Species"] <- paste0(prefix, "_Species")
  blastresults2 <- full_join(blastresults2,Diamond_maindf, by ="Transdecoder_ID")
  
}

#Domain Labels 

transdf <- read.csv(transdf_final_filtered_lax_file)
toxvsnontoxIP_df <- read.csv(toxvsnontoxIP) %>%
  dplyr::select(Rank, InterPro, Name_short)
toxvsnontoxMF_df <- read.csv(toxvsnontoxMF) %>%
  dplyr::select(Rank, GO_ID, Name)
toxvsnontoxBP_df <- read.csv(toxvsnontoxBP) %>%
  dplyr::select(Rank, GO_ID, Name)

transdf$Domain_Label <- NA
transdf$Domain_Rank <- NA
for (i in seq_len(nrow(toxvsnontoxIP_df))) {
  current_id   <- toxvsnontoxIP_df$InterPro[i]
  current_rank1 <- toxvsnontoxIP_df$Name_short[i]
  current_rank2 <- toxvsnontoxIP_df$Rank[i]
  matches1 <- grepl(current_id, transdf$InterPro_accession_Names, fixed = TRUE) & is.na(transdf$Domain_Label)
  matches2 <- grepl(current_id, transdf$InterPro_accession_Names, fixed = TRUE) & is.na(transdf$Domain_Rank)
  transdf$Domain_Label[matches1] <- current_rank1
  transdf$Domain_Rank[matches2] <- current_rank2
}


transdf$Molecular_Function <- NA
for (i in seq_len(nrow(toxvsnontoxMF_df))) {
  current_id   <- toxvsnontoxMF_df$GO_ID[i]
  current_rank1 <- toxvsnontoxMF_df$Name[i]
  matches1 <- grepl(current_id, transdf$GO_name, fixed = TRUE) & is.na(transdf$Molecular_Function)
  transdf$Molecular_Function[matches1] <- current_rank1
  
}

transdf$Biological_Process <- NA
for (i in seq_len(nrow(toxvsnontoxBP_df))) {
  current_id   <- toxvsnontoxBP_df$GO_ID[i]
  current_rank1 <- toxvsnontoxBP_df$Name[i]
  matches1 <- grepl(current_id, transdf$GO_name, fixed = TRUE) & is.na(transdf$Biological_Process)
  transdf$Biological_Process[matches1] <- current_rank1
}

transdf <- transdf %>%
  dplyr::select(Transdecoder_ID, Domain_Label, Domain_Rank, Molecular_Function, Biological_Process)             

if (!("percent_aggregates" %in% colnames(transdf_final_filtered_lax))) {
  transdf_final_filtered_lax[["percent_aggregates"]] <- NA
}

Annotationdf <- transdf_final_filtered_lax %>%
  left_join(df_to_compare, by ="Transdecoder_ID") %>%
  left_join(transdf,by ="Transdecoder_ID") %>%
  left_join(blastresults2, by = "Transdecoder_ID") 



if (!is.na(Diamondblast6) && Diamondblast6 != "NULL" && Diamondblast6 != "") {
  Annotationdf <- Annotationdf %>%
    mutate(
      Protein_Name = case_when(is.na(Annotation_Source) & !is.na(NCBInr_Name) ~ NCBInr_Name, TRUE ~ Protein_Name),
      Gene_Name = case_when(is.na(Annotation_Source) & !is.na(NCBInr_Name) ~ NCBInr_Gene_Name, TRUE ~ Gene_Name),
      Enzyme_Class = case_when(is.na(Annotation_Source) & !is.na(NCBInr_Name) ~ NCBInr_EC_class_name, TRUE ~ Enzyme_Class),
      Annotation_Source = case_when(is.na(Annotation_Source) & !is.na(NCBInr_Name) ~ Source, TRUE ~ Annotation_Source)
    ) %>%
    dplyr::select(-NCBInr_Gene_Name,-NCBInr_EC_class_name,-Source)
    
Annotationdf$NCBInr_Name <- gsub("\\[.*?\\]", "", Annotationdf$NCBInr_Name)

}

Annotationdf <- Annotationdf %>%
  dplyr::select(Transdecoder_ID,UniqueSequenceName, ORF_type,	SP,Sample_name,Species,Protein_Name,Gene_Name,Enzyme_Class,Annotation_Source,Domain_Label, Molecular_Function, Biological_Process,CysPer,Signal_Length,Mature_Length,PEP_Length,CDS_Length,tpm, est_counts,percent,percent_aggregates,cumulativepercent,Signal_Sequence,Mature_Sequence,PEP_Sequence,CDS_Sequence,InterPro_accession_Names, GO_name,Panther_ID_Name,Phobius_Name,TMHMM,Domain_Rank,everything()) %>%
distinct(UniqueSequenceName, .keep_all = TRUE)




#Final columns: Index(), UniqueSequenceName, Sample_name, Species,
#Protein_Name, Gene_Name, Molecular_Function, Biological_Function, Enzyme_Class, Annotation_Source
#PEP_Length CDS_Length, Signal_Length, mature_length, TPM, percent, cumulativepercent, CysPer, 
#Signal_Sequence, mature_sequence, PEP_Sequence, CDS_Sequence
#InterPro_accession_Names, #GO_name, #Panther_ID_Names, 
#MassSpecColumns, GenomeBlastColumns, NonToxUniProtColumns, #DiamondBlastColumns
#Categorical scores
#AnnotationScore (1,0) #ToxinDomainScore(2,1,0), #ProteomicCoverageScore(2,1,0), #KallistoExpressionScore(2,1,0), #CysPerScore(2,1,0), #BlastScore(2,1,0)
#if genomeblastcolumns are present -> GenomeSupportScore(2,1,0) else GenomeSupportScore(NA)
# 0 <- nothing, #2 <-pident ≥ 85 AND qcovs ≥ 85,  1 <- qcovs ≥ 50 AND pident ≥ 50

#CategoryColumn. 

#if genome is present, all have at least some match to the genome already 
#Category1[Toxin similarity]: AnnoationScore =1 and/or BlastScore = 2
#Category2[Domain similarity]: AnnoationScore =0 and BlastScore <2 and ToxinDomainScore =2 
#Category3[Novel]: No reviewed uniprot match, ToxinDomainScore <2, KallistoExpressionScore >2 and/or ProteomicCoverageScore>2 and/or CysPerScore >2
#4: remaining 
# ADD NA genome and proteome columns if not there 

#Pattern 1 and 2
prepare_patterns <- function(file) {
  summary_IP <- read.csv(file)
  IPInToxins <- summary_IP %>% filter(RelativeExpression > 0)
  Toxin_IPs <- unique(IPInToxins$InterPro)
  pattern1 <- if(length(Toxin_IPs) == 0) "$^" else paste0("\\b(", paste(Toxin_IPs, collapse="|"), ")\\b")
  
  IPOver <- summary_IP %>% filter(RelativeExpression > 1)
  Over_IPs <- unique(IPOver$InterPro)
  pattern2 <- if(length(Over_IPs) == 0) "$^" else paste0("\\b(", paste(Over_IPs, collapse="|"), ")\\b")
  
  return(list(pattern1 = pattern1, pattern2 = pattern2))
}

patterns <- prepare_patterns(toxvsnontoxIP)


# Add column with NA values if it doesn't exist
if (!("Top" %in% colnames(Annotationdf))) {
  Annotationdf[["Top"]] <- NA
  Annotationdf[["Coverage"]] <- NA
  Annotationdf[["Unique"]] <- NA
}

if (!("genome_qcovs" %in% colnames(Annotationdf))) {
  Annotationdf[["genome_qcovs"]] <- NA
  Annotationdf[["genome_pident"]] <- NA
}

if (!("Minimap_mapq" %in% colnames(Annotationdf))) {
  Annotationdf[["Minimap_mapq"]] <- NA
}



#Adding scoring 
library(dplyr)
library(stringr)

Annotationdf <- Annotationdf %>%
  mutate(
    AnnotationScore = case_when(
      Annotation_Source == "ToxProt" ~ 1,
      TRUE ~ 0
    )
  ) %>%
  mutate(
    ToxinDomainScore = case_when(
      str_detect(InterPro_accession_Names, patterns$pattern1) ~ 2,
      str_detect(InterPro_accession_Names, patterns$pattern2) ~ 1,
      TRUE ~ 0
    )
  ) %>%
  mutate(
    ProteomicCoverageScore = case_when(
      Top == "TRUE" & Coverage >= 50 ~ 2,
      Top == "TRUE" & Unique >= 1 ~ 1,
      is.na(Top) ~ NA_real_,
      TRUE ~ 0
    )
  ) %>%
  mutate(
    KallistoExpressionScore = case_when(
      percent >= 1 | percent_aggregates >= 1 ~ 2,
      percent > 0 | percent_aggregates > 0 ~ 1,
      TRUE ~ 0
    )
  ) %>%
  mutate(
    CysPerScore = case_when(
      CysPer >= 5 ~ 2,
      CysPer >= 1 ~ 1,
      TRUE ~ 0
    )
  ) %>%
  mutate(
    BlastScore = case_when(
      ToxProt_bitscore >= 250 ~ 2,
      ToxProt_bitscore >= 50 ~ 1,
      TRUE ~ 0
    )
  ) %>%
  mutate(
    GenomeSupportScore = case_when(
      genome_qcovs == 100 & genome_pident == 100 ~ 3,
      genome_qcovs >= 90 & genome_pident >= 90 ~ 2,
      genome_qcovs >= 80 & genome_pident >= 80 ~ 1,
      is.na(genome_qcovs) ~ NA_real_,
      TRUE ~ 0
    )
  )

# add category label 

Annotationdf <- Annotationdf %>%
  mutate(Category = case_when(
(GenomeSupportScore >= 2 |is.na(GenomeSupportScore)) &(AnnotationScore == 1 | BlastScore == 2) ~ "Category 1",
(GenomeSupportScore >= 2 |is.na(GenomeSupportScore)) & ToxinDomainScore == 2 ~ "Category 2",
    (is.na(Annotation_Source) | Annotation_Source == "NCBInr") &
      (GenomeSupportScore >= 1 |is.na(GenomeSupportScore)) &
      (CysPerScore > 0 |
         KallistoExpressionScore > 0 |
         ProteomicCoverageScore > 0 | GenomeSupportScore == 3 ) ~ "Category 3",
(Annotation_Source == "NonToxUniProt")  ~ "Putative non-toxins",
    TRUE ~ "Lower genome support transcripts"
  ))

Annotationdf$Protein_Name <- gsub("\\[.*?\\]", "", Annotationdf$Protein_Name)

Annotationdf <- Annotationdf %>%
  mutate(across(everything(), ~str_remove_all(., ","))) %>%
  select(-Code,	-Percentage_Identity,	-E_value,	-BitScore,	-Hit_species)

write.csv(Annotationdf, paste0(sample, "_all_Annotated_df.csv"), row.names = FALSE)

Annotationdf <- Annotationdf %>%
  filter(Category != "Lower genome support transcripts" & Category != "Putative non-toxins")
  
write.csv(Annotationdf, paste0(sample, "_Select_Annotated_df.csv"), row.names = FALSE)

ProtSpaceAnnoation <- Annotationdf %>%
  dplyr::select(Transdecoder_ID,UniqueSequenceName,PEP_Length,Signal_Length,Mature_Length,CysPer,tpm,Hit,Sample_name,Species,Protein_Name,Gene_Name,Enzyme_Class,Annotation_Source,Domain_Label, Molecular_Function, Biological_Process,AnnotationScore,ToxinDomainScore,ProteomicCoverageScore,KallistoExpressionScore,CysPerScore,BlastScore,GenomeSupportScore,Category,Filter, Minimap_mapq)

write.csv(ProtSpaceAnnoation, paste0(sample, "_ProtSpaceAnnotation.csv"), row.names = FALSE)

Pepseq <- AAStringSet(Annotationdf$PEP_Sequence)
names(Pepseq) <- Annotationdf$Transdecoder_ID
writeXStringSet(Pepseq, paste0(sample,"_ProtSpacePEP.fasta"))

cdsseq <- DNAStringSet(Annotationdf$CDS_Sequence)
names(cdsseq) <- Annotationdf$Transdecoder_ID
writeXStringSet(cdsseq, paste0(sample,"_ProtSpaceCDS.fasta"))


