#!/usr/bin/env Rscript

#Integrate mass spec and genomic data if available 
#generate filtered lists and venn diagrams 

#Strict
#cysper > 5%
#toxin overexpressed domain present
#Blast > 250
#ke > 1
#mass spec >50 
#Lax
#cysper > 1%
#any toxin domain present
#blast > 50 
#ke > 0 
#mass spec >0
#all distinct, distinct cds sequences, each ID appears once only, Complete ORF and signalP transcripts only.
#unfiltered: contains all columns of genome blast and mass if available, contain transcripts with overlapping regions, save everything no other criteria applied
#filtered:, distinct cds sequences

#for those with genomes, cd-hit was only run on transcriptomes not ORF, transcripts that overlap the same genomic locus are clustered only best blast match is kept
#for those without genomes, cd-hit was on both transcriptomes and ORFs, before running kallisto_transdecoder, only distinct pep will be kept there will likely be redudancy/isoforms that need to be handled downstream

library(dplyr)
library(tidyr)
library(stringr)
library(ggvenn)
library(ggplot2)
library(GenomicRanges)
library(igraph)


args <- commandArgs(trailingOnly = TRUE)
Sample_name <- args[1] #sample name 
species_name  <- args[2] #species name
Transdf_distinct_file <- args[3] #transdf_distinct
ToxnontoxIP <- args[4]
Blastn_result <- args[5] #blastn_result 
mass_spec <- args[6] #massspec csv 

#read in transdf_distinct 
transdf <- read.csv(Transdf_distinct_file, header = TRUE )
#add sample and species name to df and move to the end 
transdf$Species <- species_name
transdf$Sample_name <- Sample_name
transdf <- transdf[c("Species", "Sample_name", setdiff(names(transdf), c("Species", "Sample_name")))]
#addUniqueNameColumn for when multiple datasets might be binded downstream 
transdf$UniqueSequenceName = paste(Sample_name, transdf$Transdecoder_ID, sep = "_")

#should already have been removed but just in case  duplicates in CDS, keep one with higher BitScore to ToxinProtein (should already have been made distinct)
transdf_distinct <- transdf[order(-transdf$BitScore, transdf$E_value, -transdf$tpm ), ]
transdf_distinct <- distinct(transdf_distinct, Transdecoder_ID, .keep_all = TRUE)
transdf_distinct <- transdf_distinct %>%
  group_by(CDS_Sequence) %>%
  mutate(percent = sum(percent, na.rm = TRUE)) %>%
  mutate(tpm = sum(tpm, na.rm = TRUE)) %>%
  ungroup() %>%
  distinct(Transdecoder_ID, CDS_Sequence, .keep_all = TRUE)

transdf_unfiltered <- transdf_distinct
transdf_filtered <- transdf_distinct

#readinToxinCSV 
summary_IP <- read.csv(ToxnontoxIP)
IPInToxins <- summary_IP %>%
  filter(RelativeExpression > 0 )

Toxin_IPs <- data.frame(unique(IPInToxins$InterPro)) 
colnames(Toxin_IPs)[colnames(Toxin_IPs) == "unique.IPInToxins.InterPro."] <- "IP"
pattern <- paste0("\\b(", paste(Toxin_IPs$IP, collapse = "|"), ")\\b")
if (length(Toxin_IPs$IP) == 0) {
  pattern <- "$^"  
}
IPOverRepresentedInToxins <- summary_IP %>%
  filter(RelativeExpression > 1 )

OverExpressedIPs <- data.frame(unique(IPOverRepresentedInToxins$InterPro)) 
colnames(OverExpressedIPs)[colnames(OverExpressedIPs) == "unique.IPOverRepresentedInToxins.InterPro."] <- "IP"
pattern2 <- paste0("\\b(", paste(OverExpressedIPs$IP, collapse = "|"), ")\\b")
if (length(OverExpressedIPs$IP) == 0) {
  pattern2 <- "$^"  
}

#Mass Columns expected 
#Accession, Coverage, Top, Unique



if (
  mass_spec != "NULL" &&
  Blastn_result != "NULL"
)
{
  #read in mass spec   #rename mass spec column 
  mass_spec <- read.csv(mass_spec, header = TRUE) 
  colnames(mass_spec)[which(names(mass_spec) == "Accession")] <- "Transdecoder_ID" 
  #left_join full massspec 
  transdf_unfiltered <- left_join(transdf_unfiltered,mass_spec, by = "Transdecoder_ID")
  mass_spec <- mass_spec %>%
    dplyr::select(Transdecoder_ID,Top,Coverage,Unique)
  transdf_filtered <- left_join(transdf_filtered,mass_spec, by = "Transdecoder_ID")
  #read in Blastn   #rename blastn columns 
  Blastn_result_read <- read.table(Blastn_result,sep = "\t",header = FALSE,stringsAsFactors = FALSE)  
  colnames(Blastn_result_read) <- c("Transdecoder_ID", "genome_sseqid", "genome_pident", "genome_length", "genome_mismatch","genome_gapopen","genome_qstart","genome_qend","genome_sstart","genome_send", "genome_sstrand","genome_evalue", "genome_bitscore", "genome_qframe", "genome_qcovs")
  #order by genome_bitscore   #only keep distinct hit per transcript, keeping the hit  with higher bitscore
  Blastn_result_read_genome_bitscore <- Blastn_result_read[order(Blastn_result_read$genome_bitscore, decreasing = TRUE), ]
  Blastn_result_read_genome_bitscore_distinct <- Blastn_result_read_genome_bitscore %>% distinct(Transdecoder_ID, .keep_all = TRUE)
  #left join blastn
  transdf_unfiltered <- left_join(transdf_unfiltered,Blastn_result_read_genome_bitscore_distinct, by = "Transdecoder_ID")
  transdf_filtered <- left_join(transdf_filtered, Blastn_result_read_genome_bitscore_distinct, by = "Transdecoder_ID")
  
  #save unfiltered csv 
  write.csv(transdf_unfiltered, paste0(Sample_name, "_transdf_distinct_final_unfiltered.csv"), row.names = FALSE)
  
  transdf_blastn <- transdf_blastn %>%
    mutate(genome_sstrand = case_when(
      genome_sstrand == "plus" ~ "+",
      genome_sstrand == "minus" ~ "-"
    )) %>%
    filter(!is.na(genome_sseqid))
  
  #Filtering genome 
  #create genome range object
  gr <- GRanges(
    seqnames = trimws(as.character(transdf_blastn$genome_sseqid)),
    ranges = IRanges(
      start = pmin(transdf_blastn$genome_sstart, transdf_blastn$genome_send),
      end   = pmax(transdf_blastn$genome_sstart, transdf_blastn$genome_send)
    ),
    query_id = transdf_blastn$Transdecoder_ID,
    strand = transdf_blastn$genome_sstrand
  )
  
  #find overalps
  hits <- findOverlaps(gr, gr, ignore.strand = FALSE)
  #create dataframe of overlaps
  edges <- data.frame(
    from = queryHits(hits),
    to   = subjectHits(hits)
  )
  #identify clusters
  g <- graph_from_data_frame(edges, directed = FALSE)
  clusters <- components(g)$membership
  transdf_blastn$cluster <- clusters[seq_along(transdf_blastn$Transdecoder_ID)]
  
  transdf_blastn <- transdf_blastn%>%
    group_by(cluster) %>%
    mutate(tpm_aggregates = sum(tpm)) %>%
    relocate(tpm_aggregates, .after = CysPer)
  
  transdf_blastn <- transdf_blastn %>% arrange (
    desc(genome_qcovs),
    desc(genome_pident),
    desc(genome_bitscore)
  ) %>%
    distinct(cluster, .keep_all = TRUE) %>%
    select(-cluster)
  
  
  #only those with genome matches 
  
  #Filtering 
  #completeORF + SignalP + no Transmembrane
  Base <- transdf_blastn %>%
    filter(
      grepl("complete", ORF_type, ignore.case = TRUE),
      grepl("SP", SP, ignore.case = TRUE),
      TMHMM == FALSE
    )
  #changing criteria to numeric for filtering 
  Base[["Coverage"]] <- as.numeric(Base[["Coverage"]])
  Base[["BitScore"]] <- as.numeric(Base[["BitScore"]])
  Base[["percent"]] <- as.numeric(Base[["percent"]])
  Base[["CysPer"]] <- as.numeric(Base[["CysPer"]])
  ##VennDiagram & Filtering (Strict)
  matching_rows <- Base[
    !is.na(Base$InterPro_accession_Names) &
      grepl(pattern2, Base$InterPro_accession_Names),
  ]
  set_A<-matching_rows
  set_B <- Base[Base[["Coverage"]] >= 50 & !is.na(Base[["Coverage"]]) & Base[["Top"]] == TRUE, ]
  set_C <- Base[Base[["BitScore"]] >= 200 & !is.na(Base[["BitScore"]]), ]
  set_D <-Base[Base[["percent"]] >= 1 & !is.na(Base[["percent"]]), ]
  set_E <-Base[Base[["CysPer"]] >= 5 & !is.na(Base[["CysPer"]]), ]
  venn_list <- list(
    TD = set_A$Transdecoder_ID,
    MS = set_B$Transdecoder_ID,
    TP = set_C$Transdecoder_ID,
    KE = set_D$Transdecoder_ID,
    CP = set_E$Transdecoder_ID
    
  )
  p <-ggvenn(venn_list,c("TD", "MS","TP","KE", "CP"), fill_color = c("#E41A1C", "#377EB8", "#4DAF4A", "#EBAC4D","#782DC8" ))
  Union_ABC <- Reduce(union, list(set_A$Transdecoder_ID, set_B$Transdecoder_ID, set_C$Transdecoder_ID, set_D$Transdecoder_ID,set_E$Transdecoder_ID))
  ggsave(paste0(Sample_name,"_Venn_strict.png"), plot = p, width = 6, height = 4, dpi = 300)
  Base_union <- Base[Base$Transdecoder_ID %in% Union_ABC, ]  %>%
    dplyr::mutate(Filter = "Strict")
  transdf_final_filtered_strict <- Base_union %>%
    select(-genome_length, -genome_gapopen, -genome_qstart,-genome_qend,-genome_sstart,-genome_send,-genome_sstrand,-genome_qframe)
  write.csv(transdf_final_filtered_strict, paste0(Sample_name, "_transdf_distinct_final_filtered_strict.csv"), row.names = FALSE)
  
  ##VennDiagram & Filtering (Lax)
  matching_rows <- Base[
    !is.na(Base$InterPro_accession_Names) &
      grepl(pattern, Base$InterPro_accession_Names),
  ]
  set_A<-matching_rows
  set_B <- Base[Base[["Coverage"]] >= 0 & !is.na(Base[["Coverage"]]) & Base[["Top"]] == TRUE, ]
  set_C <- Base[Base[["BitScore"]] >= 50 & !is.na(Base[["BitScore"]]), ]
  set_D <-Base[Base[["percent"]] >= 0 & !is.na(Base[["percent"]]), ]
  set_E <-Base[Base[["CysPer"]] >= 1 & !is.na(Base[["CysPer"]]), ]
  venn_list <- list(
    TD = set_A$Transdecoder_ID,
    MS = set_B$Transdecoder_ID,
    TP = set_C$Transdecoder_ID,
    KE = set_D$Transdecoder_ID,
    CP = set_E$Transdecoder_ID
    
  )
  p <- ggvenn(venn_list,c("TD", "MS","TP","KE", "CP"), fill_color = c("#E41A1C", "#377EB8", "#4DAF4A", "#EBAC4D","#782DC8" ))
  Union_ABC <- Reduce(union, list(set_A$Transdecoder_ID, set_B$Transdecoder_ID, set_C$Transdecoder_ID, set_D$Transdecoder_ID,set_E$Transdecoder_ID))
  ggsave(paste0(Sample_name,"_Venn_lax.png"), plot = p, width = 6, height = 4, dpi = 300)
  Base_union <- Base[Base$Transdecoder_ID %in% Union_ABC, ] 
  transdf_final_filtered_lax <- transdf_final_filtered_strict %>%
    select(UniqueSequenceName, Filter) %>% 
    right_join(Base_union, by = "UniqueSequenceName") %>%
    mutate(Filter = case_when(
      Filter == "Strict" ~ "Strict",
      TRUE ~ "Lax"
    ))  %>%
    select(-genome_length, -genome_gapopen, -genome_qstart,-genome_qend,-genome_sstart,-genome_send,-genome_sstrand,-genome_qframe)
  
  write.csv(transdf_final_filtered_lax, paste0(Sample_name, "_transdf_distinct_final_filtered_lax.csv"), row.names = FALSE)
}

else if (  mass_spec == "NULL" &&
       Blastn_result == "NULL")
{
  #just save the modified transdf 
  write.csv(transdf, paste0(Sample_name, "_transdf_distinct_final_unfiltered.csv"), row.names = FALSE)
  #Filtering 
  
  #Filtering 
  transdf_distinct <- transdf[order(-transdf$BitScore, transdf$E_value, -transdf$tpm ), ]
  transdf_distinct <- distinct(transdf_distinct, PEP_Sequence, .keep_all = TRUE)
  
  #completeORF + SignalP + no Transmembrane
  Base <- transdf_distinct %>%
    filter(
      grepl("complete", ORF_type, ignore.case = TRUE),
      grepl("SP", SP, ignore.case = TRUE),
      TMHMM == FALSE
    )
  #changing criteria to numeric for filtering 
  Base[["BitScore"]] <- as.numeric(Base[["BitScore"]])
  Base[["percent"]] <- as.numeric(Base[["percent"]])
  Base[["CysPer"]] <- as.numeric(Base[["CysPer"]])
  ##VennDiagram & Filtering (Strict)
  matching_rows <- Base[
    !is.na(Base$InterPro_accession_Names) &
      grepl(pattern2, Base$InterPro_accession_Names),
  ]
  set_A<-matching_rows
  set_C <- Base[Base[["BitScore"]] >= 200 & !is.na(Base[["BitScore"]]), ]
  set_D <-Base[Base[["percent"]] >= 1 & !is.na(Base[["percent"]]), ]
  set_E <-Base[Base[["CysPer"]] >= 5 & !is.na(Base[["CysPer"]]), ]
  venn_list <- list(
    TD = set_A$Transdecoder_ID,
    TP = set_C$Transdecoder_ID,
    KE = set_D$Transdecoder_ID,
    CP = set_E$Transdecoder_ID
    
  )
  p <-ggvenn(venn_list,c("TD","TP","KE", "CP"), fill_color = c("#E41A1C", "#4DAF4A", "#EBAC4D","#782DC8" ))
  Union_ABC <- Reduce(union, list(set_A$Transdecoder_ID, set_C$Transdecoder_ID, set_D$Transdecoder_ID,set_E$Transdecoder_ID))
  ggsave(paste0(Sample_name,"_Venn_strict.png"), plot = p, width = 6, height = 4, dpi = 300)
  Base_union <- Base[Base$Transdecoder_ID %in% Union_ABC, ]  %>%
    dplyr::mutate(Filter = "Strict")
  transdf_final_filtered_strict <- Base_union %>%
    select(-genome_length, -genome_gapopen, -genome_qstart,-genome_qend,-genome_sstart,-genome_send,-genome_sstrand,-genome_qframe)
  
  
  write.csv(transdf_final_filtered_strict, paste0(Sample_name, "_transdf_distinct_final_filtered_strict.csv"), row.names = FALSE)
  
  ##VennDiagram & Filtering (Lax)
  matching_rows <- Base[
    !is.na(Base$InterPro_accession_Names) &
      grepl(pattern, Base$InterPro_accession_Names),
  ]
  set_A<-matching_rows
  set_C <- Base[Base[["BitScore"]] >= 50 & !is.na(Base[["BitScore"]]), ]
  set_D <-Base[Base[["percent"]] >= 0 & !is.na(Base[["percent"]]), ]
  set_E <-Base[Base[["CysPer"]] >= 1 & !is.na(Base[["CysPer"]]), ]
  venn_list <- list(
    TD = set_A$Transdecoder_ID,
    TP = set_C$Transdecoder_ID,
    KE = set_D$Transdecoder_ID,
    CP = set_E$Transdecoder_ID
    
  )
  p <- ggvenn(venn_list,c("TD","TP","KE", "CP"), fill_color = c("#E41A1C", "#4DAF4A", "#EBAC4D","#782DC8" ))
  Union_ABC <- Reduce(union, list(set_A$Transdecoder_ID, set_C$Transdecoder_ID, set_D$Transdecoder_ID,set_E$Transdecoder_ID))
  ggsave(paste0(Sample_name,"_Venn_lax.png"), plot = p, width = 6, height = 4, dpi = 300)
  Base_union <- Base[Base$Transdecoder_ID %in% Union_ABC, ]
  transdf_final_filtered_lax <- transdf_final_filtered_strict %>%
    select(UniqueSequenceName, Filter) %>% 
    right_join(Base_union, by = "UniqueSequenceName") %>%
    mutate(Filter = case_when(
      Filter == "Strict" ~ "Strict",
      TRUE ~ "Lax"
    ))  %>%
    select(-genome_length, -genome_gapopen, -genome_qstart,-genome_qend,-genome_sstart,-genome_send,-genome_sstrand,-genome_qframe)
  
  
  write.csv(transdf_final_filtered_lax, paste0(Sample_name, "_transdf_distinct_final_filtered_lax.csv"), row.names = FALSE)

}

else if ( mass_spec == "NULL" &&
        Blastn_result != "NULL"
) {

  #read in Blastn   #rename blastn columns 
  Blastn_result_read <- read.table(Blastn_result,sep = "\t",header = FALSE,stringsAsFactors = FALSE)  
  colnames(Blastn_result_read) <- c("Transdecoder_ID", "genome_sseqid", "genome_pident", "genome_length", "genome_mismatch","genome_gapopen","genome_qstart","genome_qend","genome_sstart","genome_send", "genome_sstrand","genome_evalue", "genome_bitscore", "genome_qframe", "genome_qcovs")
  #order by genome_bitscore   #only keep distinct hit per transcript, keeping the hit  with higher bitscore
  Blastn_result_read_genome_qcovs <- Blastn_result_read[order(Blastn_result_read$genome_bitscore, decreasing = TRUE), ]
  Blastn_result_read_genome_bitscore_distinct <- Blastn_result_read_genome_qcovs %>% distinct(Transdecoder_ID, .keep_all = TRUE)
  #left join blastn
  transdf_massspec_blastn <- left_join(transdf,Blastn_result_read_genome_bitscore_distinct, by = "Transdecoder_ID")
  transdf_blastn <- left_join(transdf_distinct, Blastn_result_read_genome_bitscore_distinct, by = "Transdecoder_ID")
  
  #save unfiltered csv 
  write.csv(transdf_massspec_blastn, paste0(Sample_name, "_transdf_distinct_final_unfiltered.csv"), row.names = FALSE)
  
  transdf_blastn <- transdf_blastn %>%
    mutate(genome_sstrand = case_when(
      genome_sstrand == "plus" ~ "+",
      genome_sstrand == "minus" ~ "-"
    )) %>%
    filter(!is.na(genome_sseqid))
  
  #Filtering genome 
  #create genome range object
  gr <- GRanges(
    seqnames = trimws(as.character(transdf_blastn$genome_sseqid)),
    ranges = IRanges(
      start = pmin(transdf_blastn$genome_sstart, transdf_blastn$genome_send),
      end   = pmax(transdf_blastn$genome_sstart, transdf_blastn$genome_send)
    ),
    query_id = transdf_blastn$Transdecoder_ID,
    strand = transdf_blastn$genome_sstrand
  )
  
  #find overalps
  hits <- findOverlaps(gr, gr, ignore.strand = FALSE)
  #create dataframe of overlaps
  edges <- data.frame(
    from = queryHits(hits),
    to   = subjectHits(hits)
  )
  #identify clusters
  g <- graph_from_data_frame(edges, directed = FALSE)
  clusters <- components(g)$membership
  transdf_blastn$cluster <- clusters[seq_along(transdf_blastn$Transdecoder_ID)]
  
  transdf_blastn <- transdf_blastn%>%
    group_by(cluster) %>%
    mutate(tpm_aggregates = sum(tpm)) %>%
    relocate(tpm_aggregates, .after = CysPer)
  
  
  transdf_blastn <- transdf_blastn %>% arrange (
    desc(genome_qcovs),
    desc(genome_pident),
    desc(genome_bitscore)
  ) %>%
    distinct(cluster, .keep_all = TRUE) %>%
    select(-cluster)
  
  
  #only those with genome matches 
  
  #Filtering 
  #completeORF + SignalP + no Transmembrane
  Base <- transdf_blastn %>%
    filter(
      grepl("complete", ORF_type, ignore.case = TRUE),
      grepl("SP", SP, ignore.case = TRUE),
      TMHMM == FALSE
    )
  #changing criteria to numeric for filtering 
  Base[["BitScore"]] <- as.numeric(Base[["BitScore"]])
  Base[["percent"]] <- as.numeric(Base[["percent"]])
  Base[["CysPer"]] <- as.numeric(Base[["CysPer"]])
  ##VennDiagram & Filtering (Strict)
  matching_rows <- Base[
    !is.na(Base$InterPro_accession_Names) &
      grepl(pattern2, Base$InterPro_accession_Names),
  ]
  set_A<-matching_rows
  set_C <- Base[Base[["BitScore"]] >= 200 & !is.na(Base[["BitScore"]]), ]
  set_D <-Base[Base[["percent"]] >= 1 & !is.na(Base[["percent"]]), ]
  set_E <-Base[Base[["CysPer"]] >= 5 & !is.na(Base[["CysPer"]]), ]
  venn_list <- list(
    TD = set_A$Transdecoder_ID,
    TP = set_C$Transdecoder_ID,
    KE = set_D$Transdecoder_ID,
    CP = set_E$Transdecoder_ID
    
  )
  p <-ggvenn(venn_list,c("TD","TP","KE", "CP"), fill_color = c("#E41A1C", "#4DAF4A", "#EBAC4D","#782DC8" ))
  Union_ABC <- Reduce(union, list(set_A$Transdecoder_ID, set_C$Transdecoder_ID, set_D$Transdecoder_ID,set_E$Transdecoder_ID))
  ggsave(paste0(Sample_name,"_Venn_strict.png"), plot = p, width = 6, height = 4, dpi = 300)
  Base_union <- Base[Base$Transdecoder_ID %in% Union_ABC, ]  %>%
    dplyr::mutate(Filter = "Strict")
  transdf_final_filtered_strict <- Base_union %>%
    select(-genome_length, -genome_gapopen, -genome_qstart,-genome_qend,-genome_sstart,-genome_send,-genome_sstrand,-genome_qframe)
  
  
  write.csv(transdf_final_filtered_strict, paste0(Sample_name, "_transdf_distinct_final_filtered_strict.csv"), row.names = FALSE)
  
  ##VennDiagram & Filtering (Lax)
  matching_rows <- Base[
    !is.na(Base$InterPro_accession_Names) &
      grepl(pattern, Base$InterPro_accession_Names),
  ]
  set_A<-matching_rows
  set_C <- Base[Base[["BitScore"]] >= 50 & !is.na(Base[["BitScore"]]), ]
  set_D <-Base[Base[["percent"]] >= 0 & !is.na(Base[["percent"]]), ]
  set_E <-Base[Base[["CysPer"]] >= 1 & !is.na(Base[["CysPer"]]), ]
  venn_list <- list(
    TD = set_A$Transdecoder_ID,
    TP = set_C$Transdecoder_ID,
    KE = set_D$Transdecoder_ID,
    CP = set_E$Transdecoder_ID
    
  )
  p <- ggvenn(venn_list,c("TD","TP","KE", "CP"), fill_color = c("#E41A1C", "#4DAF4A", "#EBAC4D","#782DC8" ))
  Union_ABC <- Reduce(union, list(set_A$Transdecoder_ID, set_C$Transdecoder_ID, set_D$Transdecoder_ID,set_E$Transdecoder_ID))
  ggsave(paste0(Sample_name,"_Venn_lax.png"), plot = p, width = 6, height = 4, dpi = 300)
  Base_union <- Base[Base$Transdecoder_ID %in% Union_ABC, ]
  
  transdf_final_filtered_lax <- transdf_final_filtered_strict %>%
    select(UniqueSequenceName, Filter) %>% 
    right_join(Base_union, by = "UniqueSequenceName") %>%
    mutate(Filter = case_when(
      Filter == "Strict" ~ "Strict",
      TRUE ~ "Lax"
    ))  %>%
    select(-genome_length, -genome_gapopen, -genome_qstart,-genome_qend,-genome_sstart,-genome_send,-genome_sstrand,-genome_qframe)
  
  write.csv(transdf_final_filtered_lax, paste0(Sample_name, "_transdf_distinct_final_filtered_lax.csv"), row.names = FALSE)
}

else if (
  mass_spec != "NULL" &&
  Blastn_result == "NULL"
) {
  #Transdf_distinct  filter
  transdf_distinct <- transdf[order(-transdf$BitScore, transdf$E_value, -transdf$tpm ), ]
  transdf_distinct <- distinct(transdf_distinct, PEP_Sequence, .keep_all = TRUE)
  
  #read in mass spec   #rename mass spec column 
  mass_spec <- read.csv(mass_spec, header = TRUE) 
  colnames(mass_spec)[which(names(mass_spec) == "Accession")] <- "Transdecoder_ID" 
  #left_join full massspec 
  transdf_massspec <- left_join(transdf,mass_spec, by = "Transdecoder_ID")
  transdf_distinct_massspec <- left_join(transdf_distinct,mass_spec, by = "Transdecoder_ID")

  #save unfiltered csv 
  write.csv(transdf_massspec, paste0(Sample_name, "_transdf_distinct_final_unfiltered.csv"), row.names = FALSE)
  


  #Filtering 
  #completeORF + SignalP + no Transmembrane
  Base <- transdf_distinct_massspec %>%
    filter(
      grepl("complete", ORF_type, ignore.case = TRUE),
      grepl("SP", SP, ignore.case = TRUE),
      TMHMM == FALSE
    )
  #changing criteria to numeric for filtering 
  Base[["Coverage"]] <- as.numeric(Base[["Coverage"]])
  Base[["BitScore"]] <- as.numeric(Base[["BitScore"]])
  Base[["percent"]] <- as.numeric(Base[["percent"]])
  Base[["CysPer"]] <- as.numeric(Base[["CysPer"]])
  ##VennDiagram & Filtering (Strict)
  matching_rows <- Base[
    !is.na(Base$InterPro_accession_Names) &
      grepl(pattern2, Base$InterPro_accession_Names),
  ]
  set_A<-matching_rows
  set_B <- Base[Base[["Coverage"]] >= 50 & !is.na(Base[["Coverage"]]) & Base[["Top"]] == TRUE, ]
  set_C <- Base[Base[["BitScore"]] >= 250 & !is.na(Base[["BitScore"]]), ]
  set_D <-Base[Base[["percent"]] >= 1 & !is.na(Base[["percent"]]), ]
  set_E <-Base[Base[["CysPer"]] >= 5 & !is.na(Base[["CysPer"]]), ]
  venn_list <- list(
    TD = set_A$Transdecoder_ID,
    MS = set_B$Transdecoder_ID,
    TP = set_C$Transdecoder_ID,
    KE = set_D$Transdecoder_ID,
    CP = set_E$Transdecoder_ID
    
  )
  p <-ggvenn(venn_list,c("TD", "MS","TP","KE", "CP"), fill_color = c("#E41A1C", "#377EB8", "#4DAF4A", "#EBAC4D","#782DC8" ))
  Union_ABC <- Reduce(union, list(set_A$Transdecoder_ID, set_B$Transdecoder_ID, set_C$Transdecoder_ID, set_D$Transdecoder_ID,set_E$Transdecoder_ID))
  ggsave(paste0(Sample_name,"_Venn_strict.png"), plot = p, width = 6, height = 4, dpi = 300)
  Base_union <- Base[Base$Transdecoder_ID %in% Union_ABC, ]  %>%
    dplyr::mutate(Filter = "Strict")
  transdf_final_filtered_strict <- Base_union %>%
    select(-genome_length, -genome_gapopen, -genome_qstart,-genome_qend,-genome_sstart,-genome_send,-genome_sstrand,-genome_qframe)
  
  write.csv(transdf_final_filtered_strict, paste0(Sample_name, "_transdf_distinct_final_filtered_strict.csv"), row.names = FALSE)
  
  ##VennDiagram & Filtering (Lax)
  matching_rows <- Base[
    !is.na(Base$InterPro_accession_Names) &
      grepl(pattern, Base$InterPro_accession_Names),
  ]
  set_A<-matching_rows
  set_B <- Base[Base[["Coverage"]] >= 0 & !is.na(Base[["Coverage"]]) & Base[["Top"]] == TRUE, ]
  set_C <- Base[Base[["BitScore"]] >= 50 & !is.na(Base[["BitScore"]]), ]
  set_D <-Base[Base[["percent"]] >= 0 & !is.na(Base[["percent"]]), ]
  set_E <-Base[Base[["CysPer"]] >= 1 & !is.na(Base[["CysPer"]]), ]
  venn_list <- list(
    TD = set_A$Transdecoder_ID,
    MS = set_B$Transdecoder_ID,
    TP = set_C$Transdecoder_ID,
    KE = set_D$Transdecoder_ID,
    CP = set_E$Transdecoder_ID
    
  )
  p <- ggvenn(venn_list,c("TD", "MS","TP","KE", "CP"), fill_color = c("#E41A1C", "#377EB8", "#4DAF4A", "#EBAC4D","#782DC8" ))
  Union_ABC <- Reduce(union, list(set_A$Transdecoder_ID, set_B$Transdecoder_ID, set_C$Transdecoder_ID, set_D$Transdecoder_ID,set_E$Transdecoder_ID))
  ggsave(paste0(Sample_name,"_Venn_lax.png"), plot = p, width = 6, height = 4, dpi = 300)
  Base_union <- Base[Base$Transdecoder_ID %in% Union_ABC, ]
  transdf_final_filtered_lax <- transdf_final_filtered_strict %>%
    select(UniqueSequenceName, Filter) %>% 
    right_join(Base_union, by = "UniqueSequenceName") %>%
    mutate(Filter = case_when(
      Filter == "Strict" ~ "Strict",
      TRUE ~ "Lax"
    ))  %>%
    select(-genome_length, -genome_gapopen, -genome_qstart,-genome_qend,-genome_sstart,-genome_send,-genome_sstrand,-genome_qframe)
  
  
  write.csv(transdf_final_filtered_lax, paste0(Sample_name, "_transdf_distinct_final_filtered_lax.csv"), row.names = FALSE)
}




