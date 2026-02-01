library(dplyr)
library(tidyr)
library(stringr)
library(ggvenn)
library(ggplot2)
args <- commandArgs(trailingOnly = TRUE)
Sample_name <- args[1] #sample name 
species_name  <- args[2] #species name
Transdf <- args[3] #transdf_distinct
ToxinDataTSV <- args[4]
Blastn_result <- args[5] #blastn_result 
mass_spec <- args[6] #massspec csv 


#read in transdf_distinct 
transdf <- read.csv(Transdf, header = TRUE )
#add sample and species name to df and move to the end 
transdf$Species <- species_name
transdf$Sample_name <- Sample_name
transdf <- transdf[c("Species", "Sample_name", setdiff(names(transdf), c("Species", "Sample_name")))]
#addUniqueNameColumn for when multiple datasets might be binded downstream 
transdf$UniqueTransdescoderID = paste(Sample_name, transdf$Transdecoder_ID, sep = "_")

#remove duplicates in CDS and PEP sequence, keep one with higher BitScore to ToxinProtein (should already have been made distinct)
transdf <- transdf[order(-transdf$BitScore, transdf$E_value ), ]
transdf <- distinct(transdf, Transdecoder_ID, .keep_all = TRUE)
transdf <- distinct(transdf, CDS_Sequence, .keep_all = TRUE)
transdf <- distinct(transdf, PEP_Sequence, .keep_all = TRUE)

#readinToxinCSV 
toxin_data <- read.delim(ToxinDataTSV, header = TRUE)
toxin_Domain_long <- toxin_data %>%
  separate_rows(InterPro, sep = ";") %>%    
  mutate(InterPro = str_trim(InterPro)) %>%  
  filter(InterPro != "") 

Toxin_IPs <- data.frame(unique(toxin_Domain_long$InterPro)) 
colnames(Toxin_IPs)[colnames(Toxin_IPs) == "unique.toxin_Domain_long.InterPro."] <- "IP"
pattern <- paste0("\\b(", paste(Toxin_IPs$IP, collapse = "|"), ")\\b")

if (
  mass_spec != "NULL" &&
  Blastn_result != "NULL"
)
  {
  #read in mass spec
  mass_spec <- read.csv(mass_spec, header = TRUE) 
  #rename mass spec column 
  colnames(mass_spec)[which(names(mass_spec) == "Accession")] <- "Transdecoder_ID" 
  #left_join full massspec 
  #left_join simplified mass spec 
  transdf_massspec <- left_join(transdf,mass_spec, by = "Transdecoder_ID")
  massspec_select <- mass_spec %>% dplyr::select(Top, Transdecoder_ID, X.10LgP, Coverage..., X.Peptides, X.Unique )
  transdf_massspec_sim <- left_join(transdf,massspec_select, by = "Transdecoder_ID")
  #read in Blastn
  Blastn_result_read <- read.table(Blastn_result,sep = "\t",header = FALSE,stringsAsFactors = FALSE)  
  #rename blastn columns 
  colnames(Blastn_result_read) <- c("Transdecoder_ID", "sseqid", "genome_pident", "genome_length", "genome_mismatch","genome_gapopen","genome_qstart","genome_qend","genome_sstart","genome_send","genome_evalue", "genome_bitscore", "genome_qframe", "genome_qcovs")
  #order by genome_qcovs
  Blastn_result_read_genome_qcovs <- Blastn_result_read[order(Blastn_result_read$genome_qcovs, decreasing = TRUE), ]
  #only keep distinct hit per transcript, keeping the hit  with higher qcovs
  Blastn_result_read_genome_qcovs_distinct <- Blastn_result_read_genome_qcovs %>% distinct(Transdecoder_ID, .keep_all = TRUE)
  transdf_massspec_genome_msfull <- left_join(transdf_massspec,Blastn_result_read_genome_qcovs_distinct, by = "Transdecoder_ID")
  transdf_massspec_genome_mssim <- left_join(transdf_massspec_sim,Blastn_result_read_genome_qcovs_distinct, by = "Transdecoder_ID")
  write.csv(transdf_massspec_genome_mssim, paste0(Sample_name, "_transdf_distinct_masspec_blastn_mssimplified.csv"), row.names = FALSE)
  write.csv(transdf_massspec_genome_msfull, paste0(Sample_name, "_transdf_distinct_masspec_blastn_msfull.csv"), row.names = FALSE)
  
  #Filtering 
  #completeORF + SignalP + no Transmembrane
  transdf_massspec_genome_mssim_filtered_base <- transdf_massspec_genome_mssim %>%
    filter(
      grepl("complete", ORF_type, ignore.case = TRUE),
      grepl("SP", SP_Prediction, ignore.case = TRUE),
      TMHMM == FALSE
    )
  #keep only those sequences which hit the genome 
  transdf_massspec_genome_mssim_filtered_base <-
    transdf_massspec_genome_mssim_filtered_base[
      !is.na(transdf_massspec_genome_mssim_filtered_base$sseqid),
    ]
  #changing criteria to numeric for filtering 
  transdf_massspec_genome_mssim_filtered_base[["Coverage..."]] <- as.numeric(transdf_massspec_genome_mssim_filtered_base[["Coverage..."]])
  transdf_massspec_genome_mssim_filtered_base[["BitScore"]] <- as.numeric(transdf_massspec_genome_mssim_filtered_base[["BitScore"]])
  transdf_massspec_genome_mssim_filtered_base[["percent"]] <- as.numeric(transdf_massspec_genome_mssim_filtered_base[["percent"]])
  ##VennDiagram & Filtering (Strict)
  matching_rows <- transdf_massspec_genome_mssim_filtered_base[
    !is.na(transdf_massspec_genome_mssim_filtered_base$InterPro_accession_Names) &
      grepl(pattern, transdf_massspec_genome_mssim_filtered_base$InterPro_accession_Names),
  ]
  set_A<-matching_rows
  set_B <- transdf_massspec_genome_mssim_filtered_base[transdf_massspec_genome_mssim_filtered_base[["Coverage..."]] > 50 & !is.na(transdf_massspec_genome_mssim_filtered_base[["Coverage..."]]), ]
  set_C <- transdf_massspec_genome_mssim_filtered_base[transdf_massspec_genome_mssim_filtered_base[["BitScore"]] > 250 & !is.na(transdf_massspec_genome_mssim_filtered_base[["BitScore"]]), ]
  set_D <-transdf_massspec_genome_mssim_filtered_base[transdf_massspec_genome_mssim_filtered_base[["percent"]] > 1 & !is.na(transdf_massspec_genome_mssim_filtered_base[["percent"]]), ]
  venn_list <- list(
    TD = set_A$Transdecoder_ID,
    MS = set_B$Transdecoder_ID,
    TP = set_C$Transdecoder_ID,
    KE = set_D$Transdecoder_ID
    
  )
  p <- ggvenn(venn_list,c("TD", "MS","TP","KE"), fill_color = c("#E41A1C", "#377EB8", "#4DAF4A", "#EBAC4D"))
  Union_ABC <- Reduce(union, list(set_A$Transdecoder_ID, set_B$Transdecoder_ID, set_C$Transdecoder_ID, set_D$Transdecoder_ID))
  ggsave(paste0(Sample_name,"_Venn_strict.png"), plot = p, width = 6, height = 4, dpi = 300)
  transdf_massspec_genome_mssim_filtered_base_union <- transdf_massspec_genome_mssim_filtered_base[transdf_massspec_genome_mssim_filtered_base$Transdecoder_ID %in% Union_ABC, ]
  write.csv(transdf_massspec_genome_mssim_filtered_base_union, paste0(Sample_name, "_Venn_Diagram_union_strict.csv"), row.names = FALSE)

  ##VennDiagram & Filtering (Lax)
  set_B <- transdf_massspec_genome_mssim_filtered_base[transdf_massspec_genome_mssim_filtered_base[["Coverage..."]] > 0 & !is.na(transdf_massspec_genome_mssim_filtered_base[["Coverage..."]]), ]
  set_C <- transdf_massspec_genome_mssim_filtered_base[transdf_massspec_genome_mssim_filtered_base[["BitScore"]] > 0 & !is.na(transdf_massspec_genome_mssim_filtered_base[["BitScore"]]), ]
  set_D <-transdf_massspec_genome_mssim_filtered_base[transdf_massspec_genome_mssim_filtered_base[["percent"]] > 0 & !is.na(transdf_massspec_genome_mssim_filtered_base[["percent"]]), ]
  venn_list <- list(
    TD = set_A$Transdecoder_ID,
    MS = set_B$Transdecoder_ID,
    TP = set_C$Transdecoder_ID,
    KE = set_D$Transdecoder_ID
    
  )
  p <- ggvenn(venn_list,c("TD", "MS","TP","KE"), fill_color = c("#E41A1C", "#377EB8", "#4DAF4A", "#EBAC4D"))
  Union_ABC <- Reduce(union, list(set_A$Transdecoder_ID, set_B$Transdecoder_ID, set_C$Transdecoder_ID, set_D$Transdecoder_ID))
  ggsave(paste0(Sample_name,"_Venn_lax.png"), plot = p, width = 6, height = 4, dpi = 300)
  transdf_massspec_genome_mssim_filtered_base_union <- transdf_massspec_genome_mssim_filtered_base[transdf_massspec_genome_mssim_filtered_base$Transdecoder_ID %in% Union_ABC, ]
  write.csv(transdf_massspec_genome_mssim_filtered_base_union, paste0(Sample_name, "_Venn_Diagram_union_lax.csv"), row.names = FALSE)
  union_lax <- transdf_massspec_genome_mssim_filtered_base_union
}

else if (  mass_spec == "NULL" &&
       Blastn_result == "NULL")
{
  #just save the modified transdf 
  write.csv(transdf, paste0(Sample_name, "_transdf_distinct_nomasspec_noblastn.csv"), row.names = FALSE)
  #Filtering 
  
  #completeORF + SignalP + no Transmembrane
  transdf_filtered_base <- transdf %>%
    filter(
      grepl("complete", ORF_type, ignore.case = TRUE),
      grepl("SP", SP_Prediction, ignore.case = TRUE),
      TMHMM == FALSE
    )
  #changing criteria to numeric for filtering 
  transdf_filtered_base[["BitScore"]] <- as.numeric(transdf_filtered_base[["BitScore"]])
  transdf_filtered_base[["percent"]] <- as.numeric(transdf_filtered_base[["percent"]])
  
  ##VennDiagram & Filtering (Strict)
  matching_rows <- transdf_filtered_base[
    !is.na(transdf_filtered_base$InterPro_accession_Names) &
      grepl(pattern, transdf_filtered_base$InterPro_accession_Names),
  ]
  set_A<-matching_rows
  set_C <- transdf_filtered_base[transdf_filtered_base[["BitScore"]] > 250 & !is.na(transdf_filtered_base[["BitScore"]]), ]
  set_D <-transdf_filtered_base[transdf_filtered_base[["percent"]] > 1 & !is.na(transdf_filtered_base[["percent"]]), ]
  venn_list <- list(
    TD = set_A$Transdecoder_ID,
    TP = set_C$Transdecoder_ID,
    KE = set_D$Transdecoder_ID
    
  )
  p <- ggvenn(venn_list,c("TD","TP","KE"), fill_color = c("#E41A1C", "#4DAF4A", "#EBAC4D"))
  Union_ABC <- Reduce(union, list(set_A$Transdecoder_ID, set_C$Transdecoder_ID, set_D$Transdecoder_ID))
  ggsave(paste0(Sample_name,"_Venn_strict.png"), plot = p, width = 6, height = 4, dpi = 300)
  transdf_filtered_base_union <- transdf_filtered_base[transdf_filtered_base$Transdecoder_ID %in% Union_ABC, ]
  write.csv(transdf_filtered_base_union, paste0(Sample_name, "_Venn_Diagram_union_strict.csv"), row.names = FALSE)
  
  ##VennDiagram & Filtering (Lax)
  set_C <- transdf_filtered_base[transdf_filtered_base[["BitScore"]] > 0 & !is.na(transdf_filtered_base[["BitScore"]]), ]
  set_D <-transdf_filtered_base[transdf_filtered_base[["percent"]] > 0 & !is.na(transdf_filtered_base[["percent"]]), ]
  venn_list <- list(
    TD = set_A$Transdecoder_ID,
    TP = set_C$Transdecoder_ID,
    KE = set_D$Transdecoder_ID
    
  )
  p <- ggvenn(venn_list,c("TD","TP","KE"), fill_color = c("#E41A1C", "#4DAF4A", "#EBAC4D"))
  Union_ABC <- Reduce(union, list(set_A$Transdecoder_ID, set_C$Transdecoder_ID, set_D$Transdecoder_ID))
  ggsave(paste0(Sample_name,"_Venn_lax.png"), plot = p, width = 6, height = 4, dpi = 300)
  transdf_filtered_base_union <- transdf_filtered_base[transdf_filtered_base$Transdecoder_ID %in% Union_ABC, ]
  write.csv(transdf_filtered_base_union, paste0(Sample_name, "_Venn_Diagram_union_lax.csv"), row.names = FALSE)
  union_lax <- transdf_filtered_base_union
  

}

else if ( mass_spec == "NULL" &&
        Blastn_result != "NULL"
) {
 #read in Blastn
  Blastn_result_read <- read.table(Blastn_result,sep = "\t",header = FALSE,stringsAsFactors = FALSE)  
  #rename blastn columns 
  colnames(Blastn_result_read) <- c("Transdecoder_ID", "sseqid", "genome_pident", "genome_length", "genome_mismatch","genome_gapopen","genome_qstart","genome_qend","genome_sstart","genome_send","genome_evalue", "genome_bitscore", "genome_qframe", "genome_qcovs")
  #order by genome_qcovs
  Blastn_result_read_genome_qcovs <- Blastn_result_read[order(Blastn_result_read$genome_qcovs, decreasing = TRUE), ]
  #only keep distinct hit per transcript, keeping the hit  with higher qcovs
  Blastn_result_read_genome_qcovs_distinct <- Blastn_result_read_genome_qcovs %>% distinct(Transdecoder_ID, .keep_all = TRUE)
  transdf_nomassspec_genome <- left_join(transdf,Blastn_result_read_genome_qcovs_distinct, by = "Transdecoder_ID")
  write.csv(transdf_nomassspec_genome, paste0(Sample_name, "_transdf_distinct_nomasspec_blastn.csv"), row.names = FALSE)
  

  #Filtering 
  
  #completeORF + SignalP + no Transmembrane
  transdf_nomassspec_genome_filtered_base <- transdf_nomassspec_genome %>%
    filter(
      grepl("complete", ORF_type, ignore.case = TRUE),
      grepl("SP", SP_Prediction, ignore.case = TRUE),
      TMHMM == FALSE
    )
  #keep only those sequences which hit the genome 
  transdf_nomassspec_genome_filtered_base <-
    transdf_nomassspec_genome_filtered_base[
      !is.na(transdf_nomassspec_genome_filtered_base$sseqid),
    ]
  #changing criteria to numeric for filtering 
  transdf_nomassspec_genome_filtered_base[["BitScore"]] <- as.numeric(transdf_nomassspec_genome_filtered_base[["BitScore"]])
  transdf_nomassspec_genome_filtered_base[["percent"]] <- as.numeric(transdf_nomassspec_genome_filtered_base[["percent"]])
  
  ##VennDiagram & Filtering (Strict)
  matching_rows <- transdf_nomassspec_genome_filtered_base[
    !is.na(transdf_nomassspec_genome_filtered_base$InterPro_accession_Names) &
      grepl(pattern, transdf_nomassspec_genome_filtered_base$InterPro_accession_Names),
  ]
  set_A<-matching_rows
  set_C <- transdf_nomassspec_genome_filtered_base[transdf_nomassspec_genome_filtered_base[["BitScore"]] > 250 & !is.na(transdf_nomassspec_genome_filtered_base[["BitScore"]]), ]
  set_D <-transdf_nomassspec_genome_filtered_base[transdf_nomassspec_genome_filtered_base[["percent"]] > 1 & !is.na(transdf_nomassspec_genome_filtered_base[["percent"]]), ]
  venn_list <- list(
    TD = set_A$Transdecoder_ID,
    TP = set_C$Transdecoder_ID,
    KE = set_D$Transdecoder_ID
    
  )
  p <- ggvenn(venn_list,c("TD","TP","KE"), fill_color = c("#E41A1C", "#4DAF4A", "#EBAC4D"))
  Union_ABC <- Reduce(union, list(set_A$Transdecoder_ID, set_C$Transdecoder_ID, set_D$Transdecoder_ID))
  ggsave(paste0(Sample_name,"_Venn_strict.png"), plot = p, width = 6, height = 4, dpi = 300)
  transdf_nomassspec_genome_filtered_base_union <- transdf_nomassspec_genome_filtered_base[transdf_nomassspec_genome_filtered_base$Transdecoder_ID %in% Union_ABC, ]
  write.csv(transdf_nomassspec_genome_filtered_base_union, paste0(Sample_name, "_Venn_Diagram_union_strict.csv"), row.names = FALSE)
  
  ##VennDiagram & Filtering (Lax)
  set_C <- transdf_nomassspec_genome_filtered_base[transdf_nomassspec_genome_filtered_base[["BitScore"]] > 0 & !is.na(transdf_nomassspec_genome_filtered_base[["BitScore"]]), ]
  set_D <-transdf_nomassspec_genome_filtered_base[transdf_nomassspec_genome_filtered_base[["percent"]] > 0 & !is.na(transdf_nomassspec_genome_filtered_base[["percent"]]), ]
  venn_list <- list(
    TD = set_A$Transdecoder_ID,
    TP = set_C$Transdecoder_ID,
    KE = set_D$Transdecoder_ID
    
  )
  p <- ggvenn(venn_list,c("TD","TP","KE"), fill_color = c("#E41A1C", "#4DAF4A", "#EBAC4D"))
  Union_ABC <- Reduce(union, list(set_A$Transdecoder_ID, set_C$Transdecoder_ID, set_D$Transdecoder_ID))
  ggsave(paste0(Sample_name,"_Venn_lax.png"), plot = p, width = 6, height = 4, dpi = 300)
  transdf_nomassspec_genome_filtered_base_union <- transdf_nomassspec_genome_filtered_base[transdf_nomassspec_genome_filtered_base$Transdecoder_ID %in% Union_ABC, ]
  write.csv(transdf_nomassspec_genome_filtered_base_union, paste0(Sample_name, "_Venn_Diagram_union_lax.csv"), row.names = FALSE)
  union_lax <- transdf_nomassspec_genome_filtered_base_union
  
  
}

else if (
  mass_spec != "NULL" &&
  Blastn_result == "NULL"
) {
  #read in mass spec
  mass_spec <- read.csv(mass_spec, header = TRUE) 
  #rename mass spec column 
  colnames(mass_spec)[which(names(mass_spec) == "Accession")] <- "Transdecoder_ID" 
  #left_join full massspec 
  #left_join simplified mass spec 
  transdf_massspec <- left_join(transdf,mass_spec, by = "Transdecoder_ID")
  massspec_select <- mass_spec %>% dplyr::select(Top, Transdecoder_ID, X.10LgP, Coverage..., X.Peptides, X.Unique )
  transdf_massspec_sim <- left_join(transdf,massspec_select, by = "Transdecoder_ID")
  write.csv(transdf_massspec_sim, paste0(Sample_name, "_transdf_distinct_masspec_noblastn_mssimplified.csv"), row.names = FALSE)
  write.csv(transdf_massspec, paste0(Sample_name, "_transdf_distinct_masspec_noblastn_msfull.csv"), row.names = FALSE)
  
  #Filtering 
  #completeORF + SignalP + no Transmembrane
  transdf_massspec_filtered_base <- transdf_massspec %>%
    filter(
      grepl("complete", ORF_type, ignore.case = TRUE),
      grepl("SP", SP_Prediction, ignore.case = TRUE),
      TMHMM == FALSE
    )

  #changing criteria to numeric for filtering 
  transdf_massspec_filtered_base[["Coverage..."]] <- as.numeric(transdf_massspec_filtered_base[["Coverage..."]])
  transdf_massspec_filtered_base[["BitScore"]] <- as.numeric(transdf_massspec_filtered_base[["BitScore"]])
  transdf_massspec_filtered_base[["percent"]] <- as.numeric(transdf_massspec_filtered_base[["percent"]])
  ##VennDiagram & Filtering (Strict)
  matching_rows <- transdf_massspec_filtered_base[
    !is.na(transdf_massspec_filtered_base$InterPro_accession_Names) &
      grepl(pattern, transdf_massspec_filtered_base$InterPro_accession_Names),
  ]
  set_A<-matching_rows
  set_B <- transdf_massspec_filtered_base[transdf_massspec_filtered_base[["Coverage..."]] > 50 & !is.na(transdf_massspec_filtered_base[["Coverage..."]]), ]
  set_C <- transdf_massspec_filtered_base[transdf_massspec_filtered_base[["BitScore"]] > 250 & !is.na(transdf_massspec_filtered_base[["BitScore"]]), ]
  set_D <-transdf_massspec_filtered_base[transdf_massspec_filtered_base[["percent"]] > 1 & !is.na(transdf_massspec_filtered_base[["percent"]]), ]
  venn_list <- list(
    TD = set_A$Transdecoder_ID,
    MS = set_B$Transdecoder_ID,
    TP = set_C$Transdecoder_ID,
    KE = set_D$Transdecoder_ID
    
  )
  p <- ggvenn(venn_list,c("TD", "MS","TP","KE"), fill_color = c("#E41A1C", "#377EB8", "#4DAF4A", "#EBAC4D"))
  Union_ABC <- Reduce(union, list(set_A$Transdecoder_ID, set_B$Transdecoder_ID, set_C$Transdecoder_ID, set_D$Transdecoder_ID))
  ggsave(paste0(Sample_name,"_Venn_strict.png"), plot = p, width = 6, height = 4, dpi = 300)
  transdf_massspec_filtered_base_union <- transdf_massspec_filtered_base[transdf_massspec_filtered_base$Transdecoder_ID %in% Union_ABC, ]
  write.csv(transdf_massspec_filtered_base_union, paste0(Sample_name, "_Venn_Diagram_union_strict.csv"), row.names = FALSE)
  
  ##VennDiagram & Filtering (Lax)
  set_B <- transdf_massspec_filtered_base[transdf_massspec_filtered_base[["Coverage..."]] > 0 & !is.na(transdf_massspec_filtered_base[["Coverage..."]]), ]
  set_C <- transdf_massspec_filtered_base[transdf_massspec_filtered_base[["BitScore"]] > 0 & !is.na(transdf_massspec_filtered_base[["BitScore"]]), ]
  set_D <-transdf_massspec_filtered_base[transdf_massspec_filtered_base[["percent"]] > 0 & !is.na(transdf_massspec_filtered_base[["percent"]]), ]
  venn_list <- list(
    TD = set_A$Transdecoder_ID,
    MS = set_B$Transdecoder_ID,
    TP = set_C$Transdecoder_ID,
    KE = set_D$Transdecoder_ID
    
  )
  p <- ggvenn(venn_list,c("TD", "MS","TP","KE"), fill_color = c("#E41A1C", "#377EB8", "#4DAF4A", "#EBAC4D"))
  Union_ABC <- Reduce(union, list(set_A$Transdecoder_ID, set_B$Transdecoder_ID, set_C$Transdecoder_ID, set_D$Transdecoder_ID))
  ggsave(paste0(Sample_name,"_Venn_lax.png"), plot = p, width = 6, height = 4, dpi = 300)
  transdf_massspec_filtered_base_union <- transdf_massspec_filtered_base[transdf_massspec_filtered_base$Transdecoder_ID %in% Union_ABC, ]
  write.csv(transdf_massspec_filtered_base_union, paste0(Sample_name, "_Venn_Diagram_union_lax.csv"), row.names = FALSE)
  union_lax <- transdf_massspec_filtered_base_union
  
}




