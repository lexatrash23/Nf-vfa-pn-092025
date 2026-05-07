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

#for those with genomes, transcripts that overlap the same genomic locus are clustered only best blast match is kept
#!/usr/bin/env Rscript

library(dplyr)
library(tidyr)
library(stringr)
library(ggvenn)
library(ggplot2)
library(GenomicRanges)
library(igraph)
library(Biostrings)
#load in command line arguments
args <- commandArgs(trailingOnly = TRUE)
Sample_name <- args[1] 
species_name  <- args[2]
Transdf_distinct_file <- args[3]
ToxnontoxIP <- args[4]
Blastn_result <- args[5]
mass_spec_file <- args[6]
setwd("/Users/praveena/Desktop/PhD_all/2025/github/oldergithubscripts/Retrial/is15/")
Sample_name <- "DP3"
species_name <- "DP"
Transdf_distinct_file <- "DP3_PD_transdf_distinct.csv"
ToxnontoxIP <- "ToxNonTox_IP.csv"
Blastn_result <- "DP3_PD.blastn.db.6.txt"
mass_spec_file <-  "NULL"
 

# Function to read in transdf and make transdf_distinct, filter for complete ORF with signal P //slightly different from other transdf distincct
read_transdf <- function(file) {
  df <- read.csv(file, header = TRUE)
  df$Species <- species_name
  df$Sample_name <- Sample_name
  df <- df[c("Species", "Sample_name", setdiff(names(df), c("Species", "Sample_name")))]
  df$UniqueSequenceName <- paste(Sample_name, df$Transdecoder_ID, sep = "_")
  
  df <- df %>% 
    arrange(desc(BitScore), E_value, desc(tpm)) %>%
    distinct(Transdecoder_ID, .keep_all = TRUE) %>%
    group_by(CDS_Sequence) %>%
    mutate(percent = sum(percent, na.rm = TRUE),
           tpm = sum(tpm, na.rm = TRUE)) %>%
    ungroup() %>%
    distinct(Transdecoder_ID, CDS_Sequence, .keep_all = TRUE)
  df <- df %>% filter(grepl("complete", ORF_type, ignore.case = TRUE),
                      grepl("SP", SP, ignore.case = TRUE),
                      TMHMM %in% c(FALSE, "FALSE"))
  for(col in c("BitScore","percent","CysPer","Coverage")) if(col %in% names(df)) df[[col]] <- as.numeric(df[[col]])

  return(df)
}

# Read in summary data for pattern generation
prepare_patterns <- function(file) {
  summary_IP <- read.csv(file)
  IPInToxins <- summary_IP %>% filter(RelativeExpression > 0)
  Toxin_IPs <- unique(IPInToxins$InterPro)
  pattern <- if(length(Toxin_IPs) == 0) "$^" else paste0("\\b(", paste(Toxin_IPs, collapse="|"), ")\\b")
  
  IPOver <- summary_IP %>% filter(RelativeExpression > 1)
  Over_IPs <- unique(IPOver$InterPro)
  pattern2 <- if(length(Over_IPs) == 0) "$^" else paste0("\\b(", paste(Over_IPs, collapse="|"), ")\\b")
  
  return(list(pattern = pattern, pattern2 = pattern2))
}

# Read mass spec if available 
read_mass_spec <- function(file) {
  mass <- read.csv(file)
  colnames(mass)[colnames(mass)=="Accession"] <- "Transdecoder_ID"
  return(mass)
}

# Read in blastn data if available
read_blast <- function(file) {
  blast <- read.table(file, sep="\t", header=FALSE, stringsAsFactors = FALSE)
  colnames(blast) <- c("Transdecoder_ID", "genome_sseqid", "genome_pident", "genome_length",
                       "genome_mismatch","genome_gapopen","genome_qstart","genome_qend",
                       "genome_sstart","genome_send", "genome_sstrand","genome_evalue", 
                       "genome_bitscore", "genome_qframe", "genome_qcovs")
  blast <- blast %>% arrange(desc(genome_bitscore)) %>% distinct(Transdecoder_ID, .keep_all = TRUE)
  return(blast)
}

# Filter and prepare Base
prepare_base <- function(df) {
  df <- df %>% filter(grepl("complete", ORF_type, ignore.case = TRUE),
                      grepl("SP", SP, ignore.case = TRUE),
                      TMHMM == FALSE)
  for(col in c("BitScore","percent","CysPer","Coverage")) if(col %in% names(df)) df[[col]] <- as.numeric(df[[col]])
  return(df)
}

# Filtering function
filter_and_venn <- function(Base, pattern, pattern2, mass = FALSE, strict = TRUE, Sample_name) {
  # patterns depending on strict vs lax
  pat <- if(strict) pattern2 else pattern
  if(strict){
    Base <- Base %>%
      filter(genome_qcovs >= 95)
  }
  matching_rows <- Base[!is.na(Base$InterPro_accession_Names) & grepl(pat, Base$InterPro_accession_Names), ]
  
  # thresholds
  if(strict) {
    sets <- list(
      TD = matching_rows$Transdecoder_ID,
      MS = if(mass) Base[Base$Coverage >= 50 & Base$Top == TRUE,]$Transdecoder_ID else NULL,
      TP = Base[Base$BitScore >= 250,]$Transdecoder_ID,
      KE = Base[Base$percent >= 1,]$Transdecoder_ID,
      CP = Base[Base$CysPer >= 5,]$Transdecoder_ID
    )
    file_suffix <- "strict"
  } else {
    sets <- list(
      TD = matching_rows$Transdecoder_ID,
      MS = if(mass) Base[Base$Coverage >= 0 & Base$Top == TRUE,]$Transdecoder_ID else NULL,
      TP = Base[Base$BitScore >= 50,]$Transdecoder_ID,
      KE = Base[Base$percent >= 0,]$Transdecoder_ID,
      CP = Base[Base$CysPer >= 1,]$Transdecoder_ID
    )
    file_suffix <- "lax"
  }
  
  #set colours
  color_map <- c(
    TD = "#E41A1C",
    MS = "#377EB8",
    TP = "#4DAF4A",
    KE = "#EBAC4D",
    CP = "#782DC8"
  )
  # remove NULL to remove MS if not available. only keep relevant colours
  sets <- sets[!sapply(sets, is.null)]
  fill_colors <- color_map[names(sets)]
  
  p <- ggvenn(sets, names(sets), fill_color = fill_colors)
  ggsave(paste0(Sample_name,"_Venn_",file_suffix,".png"), plot = p, width = 6, height = 4, dpi = 300)
  
  union_ids <- Reduce(union, sets)
  Base_union <- Base[Base$Transdecoder_ID %in% union_ids, ] %>% mutate(Filter = if(strict) "Strict" else "Lax")
  return(Base_union)
  
  Pepseq <- AAStringSet(Base_union$PEP_Sequence)
  names(Pepseq) <- Base_union$Transdecoder_ID
  writeAAStringSet(Pepseq, paste0(Sample_name,"_filtered_",file_suffix,".pep"))

}

## Loading files
transdf <- read_transdf(Transdf_distinct_file)


patterns <- prepare_patterns(ToxnontoxIP)

# handle mass spec and blast conditions
mass_spec <- if (!is.na(mass_spec_file) && mass_spec_file != "NULL") {
  read_mass_spec(mass_spec_file)
} else {
  NULL
}

blastn <- if (!is.na(Blastn_result) && Blastn_result != "NULL") {
  read_blast(Blastn_result)
} else {
  NULL
}
# Join mass spec and blast unfiltered
transdf_unfiltered <- transdf
if (!is.na(mass_spec_file) && mass_spec_file != "NULL") {
  transdf_unfiltered <- left_join(transdf_unfiltered, mass_spec, by="Transdecoder_ID")
}
if (!is.na(Blastn_result) && Blastn_result != "NULL") {
  transdf_unfiltered <- left_join(transdf_unfiltered, blastn, by="Transdecoder_ID")
}
write.csv(transdf_unfiltered, paste0(Sample_name,"_transdf_distinct_final_unfiltered.csv"), row.names = FALSE)

#Join mass spec and blast filtered
transdf_filtered <- transdf 
if (!is.na(mass_spec_file) && mass_spec_file != "NULL")
  {transdf_filtered <- left_join(transdf_filtered, mass_spec %>% select(Transdecoder_ID, Top, Coverage, Unique), by = "Transdecoder_ID")
}

if (!is.na(Blastn_result) && Blastn_result != "NULL") {
  blastn_file <- blastn %>%
    dplyr::select(Transdecoder_ID, genome_sseqid, genome_pident, genome_length,
                                   genome_mismatch,genome_sstart,genome_send, genome_sstrand,genome_evalue, 
                                   genome_bitscore, genome_qframe, genome_qcovs)
  
  transdf_filtered <- left_join(transdf_filtered, blastn_file , by = "Transdecoder_ID") %>%
    filter(!is.na(genome_sseqid)) %>%
    filter(genome_sseqid != "")
  
  transdf_filtered <- transdf_filtered %>%
    mutate(genome_sstrand = case_when(
      genome_sstrand == "plus" ~ "+",
      genome_sstrand == "minus" ~ "-"
    )) 
  
  #create genome range object
  gr <- GRanges(
    seqnames = trimws(as.character(transdf_filtered$genome_sseqid)),
    ranges = IRanges(
      start = pmin(transdf_filtered$genome_sstart, transdf_filtered$genome_send),
      end   = pmax(transdf_filtered$genome_sstart, transdf_filtered$genome_send)
    ),
    query_id = transdf_filtered$Transdecoder_ID,
    strand = transdf_filtered$genome_sstrand
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
  transdf_filtered$cluster <- clusters[seq_along(transdf_filtered$Transdecoder_ID)]
  
  transdf_filtered <- transdf_filtered%>%
    group_by(cluster) %>%
    mutate(tpm_aggregates = sum(tpm)) %>%
    relocate(tpm_aggregates, .after = CysPer)  %>%
    ungroup()
  
  transdf_filtered <- transdf_filtered %>% arrange (
    desc(genome_qcovs),
    desc(genome_pident),
    desc(genome_bitscore),
    desc(PEP_Length)
  ) %>%
    distinct(cluster, .keep_all = TRUE) %>%
    dplyr::select(-cluster)
  
  transdf_filtered <- transdf_filtered %>% 
    filter(genome_qcovs >= 50)

}

Base <- transdf_filtered
# strict filter
Base_strict <- filter_and_venn(Base, patterns$pattern, patterns$pattern2, mass = !is.null(mass_spec), strict = TRUE, Sample_name = Sample_name)
write.csv(Base_strict, paste0(Sample_name,"_transdf_distinct_final_filtered_strict.csv"), row.names = FALSE)

# lax filter
Base_lax <- filter_and_venn(Base, patterns$pattern, patterns$pattern2, mass = !is.null(mass_spec), strict = FALSE, Sample_name = Sample_name)
write.csv(Base_lax, paste0(Sample_name,"_transdf_distinct_final_filtered_lax.csv"), row.names = FALSE)