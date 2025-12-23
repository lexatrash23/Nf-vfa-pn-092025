#!/usr/bin/env Rscript


#installing and loading packages
library(Biostrings)
library(dplyr)
library(tidyr)
library(htmlwidgets)
library(plotly)


#command line arguments
args <- commandArgs(trailingOnly = TRUE)
Trans <- args[1]
Toxin_domains_file <- args[2]
fasta_file_pep <- args[3]
fasta_file_cds <- args[4]
sample <- args[4]

#Generating Interproscan plotly
#read in transdf_distinct
Trans <- read.csv(file = Trans, header = TRUE)

# Reordering columns and keeping only those of interest
keeps <- c("Transdecoder_ID", "ORF_type", "PEP_Length", "CDS_Length", "SP_Prediction", "Signal_Length", "mature_length", "percent", "cumulativepercent", "Code", "Hit", "Percentage_Identity", "E_value", "BitScore", "Hit_species","InterPro_accession_Names","GO_name","TMHMM","Signal_Sequence","mature_sequence", "PEP_Sequence","CDS_Sequence")
Trans <- Trans[keeps]


# sort by bitscore
Trans_sorted <- Trans[order(Trans$BitScore, decreasing = TRUE),]
#sort only for those complete with false tmhmm
FINAL_CSV_distinct_filtered <- Trans_sorted %>%
  filter(
    grepl("complete", ORF_type, ignore.case = TRUE),
    grepl("SP", SP_Prediction, ignore.case = TRUE),
    TMHMM == FALSE
  )
#keep only columns of interest
keeps <- c("Transdecoder_ID","PEP_Length", "InterPro_accession_Names","GO_name")
df_figures <-FINAL_CSV_distinct_filtered[keeps]



#separate the interproscan IPs by PIPs
annotation_list <- strsplit(df_figures$InterPro_accession_Names, split = "\\|")
#creates one long list of IPs and summarizes the counts for each
all_annotations <- unlist(annotation_list)
annotation_counts <- table(all_annotations)
annotation_counts_df <- as.data.frame(annotation_counts)
#sort by frequencing of IP
annotation_counts_df <- annotation_counts_df[order(-annotation_counts_df$Freq), ]

colnames(annotation_counts_df) <- c("Transdf_annotations", "Freq")

write.csv(annotation_counts_df, paste0(sample,"_annotation_counts_trandsf.csv"), row.names = FALSE)


#read in toxin domains tsv
Toxin_domains <- read.delim(file = Toxin_domains_file, header = TRUE, sep = "\t")
#same as above but different separator for this tsv
interpro_list <- strsplit(Toxin_domains$InterPro, split = ";")
all_interpro_ids <- unlist(interpro_list)
all_interpro_ids <- trimws(all_interpro_ids) #remove ws
unique_interpro_ids <- unique(all_interpro_ids) #only get unique IDs
unique_interpro_df <- data.frame(InterPro = unique_interpro_ids)
annotation_counts_toxin <- table(all_interpro_ids)
annotation_counts_toxin_df <- as.data.frame(annotation_counts_toxin)
annotation_counts_toxin_df <- annotation_counts_toxin[order(-annotation_counts_toxin$Freq), ]
colnames(annotation_counts_toxin_df) <- c("toxin_annotations", "Freq")

write.csv(annotation_counts_toxin_df, paste0(sample,"_annotation_counts_toxins.csv"), row.names = FALSE)

#iterate through the unique_interpro_df$InterPro, search df_figres_interpro_names df with this value as query, save the rows that have a match at least somewhere in  the row
matched_rows <- data.frame()
for (interpro_id in unique_interpro_df$InterPro) {
  # searches
  matches <- df_figures[apply(df_figures, 1, function(row) any(grepl(interpro_id, row, fixed = TRUE))), ]
  
  # add the whole df_Figures row row to a new dataframe
  if (nrow(matches) > 0) {
    matched_rows <- rbind(matched_rows, matches)
  }
}
#unique Transdecoder IDs only
matched_rows <- distinct(matched_rows, Transdecoder_ID, .keep_all = TRUE)

#remove any na rows
matched_rows <- matched_rows[!is.na(matched_rows$InterPro_accession_Names), ]

#Frequency list of toxin-related domains
annotation_list_matched <- strsplit(matched_rows$InterPro_accession_Names, split = "\\|")
all_annotations_matched  <- unlist(annotation_list_matched)
annotation_counts_matched <- table(all_annotations_matched)
annotation_counts_df_matched <- as.data.frame(annotation_counts_matched)
annotation_counts_df_matched <- annotation_counts_df_matched[order(-annotation_counts_df_matched$Freq), ]

colnames(annotation_counts_df_matched) <- c("matched_toxin_annotations", "Freq")

write.csv(annotation_counts_df_matched, paste0(sample,"_annotation_counts_trandf_matched_to_toxins.csv"), row.names = FALSE)




fasta <- readAAStringSet(file = fasta_file_pep)
seq_ids <- sapply(strsplit(names(fasta), " "), `[`, 1)
keep <- seq_ids %in% matched_rows$Transdecoder_ID
filtered_fasta <- fasta[keep]
writeXStringSet(filtered_fasta, file.path("filtered_sequences.fasta"))


#GROUPING TYPE 
annotation_counts_df_matched$Group <- ifelse(
  grepl("immunoglobulin", tolower(annotation_counts_df_matched$all_annotations)), 
  'Immunoglobulin-related', 
  ifelse(
    grepl("cysteine", tolower(annotation_counts_df_matched$all_annotations)),
    'Cysteine-related',
    ifelse(
      grepl("peptidase", tolower(annotation_counts_df_matched$all_annotations)),
      'Peptidase',
      ifelse(
        grepl("serine protease", tolower(annotation_counts_df_matched$all_annotations)),
        'Serine protease',
        ifelse(
          grepl("lipase", tolower(annotation_counts_df_matched$all_annotations)),
          'Lipases',
          ifelse(
            grepl("chitinase", tolower(annotation_counts_df_matched$all_annotations)),
            'Chitinases',
            ifelse(
              grepl("inhibitor", tolower(annotation_counts_df_matched$all_annotations)),
              'Inhibitors',
              ifelse(
                grepl("snake toxin-like", tolower(annotation_counts_df_matched$all_annotations)),
                'Snake toxin-like superfamily',
                ifelse(
                  grepl("kinase", tolower(annotation_counts_df_matched$all_annotations)),
                  'Kinases',
                  ifelse(
                    grepl("hydrolase", tolower(annotation_counts_df_matched$all_annotations)),
                    'Hydrolase',
                    ifelse(
                      grepl("integrin", tolower(annotation_counts_df_matched$all_annotations)),
                      'Adhesion/Integrin',
                      ifelse(
                        grepl("fibrinogen|collagen|laminin", tolower(annotation_counts_df_matched$all_annotations)),
                        'Coagulation/ECM',
                        ifelse(
                          grepl("lectin|galectin", tolower(annotation_counts_df_matched$all_annotations)),
                          'Lectin-related',
                          ifelse(
                            grepl("kunitz", tolower(annotation_counts_df_matched$all_annotations)),
                            'Protease Inhibitors',
                            ifelse(
                              grepl("egf", tolower(annotation_counts_df_matched$all_annotations)),
                              'Growth Factor Domains',
                              ifelse(
                                grepl("notch", tolower(annotation_counts_df_matched$all_annotations)),
                                'Notch Signaling',
                                ifelse(
                                  grepl("interleukin|tir", tolower(annotation_counts_df_matched$all_annotations)),
                                  'Immune Signaling',
                                  ifelse(
                                    grepl("vegf|pdgf", tolower(annotation_counts_df_matched$all_annotations)),
                                    'Growth Factor Pathways',
                                    ifelse(
                                      grepl("kringle", tolower(annotation_counts_df_matched$all_annotations)),
                                      'Coagulation',
                                      ifelse(
                                        grepl("gpcr", tolower(annotation_counts_df_matched$all_annotations)),
                                        'GPCR/Signaling',
                                        ifelse(
                                          grepl("peroxiredoxin", tolower(annotation_counts_df_matched$all_annotations)),
                                          'Antioxidant Enzymes',
                                          ifelse(
                                            grepl("phosphatase", tolower(annotation_counts_df_matched$all_annotations)),
                                            'Phosphatase',
                                            ifelse(
                                              grepl("cap", tolower(annotation_counts_df_matched$all_annotations)),
                                              'CAP',
                                              ifelse(
                                                grepl("fibronectin", tolower(annotation_counts_df_matched$all_annotations)),
                                                'Fibronectin-related',
                                                ifelse(
                                                  grepl("ef", tolower(annotation_counts_df_matched$all_annotations)),
                                                  'EF-related',
                                                  ifelse(
                                                    grepl("von willebrand factor", tolower(annotation_counts_df_matched$all_annotations)),
                                                    'von Willebrand factor',
                                                    ifelse(
                                                      grepl("extracellular matrix assembly and organization", tolower(annotation_counts_df_matched$all_annotations)),
                                                      'Extracellular Matrix Assembly and Organization',
                                                      ifelse(
                                                        grepl("thioredoxin", tolower(annotation_counts_df_matched$all_annotations)),
                                                        'Thioredoxin',
                                                        ifelse(
                                                          grepl("pan/apple", tolower(annotation_counts_df_matched$all_annotations)),
                                                          'Pan/Apple Domain',
                                                          ifelse(
                                                            grepl("chitin binding domain", tolower(annotation_counts_df_matched$all_annotations)),
                                                            'Chitin Binding Domain',
                                          ifelse(
                                            grepl("complement", tolower(annotation_counts_df_matched$all_annotations)),
                                            'Complement System',
                                            'Other'
                                          )
                                          )
                                          )
                                                      
                                          )
                                          )
                                          )
                                                )
                                              )
                                            )
                                          )
                                        )
                                      )
                                    )
                                  )
                                )
                              )
                            )
                          )
                        )
                      )
                    )
                  )
                )
              )
            )
          )
        )
      )
    )
  )
)
library(plotly)

group_totals <- annotation_counts_df_matched %>%
  group_by(Group) %>%
  summarise(TotalFreq = sum(Freq)) %>%
  mutate(GroupLabel = paste0(Group, " (", TotalFreq, ")"))
annotation_counts_df_matched <- annotation_counts_df_matched %>%
  left_join(group_totals, by = "Group")

fig <- plot_ly(annotation_counts_df_matched, 
               x = ~all_annotations,      # X-axis = annotations
               y = ~Freq,  
               color = ~GroupLabel,
               colors = "Set1",# Y-axis = frequency
               type = 'bar',                 # Bar plot
               text = ~paste('Frequency: ', Freq, 'Annotation: ', all_annotations ),  # Tooltip with frequency
               hoverinfo = 'text'        # Show tooltip on hover 
               )  # Color of bars

# Set plot title and axis labels, reduce x-axis text size
fig <- fig %>% layout(
  title = 'Frequency of InterPro Annotations found in known Toxins',
  xaxis = list(
    title = 'Annotation',
    tickangle = 45,               # Rotate the x-axis labels for better readability
    tickfont = list(size = 3)     # Reduce the sample label text size (x-axis labels)
  ),
  yaxis = list(title = 'Frequency')
)

# Show the plot
fig



# Save as static image (PNG)
htmlwidgets::saveWidget(fig, file.path("plotly_graph.html"))
