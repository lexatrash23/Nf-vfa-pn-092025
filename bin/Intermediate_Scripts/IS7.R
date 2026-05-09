#!/usr/bin/env Rscript

library(dplyr)
library(grid)  
library(ggplot2)
library(ggrepel)
library(cowplot)
library(stringr)

#Save plot and legend separate
args <- commandArgs(trailingOnly = TRUE)
Transdf_distinct_file <- args[1] 
sample <- args[2]
summary_IP_file <- args[3]
Summary_MF_file <- args[4]
Summary_BP_file <- args[5]
colourRDS <- args[6]

color_palette <- readRDS(colourRDS)
transdf_distinct <- read.csv(Transdf_distinct_file)

#filter only for complete ORFs with signal sequence 
transdf_distinct <- transdf_distinct %>%
  filter(
    grepl("complete", ORF_type, ignore.case = TRUE),
    grepl("SP", SP, ignore.case = TRUE),
    TMHMM == FALSE
  ) %>%
  filter(percent > 0) %>%
  dplyr::select(Transdecoder_ID, percent,  InterPro_accession_Names, GO_name)

#Filter only for those with IP
transdf_distinct_IP <- transdf_distinct %>%
  dplyr::select(Transdecoder_ID, percent, InterPro_accession_Names) %>%
  filter(!is.na(InterPro_accession_Names) & str_trim(InterPro_accession_Names) != "")

sum(transdf_distinct_IP$percent)
#Read in IP comparisons only keep those that are more represented in toxins than non-toxin proteins 
summary_IP <- read.csv(summary_IP_file)

#As long as it there in toxins
IPOverRepresentedInToxins <- summary_IP %>%
  filter(RelativeExpression > 0 ) %>%
  arrange(desc(Rank)) #the script below takes the last match's rank value 


#Add the highest rank domain present for each transcript
for (i in 1:nrow(IPOverRepresentedInToxins)) {
  entry_ac <- IPOverRepresentedInToxins$InterPro[i]
  matched_rows <- transdf_distinct_IP[grepl(entry_ac, transdf_distinct_IP$InterPro_accession_Names, ignore.case = TRUE), ]
  if (nrow(matched_rows) > 0) {
    best_match <- IPOverRepresentedInToxins[i, ]  # Get the current row in IPInterest (since we're looping over IPInterest)
    transdf_distinct_IP$Rank[transdf_distinct_IP$InterPro_accession_Names %in% matched_rows$InterPro_accession_Names] <- best_match$Rank
  }
}

# Only keep those that are ranked, and join with metadata, calculate combined percent for those with the same domain name
transdf_distinct_IP <- transdf_distinct_IP %>%
  filter(!is.na(Rank) & str_trim(Rank) != "") %>%
  left_join(summary_IP, by = "Rank") %>% 
  arrange(desc(percent)) %>%
  group_by(Name_short) %>%
  mutate(percent = sum(percent, na.rm = TRUE)) %>%
  ungroup() %>%
  distinct(Name_short, .keep_all = TRUE) %>%
  arrange(percent)

SelectIP <- transdf_distinct_IP %>%
  dplyr::select(Name_short,percent) %>%
  mutate(across(everything(), ~str_remove_all(., ",")))

write.csv(SelectIP, paste0(sample, "_RelativeOEDomains.csv"))

#Total percentage represented by this df 
PercentSumIPs <- sum(transdf_distinct_IP$percent)
label_text2 <- paste0("Total Expression Represented: ", round(PercentSumIPs, 2), "%")
#Threshold for labelling
threshold <- 5

#Add relativeproportion column, text labels, and positions for labels 
transdf_distinct_IP <- transdf_distinct_IP %>%
  mutate(RelativeProportion = round(((percent/PercentSumIPs)*100),1)) %>%
  dplyr::select(Transdecoder_ID,Name_short,RelativeProportion ) %>%
  mutate(label_text = ifelse( RelativeProportion >= threshold, str_wrap(paste0(Name_short, ", ", RelativeProportion, "%"), width = 15),"")) %>%
  arrange(Name_short) %>% 
  mutate(csum = rev(cumsum(rev(RelativeProportion))), 
         pos = RelativeProportion/2 + lead(csum, 1),
         pos = if_else(is.na(pos), RelativeProportion/2, pos))





#plot 
IP_plot <- ggplot(transdf_distinct_IP, aes(x = 2, y = RelativeProportion, fill = fct_inorder(Name_short))) +
  geom_bar(stat = "identity", color = "black", width = 1) +
  coord_polar(theta = "y") +
  # Create the white hole in the middle (for a donut chart effect)
  annotate("rect", xmin = 0, xmax = 1.5, ymin = 0, ymax = sum(transdf_distinct_IP$RelativeProportion),
           fill = "white", color = NA) +
  theme_void() +
  labs(fill = "Venom-associated InterPro domains") +
  theme(
    legend.position = "right",
    legend.title = element_text(),
    legend.text = element_text(size = 4, color = "#000000"),
    legend.key.size = unit(0.5, "cm"),
    legend.key.height = unit(0.3, "cm"),
    legend.key.width = unit(0.5, "cm"),
    legend.key.border = element_rect(color = "black", size = 1, linetype = "solid"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  ) +
  geom_label_repel(data = transdf_distinct_IP,
                   aes(y = pos ,x = 2, label = label_text),
                   size = 2, nudge_x = 1.5, show.legend = FALSE, box.padding = 0.5,
                   point.padding = 0.5, force =5,hjust = 0.5) +
  annotate("text", x = 0, y = 0, label = label_text2, size = 2, color = "black", fontface = "bold") +
  labs(title = "Relative expression of Toxin-associated domains") 



plot_nolegend <- IP_plot + theme(legend.position = "none")

legend <- get_legend(IP_plot)
legend_plot <- ggdraw(legend)

ggsave("Plot_1_IP.png", plot_nolegend)
ggsave("Legend_1_IP.png", legend_plot)

#Filter only for those with GO
transdf_distinct_GO <- transdf_distinct %>%
  dplyr::select(Transdecoder_ID, percent, GO_name) %>%
  filter(!is.na(GO_name) & str_trim(GO_name) != "")
transdf_distinct_MF <- transdf_distinct_GO 
#Molecular Functions
Summary_MF <- read.csv(Summary_MF_file)
MFOverRepresentedInToxins <- Summary_MF %>%
  filter(RelativeExpression > 0 ) %>%
  arrange(desc(Rank)) #the script below takes the last match's rank value 

#Add the highest MF  present for each transcript
for (i in 1:nrow(MFOverRepresentedInToxins)) {
  GO <- MFOverRepresentedInToxins$GO_ID[i]
  matched_rows <- transdf_distinct_MF[grepl(GO, transdf_distinct_MF$GO_name, ignore.case = TRUE), ]
  if (nrow(matched_rows) > 0) {
    best_match <- MFOverRepresentedInToxins[i, ]  # Get the current row in GO 
    transdf_distinct_MF$Rank[transdf_distinct_MF$GO_name %in% matched_rows$GO_name] <- best_match$Rank
  }
}


# Only keep those that are ranked, and join with metadata, calculate combined percent for those with the same domain name
transdf_distinct_MF <- transdf_distinct_MF %>%
  filter(!is.na(Rank) & str_trim(Rank) != "") %>%
  left_join(Summary_MF, by = "Rank") %>% 
  arrange(desc(percent)) %>%
  group_by(Name) %>%
  mutate(percent = sum(percent, na.rm = TRUE)) %>%
  ungroup() %>%
  distinct(Name, .keep_all = TRUE) %>%
  arrange(percent)

SelectMF <- transdf_distinct_MF %>%
  dplyr::select(Name,percent) %>%
  mutate(across(everything(), ~str_remove_all(., ",")))

write.csv(SelectMF, paste0(sample, "_RelativeOEMFs.csv"))

#Total percentage represented by this df 
PercentSumMFs <- sum(transdf_distinct_MF$percent)
label_text2 <- paste0("Total Expression Represented: ", round(PercentSumMFs, 2), "%")
#Threshold for labelling
threshold <- 5

#Add relativeproportion column, text labels, and positions for labels 
transdf_distinct_MF <- transdf_distinct_MF %>%
  mutate(RelativeProportion = round(((percent/PercentSumMFs)*100),1)) %>%
  dplyr::select(Transdecoder_ID,Name,RelativeProportion ) %>%
  mutate(label_text = ifelse( RelativeProportion >= threshold, str_wrap(paste0(Name, ", ", RelativeProportion, "%"), width = 15),"")) %>%
  arrange(Name) %>% 
  mutate(csum = rev(cumsum(rev(RelativeProportion))), 
         pos = RelativeProportion/2 + lead(csum, 1),
         pos = if_else(is.na(pos), RelativeProportion/2, pos))



MF_plot <- ggplot(transdf_distinct_MF, aes(x = 2, y = RelativeProportion, fill = fct_inorder(Name))) +
  geom_bar(stat = "identity", color = "black", width = 1) +
  coord_polar(theta = "y") +
  # Create the white hole in the middle (for a donut chart effect)
  annotate("rect", xmin = 0, xmax = 1.5, ymin = 0, ymax = sum(transdf_distinct_MF$RelativeProportion),
           fill = "white", color = NA) +
  theme_void() +
  labs(fill = "Venom-associated Molecular Functions") +
  theme(
    legend.position = "right",
    legend.title = element_text(),
    legend.text = element_text(size = 4, color = "#000000"),
    legend.key.size = unit(0.5, "cm"),
    legend.key.height = unit(0.3, "cm"),
    legend.key.width = unit(0.5, "cm"),
    legend.key.border = element_rect(color = "black", size = 1, linetype = "solid"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  ) +
  geom_label_repel(data = transdf_distinct_MF,
                   aes(y = pos ,x = 2, label = label_text),
                   size = 2, nudge_x = 1.5, show.legend = FALSE, box.padding = 0.5,
                   point.padding = 0.5, force =5,hjust = 0.5) +
  annotate("text", x = 0, y = 0, label = label_text2, size = 2, color = "black", fontface = "bold") +
  labs(title = "Relative expression of Toxin-associated Molecular Functions") 


plot_nolegend <- MF_plot + theme(legend.position = "none")
legend <- get_legend(MF_plot)
legend_plot <- ggdraw(legend)

ggsave("Plot_1_MF.png", plot_nolegend)
ggsave("Legend_1_MF.png", legend_plot)

#Biological Processess
Summary_BP <- read.csv(Summary_BP_file)
BPOverRepresentedInToxins <- Summary_BP %>%
  filter(RelativeExpression > 0 )%>%
  arrange(desc(Rank)) #the script below takes the last match's rank value 
transdf_distinct_BP <- transdf_distinct_GO 

#Add the highest MF  present for each transcript
for (i in 1:nrow(BPOverRepresentedInToxins)) {
  GO <- BPOverRepresentedInToxins$GO_ID[i]
  matched_rows <- transdf_distinct_BP[grepl(GO, transdf_distinct_BP$GO_name, ignore.case = TRUE), ]
  if (nrow(matched_rows) > 0) {
    best_match <- BPOverRepresentedInToxins[i, ]  # Get the current row in GO 
    transdf_distinct_BP$Rank[transdf_distinct_BP$GO_name %in% matched_rows$GO_name] <- best_match$Rank
  }
}


# Only keep those that are ranked, and join with metadata, calculate combined percent for those with the same domain name
transdf_distinct_BP <- transdf_distinct_BP %>%
  filter(!is.na(Rank) & str_trim(Rank) != "") %>%
  left_join(Summary_BP, by = "Rank") %>% 
  arrange(desc(percent)) %>%
  group_by(Name) %>%
  mutate(percent = sum(percent, na.rm = TRUE)) %>%
  ungroup() %>%
  distinct(Name, .keep_all = TRUE) %>%
  arrange(percent)

SelectBP <- transdf_distinct_BP %>%
  dplyr::select(Name,percent) %>%
  mutate(across(everything(), ~str_remove_all(., ",")))

write.csv(SelectBP, paste0(sample, "_RelativeOEBPs.csv"))


#Total percentage represented by this df 
PercentSumBPs <- sum(transdf_distinct_BP$percent)
label_text2 <- paste0("Total Expression Represented: ", round(PercentSumBPs, 2), "%")
#Threshold for labelling
threshold <- 5

#Add relativeproportion column, text labels, and positions for labels 
transdf_distinct_BP <- transdf_distinct_BP %>%
  mutate(RelativeProportion = round(((percent/PercentSumBPs)*100),1)) %>%
  dplyr::select(Transdecoder_ID,Name,RelativeProportion ) %>%
  mutate(label_text = ifelse( RelativeProportion >= threshold, str_wrap(paste0(Name, ", ", RelativeProportion, "%"), width = 15),"")) %>%
  arrange(Name) %>% 
  mutate(csum = rev(cumsum(rev(RelativeProportion))), 
         pos = RelativeProportion/2 + lead(csum, 1),
         pos = if_else(is.na(pos), RelativeProportion/2, pos))



BP_plot <- ggplot(transdf_distinct_BP, aes(x = 2, y = RelativeProportion, fill = fct_inorder(Name))) +
  geom_bar(stat = "identity", color = "black", width = 1) +
  coord_polar(theta = "y") +
  # Create the white hole in the middle (for a donut chart effect)
  annotate("rect", xmin = 0, xmax = 1.5, ymin = 0, ymax = sum(transdf_distinct_BP$RelativeProportion),
           fill = "white", color = NA) +
  theme_void() +
  labs(fill = "Venom-associated Molecular Functions") +
  theme(
    legend.position = "right",
    legend.title = element_text(),
    legend.text = element_text(size = 4, color = "#000000"),
    legend.key.size = unit(0.5, "cm"),
    legend.key.height = unit(0.3, "cm"),
    legend.key.width = unit(0.5, "cm"),
    legend.key.border = element_rect(color = "black", size = 1, linetype = "solid"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  ) +
  geom_label_repel(data = transdf_distinct_BP,
                   aes(y = pos ,x = 2, label = label_text),
                   size = 2, nudge_x = 1.5, show.legend = FALSE, box.padding = 0.5,
                   point.padding = 0.5, force =5,hjust = 0.5) +
  annotate("text", x = 0, y = 0, label = label_text2, size = 2, color = "black", fontface = "bold") +
  labs(title = "Relative expression of Toxin-associated Biological Processess") 

plot_nolegend <- BP_plot + theme(legend.position = "none")
legend <- get_legend(BP_plot)
legend_plot <- ggdraw(legend)

ggsave("Plot_1_BP.png", plot_nolegend)
ggsave("Legend_1_BP.png", legend_plot)
##Now just those that are overrepresented in toxins RE>1

transdf_distinct <- read.csv(Transdf_distinct_file)



#filter only for complete ORFs with signal sequence 
transdf_distinct <- transdf_distinct %>%
  filter(
    grepl("complete", ORF_type, ignore.case = TRUE),
    grepl("SP", SP, ignore.case = TRUE),
    TMHMM == FALSE
  ) %>%
  filter(percent > 0) %>%
  dplyr::select(Transdecoder_ID, percent,  InterPro_accession_Names, GO_name)

#Filter only for those with IP
transdf_distinct_IP <- transdf_distinct %>%
  dplyr::select(Transdecoder_ID, percent, InterPro_accession_Names) %>%
  filter(!is.na(InterPro_accession_Names) & str_trim(InterPro_accession_Names) != "")

#Only those that are overexpressed
#Read in IP comparisons only keep those that are more represented in toxins than non-toxin proteins 
summary_IP <- read.csv(summary_IP_file)
IPOverRepresentedInToxins <- summary_IP %>%
  filter(RelativeExpression > 1 ) %>%
  arrange(desc(Rank)) #the script below takes the last match's rank value 

#Add the highest rank domain present for each transcript
for (i in 1:nrow(IPOverRepresentedInToxins)) {
  entry_ac <- IPOverRepresentedInToxins$InterPro[i]
  matched_rows <- transdf_distinct_IP[grepl(entry_ac, transdf_distinct_IP$InterPro_accession_Names, ignore.case = TRUE), ]
  if (nrow(matched_rows) > 0) {
    best_match <- IPOverRepresentedInToxins[i, ]  # Get the current row in IPInterest (since we're looping over IPInterest)
    transdf_distinct_IP$Rank[transdf_distinct_IP$InterPro_accession_Names %in% matched_rows$InterPro_accession_Names] <- best_match$Rank
  }
}

# Only keep those that are ranked, and join with metadata, calculate combined percent for those with the same domain name
transdf_distinct_IP <- transdf_distinct_IP %>%
  filter(!is.na(Rank) & str_trim(Rank) != "") %>%
  left_join(summary_IP, by = "Rank") %>% 
  arrange(desc(percent)) %>%
  group_by(Name_short) %>%
  mutate(percent = sum(percent, na.rm = TRUE)) %>%
  ungroup() %>%
  distinct(Name_short, .keep_all = TRUE) %>%
  arrange(percent)

#Total percentage represented by this df 
PercentSumIPs <- sum(transdf_distinct_IP$percent)
label_text2 <- paste0("Total Expression Represented: ", round(PercentSumIPs, 2), "%")
#Threshold for labelling
threshold <- 5

#Add relativeproportion column, text labels, and positions for labels 
transdf_distinct_IP <- transdf_distinct_IP %>%
  mutate(RelativeProportion = round(((percent/PercentSumIPs)*100),1)) %>%
  dplyr::select(Transdecoder_ID,Name_short,RelativeProportion ) %>%
  mutate(label_text = ifelse( RelativeProportion >= threshold, str_wrap(paste0(Name_short, ", ", RelativeProportion, "%"), width = 15),"")) %>%
  arrange(Name_short) %>% 
  mutate(csum = rev(cumsum(rev(RelativeProportion))), 
         pos = RelativeProportion/2 + lead(csum, 1),
         pos = if_else(is.na(pos), RelativeProportion/2, pos))



#plot 
IP_plot <- ggplot(transdf_distinct_IP, aes(x = 2, y = RelativeProportion, fill = fct_inorder(Name_short))) +
  geom_bar(stat = "identity", color = "black", width = 1) +
  coord_polar(theta = "y") +
  # Create the white hole in the middle (for a donut chart effect)
  annotate("rect", xmin = 0, xmax = 1.5, ymin = 0, ymax = sum(transdf_distinct_IP$RelativeProportion),
           fill = "white", color = NA) +
  theme_void() +
  labs(fill = "Venom-associated InterPro domains") +
  theme(
    legend.position = "right",
    legend.title = element_text(),
    legend.text = element_text(size = 4, color = "#000000"),
    legend.key.size = unit(0.5, "cm"),
    legend.key.height = unit(0.3, "cm"),
    legend.key.width = unit(0.5, "cm"),
    legend.key.border = element_rect(color = "black", size = 1, linetype = "solid"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  ) +
  geom_label_repel(data = transdf_distinct_IP,
                   aes(y = pos ,x = 2, label = label_text),
                   size = 2, nudge_x = 1.5, show.legend = FALSE, box.padding = 0.5,
                   point.padding = 0.5, force =5,hjust = 0.5) +
  annotate("text", x = 0, y = 0, label = label_text2, size = 2, color = "black", fontface = "bold") +
  labs(title = "Relative expression of InterPro domains overrepresented in Venom") 

plot_nolegend <- IP_plot + theme(legend.position = "none")
legend <- get_legend(IP_plot)
legend_plot <- ggdraw(legend)

ggsave("Plot_2_IP.png", plot_nolegend)
ggsave("Legend_2_IP.png", legend_plot)

#Filter only for those with GO
transdf_distinct_GO <- transdf_distinct %>%
  dplyr::select(Transdecoder_ID, percent, GO_name) %>%
  filter(!is.na(GO_name) & str_trim(GO_name) != "")
transdf_distinct_MF <- transdf_distinct_GO
#Molecular Functions
Summary_MF <- read.csv(Summary_MF_file)
MFOverRepresentedInToxins <- Summary_MF %>%
  filter(RelativeExpression > 1 ) %>%
  arrange(desc(Rank)) #the script below takes the last match's rank value 

#Add the highest MF  present for each transcript
for (i in 1:nrow(MFOverRepresentedInToxins)) {
  GO <- MFOverRepresentedInToxins$GO_ID[i]
  matched_rows <- transdf_distinct_MF[grepl(GO, transdf_distinct_MF$GO_name, ignore.case = TRUE), ]
  if (nrow(matched_rows) > 0) {
    best_match <- MFOverRepresentedInToxins[i, ]  # Get the current row in GO 
    transdf_distinct_MF$Rank[transdf_distinct_MF$GO_name %in% matched_rows$GO_name] <- best_match$Rank
  }
}


# Only keep those that are ranked, and join with metadata, calculate combined percent for those with the same domain name
transdf_distinct_MF <- transdf_distinct_MF %>%
  filter(!is.na(Rank) & str_trim(Rank) != "") %>%
  left_join(Summary_MF, by = "Rank") %>% 
  arrange(desc(percent)) %>%
  group_by(Name) %>%
  mutate(percent = sum(percent, na.rm = TRUE)) %>%
  ungroup() %>%
  distinct(Name, .keep_all = TRUE) %>%
  arrange(percent)

#Total percentage represented by this df 
PercentSumMFs <- sum(transdf_distinct_MF$percent)
label_text2 <- paste0("Total Expression Represented: ", round(PercentSumMFs, 2), "%")
#Threshold for labelling
threshold <- 5

#Add relativeproportion column, text labels, and positions for labels 
transdf_distinct_MF <- transdf_distinct_MF %>%
  mutate(RelativeProportion = round(((percent/PercentSumMFs)*100),1)) %>%
  dplyr::select(Transdecoder_ID,Name,RelativeProportion ) %>%
  mutate(label_text = ifelse( RelativeProportion >= threshold, str_wrap(paste0(Name, ", ", RelativeProportion, "%"), width = 15),"")) %>%
  arrange(Name) %>% 
  mutate(csum = rev(cumsum(rev(RelativeProportion))), 
         pos = RelativeProportion/2 + lead(csum, 1),
         pos = if_else(is.na(pos), RelativeProportion/2, pos))



MF_plot <- ggplot(transdf_distinct_MF, aes(x = 2, y = RelativeProportion, fill = fct_inorder(Name))) +
  geom_bar(stat = "identity", color = "black", width = 1) +
  coord_polar(theta = "y") +
  # Create the white hole in the middle (for a donut chart effect)
  annotate("rect", xmin = 0, xmax = 1.5, ymin = 0, ymax = sum(transdf_distinct_MF$RelativeProportion),
           fill = "white", color = NA) +
  theme_void() +
  labs(fill = "Venom-associated Molecular Functions") +
  theme(
    legend.position = "right",
    legend.title = element_text(),
    legend.text = element_text(size = 4, color = "#000000"),
    legend.key.size = unit(0.5, "cm"),
    legend.key.height = unit(0.3, "cm"),
    legend.key.width = unit(0.5, "cm"),
    legend.key.border = element_rect(color = "black", size = 1, linetype = "solid"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  ) +
  geom_label_repel(data = transdf_distinct_MF,
                   aes(y = pos ,x = 2, label = label_text),
                   size = 2, nudge_x = 1.5, show.legend = FALSE, box.padding = 0.5,
                   point.padding = 0.5, force =5,hjust = 0.5) +
  annotate("text", x = 0, y = 0, label = label_text2, size = 2, color = "black", fontface = "bold") +
  labs(title = "Relative expression of Molecular Functions overrepresented in Venom") 

plot_nolegend <- MF_plot + theme(legend.position = "none")
legend <- get_legend(MF_plot)
legend_plot <- ggdraw(legend)

ggsave("Plot_2_MF.png", plot_nolegend)
ggsave("Legend_2_MF.png", legend_plot)
#Biological Processess
Summary_BP <- read.csv(Summary_BP_file)
BPOverRepresentedInToxins <- Summary_BP %>%
  filter(RelativeExpression > 1 ) %>%
  arrange(desc(Rank)) #the script below takes the last match's rank value 
transdf_distinct_BP <- transdf_distinct_GO 

#Add the highest MF  present for each transcript
for (i in 1:nrow(BPOverRepresentedInToxins)) {
  GO <- BPOverRepresentedInToxins$GO_ID[i]
  matched_rows <- transdf_distinct_BP[grepl(GO, transdf_distinct_BP$GO_name, ignore.case = TRUE), ]
  if (nrow(matched_rows) > 0) {
    best_match <- BPOverRepresentedInToxins[i, ]  # Get the current row in GO 
    transdf_distinct_BP$Rank[transdf_distinct_BP$GO_name %in% matched_rows$GO_name] <- best_match$Rank
  }
}


# Only keep those that are ranked, and join with metadata, calculate combined percent for those with the same domain name
transdf_distinct_BP <- transdf_distinct_BP %>%
  filter(!is.na(Rank) & str_trim(Rank) != "") %>%
  left_join(Summary_BP, by = "Rank") %>% 
  arrange(desc(percent)) %>%
  group_by(Name) %>%
  mutate(percent = sum(percent, na.rm = TRUE)) %>%
  ungroup() %>%
  distinct(Name, .keep_all = TRUE) %>%
  arrange(percent)

#Total percentage represented by this df 
PercentSumBPs <- sum(transdf_distinct_BP$percent)
label_text2 <- paste0("Total Expression Represented: ", round(PercentSumBPs, 2), "%")
#Threshold for labelling
threshold <- 5

#Add relativeproportion column, text labels, and positions for labels 
transdf_distinct_BP <- transdf_distinct_BP %>%
  mutate(RelativeProportion = round(((percent/PercentSumBPs)*100),1)) %>%
  dplyr::select(Transdecoder_ID,Name,RelativeProportion ) %>%
  mutate(label_text = ifelse( RelativeProportion >= threshold, str_wrap(paste0(Name, ", ", RelativeProportion, "%"), width = 15),"")) %>%
  arrange(Name) %>% 
  mutate(csum = rev(cumsum(rev(RelativeProportion))), 
         pos = RelativeProportion/2 + lead(csum, 1),
         pos = if_else(is.na(pos), RelativeProportion/2, pos))



BP_plot <- ggplot(transdf_distinct_BP, aes(x = 2, y = RelativeProportion, fill = fct_inorder(Name))) +
  geom_bar(stat = "identity", color = "black", width = 1) +
  coord_polar(theta = "y") +
  # Create the white hole in the middle (for a donut chart effect)
  annotate("rect", xmin = 0, xmax = 1.5, ymin = 0, ymax = sum(transdf_distinct_BP$RelativeProportion),
           fill = "white", color = NA) +
  theme_void() +
  labs(fill = "Venom-associated Biological Processess") +
  theme(
    legend.position = "right",
    legend.title = element_text(),
    legend.text = element_text(size = 4, color = "#000000"),
    legend.key.size = unit(0.5, "cm"),
    legend.key.height = unit(0.3, "cm"),
    legend.key.width = unit(0.5, "cm"),
    legend.key.border = element_rect(color = "black", size = 1, linetype = "solid"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  ) +
  geom_label_repel(data = transdf_distinct_BP,
                   aes(y = pos ,x = 2, label = label_text),
                   size = 2, nudge_x = 1.5, show.legend = FALSE, box.padding = 0.5,
                   point.padding = 0.5, force =5,hjust = 0.5) +
  annotate("text", x = 0, y = 0, label = label_text2, size = 2, color = "black", fontface = "bold") +
  labs(title = "Relative expression of Biological Processess overrepresented in Venom") 

plot_nolegend <- BP_plot + theme(legend.position = "none")
legend <- get_legend(BP_plot)
legend_plot <- ggdraw(legend)

ggsave("Plot_2_BP.png", plot_nolegend)
ggsave("Legend_2_BP.png", legend_plot)

# NOW ALL domains and MF and BP 
#All 
transdf_distinct <- read.csv(Transdf_distinct_file)


#filter only for complete ORFs with signal sequence 
transdf_distinct <- transdf_distinct %>%
  filter(
    grepl("complete", ORF_type, ignore.case = TRUE),
    grepl("SP", SP, ignore.case = TRUE),
    TMHMM == FALSE
  ) %>%
  filter(percent > 0) %>%
  dplyr::select(Transdecoder_ID, percent,  InterPro_accession_Names, GO_name)

#Filter only for those with IP
transdf_distinct_IP <- transdf_distinct %>%
  dplyr::select(Transdecoder_ID, percent, InterPro_accession_Names) %>%
  filter(!is.na(InterPro_accession_Names) & str_trim(InterPro_accession_Names) != "")


#Read in IP comparisons only keep those that are more represented in toxins than non-toxin proteins 
summary_IP <- read.csv(summary_IP_file)
IPOverRepresentedInToxins <- summary_IP %>%
  arrange(desc(Rank)) #the script below takes the last match's rank value 

#Add the highest rank domain present for each transcript
for (i in 1:nrow(IPOverRepresentedInToxins)) {
  entry_ac <- IPOverRepresentedInToxins$InterPro[i]
  matched_rows <- transdf_distinct_IP[grepl(entry_ac, transdf_distinct_IP$InterPro_accession_Names, ignore.case = TRUE), ]
  if (nrow(matched_rows) > 0) {
    best_match <- IPOverRepresentedInToxins[i, ]  # Get the current row in IPInterest (since we're looping over IPInterest)
    transdf_distinct_IP$Rank[transdf_distinct_IP$InterPro_accession_Names %in% matched_rows$InterPro_accession_Names] <- best_match$Rank
  }
}

# Only keep those that are ranked, and join with metadata, calculate combined percent for those with the same domain name
transdf_distinct_IP <- transdf_distinct_IP %>%
  filter(!is.na(Rank) & str_trim(Rank) != "") %>%
  left_join(summary_IP, by = "Rank") %>% 
  arrange(desc(percent)) %>%
  group_by(Name_short) %>%
  mutate(percent = sum(percent, na.rm = TRUE)) %>%
  ungroup() %>%
  distinct(Name_short, .keep_all = TRUE) %>%
  arrange(percent)

#Total percentage represented by this df 
PercentSumIPs <- sum(transdf_distinct_IP$percent)
label_text2 <- paste0("Total Expression Represented: ", round(PercentSumIPs, 2), "%")
#Threshold for labelling
threshold <- 5

#Add relativeproportion column, text labels, and positions for labels 
transdf_distinct_IP <- transdf_distinct_IP %>%
  mutate(RelativeProportion = round(((percent/PercentSumIPs)*100),1)) %>%
  dplyr::select(Transdecoder_ID,Name_short,RelativeProportion ) %>%
  mutate(label_text = ifelse( RelativeProportion >= threshold, str_wrap(paste0(Name_short, ", ", RelativeProportion, "%"), width = 15),"")) %>%
  arrange(Name_short) %>% 
  mutate(csum = rev(cumsum(rev(RelativeProportion))), 
         pos = RelativeProportion/2 + lead(csum, 1),
         pos = if_else(is.na(pos), RelativeProportion/2, pos))



#plot 
IP_plot <- ggplot(transdf_distinct_IP, aes(x = 2, y = RelativeProportion, fill = fct_inorder(Name_short))) +
  geom_bar(stat = "identity", color = "black", width = 1) +
  coord_polar(theta = "y") +
  # Create the white hole in the middle (for a donut chart effect)
  annotate("rect", xmin = 0, xmax = 1.5, ymin = 0, ymax = sum(transdf_distinct_IP$RelativeProportion),
           fill = "white", color = NA) +
  theme_void() +
  labs(fill = "InterPro domains") +
  theme(
    legend.position = "right",
    legend.title = element_text(),
    legend.text = element_text(size = 4, color = "#000000"),
    legend.key.size = unit(0.5, "cm"),
    legend.key.height = unit(0.3, "cm"),
    legend.key.width = unit(0.5, "cm"),
    legend.key.border = element_rect(color = "black", size = 1, linetype = "solid"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  ) +
  geom_label_repel(data = transdf_distinct_IP,
                   aes(y = pos ,x = 2, label = label_text),
                   size = 2, nudge_x = 1.5, show.legend = FALSE, box.padding = 0.5,
                   point.padding = 0.5, force =5,hjust = 0.5) +
  annotate("text", x = 0, y = 0, label = label_text2, size = 2, color = "black", fontface = "bold") +
  labs(title = "Relative expression of InterPro domains") 
plot_nolegend <- IP_plot + theme(legend.position = "none")
legend <- get_legend(IP_plot)
legend_plot <- ggdraw(legend)

ggsave("Plot_3_IP.png", plot_nolegend)
ggsave("Legend_3_IP.png", legend_plot, width = 20, height = 16)
#Filter only for those with GO
transdf_distinct_GO <- transdf_distinct %>%
  dplyr::select(Transdecoder_ID, percent, GO_name) %>%
  filter(!is.na(GO_name) & str_trim(GO_name) != "")
transdf_distinct_MF <- transdf_distinct_GO
#Molecular Functions
Summary_MF <- read.csv(Summary_MF_file)
MFOverRepresentedInToxins <- Summary_MF %>%
  arrange(desc(Rank)) #the script below takes the last match's rank value 

#Add the highest MF  present for each transcript
for (i in 1:nrow(MFOverRepresentedInToxins)) {
  GO <- MFOverRepresentedInToxins$GO_ID[i]
  matched_rows <- transdf_distinct_MF[grepl(GO, transdf_distinct_MF$GO_name, ignore.case = TRUE), ]
  if (nrow(matched_rows) > 0) {
    best_match <- MFOverRepresentedInToxins[i, ]  # Get the current row in GO 
    transdf_distinct_MF$Rank[transdf_distinct_MF$GO_name %in% matched_rows$GO_name] <- best_match$Rank
  }
}


# Only keep those that are ranked, and join with metadata, calculate combined percent for those with the same domain name
transdf_distinct_MF <- transdf_distinct_MF %>%
  filter(!is.na(Rank) & str_trim(Rank) != "") %>%
  left_join(Summary_MF, by = "Rank") %>% 
  arrange(desc(percent)) %>%
  group_by(Name) %>%
  mutate(percent = sum(percent, na.rm = TRUE)) %>%
  ungroup() %>%
  distinct(Name, .keep_all = TRUE) %>%
  arrange(percent)

#Total percentage represented by this df 
PercentSumMFs <- sum(transdf_distinct_MF$percent)
label_text2 <- paste0("Total Expression Represented: ", round(PercentSumMFs, 2), "%")
#Threshold for labelling
threshold <- 5

#Add relativeproportion column, text labels, and positions for labels 
transdf_distinct_MF <- transdf_distinct_MF %>%
  mutate(RelativeProportion = round(((percent/PercentSumMFs)*100),1)) %>%
  dplyr::select(Transdecoder_ID,Name,RelativeProportion ) %>%
  mutate(label_text = ifelse( RelativeProportion >= threshold, str_wrap(paste0(Name, ", ", RelativeProportion, "%"), width = 15),"")) %>%
  arrange(Name) %>% 
  mutate(csum = rev(cumsum(rev(RelativeProportion))), 
         pos = RelativeProportion/2 + lead(csum, 1),
         pos = if_else(is.na(pos), RelativeProportion/2, pos))



MF_plot <- ggplot(transdf_distinct_MF, aes(x = 2, y = RelativeProportion, fill = fct_inorder(Name))) +
  geom_bar(stat = "identity", color = "black", width = 1) +
  coord_polar(theta = "y") +
  # Create the white hole in the middle (for a donut chart effect)
  annotate("rect", xmin = 0, xmax = 1.5, ymin = 0, ymax = sum(transdf_distinct_MF$RelativeProportion),
           fill = "white", color = NA) +
  theme_void() +
  labs(fill = "Molecular Functions") +
  theme(
    legend.position = "right",
    legend.title = element_text(),
    legend.text = element_text(size = 4, color = "#000000"),
    legend.key.size = unit(0.5, "cm"),
    legend.key.height = unit(0.3, "cm"),
    legend.key.width = unit(0.5, "cm"),
    legend.key.border = element_rect(color = "black", size = 1, linetype = "solid"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  ) +
  geom_label_repel(data = transdf_distinct_MF,
                   aes(y = pos ,x = 2, label = label_text),
                   size = 2, nudge_x = 1.5, show.legend = FALSE, box.padding = 0.5,
                   point.padding = 0.5, force =5,hjust = 0.5) +
  annotate("text", x = 0, y = 0, label = label_text2, size = 2, color = "black", fontface = "bold") +
  labs(title = "Relative expression of Molecular Functions") 

plot_nolegend <- MF_plot + theme(legend.position = "none")
legend <- get_legend(MF_plot)
legend_plot <- ggdraw(legend)
ggsave("Plot_3_MF.png", plot_nolegend)
ggsave("Legend_3_MF.png", legend_plot, width = 20, height = 16)

#Biological Processess
Summary_BP <- read.csv(Summary_BP_file)
BPOverRepresentedInToxins <- Summary_BP  %>%
  arrange(desc(Rank)) 
transdf_distinct_BP <- transdf_distinct_GO  #the script below takes the last match's rank value 

#Add the highest BP  present for each transcript
for (i in 1:nrow(BPOverRepresentedInToxins)) {
  GO <- BPOverRepresentedInToxins$GO_ID[i]
  matched_rows <- transdf_distinct_BP[grepl(GO, transdf_distinct_BP$GO_name, ignore.case = TRUE), ]
  if (nrow(matched_rows) > 0) {
    best_match <- BPOverRepresentedInToxins[i, ]  # Get the current row in GO 
    transdf_distinct_BP$Rank[transdf_distinct_BP$GO_name %in% matched_rows$GO_name] <- best_match$Rank
  }
}


# Only keep those that are ranked, and join with metadata, calculate combined percent for those with the same domain name
transdf_distinct_BP <- transdf_distinct_BP %>%
  filter(!is.na(Rank) & str_trim(Rank) != "") %>%
  left_join(Summary_BP, by = "Rank") %>% 
  arrange(desc(percent)) %>%
  group_by(Name) %>%
  mutate(percent = sum(percent, na.rm = TRUE)) %>%
  ungroup() %>%
  distinct(Name, .keep_all = TRUE) %>%
  arrange(percent)

#Total percentage represented by this df 
PercentSumBPs <- sum(transdf_distinct_BP$percent)
label_text2 <- paste0("Total Expression Represented: ", round(PercentSumBPs, 2), "%")
#Threshold for labelling
threshold <- 5

#Add relativeproportion column, text labels, and positions for labels 
transdf_distinct_BP <- transdf_distinct_BP %>%
  mutate(RelativeProportion = round(((percent/PercentSumBPs)*100),1)) %>%
  dplyr::select(Transdecoder_ID,Name,RelativeProportion ) %>%
  mutate(label_text = ifelse( RelativeProportion >= threshold, str_wrap(paste0(Name, ", ", RelativeProportion, "%"), width = 15),"")) %>%
  arrange(Name) %>% 
  mutate(csum = rev(cumsum(rev(RelativeProportion))), 
         pos = RelativeProportion/2 + lead(csum, 1),
         pos = if_else(is.na(pos), RelativeProportion/2, pos))



BP_plot <- ggplot(transdf_distinct_BP, aes(x = 2, y = RelativeProportion, fill = fct_inorder(Name))) +
  geom_bar(stat = "identity", color = "black", width = 1) +
  coord_polar(theta = "y") +
  # Create the white hole in the middle (for a donut chart effect)
  annotate("rect", xmin = 0, xmax = 1.5, ymin = 0, ymax = sum(transdf_distinct_BP$RelativeProportion),
           fill = "white", color = NA) +
  theme_void() +
  labs(fill = "Biological Processess") +
  theme(
    legend.position = "right",
    legend.title = element_text(),
    legend.text = element_text(size = 4, color = "#000000"),
    legend.key.size = unit(0.5, "cm"),
    legend.key.height = unit(0.3, "cm"),
    legend.key.width = unit(0.5, "cm"),
    legend.key.border = element_rect(color = "black", size = 1, linetype = "solid"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  ) +
  geom_label_repel(data = transdf_distinct_BP,
                   aes(y = pos ,x = 2, label = label_text),
                   size = 2, nudge_x = 1.5, show.legend = FALSE, box.padding = 0.5,
                   point.padding = 0.5, force =5,hjust = 0.5) +
  annotate("text", x = 0, y = 0, label = label_text2, size = 2, color = "black", fontface = "bold") +
  labs(title = "Relative expression of Biological Processess") 

plot_nolegend <- BP_plot + theme(legend.position = "none")
legend <- get_legend(BP_plot)
legend_plot <- ggdraw(legend)

ggsave("Plot_3_BP.png", plot_nolegend)
ggsave("Legend_3_BP.png", legend_plot, width = 20, height = 16)
##ANNOTATION 
#Labelling the csv with this information - Will be used in annotate script 

transdf_distinct <- read.csv(Transdf_distinct_file)


#filter only for complete ORFs with signal sequence 
transdf_distinct <- transdf_distinct %>%
  filter(
    grepl("complete", ORF_type, ignore.case = TRUE),
    grepl("SP", SP, ignore.case = TRUE),
    TMHMM == FALSE
  ) %>%
  filter(percent > 0) %>%
  dplyr::select(Transdecoder_ID,  InterPro_accession_Names, GO_name)


summary_IP <- read.csv(summary_IP_file) %>%
  dplyr::rename(IPRank = Rank) %>%
  arrange(desc(IPRank)) #the script below takes the last match's rank value 

#Add the highest rank domain present for each transcript
for (i in 1:nrow(summary_IP)) {
  entry_ac <- summary_IP$InterPro[i]
  matched_rows <- transdf_distinct[grepl(entry_ac, transdf_distinct$InterPro_accession_Names, ignore.case = TRUE), ]
  if (nrow(matched_rows) > 0) {
    best_match <- summary_IP[i, ]  # Get the current row in IPInterest (since we're looping over IPInterest)
    transdf_distinct$IPRank[transdf_distinct$InterPro_accession_Names %in% matched_rows$InterPro_accession_Names] <- best_match$IPRank
  }
}

# Only keep those that are ranked, and join with metadata, calculate combined percent for those with the same domain name
transdf_distinct <- transdf_distinct %>%
  left_join(summary_IP, by = "IPRank") %>% 
  dplyr::select(Transdecoder_ID,GO_name, Name_short) %>%
  dplyr::rename(Domain_label = Name_short )



#Molecular Functions
Summary_MF <- read.csv(Summary_MF_file)
MFOverRepresentedInToxins <- Summary_MF %>%
  arrange(desc(Rank)) #the script below takes the last match's rank value 

#Add the highest MF  present for each transcript
for (i in 1:nrow(MFOverRepresentedInToxins)) {
  GO <- MFOverRepresentedInToxins$GO_ID[i]
  matched_rows <- transdf_distinct[grepl(GO, transdf_distinct$GO_name, ignore.case = TRUE), ]
  if (nrow(matched_rows) > 0) {
    best_match <- MFOverRepresentedInToxins[i, ]  # Get the current row in GO 
    transdf_distinct$Rank[transdf_distinct$GO_name %in% matched_rows$GO_name] <- best_match$Rank
  }
}
# Only keep those that are ranked, and join with metadata, calculate combined percent for those with the same domain name
transdf_distinct <- transdf_distinct %>%
  left_join(MFOverRepresentedInToxins, by = "Rank") %>% 
  dplyr::select(-GO_ID,-RelativeExpression, -Rank, -Ontology) %>%
  dplyr::rename(Molecular_Function = Name )


#Biological Processess
Summary_BP <- read.csv(Summary_BP_file)
BPOverRepresentedInToxins <- Summary_BP %>%
  arrange(desc(Rank)) 

#Add the highest MF  present for each transcript
for (i in 1:nrow(BPOverRepresentedInToxins)) {
  GO <- BPOverRepresentedInToxins$GO_ID[i]
  matched_rows <- transdf_distinct[grepl(GO, transdf_distinct$GO_name, ignore.case = TRUE), ]
  if (nrow(matched_rows) > 0) {
    best_match <- BPOverRepresentedInToxins[i, ]  # Get the current row in GO 
    transdf_distinct$Rank[transdf_distinct$GO_name %in% matched_rows$GO_name] <- best_match$Rank
  }
}

# Only keep those that are ranked, and join with metadata, calculate combined percent for those with the same domain name
transdf_distinct <- transdf_distinct %>%
  left_join(BPOverRepresentedInToxins, by = "Rank") %>% 
  dplyr::select(-GO_ID,-RelativeExpression, -Rank, -Ontology,-GO_name) %>%
  dplyr::rename(Biological_Process = Name )


write.csv(transdf_distinct, "Transdf_cORFSP_unfiltered_Domain_Annotation.csv")

