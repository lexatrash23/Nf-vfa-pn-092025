#!/usr/bin/env Rscript
# Figure generations

library(dplyr)
library(ggplot2)
library(ggalluvial)
library(grid)
library(stringr)
library(ggrepel)
library(forcats)
library(cowplot)

#check of current working directory. The script assumes the working directory is the
#locally
#cluster 

args <- commandArgs(trailingOnly = TRUE)
TBK <- args[1]
colours <- args[2]
#read in TBK dataframe

TBK <- read.csv(TBK, header = TRUE)
keeps <- c("Trinity_ID", "Length", "percent", "cumulativepercent", "Code", "Hit", "E_value", "BitScore")
Kallisto_Blastx_Trinity <- TBK[keeps]

# distinct Transcripts
Distinct_Transcripts <- Kallisto_Blastx_Trinity %>%
  distinct(Trinity_ID, .keep_all = TRUE)

# piechart for percentage of transcripts with no hit vs hit
total_number_of_transcripts = nrow(Distinct_Transcripts)

# create the data frame with counts
counts <- data.frame(
  Category = c("With Hit", "No Hit"),
  Count = c(sum(!is.na(Distinct_Transcripts$Hit)), sum(is.na(Distinct_Transcripts$Hit)))
)

# calcualte percentage
counts <- counts %>%
  mutate(Percentage = Count / sum(Count) * 100,
         Label = paste0("\n", Count, " (", round(Percentage, 1), "%)"))

# pie chart with labels on the outside
pie1 <- ggplot(counts, aes(x = "", y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 1, color = "black") +
  coord_polar(theta = "y", start = 0) +  # Set the starting angle for the slices
  theme_void() +  # Removes axes and background
  labs(title = "% of transcripts with a significant hit to a uniprot annotated toxin(Evalue cutoff:1e-5)") +
  geom_text(aes(label = Label),
            size = 3,  # Increased text size for better visibility
            nudge_x = 0.7,  # Adjust nudging to better position labels outside the pie
            check_overlap = TRUE,  # Avoid label overlap
            color = "#000000") + # Avoid label overlap
  scale_fill_manual(values = c("With Hit" = "#4C9E9A", "No Hit" = "#B0B0B0")) +  # Set colors
  theme(
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5, color = "black"),  # Center title with nice font
    legend.position = "right",  # Position the legend on the right
    legend.title = element_blank(), 
    legend.text = element_text(color = "black"),
    plot.margin = margin(20, 20, 20, 20)  # Add padding around the plot
  )


ggsave(filename = file.path("pie1.png"), plot = pie1, width = 8, height = 6, dpi = 600)


#alluvial graphs 
#only those with hits 
Distinct_Transcripts_hits <- Distinct_Transcripts[!is.na(Distinct_Transcripts$Hit),]
#bitscore 50 cutoff 
Distinct_Transcripts_50 <- Distinct_Transcripts_hits[(Distinct_Transcripts_hits$BitScore >= 50),]
#bitscore 250 cutoff
Distinct_Transcripts_250 <- Distinct_Transcripts_hits[(Distinct_Transcripts_hits$BitScore >= 250),]

alluvial1 <-  ggplot(data = Distinct_Transcripts_50,
                     aes(axis1 =Trinity_ID , axis2 = Hit)) +
  geom_alluvium(aes(fill = Hit)) +
  geom_stratum(aes(fill = Hit)) +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum)), size = 2, min.y = 15) +
  scale_x_discrete(limits = c("Trinity_ID", "Hit"),
                   expand = c(0.15, 0.05)) +  
  theme_void() + theme(legend.position = "right",legend.text = element_text(size = 4, color = "#000000"),  # Smaller text size and color
                       legend.key.size = unit(0.5, "cm"),  # Make the legend keys (colored boxes) smaller
                       legend.key.height = unit(0.3, "cm"),  # Adjust height of the key
                       legend.key.width = unit(0.5, "cm"), plot.title = element_text(size = 14, face = "bold", hjust = -0.5, vjust = 1)) +
  labs(title = "Most significant unitprot toxin hit per transcript")
plot_nolegend <- alluvial1 + theme(legend.position = "none")
legend <- get_legend(alluvial1)
legend_plot <- ggdraw(legend)
w <- convertWidth(sum(legend$widths), "in", valueOnly = TRUE)
h <- convertHeight(sum(legend$heights), "in", valueOnly = TRUE)

ggsave("alluvial1_plot.png", plot = plot_nolegend, width = 8, height = 6, dpi = 600)
ggsave("alluvial1_legend.png", legend_plot,width = w + 0.2,height = h + 0.2, dpi = 600))



alluvial2 <- ggplot(data = Distinct_Transcripts_250,
                    aes(axis1 =Trinity_ID , axis2 = Hit)) +
  geom_alluvium(aes(fill = Hit)) +
  geom_stratum(aes(fill = Hit)) +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum)), size = 2, min.y = 1) +
  scale_x_discrete(limits = c("Trinity_ID", "Hit"),
                   expand = c(0.15, 0.05)) +  
  theme_void() + theme(legend.position = "bottom",legend.text = element_text(size = 8, color = "#000000"),  # Smaller text size and color
                       legend.key.size = unit(0.5, "cm"),  # Make the legend keys (colored boxes) smaller
                       legend.key.height = unit(0.5, "cm"),  # Adjust height of the key
                       legend.key.width = unit(0.5, "cm"), plot.title = element_text(size = 14, face = "bold", hjust = 0.5) ) +
  labs(title = "Most significant unitprot toxin hit per transcript")
plot_nolegend <- alluvial2 + theme(legend.position = "none")
legend <- get_legend(alluvial2)
legend_plot <- ggdraw(legend)
w <- convertWidth(sum(legend$widths), "in", valueOnly = TRUE)
h <- convertHeight(sum(legend$heights), "in", valueOnly = TRUE)

ggsave("alluvial2_plot.png", plot = plot_nolegend, width = 8, height = 6, dpi = 600)
ggsave("alluvial2_legend.png", legend_plot,width = w + 0.2,height = h + 0.2, dpi = 600))

#those with expression 
sum_of_expression_of_those_with_hits = sum(Distinct_Transcripts_hits$percent)
sum_of_expression_of_those_without_hits = 100-sum_of_expression_of_those_with_hits
pie_data <- data.frame(
  Category = c("With Hits", "Without Hits"),
  Value = c(sum_of_expression_of_those_with_hits, sum_of_expression_of_those_without_hits)
)

# pie chart
pie2 <- ggplot(pie_data, aes(x = "", y = Value, fill = Category)) +
  geom_bar(stat = "identity", width = 1, color = "black") +
  coord_polar(theta = "y", start = 0) +  # Set the starting angle for the slices
  theme_void() +  # Removes axes and background
  labs(title = "% of Expression from transcripts with uniprot toxin hits(Evalue cutoff:1e-5)") +
  geom_text(aes(label = paste0(round(Value/sum(Value) * 100, 1), "%")),
            position = position_stack(vjust = 0.5), size = 3) +
  scale_fill_manual(values = c("With Hits" = "#4C9E9A", "Without Hits" = "#B0B0B0")) +
  theme(
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5, color = "black"),  # Center title with nice font
    legend.position = "right",  # Position the legend on the right
    legend.title = element_blank(), 
    legend.text = element_text(color = "black"),
    plot.margin = margin(20, 20, 20, 20)  # Add padding around the plot
  )

ggsave(filename = file.path("pie2.png"), plot = pie2, width = 8, height = 6, dpi = 600)


sum_of_expression_of_those_with_hits = sum(Distinct_Transcripts_50$percent)
sum_of_expression_of_those_without_hits = 100-sum_of_expression_of_those_with_hits
pie_data2 <- data.frame(
  Category = c("With Hits", "Without Hits"),
  Value = c(sum_of_expression_of_those_with_hits, sum_of_expression_of_those_without_hits)
)
# percentage of expression that can be explained by a transcript that has a uniprot toxin match
pie3 <- ggplot(pie_data2, aes(x = "", y = Value, fill = Category)) +
  geom_bar(stat = "identity", width = 1, color = "black") +
  coord_polar(theta = "y", start = 0) +  # Set the starting angle for the slices
  theme_void() +  # Removes axes and background
  labs(title = "% of Expression from transcripts with uniprot toxin hits") +
  geom_text(aes(label = paste0(round(Value/sum(Value) * 100, 1), "%")),
            position = position_stack(vjust = 0.5), size = 3) +
  scale_fill_manual(values = c("With Hits" = "#4C9E9A", "Without Hits" = "#B0B0B0")) +
  theme(
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5, color = "black"),  # Center title with nice font
    legend.position = "right",  # Position the legend on the right
    legend.title = element_blank(), 
    legend.text = element_text(color = "black"),
    plot.margin = margin(20, 20, 20, 20)  # Add padding around the plot
  )
ggsave(filename = file.path("pie3.png"), plot = pie3, width = 8, height = 6, dpi = 600)

#donut graph to show relative expression of each transcript with a uniprot toxin hit 
Distinct_Transcripts_50_with_kallisto <- Distinct_Transcripts_50[Distinct_Transcripts_50$percent >0 , ]
Distinct_Transcripts_50_with_kallisto$hitspercent <- (Distinct_Transcripts_50_with_kallisto$percent / sum(Distinct_Transcripts_50_with_kallisto$percent))*100
new_df2 <- Distinct_Transcripts_50_with_kallisto %>%
  group_by(Hit) %>%
  summarize(total_percentage = sum(hitspercent))
total_percent <- sum(Distinct_Transcripts_50$percent)
total_percent <- sum(Distinct_Transcripts_50_with_kallisto$percent)
new_df2 <-new_df2[order(new_df2$total_percentage, decreasing = TRUE),]
top10mostexpressed <- new_df2
write.csv(top10mostexpressed,file = file.path("Table13.csv"), row.names = FALSE)

#read in my colours and add more # i should really save this as a new palette as some point
color_palette <- readRDS(colours)
new_colors <- c(
  "#FF5733", "#33FF57", "#5733FF", "#FF33A1", "#33FFF5", "#F5FF33", "#A133FF", "#FF8333", "#33A1FF", "#A1FF33",  
  "#FF3385", "#85FF33", "#3385FF", "#FF33D1", "#D1FF33", "#33D1FF", "#FFAF33", "#33FFAF", "#AF33FF", "#FFAF85",  
  "#85FFAF", "#AF85FF", "#FFAFD1", "#D1FFAF", "#AFD1FF", "#FFD133", "#33FFD1", "#D133FF", "#FFD185", "#85FFD1",  
  "#D185FF", "#FFD1AF", "#AFFFD1", "#D1AFFD", "#FF7F50", "#6495ED", "#FFD700", "#00FA9A", "#B22222", "#20B2AA",  
  "#9370DB", "#FF69B4", "#8B4513", "#FF6347", "#2E8B57", "#4169E1", "#8A2BE2", "#00CED1", "#BA55D3", "#228B22",  
  "#E91E63", "#9C27B0", "#3F51B5", "#03A9F4", "#00BCD4", "#009688", "#4CAF50", "#8BC34A", "#CDDC39", "#FFEB3B",  
  "#FFC107", "#FF9800", "#795548", "#607D8B", "#6A5ACD", "#483D8B", "#2F4F4F", "#DAA520", "#A52A2A", "#DC143C",  
  "#800000", "#191970", "#008B8B", "#B0C4DE", "#ADD8E6", "#7FFFD4", "#48D1CC", "#87CEEB", "#F08080", "#FA8072",  
  "#D2691E", "#B8860B", "#008080", "#556B2F", "#8FBC8F", "#66CDAA", "#DB7093", "#F4A460", "#4682B4", "#EE82EE",  
  "#9400D3", "#9932CC", "#C71585", "#D8BFD8", "#FFB6C1", "#FFA07A", "#FFDAB9", "#DDA0DD", "#B0E0E6", "#5F9EA0"
)
extended_color_palette <- c(color_palette, new_colors)
sorted_df2 <- new_df2[order(-new_df2$total_percentage), ]
head(sorted_df2)
PercentSum = sum(Distinct_Transcripts_50_with_kallisto$percent)
label_text2 <- paste0("Total Expression Represented: ", round(PercentSum, 2), "%")
threshold <- 5

#Add relativeproportion column, text labels, and positions for labels
Distinct_Transcripts_50_with_kallisto_mod <- Distinct_Transcripts_50_with_kallisto %>%
  mutate(RelativeProportion = round(((percent/PercentSum)*100),1)) %>%
  dplyr::select(Transdecoder_ID,Hit,RelativeProportion ) %>%
  mutate(label_text = ifelse( RelativeProportion >= threshold, str_wrap(paste0(Hit, ", ", RelativeProportion, "%"), width = 15),"")) %>%
  arrange(Hit) %>%
  mutate(csum = rev(cumsum(rev(RelativeProportion))),
         pos = RelativeProportion/2 + lead(csum, 1),
         pos = if_else(is.na(pos), RelativeProportion/2, pos))



Plot3 <- ggplot(Distinct_Transcripts_50_with_kallisto_mod, aes(x = 2, y = RelativeProportion, fill = fct_inorder(Hit))) +
  geom_bar(stat = "identity", color = "black", width = 1) +
  coord_polar(theta = "y") +
  # Create the white hole in the middle (for a donut chart effect)
  annotate("rect", xmin = 0, xmax = 1.5, ymin = 0, ymax = sum(Distinct_Transcripts_50_with_kallisto_mod$RelativeProportion),
           fill = "white", color = NA) +
  theme_void() +
  labs(fill = "ToxProt Hits") +
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
  geom_label_repel(data = Distinct_Transcripts_50_with_kallisto_mod,
                   aes(y = pos ,x = 2, label = label_text),
                   size = 2, nudge_x = 1.5, show.legend = FALSE, box.padding = 0.5,
                   point.padding = 0.5, force =5,hjust = 0.5) +
  annotate("text", x = 0, y = 0, label = label_text2, size = 2, color = "black", fontface = "bold") +
  labs(title = "Relative expression of ToxProt Hits") +
  scale_color_manual(values = extended_color_palette)

plot_nolegend <- Plot3 + theme(legend.position = "none")
legend <- get_legend(Plot3)
legend_plot <- ggdraw(legend)
w <- convertWidth(sum(legend$widths), "in", valueOnly = TRUE)
h <- convertHeight(sum(legend$heights), "in", valueOnly = TRUE)

ggsave("pie4_plot.png", plot = plot_nolegend, width = 8, height = 6, dpi = 600)
ggsave("pie4_legend.png", legend_plot,width = w + 0.2,height = h + 0.2, dpi = 600))
