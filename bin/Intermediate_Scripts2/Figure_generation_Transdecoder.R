#!/usr/bin/env Rscript
# Figure generations
library(dplyr)
library(ggplot2)
library(ggalluvial)
library(grid)

args <- commandArgs(trailingOnly = TRUE)
transdf <- args[1]
colours <- args[2]

#check of current working directory. The script assumes the working directory is the 
#read in transdf dataframe

transdf <- read.csv(transdf, header = TRUE)
keeps <- c("Transdecoder_ID", "SP_Prediction","percent", "cumulativepercent", "Code", "Hit", "BitScore")
transdf <- transdf[keeps]

# distinct Transcripts
Distinct_Transcripts <- transdf %>%
  distinct(Transdecoder_ID, .keep_all = TRUE)

# piechart for percentage of transcripts with no hit vs hit
total_number_of_transcripts = nrow(Distinct_Transcripts)

# data frame with counts
counts <- data.frame(
  Category = c("With Hit", "No Hit"),
  Count = c(sum(!is.na(Distinct_Transcripts$Hit)), sum(is.na(Distinct_Transcripts$Hit)))
)

#  percentage
counts <- counts %>%
  mutate(Percentage = Count / sum(Count) * 100,
         Label = paste0("\n", Count, " (", round(Percentage, 1), "%)"))

# pie chart with labels on the outside
pie5 <- ggplot(counts, aes(x = "", y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 1, color = "black") +
  coord_polar(theta = "y", start = 0) +  # Set the starting angle for the slices
  theme_void() +  # Removes axes and background
  labs(title = "% of Transdecoder transcripts with a significant hit to a uniprot annotated toxin(bitscore > 50)") +
  geom_text(aes(label = Label), 
            size = 3,  # Increased text size for better visibility
            nudge_x = 0.7,  # Adjust nudging to better position labels outside the pie
            check_overlap = TRUE,  # Avoid label overlap
            color = "black") + # Avoid label overlap
  scale_fill_manual(values = c("With Hit" = "#4C9E9A", "No Hit" = "#B0B0B0")) +  # Set colours
  theme(
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5, color = "black"),  # Center title with nice font
    legend.position = "right",  # Position the legend on the right
    legend.title = element_blank(), 
    legend.text = element_text(color = "black"),
    plot.margin = margin(20, 20, 20, 20)  # Add padding around the plot
  )

ggsave(filename = file.path("pie5.png"), plot = pie5, width = 8, height = 6, dpi = 600)


#alluvial graphs 
#only those with hits 
Distinct_Transcripts_hits <- Distinct_Transcripts[!is.na(Distinct_Transcripts$Hit),]
#bitscore 50 cutoff 
Distinct_Transcripts_50 <- Distinct_Transcripts_hits[(Distinct_Transcripts_hits$BitScore > 50),]
#bitscore 300 cutoff 
Distinct_Transcripts_300 <- Distinct_Transcripts_hits[(Distinct_Transcripts_hits$BitScore > 300),]

alluvial3 <-  ggplot(data = Distinct_Transcripts_50,
                     aes(axis1 =Transdecoder_ID , axis2 = Hit)) +
  geom_alluvium(aes(fill = Hit)) +
  geom_stratum(aes(fill = Hit)) +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum)), size = 2, min.y = 15) +
  scale_x_discrete(limits = c("Transdecoder_ID", "Hit"),
                   expand = c(0.15, 0.05)) +  
  theme_void() + theme(legend.position = "right",legend.text = element_text(size = 4, color = "#000000"),  # Smaller text size and colour
                       legend.key.size = unit(0.5, "cm"),  # Make the legend keys (coloured boxes) smaller
                       legend.key.height = unit(0.3, "cm"),  # Adjust height of the key
                       legend.key.width = unit(0.5, "cm"), plot.title = element_text(size = 14, face = "bold", hjust = -0.5, vjust = 1)) +
  labs(title = "Most significant unitprot toxin hit per transcript(Bitscore >50")

ggsave(filename = file.path("alluvial3.png"), plot = alluvial3, width = 8, height = 6, dpi = 600)

alluvial4 <- ggplot(data = Distinct_Transcripts_300,
                    aes(axis1 =Transdecoder_ID , axis2 = Hit)) +
  geom_alluvium(aes(fill = Hit)) +
  geom_stratum(aes(fill = Hit)) +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum)), size = 2, min.y = 1) +
  scale_x_discrete(limits = c("Transdecoder_ID", "Hit"),
                   expand = c(0.15, 0.05)) +  
  theme_void() + theme(legend.position = "bottom",legend.text = element_text(size = 8, color = "#000000"),  # Smaller text size and color
                       legend.key.size = unit(0.5, "cm"),  # Make the legend keys (colored boxes) smaller
                       legend.key.height = unit(0.5, "cm"),  # Adjust height of the key
                       legend.key.width = unit(0.5, "cm"), plot.title = element_text(size = 14, face = "bold", hjust = 0.5) ) +
  labs(title = "Most significant unitprot toxin hit per transcript(Bitscore >300")

ggsave(filename = file.path("alluvial4.png"), plot = alluvial4, width = 8, height = 6, dpi = 600)


#those with expression 
sum_of_expression_of_those_with_hits = sum(Distinct_Transcripts_hits$percent)
sum_of_expression_of_those_without_hits = 100-sum_of_expression_of_those_with_hits
pie_data <- data.frame(
  Category = c("With Hits", "Without Hits"),
  Value = c(sum_of_expression_of_those_with_hits, sum_of_expression_of_those_without_hits)
)

#  pie chart
pie6 <- ggplot(pie_data, aes(x = "", y = Value, fill = Category)) +
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

ggsave(filename = file.path("pie6.png"), plot = pie6, width = 8, height = 6, dpi = 600)


sum_of_expression_of_those_with_hits = sum(Distinct_Transcripts_50$percent)
sum_of_expression_of_those_without_hits = 100-sum_of_expression_of_those_with_hits
pie_data2 <- data.frame(
  Category = c("With Hits", "Without Hits"),
  Value = c(sum_of_expression_of_those_with_hits, sum_of_expression_of_those_without_hits)
)
# pie chart
pie7 <- ggplot(pie_data2, aes(x = "", y = Value, fill = Category)) +
  geom_bar(stat = "identity", width = 1, color = "black") +
  coord_polar(theta = "y", start = 0) +  # Set the starting angle for the slices
  theme_void() +  # Removes axes and background
  labs(title = "% of Expression from transcripts with uniprot toxin hits(BitScore > 50 ") +
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

ggsave(filename = file.path("pie7.png"), plot = pie7, width = 8, height = 6, dpi = 600)

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
write.csv(top10mostexpressed,file = file.path("Table14.csv"), row.names = FALSE)


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
Plot4 <- ggplot(new_df2, aes(x = "", y = total_percentage, fill = Hit)) +
  geom_bar(stat = "identity", width = 1, color = "black") +  # Add border to the slices
  coord_polar(theta = "y") +
  geom_point(data = data.frame(x = 0, y = 0), aes(x = x, y = y), size = 75, color = "white", fill = "white") +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 4, color = "#000000"),  # Smaller text size and colour
    legend.key.size = unit(0.5, "cm"),  # Make the legend keys (coloured boxes) smaller
    legend.key.height = unit(0.3, "cm"),  # Adjust height of the key
    legend.key.width = unit(0.5, "cm"), 
    legend.key.border = element_rect(color = "black", size = 1.5, linetype = "solid"),  # Bold border around legend keys
    plot.title = element_text(size = 14, face = "bold", hjust = -0.1, vjust = 1)
  ) +
  labs(title = "Relative expression of transcripts with significant uniprot toxin (BitScore >50)")
ggsave(filename = file.path("pie8.png"), plot = Plot4, width = 8, height = 6, dpi = 600)




