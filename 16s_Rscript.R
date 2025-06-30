rm(list = ls())
library(ggplot2)
library(tidyverse)
library(rstatix)
library(qiime2R)
library(dplyr)
library(PERMANOVA)
library(vegan)
library(patchwork)

#Location - /home/shrinath/Aktipis lab/new-lau/Qiime-2 Analysis (Nov-2019)/NOD Kombucha-(Nov 2019)/16s/

metadata<-read_q2metadata("/home/shrinath/Aktipis lab/new-lau/Qiime-2 Analysis (Nov-2019)/NOD Kombucha-(Nov 2019)/16s/metadata.NOD_16s.tsv")
metadata <- subset(metadata,collection_day != "0")
##########################################################
#PCoA plots
#jaccord-braycurtis combined
jaccard_16s <-read_qza("/home/shrinath/Aktipis lab/new-lau/Qiime-2 Analysis (Nov-2019)/NOD Kombucha-(Nov 2019)/16s/core-metrics-results_16s_5910/jaccard_pcoa_results.qza")$data$Vectors
braycurtis_16s <- read_qza("/home/shrinath/Aktipis lab/new-lau/Qiime-2 Analysis (Nov-2019)/NOD Kombucha-(Nov 2019)/16s/core-metrics-results_16s_5910/bray_curtis_pcoa_results.qza")$data$Vectors

jaccard_ITS <- read_qza ("/home/shrinath/Aktipis lab/new-lau/Qiime-2 Analysis (Nov-2019)/NOD Kombucha-(Nov 2019)/ITS/core-metrics-results-filtered-2/jaccard_pcoa_results.qza")$data$Vectors
braycurtis_ITS <- read_qza("/home/shrinath/Aktipis lab/new-lau/Qiime-2 Analysis (Nov-2019)/NOD Kombucha-(Nov 2019)/ITS/core-metrics-results-filtered-2/bray_curtis_pcoa_results.qza")$data$Vectors



# Prepare 16S (bacterial) data
bacteria_jaccard <- dplyr::select(jaccard_16s, SampleID, PC1, PC2)  %>%
  mutate(Method = "Jaccard", Kingdom = "Bacteria")

bacteria_braycurtis <- dplyr::select(braycurtis_16s, SampleID, PC1, PC2)  %>%
  mutate(Method = "Braycurtis", Kingdom = "Bacteria")

# Prepare ITS (fungal) data
fungi_jaccard <- dplyr::select(jaccard_ITS, SampleID, PC1, PC2)  %>%
  mutate(Method = "Jaccard", Kingdom = "Fungi")

fungi_braycurtis <- dplyr::select(braycurtis_ITS, SampleID, PC1, PC2)  %>%
  mutate(Method = "Braycurtis", Kingdom = "Fungi")

# Combine all and join metadata
combined_df <- bind_rows(
  bacteria_jaccard,
  bacteria_braycurtis,
  fungi_jaccard,
  fungi_braycurtis
) %>%
  left_join(metadata, by = "SampleID")

ggplot(combined_df, aes(x = PC1, y = PC2, shape = Type, color = Type)) +
  geom_point(size = 5, alpha = 0.9) +
  facet_grid(Kingdom ~ Method) +  # Facet by Kingdom and Method if applicable
  theme_q2r() +
  scale_color_manual(
    name = "Sample Type",
    values = c(
      "Tea-Day1" = "#1f78b4",
      "Tea-Day30" = "#a6cee3",
      "Kombucha-Day1" = "#67001f",
      "Kombucha-Day30" = "#ef8a62"
    )
  ) +
  scale_shape_manual(
    name = "Sample Type",
    values = c(
      "Tea-Day1" = 15,       # e.g., square
      "Tea-Day30" = 16,      # e.g., filled circle
      "Kombucha-Day1" = 17,    # e.g., triangle
      "Kombucha-Day30" = 18    # e.g., diamond
    )
  ) +
  ggtitle("Beta Diversity") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
    strip.text = element_text(size = 18, face = "bold"),
    text = element_text(size = 16),
    legend.position = "right",
    legend.key.size = unit(1.5, "lines")
  ) +
  guides(color = guide_legend(override.aes = list(size = 5)))



ggsave("/home/shrinath/Aktipis lab/new-lau/Qiime-2 Analysis (Nov-2019)/NOD Kombucha-(Nov 2019)/PCoA_braycurtis vs jaccard_combined.pdf", height=5, width=10, device="pdf") # save a PDF 3 inches by 4 inches



#Weighted unweighted 16s

# Read in weighted and unweighted PCoA results (update file paths as needed)
weighted_16s <- read_qza("/home/shrinath/Aktipis lab/new-lau/Qiime-2 Analysis (Nov-2019)/NOD Kombucha-(Nov 2019)/16s/core-metrics-results_16s_1343/weighted_unifrac_pcoa_results.qza")$data$Vectors
unweighted_16s <- read_qza("/home/shrinath/Aktipis lab/new-lau/Qiime-2 Analysis (Nov-2019)/NOD Kombucha-(Nov 2019)/16s/core-metrics-results_16s_1343/unweighted_unifrac_pcoa_results.qza")$data$Vectors

# Convert to data frames and select columns; add method and Kingdom
weighted_df <- as.data.frame(weighted_16s) %>%
  dplyr::select(SampleID, PC1, PC2) %>%
  mutate(Method = "Weighted", Kingdom = "Bacteria")

unweighted_df <- as.data.frame(unweighted_16s) %>%
  dplyr::select(SampleID, PC1, PC2) %>%
  mutate(Method = "Unweighted", Kingdom = "Bacteria")

# Combine the data frames
combined_df_16s <- bind_rows(weighted_df, unweighted_df)

# Join with metadata (ensure metadata has a column "SampleID" that matches)
# Replace 'metadata' with your metadata data frame variable
combined_df_16s <- left_join(combined_df_16s, metadata, by = "SampleID")

# Plot using ggplot2; facet by Method to show Weighted (left) and Unweighted (right)
ggplot(combined_df, aes(x = PC1, y = PC2, shape = Type, color = Type)) +
  geom_point(size = 5, alpha = 0.9) +
  facet_grid(. ~ Method) +  # Facet into columns by Method
  theme_q2r() +
  scale_color_manual(
    name = "Sample Type",
    values = c(
      "Tea-Day1" = "#1f78b4",
      "Tea-Day30" = "#a6cee3",
      "Kombucha-Day1" = "#67001f",
      "Kombucha-Day30" = "#ef8a62"
    )
  ) +
  scale_shape_manual(
    name = "Sample Type",
    values = c(
      "Tea-Day1" = 15,  # square
      "Tea-Day30" = 16, # circle
      "Kombucha-Day1" = 17, # triangle
      "Kombucha-Day30" = 18  # diamond
    )
  ) +
  ggtitle("Bacteria Beta Diversity UniFrac") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
    strip.text = element_text(size = 18, face = "bold"),
    text = element_text(size = 16),
    legend.position = "right",
    legend.key.size = unit(1.5, "lines")
  ) +
  guides(color = guide_legend(override.aes = list(size = 5)))

# Optionally save the plot
ggsave("/home/shrinath/Aktipis lab/new-lau/Qiime-2 Analysis (Nov-2019)/NOD Kombucha-(Nov 2019)/16s_Beta_Unifrac_Weighted_vs_Unweighted.pdf", height = 5, width = 10, device = "pdf")


#Braycurtis jaccard ITS


# Load ITS PCoA results
jaccard_ITS <- read_qza("/home/shrinath/Aktipis lab/new-lau/Qiime-2 Analysis (Nov-2019)/NOD Kombucha-(Nov 2019)/ITS/core-metrics-results-filtered-2/jaccard_pcoa_results.qza")$data$Vectors
braycurtis_ITS <- read_qza("/home/shrinath/Aktipis lab/new-lau/Qiime-2 Analysis (Nov-2019)/NOD Kombucha-(Nov 2019)/ITS/core-metrics-results-filtered-2/bray_curtis_pcoa_results.qza")$data$Vectors

# Prepare ITS (fungal) data
fungi_jaccard <- as.data.frame(jaccard_ITS) %>%
  dplyr::select(SampleID, PC1, PC2) %>%
  mutate(Method = "Jaccard", Kingdom = "Fungi")

fungi_braycurtis <- as.data.frame(braycurtis_ITS) %>%
  dplyr::select(SampleID, PC1, PC2) %>%
  mutate(Method = "Braycurtis", Kingdom = "Fungi")

# Combine fungal data
combined_ITS_df <- bind_rows(fungi_jaccard, fungi_braycurtis)

# Join with metadata (ensure that your metadata object 'metadata' has a column named SampleID)
combined_ITS_df <- left_join(combined_ITS_df, metadata, by = "SampleID")

# Create a ggplot for ITS data: facet by Method (Bray–Curtis vs. Jaccard)
p <- ggplot(combined_ITS_df, aes(x = PC1, y = PC2, shape = Type, color = Type)) +
  geom_point(size = 5, alpha = 0.9) +
  facet_grid(Kingdom ~ Method) +  # Only Fungi will show, faceted by Method
  theme_q2r() +
  scale_color_manual(
    name = "Sample Type",
    values = c(
      "Tea-Day1" = "#1f78b4",       # Medium blue
      "Tea-Day30" = "#a6cee3",      # Lighter blue
      "Kombucha-Day1" = "#67001f",   # Dark red
      "Kombucha-Day30" = "#ef8a62"   # Red-orange
    )
  ) +
  scale_shape_manual(
    name = "Sample Type",
    values = c(
      "Tea-Day1" = 15,       # Square
      "Tea-Day30" = 16,      # Circle
      "Kombucha-Day1" = 17,   # Triangle
      "Kombucha-Day30" = 18   # Diamond
    )
  ) +
  ggtitle("Fungal Beta Diversity") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
    strip.text = element_text(size = 18, face = "bold"),
    text = element_text(size = 16),
    legend.position = "right",
    legend.key.size = unit(1.5, "lines")
  ) +
  guides(color = guide_legend(override.aes = list(size = 5)))

# Save the plot as a PDF; adjust the path as needed.
ggsave("/home/shrinath/Aktipis lab/new-lau/Qiime-2 Analysis (Nov-2019)/NOD Kombucha-(Nov 2019)/PCoA_ITS_Braycurtis_vs_Jaccard_combined.pdf", plot = p, height = 5, width = 10, device = "pdf")


# Combination of Weighted-Unweighted and Jaccard-Braycutris for entire dataset


# Example: Assume you already have combined_df_16s and combined_df_ITS prepared.
# Here, combined_df_16s should have rows for weighted/unweighted UniFrac of 16S with Kingdom = "Bacteria"
# and combined_df_ITS should have rows for Bray–Curtis/Jaccard of ITS with Kingdom = "Fungi".

# Define your custom palettes (update these if needed)
custom_colors <- c(
  "Tea-Day1" = "#1f78b4",       # Medium blue
  "Tea-Day30" = "#a6cee3",      # Lighter blue
  "Kombucha-Day1" = "#67001f",   # Dark red
  "Kombucha-Day30" = "#ef8a62"   # Red-orange
)
custom_shapes <- c(
  "Tea-Day1" = 15,       # Square
  "Tea-Day30" = 16,      # Circle
  "Kombucha-Day1" = 17,   # Triangle
  "Kombucha-Day30" = 18   # Diamond
)

combined_df_16s$Method <- factor(combined_df_16s$Method ,levels = c("Weighted","Unweighted"))
combined_ITS_df$Method <- factor(combined_ITS_df$Method,levels = c("Braycurtis","Jaccard"))

## 16S Plot: Bacterial Beta Diversity (Weighted vs. Unweighted UniFrac)
p16s <- ggplot(combined_df_16s, aes(x = PC1, y = PC2, color = Type, shape = Type)) +
  geom_point(size = 5, alpha = 0.9) +
  facet_grid(Kingdom ~ Method) +
  theme_q2r() +
  scale_color_manual(name = "Sample Type", values = custom_colors) +
  scale_shape_manual(name = "Sample Type", values = custom_shapes) +
  ggtitle("") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
    strip.text = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  )

pITS <- ggplot(combined_ITS_df, aes(x = PC1, y = PC2, color = Type, shape = Type)) +
  geom_point(size = 5, alpha = 0.9) +
  facet_grid(Kingdom ~ Method, switch = "x") +
  theme_q2r() +
  scale_color_manual(name = "Sample Type", values = custom_colors) +
  scale_shape_manual(name = "Sample Type", values = custom_shapes) +
  ggtitle("") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
    strip.text = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  )

# Combine the two plots vertically with a common title and shared legend
combined_plot <- (p16s / pITS) +
  plot_annotation(title = "Beta Diversity",
                  theme = theme(plot.title = element_text(hjust = 0.35))) &
  theme(legend.position = "right",
        legend.text = element_text(size = 12),      # Increase legend text size
        legend.title = element_text(size = 12))

# Collect the legends
combined_plot <- combined_plot + plot_layout(guides = "collect")

# Print the combined plot
print(combined_plot)

# Save the final plot as a PDF (adjust the path as needed)
ggsave("/home/shrinath/Aktipis lab/new-lau/Qiime-2 Analysis (Nov-2019)/NOD Kombucha-(Nov 2019)/Combined_Beta_Diversity.pdf", 
       plot = combined_plot, height = 8, width = 9, device = "pdf")




###########################################################################################################################
#Beta diversity statistical tests - 16s

#Unweighted tests

unweighted_dist_16s <- read_qza("/home/shrinath/Aktipis lab/new-lau/Qiime-2 Analysis (Nov-2019)/NOD Kombucha-(Nov 2019)/16s/core-metrics-results_16s_5910/unweighted_unifrac_distance_matrix.qza")$data
weighted_dist_16s <- read_qza("/home/shrinath/Aktipis lab/new-lau/Qiime-2 Analysis (Nov-2019)/NOD Kombucha-(Nov 2019)/16s/core-metrics-results_16s_5910/weighted_unifrac_distance_matrix.qza")$data


# Extract sample names from the distance matrix
sample_names <- rownames(as.matrix(unweighted_dist_16s))

# Filter and reorder metadata
metadata_matched <- metadata %>%
  filter(SampleID %in% sample_names) %>%
  arrange(match(SampleID, sample_names))

# Set rownames
rownames(metadata_matched) <- metadata_matched$SampleID

all(rownames(metadata_matched) == sample_names)  # should return TRUE

# Filter metadata
kombucha_meta <- metadata_matched %>%
  filter( isol_growth_condt == "Kombucha" & collection_day %in% c(1, 30))

# Subset distance matrix
kombucha_samples <- kombucha_meta$SampleID
kombucha_dist <- as.matrix(unweighted_dist_16s)[kombucha_samples, kombucha_samples]

# Convert to dist object
kombucha_dist <- as.dist(kombucha_dist)

# Ensure rownames match
rownames(kombucha_meta) <- kombucha_meta$SampleID

# Run PERMANOVA
adonis_kombucha <- adonis2(kombucha_dist ~ collection_day, data = kombucha_meta, permutations = 999)
print(adonis_kombucha)


# Convert to matrix (important!)
dist_matrix <- as.matrix(unweighted_dist_16s)

# Get tea sample IDs
tea_samples <- metadata_matched %>%
  filter(isol_growth_condt == "Tea", collection_day %in% c(1, 30)) %>%
  pull(SampleID)

# Keep only samples that exist in the distance matrix
tea_samples <- intersect(tea_samples, rownames(dist_matrix))

# Subset metadata
tea_meta <- metadata_matched %>%
  filter(SampleID %in% tea_samples) %>%
  arrange(match(SampleID, tea_samples))
rownames(tea_meta) <- tea_meta$SampleID

# Subset the matrix to tea samples
tea_dist <- as.dist(dist_matrix[tea_samples, tea_samples])

# Run PERMANOVA
adonis_tea <- adonis2(tea_dist ~ collection_day, data = tea_meta, permutations = 999)
print(adonis_tea)



sink("/home/shrinath/Aktipis lab/new-lau/Qiime-2 Analysis (Nov-2019)/NOD Kombucha-(Nov 2019)/permanova_unweighted_16s_results_tea_kombucha_day1vsday30.txt")
cat("Tea - Unweighted UniFrac Day1 vs Day30 PERMANOVA Results\n")
print(adonis_tea)
cat("\n-----------------------------\n\n")

cat("Kombucha - Unweighted UniFrac Day1 vs Day30 PERMANOVA Results\n")
print(adonis_kombucha)

sink()

------------------------------------------------------------------------------------------------------------------------------
#Weighted unifrac significance 16s  beta

weighted_dist_16s <- read_qza("/home/shrinath/Aktipis lab/new-lau/Qiime-2 Analysis (Nov-2019)/NOD Kombucha-(Nov 2019)/16s/core-metrics-results_16s_5910/weighted_unifrac_distance_matrix.qza")$data

metadata<-read_q2metadata("/home/shrinath/Aktipis lab/new-lau/Qiime-2 Analysis (Nov-2019)/NOD Kombucha-(Nov 2019)/16s/metadata.NOD_16s.tsv")


weighted_matrix <- as.matrix(weighted_dist_16s)


# Tea sample IDs
tea_samples_w <- metadata %>%
  filter(isol_growth_condt == "Tea", collection_day %in% c(1, 30)) %>%
  pull(SampleID)

# Ensure samples are present
tea_samples_w <- intersect(tea_samples_w, rownames(weighted_matrix))

# Subset metadata
tea_meta_w <- metadata %>%
  filter(SampleID %in% tea_samples_w) %>%
  arrange(match(SampleID, tea_samples_w))
rownames(tea_meta_w) <- tea_meta_w$SampleID

# Subset distance matrix
tea_dist_weighted <- as.dist(weighted_matrix[tea_samples_w, tea_samples_w])

# Run PERMANOVA
adonis_tea_weighted <- adonis2(tea_dist_weighted ~ collection_day, data = tea_meta_w, permutations = 999)
print(adonis_tea_weighted)



# Kombucha sample IDs
kombucha_samples_w <- metadata %>%
  filter(isol_growth_condt == "Kombucha", collection_day %in% c(1, 30)) %>%
  pull(SampleID)

# Ensure samples are present
kombucha_samples_w <- intersect(kombucha_samples_w, rownames(weighted_matrix))

# Subset metadata
kombucha_meta_w <- metadata %>%
  filter(SampleID %in% kombucha_samples_w) %>%
  arrange(match(SampleID, kombucha_samples_w))
rownames(kombucha_meta_w) <- kombucha_meta_w$SampleID

# Subset distance matrix
kombucha_dist_weighted <- as.dist(weighted_matrix[kombucha_samples_w, kombucha_samples_w])

# Run PERMANOVA
adonis_kombucha_weighted <- adonis2(kombucha_dist_weighted ~ collection_day, data = kombucha_meta_w, permutations = 999)
print(adonis_kombucha_weighted)

sink("/home/shrinath/Aktipis lab/new-lau/Qiime-2 Analysis (Nov-2019)/NOD Kombucha-(Nov 2019)/permanova_weighted_16s_results_tea_kombucha_day1vsday30.txt")
cat("Tea - Weighted UniFrac Day1 vs Day30 PERMANOVA Results\n")
print(adonis_tea)
cat("\n-----------------------------\n\n")

cat("Kombucha - Weighted UniFrac Day1 vs Day30 PERMANOVA Results\n")
print(adonis_kombucha)

sink()
###########################################################################################################################################
# Load required packages
library(vegan)
library(dplyr)


# Read in the weighted UniFrac distance matrix (exported as TSV from QIIME 2)
# The file should have sample IDs as row names and column names.
weighted_df <- read_qza("/home/shrinath/Aktipis lab/new-lau/Qiime-2 Analysis (Nov-2019)/NOD Kombucha-(Nov 2019)/16s/core-metrics-results_16s_5910/weighted_unifrac_distance_matrix.qza")$data
weighted_dist <- as.dist(weighted_df)

# Read in the unweighted UniFrac distance matrix
unweighted_df <- read_qza("/home/shrinath/Aktipis lab/new-lau/Qiime-2 Analysis (Nov-2019)/NOD Kombucha-(Nov 2019)/16s/core-metrics-results_16s_5910/unweighted_unifrac_distance_matrix.qza")$data
unweighted_dist <- as.dist(unweighted_df)


missing_sample <- setdiff(metadata$SampleID,rownames(as.matrix(weighted_dist)))
metadata <-subset(metadata,Description != missing_sample) # lost at depth 5910

# Subset metadata for Tea and Kombucha
metadata_Tea_Day1 <- metadata %>% filter(isol_growth_condt == "Tea",collection_day =="1")
metadata_Kombucha_Day1 <- metadata %>% filter(isol_growth_condt == "Kombucha",collection_day =="1")
metadata_Tea_Day30 <- metadata %>% filter(isol_growth_condt == "Tea",collection_day =="30")
metadata_Kombucha_Day30 <- metadata %>% filter(isol_growth_condt == "Kombucha",collection_day =="30")

# Subset the distance matrices to the samples present in each subset
all_samples_weighted <- rownames(as.matrix(weighted_dist))
all_samples_unweighted <- rownames(as.matrix(unweighted_dist))

# For Tea:
common_samples_Tea_Day1 <- intersect(all_samples_weighted, metadata_Tea_Day1$SampleID)
common_samples_Tea_Day30 <- intersect(all_samples_weighted, metadata_Tea_Day30$SampleID)

weighted_Tea_day1 <- as.dist(as.matrix(weighted_dist)[common_samples_Tea_Day1, common_samples_Tea_Day1])
weighted_Tea_day30 <- as.dist(as.matrix(weighted_dist)[common_samples_Tea_Day30, common_samples_Tea_Day30])

unweighted_Tea_day1 <- as.dist(as.matrix(unweighted_dist)[common_samples_Tea_Day1, common_samples_Tea_Day1])
unweighted_Tea_day30 <- as.dist(as.matrix(weighted_dist)[common_samples_Tea_Day30, common_samples_Tea_Day30])

# For Kombucha:

common_samples_Kombucha_Day1 <- intersect(all_samples_weighted, metadata_Kombucha_Day1$SampleID)
common_samples_Kombucha_Day30 <- intersect(all_samples_weighted, metadata_Kombucha_Day30$SampleID)

weighted_Kombucha_day1 <- as.dist(as.matrix(weighted_dist)[common_samples_Kombucha_Day1, common_samples_Kombucha_Day1])
weighted_Kombucha_day30 <- as.dist(as.matrix(weighted_dist)[common_samples_Kombucha_Day30, common_samples_Kombucha_Day30])

unweighted_Kombucha_day1 <- as.dist(as.matrix(unweighted_dist)[common_samples_Kombucha_Day1, common_samples_Kombucha_Day1])
unweighted_Kombucha_day30 <- as.dist(as.matrix(weighted_dist)[common_samples_Kombucha_Day30, common_samples_Kombucha_Day30])





# Run PERMANOVA for weighted UniFrac:
adonis_weighted_Tea <- adonis2(weighted_Tea ~ collection_day, data = metadata_Tea, permutations = 999)
adonis_weighted_Kombucha <- adonis2(weighted_Kombucha ~ collection_day, data = metadata_Kombucha, permutations = 999)

# Run PERMANOVA for unweighted UniFrac:
adonis_unweighted_Tea <- adonis2(unweighted_Tea ~ collection_day, data = metadata_Tea, permutations = 999)
adonis_unweighted_Kombucha <- adonis2(unweighted_Kombucha ~ collection_day, data = metadata_Kombucha, permutations = 999)

# Extract overall p-values from the adonis2 output (the first row corresponds to the model term)
pvals <- c(
  weighted_Tea = adonis_weighted_Tea$`Pr(>F)`[1],
  weighted_Kombucha = adonis_weighted_Kombucha$`Pr(>F)`[1],
  unweighted_Tea = adonis_unweighted_Tea$`Pr(>F)`[1],
  unweighted_Kombucha = adonis_unweighted_Kombucha$`Pr(>F)`[1]
)

# Run PERMANOVA for weighted UniFrac:
adonis_weighted <- adonis2(weighted_dist ~ Type, data = metadata, permutations = 999)

# Run PERMANOVA for unweighted UniFrac:
adonis_unweighted <- adonis2(unweighted_dist ~ Type, data = metadata, permutations = 999)

# Extract overall p-values from the adonis2 output (the first row corresponds to the model term)
pvals <- c(
  weighted = adonis_weighted$`Pr(>F)`[1],
  unweighted = adonis_unweighted$`Pr(>F)`[1]
)


# Apply Bonferroni adjustment
adj_pvals <- p.adjust(pvals, method = "bonferroni")

# Output the raw and Bonferroni-adjusted p-values
print("Raw PERMANOVA p-values:")
print(pvals)
print("Bonferroni-adjusted p-values:")
print(adj_pvals)





##################################################################################################################################
#Alpha diversity statistical tests - Shannon-1343
sink(file = "/home/shrinath/Aktipis lab/new-lau/Qiime-2 Analysis (Nov-2019)/NOD Kombucha-(Nov 2019)/Alpha diversity_16s .txt")
print("Alpha diversity statistical tests - Shannon-1343")
print("/n")
shannon<-read_qza("/home/shrinath/Aktipis lab/new-lau/Qiime-2 Analysis (Nov-2019)/NOD Kombucha-(Nov 2019)/16s/core-metrics-results_16s_1343/shannon_vector.qza")$data
shannon$Type <- metadata$Type
class(shannon$shannon_entropy) = "Numeric"
storage.mode(shannon$shannon_entropy) = "double"



res.kruskal_shannon <- kruskal_test(shannon,shannon_entropy ~ Type)

res.kruskal_shannon

pairwise.wilcox.test(shannon$shannon_entropy,shannon$Type,p.adjust.method = "fdr")
dunn_test(shannon, shannon_entropy ~ Type, p.adjust.method = "fdr")


##############
#Alpha diversity statistical tests - Observed features
print("Observed features")
ASVs<-read_qza("/home/shrinath/Aktipis lab/new-lau/Qiime-2 Analysis (Nov-2019)/NOD Kombucha-(Nov 2019)/16s/core-metrics-results_16s_1343/observed_features_vector.qza")$data %>% rownames_to_column("SampleID") 

ASVs_1 <- merge(ASVs,metadata,by = "SampleID")

class(ASVs_1$observed_features) = "Numeric"
storage.mode(ASVs_1$observed_features) = "double"

res.kruskal_ASVs <- kruskal_test(ASVs_1,observed_features ~ Type)

res.kruskal_ASVs
pairwise.wilcox.test(ASVs_1$observed_features,ASVs_1$Type,p.adjust.method = "fdr")
dunn_test(ASVs_1,observed_features ~ Type, p.adjust.method = "fdr")





#Alpha diversity statistical tests - Shannon-5910
print("Alpha diversity statistical tests - Shannon-5910")
print("/n")
shannon<-read_qza("/home/shrinath/Aktipis lab/new-lau/Qiime-2 Analysis (Nov-2019)/NOD Kombucha-(Nov 2019)/16s/core-metrics-results_16s_5910/shannon_vector.qza")$data
shannon$SampleID <- rownames(shannon)
shannon <- left_join(shannon, metadata %>% dplyr:: select(SampleID, Type), by = "SampleID")


class(shannon$shannon_entropy) = "Numeric"
storage.mode(shannon$shannon_entropy) = "double"



res.kruskal_shannon <- kruskal_test(shannon,shannon_entropy ~ Type)

res.kruskal_shannon
pairwise.wilcox.test(shannon$shannon_entropy,shannon$Type,p.adjust.method = "fdr")
dunn_test(shannon, shannon_entropy ~ Type, p.adjust.method = "fdr")



##############
#Alpha diversity statistical tests - Observed features
print("Observed features")
ASVs<-read_qza("/home/shrinath/Aktipis lab/new-lau/Qiime-2 Analysis (Nov-2019)/NOD Kombucha-(Nov 2019)/16s/core-metrics-results_16s_5910/observed_features_vector.qza")$data %>% rownames_to_column("SampleID") 

ASVs_1 <- merge(ASVs,metadata,by = "SampleID")

class(ASVs_1$observed_features) = "Numeric"
storage.mode(ASVs_1$observed_features) = "double"

res.kruskal_ASVs <- kruskal_test(ASVs_1,observed_features ~ Type)
summary(res.kruskal_ASVs)  
res.kruskal_ASVs
pairwise.wilcox.test(ASVs_1$observed_features,ASVs_1$Type,p.adjust.method = "fdr")
dunn_test(ASVs_1,observed_features ~ Type, p.adjust.method = "fdr")


#Alpha diversity statistical tests - Shannon-15204
print("Alpha diversity statistical tests - Shannon-15204")
print("/n")
shannon<-read_qza("/home/shrinath/Aktipis lab/new-lau/Qiime-2 Analysis (Nov-2019)/NOD Kombucha-(Nov 2019)/16s/core-metrics-results_16s_15204/shannon_vector.qza")$data
shannon$SampleID <- rownames(shannon)
shannon <- left_join(shannon, metadata %>% dplyr:: select(SampleID, Type), by = "SampleID")


class(shannon$shannon_entropy) = "Numeric"
storage.mode(shannon$shannon_entropy) = "double"



res.kruskal_shannon <- kruskal_test(shannon,shannon_entropy ~ Type)
res.kruskal_shannon
pairwise.wilcox.test(shannon$shannon_entropy,shannon$Type,p.adjust.method = "fdr")
dunn_test(shannon, shannon_entropy ~ Type, p.adjust.method = "fdr")



##############
#Alpha diversity statistical tests - Observed features
print("Observed features")
ASVs<-read_qza("/home/shrinath/Aktipis lab/new-lau/Qiime-2 Analysis (Nov-2019)/NOD Kombucha-(Nov 2019)/16s/core-metrics-results_16s_15204/observed_features_vector.qza")$data %>% rownames_to_column("SampleID") 

ASVs_1 <- merge(ASVs,metadata,by = "SampleID")

class(ASVs_1$observed_features) = "Numeric"
storage.mode(ASVs_1$observed_features) = "double"

res.kruskal_ASVs <- kruskal_test(ASVs_1,observed_features ~ Type)

res.kruskal_ASVs
pairwise.wilcox.test(shannon$shannon_entropy,shannon$Type,p.adjust.method = "fdr")
dunn_test(shannon, shannon_entropy ~ Type, p.adjust.method = "fdr")

sink()
########################################################################

#Alpha diversity statistical tests - Shannon-1221- ITS
sink(file = "/home/shrinath/Aktipis lab/new-lau/Qiime-2 Analysis (Nov-2019)/NOD Kombucha-(Nov 2019)/Alpha diversity_ITS .txt")
print("Alpha diversity statistical tests - Shannon-1221")
print("/n")
shannon<-read_qza("/home/shrinath/Aktipis lab/new-lau/Qiime-2 Analysis (Nov-2019)/NOD Kombucha-(Nov 2019)/ITS/ITS_core_metrics_1221/shannon_vector.qza")$data
shannon$Type <- metadata$Type
class(shannon$shannon_entropy) = "Numeric"
storage.mode(shannon$shannon_entropy) = "double"



res.kruskal_shannon <- kruskal_test(shannon,shannon_entropy ~ Type)

res.kruskal_shannon

pairwise.wilcox.test(ASVs_1$observed_features,ASVs_1$Type,p.adjust.method = "fdr")
dunn_test(ASVs_1,observed_features ~ Type, p.adjust.method = "fdr")


##############
#Alpha diversity statistical tests - Observed features
print("Observed features")
ASVs<-read_qza("/home/shrinath/Aktipis lab/new-lau/Qiime-2 Analysis (Nov-2019)/NOD Kombucha-(Nov 2019)/ITS/ITS_core_metrics_1221/observed_features_vector.qza")$data %>% rownames_to_column("SampleID") 

ASVs_1 <- merge(ASVs,metadata,by = "SampleID")

class(ASVs_1$observed_features) = "Numeric"
storage.mode(ASVs_1$observed_features) = "double"

res.kruskal_ASVs <- kruskal_test(ASVs_1,observed_features ~ Type)

res.kruskal_ASVs
pairwise.wilcox.test(ASVs_1$observed_features,ASVs_1$Type,p.adjust.method = "fdr")
dunn_test(ASVs_1,observed_features ~ Type, p.adjust.method = "fdr")



sink()




#########################################################################################################################################################################################





# Get sample names from the distance matrix
dist_samples <- rownames(as.matrix(unweighted_dist_16s))

# Find missing samples
missing_tea_samples <- setdiff(tea_samples, dist_samples)
print(missing_tea_samples)






