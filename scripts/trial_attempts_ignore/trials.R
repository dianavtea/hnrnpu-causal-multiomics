######transcriptomics_alltimepoints_prep.R#####

#ASD vs CTRL at D5 relative to D0
#resD5 <- results(dds, name = "DayD5.conditionHNRNPUdel")
#ASD vs CTRL at D28 relative to D0
#resD28 <- results(dds, name = "DayD28.conditionHNRNPUdel")

###### VST ######
vst_all <- vst(dds, blind = FALSE)
vstD28only <- vst_all[, colData(dds)$Day == "D28"]
vstmatD28only <- assay(vstD28only) #127572 elements
dim(vstmatD28only) #21262 transcripts x 6 samples

plotPCA(vstD28only, intgroup = "condition")
#82% variance explained in PC1, 10% in PC2, distinction between ASD and CTRL groups

library(ggplot2)
# PCA plot
pca_D28only <- plotPCA(vstD28only, intgroup = "condition") +
  theme_bw() +
  theme(
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10)
  ) +
  labs(title = "PCA of VST-transformed counts")

# colours 
# "#F8766D" orange for CTRL 1
# "#00BFC4" turq for HNRNPUdel

ggsave(
  filename = "PCA_D28only.png",      # use PDF or PNG
  plot = pca_D28only,           
  width = 6, height = 5,          # size in inches
  dpi = 300                       # resolution for PNG (ignored for PDF)
)

#clean vstmatD28only
# Subset the matrix to keep ONLY the Ensembl IDs that survived the "Best ID" filter above
# This ensures your matrix has exactly the same genes as your stats table
vst_D28only_clean <- vstmatD28only[ rownames(vstmatD28only) %in% res_D28only_clean$ensembl, ]

# Now we need to swap the rownames from Ensembl -> Symbol.
# We must match the order perfectly.
matcher <- match(rownames(vst_D28only_clean), res_D28only_clean$ensembl)
rownames(vst_D28only_clean) <- res_D28only_clean$hgnc_symbol[matcher]

any(rownames(vst_D28only_clean) == "")       # should be FALSE
any(duplicated(rownames(vst_D28only_clean))) # Must be FALSE
dim(vst_D28only_clean)                       # Should match rows in res_clean
head(vst_D28only_clean) #102102
# Save for PROGENy
write.csv(vst_D28only_clean, "~/hnrnpu-causal-multiomics/processeddata/vst_progD28only_clean.csv", row.names = TRUE)
vst_D28only_clean <- read.csv("~/hnrnpu-causal-multiomics/processeddata/vst_progD28only_clean.csv")
vstD28only_annot <- cbind(res_D28only_clean[, c("ensembl", "entrez", "hgnc_symbol")], vst_D28only_clean)
head(vstD28only_annot)
write.csv(vstD28only_annot, "VSTcounts_annot_D28only.csv", row.names = FALSE)
vstD28only_annot <- read.csv("VSTcounts_annot_D28only.csv")

#ADDITION 28.12.25
# Remove rows where symbol is NA or blank
#res_D5only_clean <- res_D5_only.annot[ !is.na(res_D5_only.annot$hgnc_symbol) & res_D5_only.annot$hgnc_symbol != "", ]
#17017
# Sort by baseMean (highest expression first)
#res_D5only_clean <- res_D5only_clean[order(res_D5only_clean$baseMean, decreasing = TRUE), ]
# Remove duplicates
# !duplicated() keeps the FIRST instance (which is now the highest expressed one)
#res_D5only_clean <- res_D5only_clean[!duplicated(res_D5only_clean$hgnc_symbol), ]

###### VST ######
vst_all <- vst(dds, blind = FALSE)
vstD5only <- vst_all[, colData(dds)$Day == "D5"]
vstmatD5only <- assay(vstD5only) #127572 elements
dim(vstmatD5only) #21262 transcripts x 6 samples

plotPCA(vstD5only, intgroup = "condition")
#74% variance explained in PC1, 12% in PC2, distinction between ASD and CTRL groups

library(ggplot2)
# PCA plot
pca_D5only <- plotPCA(vstD5only, intgroup = "condition") +
  theme_bw() +
  theme(
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10)
  ) +
  labs(title = "PCA of VST-transformed counts")

# colours 
# "#F8766D" orange for CTRL 1
# "#00BFC4" turq for HNRNPUdel

ggsave(
  filename = "PCA_D5only.png",      # use PDF or PNG
  plot = pca_D5only,           
  width = 6, height = 5,          # size in inches
  dpi = 300                       # resolution for PNG (ignored for PDF)
)

#clean vstmatD5only
# Subset the matrix to keep ONLY the Ensembl IDs that survived the "Best ID" filter above
# This ensures your matrix has exactly the same genes as your stats table
vst_D5only_clean <- vstmatD5only[ rownames(vstmatD5only) %in% res_D5only_clean$ensembl, ]

# Now we need to swap the rownames from Ensembl -> Symbol.
# We must match the order perfectly.
matcher <- match(rownames(vst_D5only_clean), res_D5only_clean$ensembl)
rownames(vst_D5only_clean) <- res_D5only_clean$hgnc_symbol[matcher]

any(rownames(vst_D5only_clean) == "")       # should be FALSE
any(duplicated(rownames(vst_D5only_clean))) # Must be FALSE
dim(vst_D5only_clean)                       # Should match rows in res_clean
head(vst_D5only_clean) #102102
# Save for PROGENy
write.csv(vst_D5only_clean, "~/hnrnpu-causal-multiomics/processeddata/vst_progD5only_clean.csv", row.names = TRUE)
vst_D5only_clean <- read.csv("~/hnrnpu-causal-multiomics/processeddata/vst_progD5only_clean.csv")

vstD5only_annot <- cbind(res_D5only_clean[, c("ensembl", "entrez", "hgnc_symbol")], vst_D5only_clean)
head(vstD5only_annot)
write.csv(vstD5only_annot, "VSTcounts_annot_D5only.csv", row.names = FALSE)
vstD5only_annot <- read.csv("VSTcounts_annot_D5only.csv")

###### VST ######
vst_all <- vst(dds, blind = FALSE)
vstD0 <- vst_all[, colData(dds)$Day == "D0"]
vstmatD0 <- assay(vstD0) #21620 elements
dim(vstmatD0) #21262 transcripts x 10 samples

plotPCA(vstD0, intgroup = "condition")
#73% variance explained in PC1, 18% in PC2, distinction between ASD and CTRL groups

library(ggplot2)
# PCA plot
pca_D0 <- plotPCA(vstD0 , intgroup = "condition") +
  theme_bw() +
  theme(
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10)
  ) +
  labs(title = "PCA of VST-transformed counts")

# colours 
# "#F8766D" orange for CTRL 1
# "#00BFC4" turq for HNRNPUdel

ggsave(
  filename = "PCA_D0.png",      # use PDF or PNG
  plot = pca_D0,           
  width = 6, height = 5,          # size in inches
  dpi = 300                       # resolution for PNG (ignored for PDF)
)

#clean vstmatD0
# Subset the matrix to keep ONLY the Ensembl IDs that survived the "Best ID" filter above
# This ensures your matrix has exactly the same genes as your stats table
vst_D0_clean <- vstmatD0[ rownames(vstmatD0) %in% resD0_clean$ensembl, ]

# Now we need to swap the rownames from Ensembl -> Symbol.
# We must match the order perfectly.
matcher <- match(rownames(vst_D0_clean), resD0_clean$ensembl)
rownames(vst_D0_clean) <- resD0_clean$hgnc_symbol[matcher]

any(rownames(vst_D0_clean) == "")       # should be FALSE
any(duplicated(rownames(vst_D0_clean))) # Must be FALSE
dim(vst_D0_clean)                       # Should match rows in res_clean
head(vst_D0_clean) #102102
# Save for PROGENy
write.csv(vst_D0_clean, "~/hnrnpu-causal-multiomics/processeddata/vst_progD0_clean.csv", row.names = TRUE)
vst_D0_clean <- read.csv("~/hnrnpu-causal-multiomics/processeddata/vst_progD0_clean.csv")

vstD0_annot <- cbind(resD0_clean[, c("ensembl", "entrez", "hgnc_symbol")], vst_D0_clean)
head(vstD0_annot)
write.csv(vstD0_annot, "VSTcounts_annot_D0.csv", row.names = FALSE)
vstD0_annot <- read.csv("VSTcounts_annot_D0.csv")

#####pathwayactivities.R#####
#top 100 genes
PathwayActivity_countsD0only <- progeny(nc_D0only_mat, scale = TRUE,
                                        organism = "Human", top = 100)

#define heatmap colour breaks
Activity_countsD0only <- as.vector(PathwayActivity_countsD0only) 

paletteLength <- 100
myColor <- 
  colorRampPalette(c("darkblue", "whitesmoke","indianred"))(paletteLength)

progenyBreaksD0only <- c(seq(min(Activity_countsD0only), 0, 
                             length.out=ceiling(paletteLength/2) + 1),
                         seq(max(Activity_countsD0only)/paletteLength, 
                             max(Activity_countsD0only), 
                             length.out=floor(paletteLength/2)))

annotation_col <- metafiltD0 %>%
  dplyr::select(Sample_ID, condition) %>%
  as.data.frame()

rownames(annotation_col) <- annotation_col$Sample_ID
annotation_col$Sample_ID <- NULL

# force condition order
annotation_col$condition <- factor(
  annotation_col$condition,
  levels = c("CTRL", "HNRNPUdel")
)

# order samples
sample_order <- rownames(annotation_col)[order(annotation_col$condition)]

# reorder matrix
matD0only <- t(PathwayActivity_countsD0only)[, sample_order, drop = FALSE]

# reorder annotation (IMPORTANT: keep as data.frame)
annotation_col2 <- annotation_col[sample_order, , drop = FALSE]

# plot
progeny_hmapD0only <- pheatmap(
  matD0only,
  fontsize = 14,
  fontsize_row = 10,
  fontsize_col = 10,
  color = myColor,
  breaks = progenyBreaksD0only,
  main = "PROGENy (100)",
  angle_col = 45,
  treeheight_col = 0,
  border_color = NA,
  annotation_col = annotation_col2,
  cluster_cols = FALSE,
  gaps_col = sum(annotation_col2$condition == "CTRL")
)


ggsave(
  filename = "PROGENY_heatmapD0only.png",
  plot = progeny_hmapD0only,
  width = 8,
  height = 12,
  dpi = 300
)

#enrichment analysis using competitive permutation approach to assess
#significance of pathway activity to end with NES for each pathway
PathwayActivity_zscoreD0only <- progeny(DED0only_mat, 
                                        scale=TRUE, organism="Human", top = 100, perm = 10000, z_scores = TRUE) %>%
  t()
colnames(PathwayActivity_zscoreD0only) <- "NES"

PathwayActivity_zscoreD0only_df <- as.data.frame(PathwayActivity_zscoreD0only) %>% 
  rownames_to_column(var = "Pathway") %>%
  dplyr::arrange(NES) %>%
  dplyr::mutate(Pathway = factor(Pathway))

NES_pathwayD0only <- ggplot(PathwayActivity_zscoreD0only_df,aes(x = reorder(Pathway, NES), y = NES)) + 
  geom_bar(aes(fill = NES), stat = "identity") +
  scale_fill_gradient2(low = "darkblue", high = "indianred", 
                       mid = "whitesmoke", midpoint = 0) + 
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = 
          element_text(angle = 45, hjust = 1, size =10, face= "bold"),
        axis.text.y = element_text(size =10, face= "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("Pathways")

ggsave(
  filename = "PROGENY_NES_pathwaysD0only.png",
  plot = NES_pathwayD0only,
  width = 10,
  height = 8,
  dpi = 300
)

# pathway most active - visualise most responsive genes (progeny weights)
#along with their statistic values to interpret results
#scatterplot shows which genes contributing most to the pathway enrichment

prog_matrixD0only <- getModel("Human", top=100) %>% 
  as.data.frame()  %>%
  tibble::rownames_to_column("GeneID")

DED0only_df <- DED0only_mat %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("GeneID")

sum(is.na(DED0only_df$statistic))
DED0only_df <- DED0only_df %>%
  dplyr::filter(!is.na(statistic))

scat_plotsD0 <- progeny::progenyScatter(df = DED0only_df, 
                                        weight_matrix = prog_matrixD0only, 
                                        statName = "statistic", verbose = FALSE)
plot(scat_plotsD0[[1]]$`TGFb`) 

#top 100 genes
PathwayActivity_countsD28only <- progeny(nc_D28only_mat, scale = TRUE,
                                         organism = "Human", top = 100)

#define heatmap colour breaks
Activity_countsD28only <- as.vector(PathwayActivity_countsD28only) 

paletteLength <- 100
myColor <- 
  colorRampPalette(c("darkblue", "whitesmoke","indianred"))(paletteLength)

progenyBreaksD28only <- c(seq(min(Activity_countsD28only), 0, 
                              length.out=ceiling(paletteLength/2) + 1),
                          seq(max(Activity_countsD28only)/paletteLength, 
                              max(Activity_countsD28only), 
                              length.out=floor(paletteLength/2)))

annotation_col <- metafiltD28 %>%
  dplyr::select(Sample_ID, condition) %>%
  as.data.frame()

rownames(annotation_col) <- annotation_col$Sample_ID
annotation_col$Sample_ID <- NULL

# force condition order
annotation_col$condition <- factor(
  annotation_col$condition,
  levels = c("CTRL", "HNRNPUdel")
)

# order samples
sample_order <- rownames(annotation_col)[order(annotation_col$condition)]

# reorder matrix
matD28only <- t(PathwayActivity_countsD28only)[, sample_order, drop = FALSE]

# reorder annotation (IMPORTANT: keep as data.frame)
annotation_col2 <- annotation_col[sample_order, , drop = FALSE]

# plot
progeny_hmapD28only <- pheatmap(
  matD28only,
  fontsize = 14,
  fontsize_row = 10,
  fontsize_col = 10,
  color = myColor,
  breaks = progenyBreaksD28only,
  main = "PROGENy (100)",
  angle_col = 45,
  treeheight_col = 0,
  border_color = NA,
  annotation_col = annotation_col2,
  cluster_cols = FALSE,
  gaps_col = sum(annotation_col2$condition == "CTRL")
)


ggsave(
  filename = "PROGENY_heatmapD28only.png",
  plot = progeny_hmapD28only,
  width = 8,
  height = 12,
  dpi = 300
)

#enrichment analysis using competitive permutation approach to assess
#significance of pathway activity to end with NES for each pathway
PathwayActivity_zscoreD28only <- progeny(DED28only_mat, 
                                         scale=TRUE, organism="Human", top = 100, perm = 10000, z_scores = TRUE) %>%
  t()
colnames(PathwayActivity_zscoreD28only) <- "NES"

PathwayActivity_zscoreD28only_df <- as.data.frame(PathwayActivity_zscoreD28only) %>% 
  rownames_to_column(var = "Pathway") %>%
  dplyr::arrange(NES) %>%
  dplyr::mutate(Pathway = factor(Pathway))

NES_pathwayD28only <- ggplot(PathwayActivity_zscoreD28only_df,aes(x = reorder(Pathway, NES), y = NES)) + 
  geom_bar(aes(fill = NES), stat = "identity") +
  scale_fill_gradient2(low = "darkblue", high = "indianred", 
                       mid = "whitesmoke", midpoint = 0) + 
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = 
          element_text(angle = 45, hjust = 1, size =10, face= "bold"),
        axis.text.y = element_text(size =10, face= "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("Pathways")

ggsave(
  filename = "PROGENY_NES_pathwaysD28only.png",
  plot = NES_pathwayD28only,
  width = 10,
  height = 8,
  dpi = 300
)


# pathway most active - visualise most responsive genes (progeny weights)
#along with their statistic values to interpret results
#scatterplot shows which genes contributing most to the pathway enrichment

prog_matrixD28only <- getModel("Human", top=100) %>% 
  as.data.frame()  %>%
  tibble::rownames_to_column("GeneID")

DED28only_df <- DED28only_mat %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("GeneID")

sum(is.na(DED28only_df$statistic))
DED28only_df <- DED28only_df %>%
  dplyr::filter(!is.na(statistic))

scat_plotsD28 <- progeny::progenyScatter(df = DED28only_df, 
                                         weight_matrix = prog_matrixD28only, 
                                         statName = "statistic", verbose = FALSE)
plot(scat_plotsD28[[1]]$`TGFb`) 

#top 100 genes
PathwayActivity_countsD5only <- progeny(nc_D5only_mat, scale = TRUE,
                                        organism = "Human", top = 100)

#define heatmap colour breaks
Activity_countsD5only <- as.vector(PathwayActivity_countsD5only) 

paletteLength <- 100
myColor <- 
  colorRampPalette(c("darkblue", "whitesmoke","indianred"))(paletteLength)

progenyBreaksD5only <- c(seq(min(Activity_countsD5only), 0, 
                             length.out=ceiling(paletteLength/2) + 1),
                         seq(max(Activity_countsD5only)/paletteLength, 
                             max(Activity_countsD5only), 
                             length.out=floor(paletteLength/2)))

annotation_col <- metafiltD5 %>%
  dplyr::select(Sample_ID, condition) %>%
  as.data.frame()

rownames(annotation_col) <- annotation_col$Sample_ID
annotation_col$Sample_ID <- NULL

# force condition order
annotation_col$condition <- factor(
  annotation_col$condition,
  levels = c("CTRL", "HNRNPUdel")
)

# order samples
sample_order <- rownames(annotation_col)[order(annotation_col$condition)]

# reorder matrix
matD5only <- t(PathwayActivity_countsD5only)[, sample_order, drop = FALSE]

# reorder annotation (IMPORTANT: keep as data.frame)
annotation_col2 <- annotation_col[sample_order, , drop = FALSE]

# plot
progeny_hmapD5only <- pheatmap(
  matD5only,
  fontsize = 14,
  fontsize_row = 10,
  fontsize_col = 10,
  color = myColor,
  breaks = progenyBreaksD5only,
  main = "PROGENy (100)",
  angle_col = 45,
  treeheight_col = 0,
  border_color = NA,
  annotation_col = annotation_col2,
  cluster_cols = FALSE,
  gaps_col = sum(annotation_col2$condition == "CTRL")
)


ggsave(
  filename = "PROGENY_heatmapD5only.png",
  plot = progeny_hmapD5only,
  width = 8,
  height = 12,
  dpi = 300
)

#enrichment analysis using competitive permutation approach to assess
#significance of pathway activity to end with NES for each pathway
PathwayActivity_zscoreD5only <- progeny(DED5only_mat, 
                                        scale=TRUE, organism="Human", top = 100, perm = 10000, z_scores = TRUE) %>%
  t()
colnames(PathwayActivity_zscoreD5only) <- "NES"

PathwayActivity_zscoreD5only_df <- as.data.frame(PathwayActivity_zscoreD5only) %>% 
  rownames_to_column(var = "Pathway") %>%
  dplyr::arrange(NES) %>%
  dplyr::mutate(Pathway = factor(Pathway))

NES_pathwayD5only <- ggplot(PathwayActivity_zscoreD5only_df,aes(x = reorder(Pathway, NES), y = NES)) + 
  geom_bar(aes(fill = NES), stat = "identity") +
  scale_fill_gradient2(low = "darkblue", high = "indianred", 
                       mid = "whitesmoke", midpoint = 0) + 
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = 
          element_text(angle = 45, hjust = 1, size =10, face= "bold"),
        axis.text.y = element_text(size =10, face= "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("Pathways")

ggsave(
  filename = "PROGENY_NES_pathwaysD5only.png",
  plot = NES_pathwayD5only,
  width = 10,
  height = 8,
  dpi = 300
)

# pathway most active - visualise most responsive genes (progeny weights)
#along with their statistic values to interpret results
#scatterplot shows which genes contributing most to the pathway enrichment

prog_matrixD5only <- getModel("Human", top=100) %>% 
  as.data.frame()  %>%
  tibble::rownames_to_column("GeneID")

DED5only_df <- DED5only_mat %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("GeneID")

sum(is.na(DED5only_df$statistic))
DED5only_df <- DED5only_df %>%
  dplyr::filter(!is.na(statistic))

scat_plotsD5 <- progeny::progenyScatter(df = DED5only_df, 
                                        weight_matrix = prog_matrixD5only, 
                                        statName = "statistic", verbose = FALSE)
plot(scat_plotsD5[[1]]$`TGFb`) 

############## ignore
# 1. Setup the list of raw data
raw_data_list <- list(
  D0 = PathwayActivity_countsD0only, 
  D5 = PathwayActivity_countsD5only, 
  D28 = PathwayActivity_countsD28only
)

# Fix JAK-STAT in clean_matrix
rownames(clean_matrix) <- gsub("JAK\\.STAT", "JAK-STAT", rownames(clean_matrix))

star_matrix <- matrix("", nrow = nrow(clean_matrix), ncol = ncol(clean_matrix))
colnames(star_matrix) <- c("D0", "D5", "D28")
rownames(star_matrix) <- rownames(clean_matrix)

# 2. Run the Loop
for (day in c("D0", "D5", "D28")) {
  
  # Get data & Match metadata (Same as before)
  current_counts <- as.data.frame(raw_data_list[[day]])
  meta_indices <- match(rownames(current_counts), meta_filtered$Sample_ID)
  conditions <- meta_filtered$condition[meta_indices]
  
  for (pathway in rownames(clean_matrix)) {
    if (pathway %in% colnames(current_counts)) {
      scores <- current_counts[[pathway]]
      valid_conditions <- na.omit(conditions)
      
      if (length(unique(valid_conditions)) >= 2) {
        tryCatch({
          test_res <- t.test(scores ~ conditions)
          pval <- test_res$p.value
          
          # --- NEW LOGIC FOR STARS ---
          if (pval < 0.001) {
            star_matrix[pathway, day] <- "***"  # Very strong
          } else if (pval < 0.01) {
            star_matrix[pathway, day] <- "**"   # Strong
          } else if (pval < 0.05) {
            star_matrix[pathway, day] <- "*"    # Significant
          }
          # ---------------------------
          
        }, error = function(e) { return(NULL) })
      }
    }
  }
}

# Check the upgrade
print(star_matrix)

library(pheatmap)
library(RColorBrewer)

pheatmap(clean_matrix,
         # Clustering
         cluster_cols = FALSE,      
         cluster_rows = TRUE,       
         
         # Colors
         color = my_colors,
         breaks = my_breaks,
         border_color = "white",    
         
         # --- STARS ---
         display_numbers = star_matrix, 
         fontsize_number = 12,          # Slightly smaller to fit *** if needed
         number_color = "black",        
         
         # Layout
         fontsize = 12,             
         cellwidth = 45,                # Wider cells for 3 stars
         cellheight = 20,
         
         main = "PROGENy",
         
         filename = "Figure1_Final_Significance.pdf", 
         width = 6, height = 7
)

# 1. Create a matrix to store numeric P-values (instead of stars)
numeric_pval_matrix <- matrix(NA, nrow = nrow(clean_matrix), ncol = 3)
colnames(numeric_pval_matrix) <- c("D0", "D5", "D28")
rownames(numeric_pval_matrix) <- rownames(clean_matrix)

# 2. Run the Loop (Same logic, just saving numbers)
for (day in c("D0", "D5", "D28")) {
  
  # Get data & Match metadata
  current_counts <- as.data.frame(raw_data_list[[day]])
  meta_indices <- match(rownames(current_counts), meta_filtered$Sample_ID)
  conditions <- meta_filtered$condition[meta_indices]
  
  for (pathway in rownames(clean_matrix)) {
    if (pathway %in% colnames(current_counts)) {
      scores <- current_counts[[pathway]]
      valid_conditions <- na.omit(conditions)
      
      if (length(unique(valid_conditions)) >= 2) {
        tryCatch({
          # Run Test
          test_res <- t.test(scores ~ conditions)
          
          # SAVE THE NUMBER
          numeric_pval_matrix[pathway, day] <- test_res$p.value
          
        }, error = function(e) { return(NULL) })
      }
    }
  }
}

# 3. View the table (Rounded for easier reading)
print("Raw P-Values:")
print(round(numeric_pval_matrix, 5)) # Shows 5 decimal places

# 4. (Optional) Save to CSV for your thesis/paper
write.csv(numeric_pval_matrix, "Table_S1_Progeny_PValues.csv")

# Load necessary libraries
library(tidyverse)
install.packages("ggvenn")
library(ggvenn) # Great for simple, pretty Venn diagrams

#####TFactivities.R#####
n_tfs <- 25
## top 25 TF in general (abs)
# Filter for top 25 most active TFs (by absolute score)
top_tfsD0 <- TFscores_D0 %>%
  arrange(desc(abs(score))) %>%
  slice_head(n = 25)

# Plot them directly 
top25allD0 <- ggplot(top_tfsD0, aes(x = reorder(source, score), y = score, fill = score > 0)) +
  geom_col() +
  coord_flip() + # Makes it horizontal
  theme_minimal() +
  labs(title = "Top 25 TF Activities (ASD vs CTRL)", 
       y = "Normalized Enrichment Score (NES)", 
       x = "Transcription Factor") +
  scale_fill_manual(values = c("TRUE" = "#E41A1C", "FALSE" = "#377EB8"), 
                    labels = c("Activated", "Inhibited"),
                    name = "Status")

ggsave(
  filename = "top25abs_D0.png",      # use PDF or PNG
  plot = top25allD0,           
  width = 8, height = 10,          # size in inches
  dpi = 300                       # resolution for PNG (ignored for PDF)
)

##top 15 for activated and top 15 for inhibited TFs
# Select Top 15 Activated (Highest Positive Scores)
top_up <- TFscores_D0 %>%
  arrange(desc(score)) %>%
  slice_head(n = 15)

# Select Top 15 Inhibited (Lowest Negative Scores)
top_down <- TFscores_D0 %>%
  arrange(score) %>%
  slice_head(n = 15)

# Combine them into one dataframe
top_30_balanced <- bind_rows(top_up, top_down)

# We use reorder() to ensure they appear sorted by score, not alphabetically
top30D0 <- ggplot(top_30_balanced, aes(x = reorder(source, score), y = score, fill = score > 0)) +
  geom_col() +
  coord_flip() + # Flips to horizontal bars for readability
  theme_minimal() +
  labs(title = "Top 15 Activated & Inhibited TFs",
       subtitle = "Differential Analysis (ASD vs CTRL)",
       x = "Transcription Factor",
       y = "NES (Activity Score)",
       fill = "Regulation") +
  scale_fill_manual(values = c("TRUE" = "#E41A1C", "FALSE" = "#377EB8"), 
                    labels = c("Inhibited", "Activated")) +
  theme(
    axis.text.y = element_text(size = 10), # Adjust text size if needed
    legend.position = "top"
  )

ggsave(
  filename = "top30TFs_D0.png",      # use PDF or PNG
  plot = top30D0,           
  width = 10, height = 8,          # size in inches
  dpi = 300                       # resolution for PNG (ignored for PDF)
)

## For the volcano plot (related to support functions)
library(ggrepel)
## We also load the support functions
source("support_functions.R")

#To interpret the results, we can look at the expression of targets of one of the most deregulated TFs
targets_ETS1 <- net$target[net$source == "ETS1"]
DEG_ETS1 <- DED0only[DED0only$hgnc_symbol %in% targets_ETS1, ]
volcano_nice(as.data.frame(DEG_ETS1), 
             FCIndex = 3,     # log2FoldChange column
             pValIndex = 7,   # padj column
             IDIndex = 9,     # gene ID column
             nlabels = 20, 
             label = TRUE, 
             straight = FALSE)

n_tfs <- 25

library(dplyr)
## top 25 TF in general (abs)
# Filter for top 25 most active TFs (by absolute score)
top_tfsD5 <- TFscores_D5 %>%
  arrange(desc(abs(score))) %>%
  slice_head(n = 25)

# Plot them directly 
top25allD5 <- ggplot(top_tfsD5, aes(x = reorder(source, score), y = score, fill = score > 0)) +
  geom_col() +
  coord_flip() + # Makes it horizontal
  theme_minimal() +
  labs(title = "Top 25 TF Activities (ASD vs CTRL)", 
       y = "Normalized Enrichment Score (NES)", 
       x = "Transcription Factor") +
  scale_fill_manual(values = c("TRUE" = "#E41A1C", "FALSE" = "#377EB8"), 
                    labels = c("Activated", "Inhibited"),
                    name = "Status")

ggsave(
  filename = "top25abs_D5.png",      # use PDF or PNG
  plot = top25allD5,           
  width = 8, height = 10,          # size in inches
  dpi = 300                       # resolution for PNG (ignored for PDF)
)

##top 15 for activated and top 15 for inhibited TFs
# Select Top 15 Activated (Highest Positive Scores)
top_up <- TFscores_D5 %>%
  arrange(desc(score)) %>%
  slice_head(n = 15)

# Select Top 15 Inhibited (Lowest Negative Scores)
top_down <- TFscores_D5 %>%
  arrange(score) %>%
  slice_head(n = 15)

# Combine them into one dataframe
top_30_balanced <- bind_rows(top_up, top_down)

# We use reorder() to ensure they appear sorted by score, not alphabetically
top30D5 <- ggplot(top_30_balanced, aes(x = reorder(source, score), y = score, fill = score > 0)) +
  geom_col() +
  coord_flip() + # Flips to horizontal bars for readability
  theme_minimal() +
  labs(title = "Top 15 Activated & Inhibited TFs",
       subtitle = "Differential Analysis (ASD vs CTRL)",
       x = "Transcription Factor",
       y = "NES (Activity Score)",
       fill = "Regulation") +
  scale_fill_manual(values = c("TRUE" = "#E41A1C", "FALSE" = "#377EB8"), 
                    labels = c("Inhibited", "Activated")) +
  theme(
    axis.text.y = element_text(size = 10), # Adjust text size if needed
    legend.position = "top"
  )

ggsave(
  filename = "top30TFs_D5.png",      # use PDF or PNG
  plot = top30D5,           
  width = 10, height = 8,          # size in inches
  dpi = 300                       # resolution for PNG (ignored for PDF)
)

## For the volcano plot (related to support functions)
library(ggrepel)
## We also load the support functions
source("support_functions.R")

#To interpret the results, we can look at the expression of targets of one of the most deregulated TFs
targets_SP1 <- net$target[net$source == "SP1"]
DEG_SP1 <- DED5only[DED5only$hgnc_symbol %in% targets_SP1, ]
volcano_nice(as.data.frame(DEG_SP1), 
             FCIndex = 3,     # log2FoldChange column
             pValIndex = 7,   # padj column
             IDIndex = 9,     # gene ID column
             nlabels = 20, 
             label = TRUE, 
             straight = FALSE)


n_tfs <- 25

library(dplyr)
## top 25 TF in general (abs)
# Filter for top 25 most active TFs (by absolute score)
top_tfsD28 <- TFscores_D28 %>%
  arrange(desc(abs(score))) %>%
  slice_head(n = 25)

# Plot them directly 
top25allD28 <- ggplot(top_tfsD28, aes(x = reorder(source, score), y = score, fill = score > 0)) +
  geom_col() +
  coord_flip() + # Makes it horizontal
  theme_minimal() +
  labs(title = "Top 25 TF Activities (ASD vs CTRL)", 
       y = "Normalized Enrichment Score (NES)", 
       x = "Transcription Factor") +
  scale_fill_manual(values = c("TRUE" = "#E41A1C", "FALSE" = "#377EB8"), 
                    labels = c("Activated", "Inhibited"),
                    name = "Status")

ggsave(
  filename = "top25abs_D28.png",      # use PDF or PNG
  plot = top25allD28,           
  width = 8, height = 10,          # size in inches
  dpi = 300                       # resolution for PNG (ignored for PDF)
)

##top 15 for activated and top 15 for inhibited TFs
# Select Top 15 Activated (Highest Positive Scores)
top_up <- TFscores_D28 %>%
  arrange(desc(score)) %>%
  slice_head(n = 15)

# Select Top 15 Inhibited (Lowest Negative Scores)
top_down <- TFscores_D28 %>%
  arrange(score) %>%
  slice_head(n = 15)

# Combine them into one dataframe
top_30_balanced <- bind_rows(top_up, top_down)

# We use reorder() to ensure they appear sorted by score, not alphabetically
top30D28 <- ggplot(top_30_balanced, aes(x = reorder(source, score), y = score, fill = score > 0)) +
  geom_col() +
  coord_flip() + # Flips to horizontal bars for readability
  theme_minimal() +
  labs(title = "Top 15 Activated & Inhibited TFs",
       subtitle = "Differential Analysis (ASD vs CTRL)",
       x = "Transcription Factor",
       y = "NES (Activity Score)",
       fill = "Regulation") +
  scale_fill_manual(values = c("TRUE" = "#E41A1C", "FALSE" = "#377EB8"), 
                    labels = c("Inhibited", "Activated")) +
  theme(
    axis.text.y = element_text(size = 10), # Adjust text size if needed
    legend.position = "top"
  )

ggsave(
  filename = "top30TFs_D28.png",      # use PDF or PNG
  plot = top30D28,           
  width = 10, height = 8,          # size in inches
  dpi = 300                       # resolution for PNG (ignored for PDF)
)

## For the volcano plot (related to support functions)
library(ggrepel)
## We also load the support functions
source("support_functions.R")

#To interpret the results, we can look at the expression of targets of one of the most deregulated TFs
targets_SP1 <- net$target[net$source == "SP1"]
DEG_SP1 <- DED28only[DED28only$hgnc_symbol %in% targets_SP1, ]
volcano_nice(as.data.frame(DEG_SP1), 
             FCIndex = 3,     # log2FoldChange column
             pValIndex = 7,   # padj column
             IDIndex = 9,     # gene ID column
             nlabels = 20, 
             label = TRUE, 
             straight = FALSE)


###### OPTIONAL #######
#ALT way of representing TFs
tf_activities_statf_top25 <- TFscores_D0 %>%
  as.data.frame() %>% 
  dplyr::filter(p_value < 0.05) %>%  # optional: only significant TFs
  dplyr::top_n(25, wt = abs(score)) %>%
  dplyr::arrange(score) %>% 
  dplyr::mutate(source = factor(source, levels = unique(source)))

contrast_tftop25f <- ggplot(tf_activities_statf_top25,aes(x = reorder(source, score), y = score)) + 
  geom_bar(aes(fill = score), stat = "identity") +
  scale_fill_gradient2(low = "darkblue", high = "indianred", 
                       mid = "whitesmoke", midpoint = 0) + 
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = 
          element_text(angle = 45, hjust = 1, size =10, face= "bold"),
        axis.text.y = element_text(size =10, face= "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("Transcription Factors")

ggsave(
  filename = "barplot_top25tfscores.png",
  plot = contrast_tftop25f,
  width = 12,
  height = 8,
  dpi = 300
)

##### MOON #####
##### MOON #####
library(cosmosR)

#TF scores for MOON
TFscoresMOON <- TFscores_D0[, c('source', 'score')]
names(TFscoresMOON) <- c("TF", "score")
head(TFscoresMOON)
TFscoresMOON <- as.data.frame(TFscoresMOON)
row.names(TFscoresMOON) <- TFscoresMOON$TF
head(TFscoresMOON)

clean_PKN

#DEA for MOON
DEA_MOON <- resD0.annot[, c("hgnc_symbol", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
head(DEA_MOON)
#data frame, 13821 transcripts

# Remove rows with blank hgnc_symbol
DEA_MOON <- DEA_MOON[DEA_MOON$hgnc_symbol != "", ]
# Remove duplicated hgnc_symbol rows (keeping the first one)
DEA_MOON <- DEA_MOON[!duplicated(DEA_MOON$hgnc_symbol), ]
row.names(DEA_MOON) <- DEA_MOON$hgnc_symbol
#12878 transcripts

names(DEA_MOON)[1] <- "gene"
row.names(DEA_MOON) <- DEA_MOON$gene
head(DEA_MOON)
DEA_clean_MOON <- DEA_MOON[, "stat", drop = FALSE] 

# Rename the column to represent comparison
names(DEA_clean_MOON) <- "ASD_vs_CTRL"
DEA_moon <- DEA_clean_MOON
head(DEA_moon)
conditions_of_interest <- "ASD_vs_CTRL"

head(DEA_moon)

PKND0 <- cosmosR:::filter_pkn_expressed_genes(row.names(DEA_moon), meta_pkn = clean_PKN)

#####
DEA_moon 
obsgenes <- TFscoresMOON$TF
PKND0

head(names(TFscoresMOON))
head(obsgenes)

if(!"source" %in% names(PKND0)) { names(PKND0)[1] <- "source" }
if(!"target" %in% names(PKND0)) { names(PKND0)[2] <- "target" }

network_nodesD0 <- unique(c(PKND0$source, PKND0$target))
overlap_countD0 <- sum(obsgenes %in% network_nodesD0)
print(paste("Original genes to observe:", length(obsgenes)))
print(paste("Genes found in the network:", overlap_countD0))

obsgenes_D0 <- obsgenes[obsgenes %in% network_nodesD0]
# Prune the Network (Keep only nodes upstream of your TFs)
n_steps <- 6
PKNfiltD0 <- cosmosR:::keep_observable_neighbours(PKND0, n_steps, obsgenes_D0)

# Compress the network (Simplify redundant paths)
meta_network_compressed_listD0 <- compress_same_children(PKNfiltD0, sig_input = NULL, metab_input = obsgenes_D0)
meta_network_compressed_D0 <- meta_network_compressed_listD0$compressed_network

node_signaturesD0 <- meta_network_compressed_listD0$node_signaturesD0
duplicated_parentsD0 <- meta_network_compressed_listD0$duplicated_signatures
meta_network_compressed_D0 <- meta_network_cleanup(meta_network_compressed_D0)


# 7. The Main Loop (The "MOON" Algorithm)
# This iterates to find the best signaling path explaining your data
meta_network_TFD0 <- meta_network_compressed_D0
i <- 1
before <- 1
after <- 0

load("collectri_regulon_R.RData")

TFscoresMOON_vec <- TFscoresMOON$score
names(TFscoresMOON_vec) <- TFscoresMOON$TF

print("Starting COSMOS optimization...")
while (before != after & i < 10) {
  before <- length(meta_network_TFD0[,1])
  
  # Run MOON: Upstream=NULL because we don't know the drug target
  moon_resD0 <- moon(upstream_input = NULL, 
                     downstream_input = TFscoresMOON_vec, 
                     meta_network = meta_network_TFD0, 
                     n_layers = n_steps, 
                     statistic = "ulm") 
  
  # Prune incoherent links (contradictions between Kinase -> TF -> Gene)
  meta_network_TFD0 <- filter_incohrent_TF_target(moon_resD0, net, meta_network_TFD0, DEA_moon)
  
  after <- length(meta_network_TFD0[,1])
  i <- i + 1
  print(paste("Iteration:", i, "- Network size:", after))
}

getAnywhere("filter_incohrent_TF_target")

nrow(meta_network_TFD0)
head(meta_network_TFD0)


# Look at the top regulators (Kinases/Signaling Proteins)
top_regulators <- moon_resD0 %>%
  dplyr::arrange(desc(abs(score))) # Sort by magnitude of score

# Display the top 10 potential drivers of your disease
print(head(top_regulators, 10))

# Visualize it nicely
topregD0 <- ggplot(head(top_regulators, 20), aes(x = reorder(source, score), y = score, fill = score > 0)) +
  geom_col() +
  coord_flip() +
  theme_minimal() +
  labs(title = "Top Upstream Regulators (Disease vs Control)",
       x = "Signaling Protein",
       y = "Predicted Activity Score") +
  scale_fill_manual(values = c("TRUE" = "firebrick", "FALSE" = "steelblue"), 
                    labels = c("Activated", "Inhibited"))

ggsave(
  filename = "top_upstream_regulatorsD0.png",
  plot = topregD0,
  width = 8,
  height = 6,
  dpi = 300
)

moonMAP2K2 <- get_moon_scoring_network("MAP2K2", meta_network_TFD0, moon_resD0)
library(igraph)


# Create igraph object
g_moon <- graph_from_data_frame(
  d = moonMAP2K2$SIF,
  directed = TRUE,
  vertices = moonMAP2K2$ATT  # optional, adds attributes
)

plot(
  g_moon,
  vertex.size = 6,
  vertex.label.cex = 0.7,
  vertex.color = "steelblue",
  edge.arrow.size = 0.3
)

library(ggraph)
library(ggplot2)

# Optionally, color nodes by MOON score
V(g_moon)$score <- moonMAP2K2$ATT$score[match(V(g_moon)$name, moonMAP2K2$ATT$source)]

ggraph(g_moon, layout = "fr") +
  geom_edge_link(aes(edge_alpha = 0.5, edge_color = factor(E(g_moon)$interaction))) +
  geom_node_point(aes(color = score), size = 5) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  scale_color_gradient2(low = "steelblue", mid = "white", high = "firebrick", midpoint = 0) +
  scale_edge_color_manual(values = c("1" = "darkgreen", "-1" = "red")) +
  theme_void() +
  labs(title = "MOON Network for MAP2K2")


##### D5 #####
library(cosmosR)

#TF scores for MOON
TFscoresMOOND5 <- TFscores_D5[, c('source', 'score')]
names(TFscoresMOOND5) <- c("TF", "score")
head(TFscoresMOOND5)
TFscoresMOOND5 <- as.data.frame(TFscoresMOOND5)
row.names(TFscoresMOOND5) <- TFscoresMOOND5$TF
head(TFscoresMOOND5)

clean_PKN

#DEA for MOON
DEA_MOOND5 <- res_D5_only.annot[, c("hgnc_symbol", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
head(DEA_MOOND5)


# Remove rows with blank hgnc_symbol
DEA_MOOND5 <- DEA_MOOND5[DEA_MOOND5$hgnc_symbol != "", ]
# Remove duplicated hgnc_symbol rows (keeping the first one)
DEA_MOOND5 <- DEA_MOOND5[!duplicated(DEA_MOOND5$hgnc_symbol), ]
row.names(DEA_MOOND5) <- DEA_MOOND5$hgnc_symbol
#12878 transcripts

names(DEA_MOOND5)[1] <- "gene"
row.names(DEA_MOOND5) <- DEA_MOOND5$gene
head(DEA_MOOND5)
DEA_clean_MOOND5 <- DEA_MOOND5[, "stat", drop = FALSE] 

# Rename the column to represent comparison
names(DEA_clean_MOOND5) <- "ASD_vs_CTRL"
DEA_moonD5 <- DEA_clean_MOOND5
head(DEA_moonD5)
conditions_of_interest <- "ASD_vs_CTRL"

head(DEA_moonD5)

PKND5 <- cosmosR:::filter_pkn_expressed_genes(row.names(DEA_moonD5), meta_pkn = clean_PKN)

#####
DEA_moonD5
obsgenesD5 <- TFscoresMOOND5$TF
PKND5

head(names(TFscoresMOOND5))
head(obsgenesD5)

if(!"source" %in% names(PKND5)) { names(PKND5)[1] <- "source" }
if(!"target" %in% names(PKND5)) { names(PKND5)[2] <- "target" }

network_nodesD5 <- unique(c(PKND5$source, PKND5$target))
overlap_countD5 <- sum(obsgenesD5 %in% network_nodesD5)
print(paste("Original genes to observe:", length(obsgenesD5)))
print(paste("Genes found in the network:", overlap_countD5))

obsgenes_D5 <- obsgenesD5[obsgenesD5 %in% network_nodesD5]
# Prune the Network (Keep only nodes upstream of your TFs)
n_steps <- 6
PKNfiltD5 <- cosmosR:::keep_observable_neighbours(PKND5, n_steps, obsgenes_D5)

# Compress the network (Simplify redundant paths)
meta_network_compressed_listD5 <- compress_same_children(PKNfiltD5, sig_input = NULL, metab_input = obsgenes_D5)
meta_network_compressed_D5 <- meta_network_compressed_listD5$compressed_network

node_signaturesD5 <- meta_network_compressed_listD5$node_signaturesD5
duplicated_parentsD5 <- meta_network_compressed_listD5$duplicated_signatures
meta_network_compressed_D5 <- meta_network_cleanup(meta_network_compressed_D5)


# 7. The Main Loop (The "MOON" Algorithm)
# This iterates to find the best signaling path explaining your data
meta_network_TFD5 <- meta_network_compressed_D5
i <- 1
before <- 1
after <- 0

load("collectri_regulon_R.RData")

TFscoresMOON_vecD5 <- TFscoresMOOND5$score
names(TFscoresMOON_vecD5) <- TFscoresMOOND5$TF

print("Starting COSMOS optimization...")
while (before != after & i < 10) {
  before <- length(meta_network_TFD5[,1])
  
  # Run MOON: Upstream=NULL because we don't know the drug target
  moon_resD5 <- moon(upstream_input = NULL, 
                     downstream_input = TFscoresMOON_vecD5, 
                     meta_network = meta_network_TFD5, 
                     n_layers = n_steps, 
                     statistic = "ulm") 
  
  # Prune incoherent links (contradictions between Kinase -> TF -> Gene)
  meta_network_TFD5 <- filter_incohrent_TF_target(moon_resD5, net, meta_network_TFD5, DEA_moonD5)
  
  after <- length(meta_network_TFD5[,1])
  i <- i + 1
  print(paste("Iteration:", i, "- Network size:", after))
}

getAnywhere("filter_incohrent_TF_target")

nrow(meta_network_TFD5)
head(meta_network_TFD5)


# Look at the top regulators (Kinases/Signaling Proteins)
top_regulatorsD5 <- moon_resD5 %>%
  dplyr::arrange(desc(abs(score))) # Sort by magnitude of score

# Display the top 10 potential drivers of your disease
print(head(top_regulatorsD5, 10))

# Visualize it nicely
topregD5 <- ggplot(head(top_regulatorsD5, 20), aes(x = reorder(source, score), y = score, fill = score > 0)) +
  geom_col() +
  coord_flip() +
  theme_minimal() +
  labs(title = "Top Upstream Regulators (Disease vs Control)",
       x = "Signaling Protein",
       y = "Predicted Activity Score") +
  scale_fill_manual(values = c("TRUE" = "firebrick", "FALSE" = "steelblue"), 
                    labels = c("Activated", "Inhibited"))

ggsave(
  filename = "top_upstream_regulatorsD5.png",
  plot = topregD5,
  width = 8,
  height = 6,
  dpi = 300
)
moonMAP2K2 <- get_moon_scoring_network("MAP2K2", meta_network_TFD5, moon_resD5)
library(igraph)


# Create igraph object
g_moon <- graph_from_data_frame(
  d = moonMAP2K2$SIF,
  directed = TRUE,
  vertices = moonMAP2K2$ATT  # optional, adds attributes
)

plot(
  g_moon,
  vertex.size = 6,
  vertex.label.cex = 0.7,
  vertex.color = "steelblue",
  edge.arrow.size = 0.3
)

library(ggraph)
library(ggplot2)

# Optionally, color nodes by MOON score
V(g_moon)$score <- moonTUB$ATT$score[match(V(g_moon)$name, moonTUB$ATT$source)]

ggraph(g_moon, layout = "fr") +
  geom_edge_link(aes(edge_alpha = 0.5, edge_color = factor(E(g_moon)$interaction))) +
  geom_node_point(aes(color = score), size = 5) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  scale_color_gradient2(low = "steelblue", mid = "white", high = "firebrick", midpoint = 0) +
  scale_edge_color_manual(values = c("1" = "darkgreen", "-1" = "red")) +
  theme_void() +
  labs(title = "MOON Network for MAP2K2")

##### D28 #####
library(cosmosR)

#TF scores for MOON
TFscoresMOOND28 <- TFscores_D28[, c('source', 'score')]
names(TFscoresMOOND28) <- c("TF", "score")
head(TFscoresMOOND28)
TFscoresMOOND28 <- as.data.frame(TFscoresMOOND28)
row.names(TFscoresMOOND28) <- TFscoresMOOND28$TF
head(TFscoresMOOND28)

clean_PKN

#DEA for MOON
DEA_MOOND28 <- res_D28_only.annot[, c("hgnc_symbol", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
head(DEA_MOOND28)


# Remove rows with blank hgnc_symbol
DEA_MOOND28 <- DEA_MOOND28[DEA_MOOND28$hgnc_symbol != "", ]
# Remove duplicated hgnc_symbol rows (keeping the first one)
DEA_MOOND28 <- DEA_MOOND28[!duplicated(DEA_MOOND28$hgnc_symbol), ]
row.names(DEA_MOOND28) <- DEA_MOOND28$hgnc_symbol
#12878 transcripts

names(DEA_MOOND28)[1] <- "gene"
row.names(DEA_MOOND28) <- DEA_MOOND28$gene
head(DEA_MOOND28)
DEA_clean_MOOND28 <- DEA_MOOND28[, "stat", drop = FALSE] 

# Rename the column to represent comparison
names(DEA_clean_MOOND28) <- "ASD_vs_CTRL"
DEA_moonD28 <- DEA_clean_MOOND28
head(DEA_moonD28)
conditions_of_interest <- "ASD_vs_CTRL"

head(DEA_moonD28)

PKND28 <- cosmosR:::filter_pkn_expressed_genes(row.names(DEA_moonD28), meta_pkn = clean_PKN)

#####
DEA_moonD28
obsgenesD28 <- TFscoresMOOND28$TF
PKND28

head(names(TFscoresMOOND28))
head(obsgenesD28)

if(!"source" %in% names(PKND28)) { names(PKND28)[1] <- "source" }
if(!"target" %in% names(PKND28)) { names(PKND28)[2] <- "target" }

network_nodesD28 <- unique(c(PKND28$source, PKND28$target))
overlap_countD28 <- sum(obsgenesD28 %in% network_nodesD28)
print(paste("Original genes to observe:", length(obsgenesD28)))
print(paste("Genes found in the network:", overlap_countD28))

obsgenes_D28 <- obsgenesD28[obsgenesD28 %in% network_nodesD28]
# Prune the Network (Keep only nodes upstream of your TFs)
n_steps <- 6
PKNfiltD28 <- cosmosR:::keep_observable_neighbours(PKND28, n_steps, obsgenes_D28)

# Compress the network (Simplify redundant paths)
meta_network_compressed_listD28 <- compress_same_children(PKNfiltD28, sig_input = NULL, metab_input = obsgenes_D28)
meta_network_compressed_D28 <- meta_network_compressed_listD28$compressed_network

node_signaturesD28 <- meta_network_compressed_listD28$node_signaturesD28
duplicated_parentsD28 <- meta_network_compressed_listD28$duplicated_signatures
meta_network_compressed_D28 <- meta_network_cleanup(meta_network_compressed_D28)


# 7. The Main Loop (The "MOON" Algorithm)
# This iterates to find the best signaling path explaining your data
meta_network_TFD28 <- meta_network_compressed_D28
i <- 1
before <- 1
after <- 0

load("collectri_regulon_R.RData")

TFscoresMOON_vecD28 <- TFscoresMOOND28$score
names(TFscoresMOON_vecD28) <- TFscoresMOOND28$TF

print("Starting COSMOS optimization...")
while (before != after & i < 10) {
  before <- length(meta_network_TFD28[,1])
  
  # Run MOON: Upstream=NULL because we don't know the drug target
  moon_resD28 <- moon(upstream_input = NULL, 
                      downstream_input = TFscoresMOON_vecD28, 
                      meta_network = meta_network_TFD28, 
                      n_layers = n_steps, 
                      statistic = "ulm") 
  
  # Prune incoherent links (contradictions between Kinase -> TF -> Gene)
  meta_network_TFD28 <- filter_incohrent_TF_target(moon_resD28, net, meta_network_TFD28, DEA_moonD28)
  
  after <- length(meta_network_TFD28[,1])
  i <- i + 1
  print(paste("Iteration:", i, "- Network size:", after))
}

getAnywhere("filter_incohrent_TF_target")

nrow(meta_network_TFD28)
head(meta_network_TFD28)


# Look at the top regulators (Kinases/Signaling Proteins)
top_regulatorsD28 <- moon_resD28 %>%
  dplyr::arrange(desc(abs(score))) # Sort by magnitude of score

# Display the top 10 potential drivers of your disease
print(head(top_regulatorsD28, 10))

# Visualize it nicely
topregD28 <- ggplot(head(top_regulatorsD28, 20), aes(x = reorder(source, score), y = score, fill = score > 0)) +
  geom_col() +
  coord_flip() +
  theme_minimal() +
  labs(title = "Top Upstream Regulators (Disease vs Control)",
       x = "Signaling Protein",
       y = "Predicted Activity Score") +
  scale_fill_manual(values = c("TRUE" = "firebrick", "FALSE" = "steelblue"), 
                    labels = c("Activated", "Inhibited"))
ggsave(
  filename = "top_upstream_regulatorsD28.png",
  plot = topregD28,
  width = 8,
  height = 6,
  dpi = 300
)

moonMAP2K2 <- get_moon_scoring_network("MAP2K2", meta_network_TFD28, moon_resD28)
library(igraph)


# Create igraph object
g_moon <- graph_from_data_frame(
  d = moonMAP2K2$SIF,
  directed = TRUE,
  vertices = moonMAP2K2$ATT  # optional, adds attributes
)

plot(
  g_moon,
  vertex.size = 6,
  vertex.label.cex = 0.7,
  vertex.color = "steelblue",
  edge.arrow.size = 0.3
)

library(ggraph)
library(ggplot2)

# Optionally, color nodes by MOON score
V(g_moon)$score <- moonTUB$ATT$score[match(V(g_moon)$name, moonTUB$ATT$source)]

ggraph(g_moon, layout = "fr") +
  geom_edge_link(aes(edge_alpha = 0.5, edge_color = factor(E(g_moon)$interaction))) +
  geom_node_point(aes(color = score), size = 5) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  scale_color_gradient2(low = "steelblue", mid = "white", high = "firebrick", midpoint = 0) +
  scale_edge_color_manual(values = c("1" = "darkgreen", "-1" = "red")) +
  theme_void() +
  labs(title = "MOON Network for MAP2K2")

#####CARNIVALanalysis.R#####
full_pkn <- as.data.frame(omnipath_interactions()) #85217 interactions
full_pkn <- full_pkn[!is.na(full_pkn$references),]
head(full_pkn)
clean_PKN <- full_pkn[full_pkn$consensus_stimulation == 1 | full_pkn$consensus_inhibition == 1,]

clean_PKN$sign <- clean_PKN$consensus_stimulation - clean_PKN$consensus_inhibition

clean_PKN <- clean_PKN[,c(3,4,16)]
head(clean_PKN)

clean_PKN_supp <- clean_PKN[clean_PKN$sign == 0,]
clean_PKN_supp$sign <- -1
clean_PKN[clean_PKN$sign == 0,"sign"] <- 1

clean_PKN <- as.data.frame(rbind(clean_PKN, clean_PKN_supp))
#71273 interactions

names(clean_PKN) <- c("source","target","interaction")
write_csv(x = clean_PKN, file = "clean_omnipath_PKN.csv")



head(TFscoresD0)
TFscoresD0_CARNIVAL <- TFscoresD0[, "score"]
class(TFscoresD0_CARNIVAL)
# numeric
names(TFscoresD0_CARNIVAL)
# "AHR" "AP1" "APEX1" ...

head(PathwayActivity_CARNIVALinputD0only)
Pathway_CARNIVALD0 <- PathwayActivity_CARNIVALinputD0only$score
names(Pathway_CARNIVALD0) <- PathwayActivity_CARNIVALinputD0only$Pathway
head(Pathway_CARNIVALD0)
class(Pathway_CARNIVALD0)
# numeric
names(Pathway_CARNIVALD0)
# "Androgen" "EGFR" "Estrogen" "Hypoxia" "JAK.STAT" "MAPK"


##### CARNIVAL #####
sif_clean <- read.csv("clean_omnipath_PKN.csv")



?runCARNIVAL
carnival_D0 <- runCARNIVAL(
  inputObj  = NULL,             
  measObj   = TFscoresD0_CARNIVAL,             # TF activity measurements
  netObj    = clean_PKN,                    # prior knowledge network
  weightObj = Pathway_CARNIVALD0,         # pathway activity scores (or progeny_matrix)
  solverPath = "C:/Program Files/IBM/ILOG/CPLEX_Studio2211/cplex/bin/x64_win64/cplex.exe", 
  solver     = "cplex",
  timelimit  = 7200, 
  mipGAP     = 0.05,
  poolrelGAP = 0
)
saveRDS(carnival_D0,"~/hnrnpu-causal-multiomics/processeddata/carnival_D0.rds")
carnival_D0 <- readRDS("~/hnrnpu-causal-multiomics/processeddata/carnival_D0.rds")
head(carnival_D0)

#inputObj: data frame of list fir target of perturbation, OPTIONAL or default set to NULL to run invCARNIVAL when inputs are not known

carnival_D0 <- runCARNIVAL(
  inputObj  = iniciators,             
  measObj   = TFscoresD0_CARNIVAL,             # TF activity measurements
  netObj    = sif,                    # prior knowledge network
  weightObj = Pathway_CARNIVALD0,         # pathway activity scores (or progeny_matrix)
  solverPath = "C:/Program Files/IBM/ILOG/CPLEX_Studio2211/cplex/bin/x64_win64/cplex.exe", 
  solver     = "cplex",
  timelimit  = 7200, 
  mipGAP     = 0.05,
  poolrelGAP = 0
)
saveRDS(carnival_D0,"~/hnrnpu-causal-multiomics/processeddata/carnival_D0.rds")

clean_edgesD0 <- carnival_D0$weightedSIF %>%
  filter(Weight > 0)

active_node_namesD0 <- unique(c(clean_edgesD0$Node1, clean_edgesD0$Node2))

clean_nodesD0 <- carnival_D0$nodesAttributes %>%
  filter(Node %in% active_node_namesD0)

###### visual ######
library(RCy3)
library(igraph)
library(dplyr)
write.csv(clean_edgesD0, "Clean_CARNIVAL_EdgesD0.csv", row.names = FALSE)
write.csv(clean_nodesD0, "Clean_CARNIVAL_NodesD0.csv", row.names = FALSE)

cytoscapePing() 

# Rename columns to standard RCy3 formats to be safe
# (RCy3 prefers 'source', 'target', 'interaction' but handles others too)
edges_for_cyD0 <- clean_edgesD0
colnames(edges_for_cyD0) <- c("source", "sign", "target", "weight") # Adjust based on your actual column order!

nodes_for_cyD0 <- clean_nodesD0
colnames(nodes_for_cyD0)[1] <- "id" # Ensure first column is 'id'

# Send to Cytoscape
createNetworkFromDataFrames(nodes = nodes_for_cyD0, 
                            edges = edges_for_cyD0, 
                            title = "CARNIVAL_Automated", 
                            collection = "CARNIVAL_Analysis")

setNodeColorMapping(table.column = "AvgAct", 
                    table.column.values = c(-100, 0, 100), 
                    colors = c("#0000FF", "#FFFFFF", "#FF0000"), 
                    mapping.type = "continuous")

## top 30
top30tfs_D0 <- TFscores_D0 %>%
  arrange(desc(abs(score))) %>%
  slice_head(n = 30)
head(top30tfs_D0)

top30tfs_D0 <- top30tfs_D0 %>%
  dplyr::select(source, score) %>%
  deframe()

carnival_D0top30 <- runCARNIVAL(
  inputObj  = iniciators,             
  measObj   = top30tfs_D0,             # TF activity measurements
  netObj    = sif,                    # prior knowledge network
  weightObj = Pathway_CARNIVALD0,         # pathway activity scores (or progeny_matrix)
  solverPath = "C:/Program Files/IBM/ILOG/CPLEX_Studio2211/cplex/bin/x64_win64/cplex.exe", 
  solver     = "cplex",
  timelimit  = 900, 
  mipGAP     = 0.05,
  poolrelGAP = 0
)
saveRDS(carnival_D0top30,"~/hnrnpu-causal-multiomics/processeddata/carnival_D0top30.rds")
carnival_D0top30 <- readRDS("~/hnrnpu-causal-multiomics/processeddata/carnival_D0top30.rds")

clean_edgesD0_30 <- carnival_D0top30$weightedSIF %>%
  filter(Weight > 0)

active_node_namesD0_30 <- unique(c(clean_edgesD0_30$Node1, clean_edgesD0_30$Node2))

clean_nodesD0_30 <- carnival_D0top30$nodesAttributes %>%
  filter(Node %in% active_node_namesD0_30)

write.csv(clean_edgesD0_30, "Clean_CARNIVAL_EdgesD28_30.csv", row.names = FALSE)
write.csv(clean_nodesD0_30, "Clean_CARNIVAL_NodesD28_30.csv", row.names = FALSE)

cytoscapePing() 

# Rename columns to standard RCy3 formats to be safe
# (RCy3 prefers 'source', 'target', 'interaction' but handles others too)
edges_for_cyD0_30 <- clean_edgesD0_30
colnames(edges_for_cyD0_30) <- c("source", "sign", "target", "weight") # Adjust based on your actual column order!

nodes_for_cyD0_30 <- clean_nodesD0_30
colnames(nodes_for_cyD0_30)[1] <- "id" # Ensure first column is 'id'

# Send to Cytoscape
createNetworkFromDataFrames(nodes = nodes_for_cyD0_30, 
                            edges = edges_for_cyD0_30, 
                            title = "CARNIVAL_D0_top30", 
                            collection = "CARNIVAL_Analysis_D0_top30")

carnival_D5 <- runCARNIVAL(
  inputObj  = iniciators,             
  measObj   = TFscoresD5_CARNIVAL,             # TF activity measurements
  netObj    = sif,                    # prior knowledge network
  weightObj = Pathway_CARNIVALD5,         # pathway activity scores (or progeny_matrix)
  solverPath = "C:/Program Files/IBM/ILOG/CPLEX_Studio2211/cplex/bin/x64_win64/cplex.exe", 
  solver     = "cplex",
  timelimit  = 7200, 
  mipGAP     = 0.05,
  poolrelGAP = 0
)
saveRDS(carnival_D5, "~/hnrnpu-causal-multiomics/processeddata/carnival_D5.rds")
carnival_D5 <- readRDS("~/hnrnpu-causal-multiomics/processeddata/carnival_D5.rds")
clean_edgesD5 <- carnival_D5$weightedSIF %>%
  filter(Weight > 0)

active_node_namesD5 <- unique(c(clean_edgesD5$Node1, clean_edgesD5$Node2))

clean_nodesD5 <- carnival_D5$nodesAttributes %>%
  filter(Node %in% active_node_namesD5)

###### visual ######
library(RCy3)
library(igraph)
library(dplyr)
write.csv(clean_edgesD5, "Clean_CARNIVAL_EdgesD5.csv", row.names = FALSE)
write.csv(clean_nodesD5, "Clean_CARNIVAL_NodesD5.csv", row.names = FALSE)

cytoscapePing() 

# Rename columns to standard RCy3 formats to be safe
# (RCy3 prefers 'source', 'target', 'interaction' but handles others too)
edges_for_cyD5 <- clean_edgesD5
colnames(edges_for_cyD5) <- c("source", "sign", "target", "weight") # Adjust based on your actual column order!

nodes_for_cyD5 <- clean_nodesD5
colnames(nodes_for_cyD5)[1] <- "id" # Ensure first column is 'id'

# Send to Cytoscape
createNetworkFromDataFrames(nodes = nodes_for_cyD5, 
                            edges = edges_for_cyD5, 
                            title = "CARNIVAL_D5", 
                            collection = "CARNIVAL_AnalysisD5")

setNodeColorMapping(table.column = "AvgAct", 
                    table.column.values = c(-100, 0, 100), 
                    colors = c("#0000FF", "#FFFFFF", "#FF0000"), 
                    mapping.type = "continuous")

## top 30
top30tfs_D5 <- TFscores_D5 %>%
  arrange(desc(abs(score))) %>%
  slice_head(n = 30)
head(top30tfs_D5)

top30tfs_D5 <- top30tfs_D5 %>%
  dplyr::select(source, score) %>%
  deframe()

carnival_D5top30 <- runCARNIVAL(
  inputObj  = iniciators,             
  measObj   = top30tfs_D5,             # TF activity measurements
  netObj    = sif,                    # prior knowledge network
  weightObj = Pathway_CARNIVALD5,         # pathway activity scores (or progeny_matrix)
  solverPath = "C:/Program Files/IBM/ILOG/CPLEX_Studio2211/cplex/bin/x64_win64/cplex.exe", 
  solver     = "cplex",
  timelimit  = 900, 
  mipGAP     = 0.05,
  poolrelGAP = 0
)
saveRDS(carnival_D5top30,"~/hnrnpu-causal-multiomics/processeddata/carnival_D5top30.rds")
carnival_D5top30 <- readRDS("~/hnrnpu-causal-multiomics/processeddata/carnival_D5top30.rds")

clean_edgesD5_30 <- carnival_D5top30$weightedSIF %>%
  filter(Weight > 0)

active_node_namesD5_30 <- unique(c(clean_edgesD5_30$Node1, clean_edgesD5_30$Node2))

clean_nodesD5_30 <- carnival_D5top30$nodesAttributes %>%
  filter(Node %in% active_node_namesD5_30)

write.csv(clean_edgesD5_30, "Clean_CARNIVAL_EdgesD28_30.csv", row.names = FALSE)
write.csv(clean_nodesD5_30, "Clean_CARNIVAL_NodesD28_30.csv", row.names = FALSE)

library(visNetwork)
source("~/hnrnpu-causal-multiomics/carnival_visNet.R")
carnival_visNet(evis = carnival_D5top30$weightedSIF,
                nvis = carnival_D5top30$nodesAttributes)

cytoscapePing() 

# Rename columns to standard RCy3 formats to be safe
# (RCy3 prefers 'source', 'target', 'interaction' but handles others too)
edges_for_cyD5_30 <- clean_edgesD5_30
colnames(edges_for_cyD5_30) <- c("source", "sign", "target", "weight") # Adjust based on your actual column order!

nodes_for_cyD5_30 <- clean_nodesD5_30
colnames(nodes_for_cyD5_30)[1] <- "id" # Ensure first column is 'id'

# Send to Cytoscape
createNetworkFromDataFrames(nodes = nodes_for_cyD5_30, 
                            edges = edges_for_cyD5_30, 
                            title = "CARNIVAL_D5_top30", 
                            collection = "CARNIVAL_Analysis_D5_top30")

carnival_D28 <- runCARNIVAL(
  inputObj  = iniciators,             
  measObj   = TFscoresD28_CARNIVAL,             # TF activity measurements
  netObj    = sif,                    # prior knowledge network
  weightObj = Pathway_CARNIVALD28,         # pathway activity scores (or progeny_matrix)
  solverPath = "C:/Program Files/IBM/ILOG/CPLEX_Studio2211/cplex/bin/x64_win64/cplex.exe", 
  solver     = "cplex",
  timelimit  = 7200, 
  mipGAP     = 0.05,
  poolrelGAP = 0
)

clean_edgesD28 <- carnival_D28$weightedSIF %>%
  filter(Weight > 0)

active_node_namesD28 <- unique(c(clean_edgesD28$Node1, clean_edgesD28$Node2))

clean_nodesD28 <- carnival_D28$nodesAttributes %>%
  filter(Node %in% active_node_namesD28)

###### visual ######
library(RCy3)
library(igraph)
library(dplyr)
write.csv(clean_edgesD28, "Clean_CARNIVAL_EdgesD28.csv", row.names = FALSE)
write.csv(clean_nodesD28, "Clean_CARNIVAL_NodesD28.csv", row.names = FALSE)

cytoscapePing() 

# Rename columns to standard RCy3 formats to be safe
# (RCy3 prefers 'source', 'target', 'interaction' but handles others too)
edges_for_cyD28 <- clean_edgesD28
colnames(edges_for_cyD28) <- c("source", "sign", "target", "weight") # Adjust based on your actual column order!

nodes_for_cyD28 <- clean_nodesD28
colnames(nodes_for_cyD28)[1] <- "id" # Ensure first column is 'id'

# Send to Cytoscape
createNetworkFromDataFrames(nodes = nodes_for_cyD28, 
                            edges = edges_for_cyD28, 
                            title = "CARNIVAL_Automated", 
                            collection = "CARNIVAL_Analysis")

setNodeColorMapping(table.column = "AvgAct", 
                    table.column.values = c(-100, 0, 100), 
                    colors = c("#0000FF", "#FFFFFF", "#FF0000"), 
                    mapping.type = "continuous")

##### top 30 TFs CARNIVAL #####
top30tfs_D28 <- TFscores_D28 %>%
  arrange(desc(abs(score))) %>%
  slice_head(n = 30)
head(top30tfs_D28)

top30tfs_D28 <- top30tfs_D28 %>%
  dplyr::select(source, score) %>%
  deframe()

carnival_D28top30 <- runCARNIVAL(
  inputObj  = iniciators,             
  measObj   = top30tfs_D28,             # TF activity measurements
  netObj    = sif,                    # prior knowledge network
  weightObj = Pathway_CARNIVALD28,         # pathway activity scores (or progeny_matrix)
  solverPath = "C:/Program Files/IBM/ILOG/CPLEX_Studio2211/cplex/bin/x64_win64/cplex.exe", 
  solver     = "cplex",
  timelimit  = 900, 
  mipGAP     = 0.05,
  poolrelGAP = 0
)

#accidentally saved as carnival_D28 
saveRDS(carnival_D28, "~/hnrnpu-causal-multiomics/processeddata/carnival_D28top30.rds")
carnival_D28top30 <- readRDS("~/hnrnpu-causal-multiomics/processeddata/carnival_D28top30.rds")
clean_edgesD28_30 <- carnival_D28top30$weightedSIF %>%
  filter(Weight > 0)

active_node_namesD28_30 <- unique(c(clean_edgesD28_30$Node1, clean_edgesD28_30$Node2))

clean_nodesD28_30 <- carnival_D28top30$nodesAttributes %>%
  filter(Node %in% active_node_namesD28_30)

write.csv(clean_edgesD28_30, "Clean_CARNIVAL_EdgesD28_30.csv", row.names = FALSE)
write.csv(clean_nodesD28_30, "Clean_CARNIVAL_NodesD28_30.csv", row.names = FALSE)

cytoscapePing() 

# Rename columns to standard RCy3 formats to be safe
# (RCy3 prefers 'source', 'target', 'interaction' but handles others too)
edges_for_cyD28_30 <- clean_edgesD28_30
colnames(edges_for_cyD28_30) <- c("source", "sign", "target", "weight") # Adjust based on your actual column order!

nodes_for_cyD28_30 <- clean_nodesD28_30
colnames(nodes_for_cyD28_30)[1] <- "id" # Ensure first column is 'id'

# Send to Cytoscape
createNetworkFromDataFrames(nodes = nodes_for_cyD28_30, 
                            edges = edges_for_cyD28_30, 
                            title = "CARNIVAL_D28_top30", 
                            collection = "CARNIVAL_Analysis_D28_top30")

#####visualattempts#####
#CARNIVAL visualisation trial

#gemini approach
BiocManager::install("RCy3")

library(RCy3)
library(igraph)

library(dplyr)

clean_edges <- carnival_D0$weightedSIF %>%
  filter(Weight > 0)
head(clean_edges)
# 3. Filter Nodes to match
# We grab 'nodesAttributes' and keep only nodes present in our active edges
active_node_names <- unique(c(clean_edges$Node1, clean_edges$Node2))

clean_nodes <- carnival_D0$nodesAttributes %>%
  filter(Node %in% active_node_names)

# 4. Verify columns exist (You should see "Sign" here)
print(colnames(clean_edges))
print(paste("Final Edge Count:", nrow(clean_edges)))
write.csv(clean_edges, "Clean_CARNIVAL_Edges.csv", row.names = FALSE)
write.csv(clean_nodes, "Clean_CARNIVAL_Nodes.csv", row.names = FALSE)

#rcy3
library(RCy3)

# ---------------------------------------------------------
# STEP 1: Connect and Create Network
# ---------------------------------------------------------
# Make sure Cytoscape is open!
cytoscapePing() 

# Rename columns to standard RCy3 formats to be safe
# (RCy3 prefers 'source', 'target', 'interaction' but handles others too)
edges_for_cy <- clean_edges
colnames(edges_for_cy) <- c("source", "sign", "target", "weight") # Adjust based on your actual column order!

nodes_for_cy <- clean_nodes
colnames(nodes_for_cy)[1] <- "id" # Ensure first column is 'id'

# Send to Cytoscape
createNetworkFromDataFrames(nodes = nodes_for_cy, 
                            edges = edges_for_cy, 
                            title = "CARNIVAL_Automated", 
                            collection = "CARNIVAL_Analysis")

# ---------------------------------------------------------
# STEP 2: Style the Nodes (Activity Color)
# ---------------------------------------------------------
# Map 'activity_D0' to Fill Color: Blue (-100) -> White (0) -> Red (100)
setNodeColorMapping(table.column = "AvgAct", 
                    table.column.values = c(-100, 0, 100), 
                    colors = c("#0000FF", "#FFFFFF", "#FF0000"), 
                    mapping.type = "continuous")

# Lock node size (optional, makes it cleaner)
lockNodeDimensions(TRUE)
setNodeSizeDefault(40)

# ---------------------------------------------------------
# STEP 3: Style the Edges (Signs & Weights)
# ---------------------------------------------------------
# 1. Map 'sign' to Arrow Shape
# 1 = Delta (Arrow), -1 = T (Inhibitor)
setEdgeTargetArrowShapeMapping(table.column = "sign", 
                               table.column.values = c(1, -1), 
                               shapes = c("DELTA", "T"), 
                               mapping.type = "discrete")

# 2. Map 'sign' to Color (Blue=Activate, Red=Inhibit)
setEdgeColorMapping(table.column = "sign", 
                    table.column.values = c(1, -1), 
                    colors = c("#3399FF", "#FF3333"), 
                    mapping.type = "discrete")

# 3. Map 'weight' to Thickness
setEdgeLineWidthMapping(table.column = "weight", 
                        table.column.values = c(1, 100), 
                        widths = c(1, 5), 
                        mapping.type = "continuous")

# ---------------------------------------------------------
# STEP 4: Apply Layout
# ---------------------------------------------------------
# Force-directed layout to untangle the hairball
layoutNetwork(layout.name = "force-directed")

# Optional: If you have yFiles installed
layoutNetwork(layout.name = "yfiles-organic")

##### trial again 

library(dplyr)
library(stringr)

edges_D0 <- carnival_D0$weightedSIF %>%
  mutate(
    Node1 = str_trim(Node1),
    Node2 = str_trim(Node2)
  ) %>%
  filter(abs(Weight) >= 0.3)   # try 0.3–0.7 if needed

edges_D0 <- edges_D0 %>%
  filter(Node1 != "Perturbation" & Node2 != "Perturbation")
#### d28 

carnival_D28 <- readRDS("~/hnrnpu-causal-multiomics/processeddata/carnival_D28.rds")
carnival_res1 <- readRDS("~/hnrnpu-causal-multiomics/processeddata/Previous Trials/carnival_result1.rds")
edges <- carnival_res1$weightedSIF %>%
  mutate(
    Node1 = str_trim(Node1),
    Node2 = str_trim(Node2)
  ) %>%
  filter(abs(Weight) >= 0.3) 

#CARNIVAL visualisation trial

#gemini approach
BiocManager::install("RCy3")

library(RCy3)
library(igraph)

library(dplyr)

clean_edgesD28 <- carnival_D28$weightedSIF %>%
  filter(Weight > 0)
head(clean_edgesD28)
# 3. Filter Nodes to match
# We grab 'nodesAttributes' and keep only nodes present in our active edges
active_node_names <- unique(c(clean_edges$Node1, clean_edges$Node2))

clean_nodes <- carnival_D0$nodesAttributes %>%
  filter(Node %in% active_node_names)

# 4. Verify columns exist (You should see "Sign" here)
print(colnames(clean_edges))
print(paste("Final Edge Count:", nrow(clean_edges)))
write.csv(clean_edges, "Clean_CARNIVAL_Edges.csv", row.names = FALSE)
write.csv(clean_nodes, "Clean_CARNIVAL_Nodes.csv", row.names = FALSE)

#rcy3
library(RCy3)

# ---------------------------------------------------------
# STEP 1: Connect and Create Network
# ---------------------------------------------------------
# Make sure Cytoscape is open!
cytoscapePing() 

# Rename columns to standard RCy3 formats to be safe
# (RCy3 prefers 'source', 'target', 'interaction' but handles others too)
edges_for_cy <- clean_edges
colnames(edges_for_cy) <- c("source", "sign", "target", "weight") # Adjust based on your actual column order!

nodes_for_cy <- clean_nodes
colnames(nodes_for_cy)[1] <- "id" # Ensure first column is 'id'

# Send to Cytoscape
createNetworkFromDataFrames(nodes = nodes_for_cy, 
                            edges = edges_for_cy, 
                            title = "CARNIVAL_Automated", 
                            collection = "CARNIVAL_Analysis")

# ---------------------------------------------------------
# STEP 2: Style the Nodes (Activity Color)
# ---------------------------------------------------------
# Map 'activity_D0' to Fill Color: Blue (-100) -> White (0) -> Red (100)
setNodeColorMapping(table.column = "AvgAct", 
                    table.column.values = c(-100, 0, 100), 
                    colors = c("#0000FF", "#FFFFFF", "#FF0000"), 
                    mapping.type = "continuous")

# Lock node size (optional, makes it cleaner)
lockNodeDimensions(TRUE)
setNodeSizeDefault(40)

# ---------------------------------------------------------
# STEP 3: Style the Edges (Signs & Weights)
# ---------------------------------------------------------
# 1. Map 'sign' to Arrow Shape
# 1 = Delta (Arrow), -1 = T (Inhibitor)
setEdgeTargetArrowShapeMapping(table.column = "sign", 
                               table.column.values = c(1, -1), 
                               shapes = c("DELTA", "T"), 
                               mapping.type = "discrete")

# 2. Map 'sign' to Color (Blue=Activate, Red=Inhibit)
setEdgeColorMapping(table.column = "sign", 
                    table.column.values = c(1, -1), 
                    colors = c("#3399FF", "#FF3333"), 
                    mapping.type = "discrete")

# 3. Map 'weight' to Thickness
setEdgeLineWidthMapping(table.column = "weight", 
                        table.column.values = c(1, 100), 
                        widths = c(1, 5), 
                        mapping.type = "continuous")

# ---------------------------------------------------------
# STEP 4: Apply Layout
# ---------------------------------------------------------
# Force-directed layout to untangle the hairball
layoutNetwork(layout.name = "force-directed")

# Optional: If you have yFiles installed
layoutNetwork(layout.name = "yfiles-organic")

##### CARNIVAL ATTEMPTS 07.01.26 #####
# Extract the Node Attributes
nodesd5top30 <- carnival_D5top30$nodesAttributes

# Filter for nodes that act as initiators (Activity inferred as 1 or -1)
# Note: In your inputs, initiators were nodes with no upstream parents in the network
inferred_driversd5 <- nodesd5top30 %>%
  filter(Node %in% colnames(iniciators)) %>% # Filter only the nodes you provided as potential starts
  filter(AvgAct != 0) %>% # Keep only those CARNIVAL found to be active
  arrange(desc(abs(AvgAct))) # Sort by activity strength

print(inferred_driversd5)

library(igraph)

# 1. Create an igraph object from your clean edges
# (Using the clean_edgesD5_30 you created earlier)
net_igraphd5 <- graph_from_data_frame(d = clean_edgesD5_30, directed = TRUE)

# 2. Calculate "Betweenness Centrality" (Who controls the flow?)
node_betweennessd5 <- betweenness(net_igraphd5, directed = TRUE)
top_hubsd5 <- sort(node_betweennessd5, decreasing = TRUE)[1:10]

# 3. Calculate "Degree" (Who has the most connections?)
node_degreed5 <- degree(net_igraphd5, mode = "all")
top_connectedd5 <- sort(node_degreed5, decreasing = TRUE)[1:10]

print(top_hubsd5)

BiocManager::install("clusterProfiler")
library(clusterProfiler)
library(org.Hs.eg.db) # Assuming human data

# Convert symbols to Entrez IDs for the tool
gene_list <- clean_nodesD5_30$Node
gene.df <- bitr(gene_list, fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = org.Hs.eg.db)

# Run KEGG Enrichment
kegg_res <- enrichKEGG(gene = gene.df$ENTREZID,
                       organism = 'hsa',
                       pvalueCutoff = 0.05)

# Visualize
dotplot(kegg_res, showCategory=15)

