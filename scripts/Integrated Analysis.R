##### integrated analysis trials #####

library(pheatmap)

#Consolidate your PROGENy scores
all_pathways <- unique(c(names(PathwayActivity_CARNIVALinputD0only), names(PathwayActivity_CARNIVALinputD5only), names(PathwayActivity_CARNIVALinputD28only)))

# Create a matrix
prog_matrix <- data.frame(
  D0 = PathwayActivity_CARNIVALinputD0only[all_pathways],
  D5 = PathwayActivity_CARNIVALinputD5only[all_pathways],   # Assuming you have this
  D28 = PathwayActivity_CARNIVALinputD28only[all_pathways]
)


library(pheatmap)
library(RColorBrewer)


# 1. Clean the Data
# Extract only the score columns (2, 4, 6)
clean_matrix <- prog_matrix[, c(2, 4, 6)]

# Rename the columns to simply "D0", "D5", "D28"
colnames(clean_matrix) <- c("D0", "D5", "D28")

# Set the row names using the first column (Pathway names)
rownames(clean_matrix) <- prog_matrix[, 1]

# Ensure it is a matrix (not a dataframe) for pheatmap
clean_matrix <- as.matrix(clean_matrix)

# 1. Define a professional color palette
# Rev(RdBu) gives: Blue (Low) -> White (Zero) -> Red (High)
my_colors <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)

# 2. Define Breaks to ensure 0 is White
# This locks the scale from -1 to 1, even if your data only goes to 0.9
my_breaks <- seq(-1, 1, length.out = 101)

# 3. Create the Plot
pheatmap(clean_matrix,
         # Clustering
         cluster_cols = FALSE,      # Keep Time ordered
         cluster_rows = TRUE,       # Cluster pathways to show patterns
         treeheight_row = 20,       # Make the dendrogram subtle
         
         # Colors & Scale
         color = my_colors,
         breaks = my_breaks,
         border_color = "white",    # Thin white border separates cells cleanly
         
         # Text & Labels
         display_numbers = FALSE,   # Turn off numbers for a cleaner look (optional)
         fontsize = 12,             # Base font size
         fontsize_row = 10,         # Row label size
         fontsize_col = 12,         # Column label size
         angle_col = 0,             # Keep column labels horizontal
         
         # Dimensions (Adjust these if cells look squashed)
         cellwidth = 40,
         cellheight = 15,
         
         # Titles
         main = "PROGENy Pathway Activity",
         
         # Save directly to high-res PDF (Optional)
  filename = "Figure1_Pathway_Heatmap.pdf", 
  width = 5, height = 6
)

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

# ==========================================
# 1. PREPARE THE DATA
# ==========================================
# extract the list of Node IDs from your dataframes
# Note: Check if your column is named "nodes", "Node", or "id". 
# I am assuming it is "nodes" based on typical CARNIVAL output.
list_d0  <- unique(clean_nodesD0_100$Node)
list_d5  <- unique(clean_nodesD5_100$Node)
list_d28 <- unique(clean_nodesD28_100$Node)

# Create a named list for the plotting function
venn_data <- list(
  "Day 0" = list_d0,
  "Day 5" = list_d5,
  "Day 28" = list_d28
)

# ==========================================
# 2. IDENTIFY UNIQUE & SHARED NODES
# ==========================================

# A. The "Core Machinery" (Present in ALL 3 days)
common_all <- Reduce(intersect, venn_data)

# B. Unique to Each Day (The "Stage-Specific" Drivers)
unique_d0  <- setdiff(list_d0,  union(list_d5, list_d28))
unique_d5  <- setdiff(list_d5,  union(list_d0, list_d28))
unique_d28 <- setdiff(list_d28, union(list_d0, list_d5))

# C. Present in D5 and D28 but NOT D0 (The "Differentiation" module)
shared_d5_d28 <- setdiff(intersect(list_d5, list_d28), list_d0)

# ==========================================
# 3. PRINT RESULTS TO CONSOLE
# ==========================================
cat("\n=== NETWORK OVERLAP SUMMARY ===\n")
cat("Total Nodes in Day 0: ", length(list_d0), "\n")
cat("Total Nodes in Day 5: ", length(list_d5), "\n")
cat("Total Nodes in Day 28:", length(list_d28), "\n\n")

cat("--- CORE NODES (Preserved across all stages) ---\n")
print(common_all)

cat("\n--- UNIQUE NODES (Specific to Day 0) ---\n")
print(unique_d0)

cat("\n--- UNIQUE NODES (Specific to Day 5) ---\n")
print(unique_d5)

cat("\n--- UNIQUE NODES (Specific to Day 28) ---\n")
print(unique_d28)

# ==========================================
# 4. VISUALIZE (Venn Diagram)
# ==========================================
# This creates a nice plot showing the overlaps
ggvenn(
  venn_data, 
  fill_color = c("#E41A1C", "#377EB8", "#4DAF4A"), # Red, Blue, Green
  stroke_size = 0.5, 
  set_name_size = 5
) +
  ggtitle("Overlap of Signaling Nodes Across Differentiation") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))

# Save the plot
ggsave("Node_Overlap_Venn.png", width = 8, height = 8, dpi = 300)

# ==========================================
# 5. EXPORT LISTS (Optional)
# ==========================================
# If you want to save these lists to a CSV to browse later:
max_len <- max(length(common_all), length(unique_d0), length(unique_d5), length(unique_d28))

# Helper function to pad lists with NA
pad_na <- function(x, n) { c(x, rep(NA, n - length(x))) }

export_df <- data.frame(
  Common_Nodes = pad_na(common_all, max_len),
  Unique_D0    = pad_na(unique_d0, max_len),
  Unique_D5    = pad_na(unique_d5, max_len),
  Unique_D28   = pad_na(unique_d28, max_len),
  Shared_d5_d28 = pad_na(shared_d5_d28, max_len)
)

write.csv(export_df, "Node_Comparison_Lists.csv", row.names = FALSE)

####### tf activities plot #####
library(tidyverse)

# 1. PREPARE THE DATA
# -------------------
# I assume you have 3 dataframes: top100tf_D0, top100tf_D5, top100tf_D28
# Each has columns: 'TF', 'score', 'p_value'

top10_D0 <- TFscores_D0 %>%
  arrange(desc(abs(score))) %>%
  slice_head(n = 10)
head(top10_D0)

top10_D0 <- top10_D0 %>%
  dplyr::select(TF = source, score, p_value) %>%
  dplyr::mutate(timepoint = "D0")

top10_D5 <- TFscores_D5 %>%
  arrange(desc(abs(score))) %>%
  slice_head(n = 10)
head(top10_D5)

top10_D5 <- top10_D5 %>%
  dplyr::select(TF = source, score, p_value) %>%
  dplyr::mutate(timepoint = "D5")

top10_D28 <- TFscores_D28 %>%
  arrange(desc(abs(score))) %>%
  slice_head(n = 10)
head(top10_D28)

top10_D28 <- top10_D28 %>%
  dplyr::select(TF = source, score, p_value) %>%
  dplyr::mutate(timepoint = "D28")

# Combine them into one large dataframe
df_all <- bind_rows(top10_D0, top10_D5, top10_D28)

# 2. PROCESS FOR PLOTTING
# -----------------------
df_plot <- df_all %>%
  # Create the inverse P-value score for dot size
  # We use -log10 so that smaller P-values become larger numbers
  dplyr::mutate(log_p = -log10(p_value)) %>%
  
  # Ensure the timepoints appear in the correct chronological order (not alphabetical)
  dplyr::mutate(timepoint = factor(timepoint, levels = c("D0", "D5", "D28")))

# Optional: If plotting 100 TFs is too crowded, you can filter for the top variable ones here.
# For now, we plot all of them.
limit <- max(abs(df_plot$score)) * c(-1, 1)
# 3. GENERATE THE BALLOON PLOT
# ----------------------------
TFplot <- ggplot(df_plot, aes(x = timepoint, y = reorder(TF, score))) +
  geom_point(aes(color = score, size = log_p)) +
  
  # COLOR: Diverging scale (Blue = Repressed, Red = Activated)
  scale_color_gradientn(
    colors = my_colors,
    limits = limit,    # Forces 0 to be the exact center
    name = "Activity\nScore"
  ) +
  
  # SIZE: Significance
  scale_size_continuous(
    name = "-log10(P-val)",
    range = c(5, 10) # Adjust these numbers to make dots bigger/smaller
  ) +
  
  # THEME & LABELS
  theme_minimal(base_size = 12) +
  labs(
    title = "Differential TF Activity (decoupleR)",
    x = "Timepoint",
    y = "Transcription Factor"
  ) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 10), # Make text smaller if you have 100 TFs
    panel.grid.major.x = element_blank()  # Clean up vertical grid lines
  )

ggsave("TFdotplot.png", width = 8, height = 10, dpi = 300)