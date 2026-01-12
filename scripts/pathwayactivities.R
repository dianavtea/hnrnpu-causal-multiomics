#inferring pathway activities using PROGENy for day 0, 5 and 28

##### libraries #####
library(progeny)
library(dorothea)
library(tibble)
library(tidyr)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(readr)
library(RColorBrewer)

##### PROGENy #####
##### D0 #####
## We read the normalised counts and the experimental design 
nc_D0only <- read_csv("~/hnrnpu-causal-multiomics/processeddata/VST_progD0_clean.csv") 
#vst_progD0only
head(nc_D0only)

metafiltD0 <- meta_filtered[c(1:5, 9:13),]

## We read the results from the differential analysis. 
DED0only <- read_csv("~/hnrnpu-causal-multiomics/DESeq2_res_D0only_unique.csv") #17017 genes
#res_D0_only.annot

nc_D0only_mat <- as.matrix(
  data.frame(nc_D0only, row.names = nc_D0only$...1)[, -1]
)
head(nc_D0only_mat)
# literally same as just vst_progD0only

#clean so hgnc symbols all recognised properly
DED0only <- DED0only %>%
  mutate(hgnc_symbol = as.character(hgnc_symbol)) %>%
  mutate(hgnc_symbol = trimws(hgnc_symbol))

#format for progeny
DED0only_mat <- as.matrix(DED0only$stat)
rownames(DED0only_mat) <- DED0only$hgnc_symbol
colnames(DED0only_mat) <- "statistic"
head(DED0only_mat)

#progeny results as input for CARNIVAL
#CARNIVAL sets weights based on PROGENy scores in each pathway-related node
#to find more relevant solutions
#set progeny zscores false to return pathway activity values between 1 and -1

PathwayActivity_CARNIVALinputD0only <- progeny(DED0only_mat, 
                                               scale=TRUE, organism="Human", top = 100, perm = 10000, z_scores = FALSE) %>%
  t () %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "Pathway") 
colnames(PathwayActivity_CARNIVALinputD0only)[2] <- "score"
head(PathwayActivity_CARNIVALinputD0only)
write_csv(PathwayActivity_CARNIVALinputD0only, 
          "~/hnrnpu-causal-multiomics/processeddata/PathwayActivity_CARNIVALinputD0.csv")


##### D5 only #####
## We read the normalised counts and the experimental design 
nc_D5only <- read_csv("~/hnrnpu-causal-multiomics/processeddata/VST_progD5only_clean.csv") 
#vst_D5only_clean
head(nc_D5only)

metafiltD5 <- meta_filtered[c(6:8, 14:20), ]

## We read the results from the differential analysis. 
DED5only <- read_csv("~/hnrnpu-causal-multiomics/DESeq2_res_D5only_unique.csv") #21262 genes
#res_D5_only.annot

nc_D5only_mat <- as.matrix(
  data.frame(nc_D5only, row.names = nc_D5only$...1)[, -1]
)
head(nc_D5only_mat)
# literally same as just vst_progD5only

#clean so hgnc symbols all recognised properly
DED5only <- DED5only %>%
  mutate(hgnc_symbol = as.character(hgnc_symbol)) %>%
  mutate(hgnc_symbol = trimws(hgnc_symbol))

#format for progeny
DED5only_mat <- as.matrix(DED5only$stat)
rownames(DED5only_mat) <- DED5only$hgnc_symbol
colnames(DED5only_mat) <- "statistic"
head(DED5only_mat)

#progeny results as input for CARNIVAL
#CARNIVAL sets weights based on PROGENy scores in each pathway-related node
#to find more relevant solutions
#set progeny zscores false to return pathway activity values between 1 and -1

PathwayActivity_CARNIVALinputD5only <- progeny(DED5only_mat, 
                                                scale=TRUE, organism="Human", top = 100, perm = 10000, z_scores = FALSE) %>%
  t () %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "Pathway") 
colnames(PathwayActivity_CARNIVALinputD5only)[2] <- "score"
head(PathwayActivity_CARNIVALinputD5only)

write_csv(PathwayActivity_CARNIVALinputD5only, 
          "~/hnrnpu-causal-multiomics/processeddata/PathwayActivity_CARNIVALinputD5.csv")

##### D28 only #####
## We read the normalised counts and the experimental design 
nc_D28only <- read_csv("~/hnrnpu-causal-multiomics/processeddata/VST_progD28only.csv") 
#vst_progD28only
head(nc_D28only)

metafiltD28 <- meta_filtered[21:26,]

## We read the results from the differential analysis. 
#DED28only <- read_csv("~/hnrnpu-causal-multiomics/files/DESeq2_res_D28_only.csv") #21262 
#res_D28_only.annot

### ERROR: wrong file chosen above, do this one instead with 17017:
DED28only <- read_csv("~/hnrnpu-causal-multiomics/DESeq2_res_D28only_unique.csv")

nc_D28only_mat <- as.matrix(
  data.frame(nc_D28only, row.names = nc_D28only$...1)[, -1]
)
head(nc_D28only_mat)
#vst_progD28only

#clean so hgnc symbols all recognised properly
DED28only <- DED28only %>%
  mutate(hgnc_symbol = as.character(hgnc_symbol)) %>%
  mutate(hgnc_symbol = trimws(hgnc_symbol))

#format for progeny
DED28only_mat <- as.matrix(DED28only$stat)
rownames(DED28only_mat) <- DED28only$hgnc_symbol
colnames(DED28only_mat) <- "statistic"
head(DED28only_mat)

#progeny results as input for CARNIVAL
#CARNIVAL sets weights based on PROGENy scores in each pathway-related node
#to find more relevant solutions
#set progeny zscores false to return pathway activity values between 1 and -1

PathwayActivity_CARNIVALinputD28only <- progeny(DED28only_mat, 
                                                scale=TRUE, organism="Human", top = 100, perm = 10000, z_scores = FALSE) %>%
  t () %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "Pathway") 
colnames(PathwayActivity_CARNIVALinputD28only)[2] <- "score"
head(PathwayActivity_CARNIVALinputD28only)

write_csv(PathwayActivity_CARNIVALinputD28only, 
          "~/hnrnpu-causal-multiomics/processeddata/PathwayActivity_CARNIVALinputD28.csv")

##### Visualisation integrated #####

#PROGENy scores compilation
all_pathways <- unique(c(names(PathwayActivity_CARNIVALinputD0only), names(PathwayActivity_CARNIVALinputD5only), names(PathwayActivity_CARNIVALinputD28only)))

# Create a matrix
prog_matrix <- data.frame(
  D0 = PathwayActivity_CARNIVALinputD0only[all_pathways],
  D5 = PathwayActivity_CARNIVALinputD5only[all_pathways],   
  D28 = PathwayActivity_CARNIVALinputD28only[all_pathways]
)

#Clean the Data
clean_matrix <- prog_matrix[, c(2, 4, 6)]

# Rename the columns to simply "D0", "D5", "D28"
colnames(clean_matrix) <- c("D0", "D5", "D28")
# Set the row names using the first column (Pathway names)
rownames(clean_matrix) <- prog_matrix[, 1]

# Ensure it is a matrix (not a dataframe) for pheatmap
clean_matrix <- as.matrix(clean_matrix)

# Rev(RdBu) gives: Blue (Low) -> White (Zero) -> Red (High)
my_colors <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)
my_breaks <- seq(-1, 1, length.out = 101)

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
