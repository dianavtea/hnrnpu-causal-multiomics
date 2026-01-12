library(OmnipathR)
interaction_datasets()

library(readr)
library(ggplot2)
library(reshape)
library(pheatmap)
library(gridExtra)
library(grid)
library(cowplot)
library(ggrepel)
library(hexbin)
library(DESeq2)
library(ggpubr)
library(biomaRt)
library(RSQLite)
library(org.Hs.eg.db)
library(readxl)
library(dplyr)

source("support_functions.R") 
merged_gene_counts <- as.data.frame(
  read.delim("~/hnrnpu-causal-multiomics/rawdata/merged_gene_counts 1.txt", header=TRUE, stringsAsFactors = FALSE,row.names=1))
#63677 transcripts

#simple summary of experimental design
targets <- as.data.frame(matrix(NA, length(names(merged_gene_counts)),2))
names(targets) <- c("sample", "condition")
targets$sample <- names(merged_gene_counts)
targets$condition <- c("CTRL", "CTRL", "CTRL", "ASD", "ASD", "ASD")
head(targets)

#merged gene counts without NA
# Keep genes expressed (>50 counts) per group (CTRL vs ASD)
gene_counts_CTRL <- rowSums(merged_gene_counts[,1:3] > 50) >= 3
gene_counts_ASD <- rowSums(merged_gene_counts[,4:6] > 50) >= 3
merged_gene_counts_filtered <- merged_gene_counts[gene_counts_CTRL | gene_counts_ASD, ] #13821 genes
any(is.na(merged_gene_counts_filtered)) #FALSE, no NA

dds <- DESeqDataSetFromMatrix(countData = merged_gene_counts_filtered,
                               colData = targets,
                               design= ~ condition)

dds$condition<- relevel( dds$condition, "CTRL" )
dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds)

res.mat <- cbind(counts(dds, normalized=TRUE), res$log2FoldChange, res$padj)
mcols(res, use.names = TRUE)
saveRDS(dds, file="dds_DEG_ASDvsCTRL_D28.rds")
head(res)

##### annotations #####
#Connect to GRCh37 Ensembl
# Use rownames of res directly
ensembl_ids <- rownames(res)

mart <- useEnsembl(
  biomart = "ensembl",
  dataset = "hsapiens_gene_ensembl",
  GRCh = 37
)

annot <- getBM(
  attributes = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = ensembl_ids,
  mart = mart
)

res.df <- as.data.frame(res)
res.df$ensembl <- ensembl_ids

res.annot <- merge(
  res.df,
  annot,
  by.x = "ensembl",
  by.y = "ensembl_gene_id",
  all.x = TRUE,
  sort = FALSE
)

res.annot <- res.annot[match(ensembl_ids, res.annot$ensembl), ]

colnames(res.annot)[colnames(res.annot) == "entrezgene_id"] <- "entrez"
colnames(res.annot)[colnames(res.annot) == "hgnc_symbol"]   <- "hgnc_symbol"
head(res.annot)
write.csv(res.annot, "DESeq2_results50.csv", row.names = FALSE)

##### VST #####
vst_dds <- vst(dds, blind = FALSE)
vst_mat <- assay(vst_dds)
#82926 elements
dim(vst_mat) #13821 transcripts x 6 samples

# Ensure rows of res.annot match vst_mat
res.annot <- res.annot[match(rownames(vst_mat), res.annot$ensembl), ]

vst_annot <- cbind(res.annot[, c("ensembl", "entrez", "hgnc_symbol")], vst_mat)

write.csv(vst_annot, "VST_counts_with_annot50.csv", row.names = FALSE)

#vst for progeny
vst_progeny <- vst_mat
rownames(vst_progeny) <- res.annot$hgnc_symbol
head(vst_progeny)
any(is.na(vst_progeny)) #false

rownames(vst_progeny)[duplicated(rownames(vst_progeny))] 
vst_progeny <- vst_progeny[rownames(vst_progeny) != "", ] #remove blanks
vst_progeny <- vst_progeny[!duplicated(rownames(vst_progeny)), ] #remove dups

any(rownames(vst_progeny) == "")       # should be FALSE
any(duplicated(rownames(vst_progeny))) # should be FALSE

vst_progenydf <- as.data.frame(vst_progeny)
vst_progenydf$ID <- row.names(vst_progenydf)
num_samples <- ncol(vst_progenydf) - 1 

# Select the last column (the ID column) first, followed by the original sample columns (1 to num_samples)
vst_progenydf <- vst_progenydf[, c(ncol(vst_progenydf), 1:num_samples)]

write.csv(vst_progenydf, "VST_progeny50.csv", row.names = TRUE)
write.csv(targets, "targets.csv", row.names = FALSE)

plotPCA(vst_dds, intgroup = "condition")
#83% variance explained in PC1, 9% in PC2, distinction between ASD and CTRL groups

# PCA plot
pca_plot <- plotPCA(vst_dds, intgroup = "condition") +
  theme_bw() +
  theme(
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10)
  ) +
  labs(title = "PCA of VST-transformed counts")



# Save as high-resolution image for publication
ggsave(
  filename = "PCA_plot_VST.pdf",      # use PDF or PNG
  plot = pca_plot,           
  width = 6, height = 5,          # size in inches
  dpi = 300                       # resolution for PNG (ignored for PDF)
)

##### prep for inferring TF activities #####

##### collectri #####
library(progeny)
library(dorothea)
library(tibble)
library(tidyr)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(readr)
library(viper)
library(ggrepel)
source("support_functions.R")

collectri_grn <- collectri() #64516 interactions
collectri_grn

##### OmniPath tutorials#####
import_omnipath_interactions # interactions
import_pathwayextra_interactions #directed and signed interactions
#gene regulatory interactions
import_transcriptional_interactions
import_dorothea_interactions #with confidence levels
import_tf_target_interactions #literature curated

gri <- import_transcriptional_interactions()
gri
gr_graph <- interaction_graph(gri)
gr_graph
paths <- find_all_paths(
  graph = gr_graph,
  start = c('EGFR', 'STAT3'),
  end = c('AKT1', 'ULK1'),
  attr = 'name'
)

net <- decoupleR::get_collectri(organism='human', split_complexes=FALSE)
net

head(res.annot)

ttop_matrix_df <- res.annot[, c("hgnc_symbol", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
head(ttop_matrix_df)

# Remove rows with blank hgnc_symbol
ttop_matrix_df <- ttop_matrix_df[ttop_matrix_df$hgnc_symbol != "", ]
# Remove duplicated hgnc_symbol rows (keeping the first one)
ttop_matrix_df <- ttop_matrix_df[!duplicated(ttop_matrix_df$hgnc_symbol), ]
row.names(ttop_matrix_df) <- ttop_matrix_df$hgnc_symbol

names(ttop_matrix_df)[1] <- "ID"
head(ttop_matrix_df)
ttop_matrix_df <- unique(ttop_matrix_df)

eset <- ttop_matrix_df$stat
names(eset) <- ttop_matrix_df$ID


viperRes_decoupleR <- decoupleR::run_viper(
  mat = eset, # Your statistic matrix (equivalent to eset)
  net = net,                 # <-- Uses the unconverted tidy data frame
  .source = 'source', 
  .target = 'target', 
  .mor = 'mor', 
  minsize = 25
)
viperRes <- viperRes_decoupleR[, c('source', 'score')]
names(viperRes) <- c("ID", "NES")

viperRes <- as.data.frame(viperRes)
row.names(viperRes) <- viperRes$ID
head(viperRes)
write.csv(viperRes, "TF_scores.csv", row.names = TRUE) #217

n_tfs <- 25

# Transform to wide matrix
viperRes_mat <- viperRes_decoupleR %>%
  tidyr::pivot_wider(id_cols = 'condition', 
                     names_from = 'source',
                     values_from = 'score') %>%
  tibble::column_to_rownames('condition') %>%
  as.matrix()

# Get top tfs with more variable means across clusters
tfs <- sample_acts %>%
  dplyr::group_by(source) %>%
  dplyr::summarise(std = stats::sd(score)) %>%
  dplyr::arrange(-abs(std)) %>%
  head(n_tfs) %>%
  dplyr::pull(source)

sample_acts_mat <- sample_acts_mat[,tfs]

# Scale per sample
sample_acts_mat <- scale(sample_acts_mat)

# Choose color palette
colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))
colors.use <- grDevices::colorRampPalette(colors = colors)(100)

my_breaks <- c(seq(-2, 0, length.out = ceiling(100 / 2) + 1),
               seq(0.05, 2, length.out = floor(100 / 2)))

# Plot
pheatmap::pheatmap(mat = sample_acts_mat,
                   color = colors.use,
                   border_color = "white",
                   breaks = my_breaks,
                   cellwidth = 15,
                   cellheight = 15,
                   treeheight_row = 20,
                   treeheight_col = 20)

##### PROGENy #####
library(progeny)
library(dorothea)
library(tibble)
library(tidyr)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(readr)
## We read the normalised counts and the experimental design 
Normalised_counts <- read_csv("~/hnrnpu-causal-multiomics/VST_progeny50.csv") #12878 obs
# OR vst_progenydf
normcounts <- vst_progenydf
Experimental_design <- read_csv("~/hnrnpu-causal-multiomics/targets.csv")
targets

## We read the results from the differential analysis. 
DEGASDvsCTRL <- read_csv("~/hnrnpu-causal-multiomics/DESeq2_results50.csv") #13821 genes

head(normcounts)
colnames(normcounts)[1] <- "gene"
normcounts_mat <- as.matrix(normcounts[,-1])  # remove gene column
rownames(normcounts_mat) <- normcounts$gene     # set rownames

head(normcounts_mat)

DEGASDvsCTRL_mat <- as.matrix(DEGASDvsCTRL$stat)
rownames(DEGASDvsCTRL_mat) <- DEGASDvsCTRL$hgnc_symbol
colnames(DEGASDvsCTRL_mat) <- "stat"
head(DEGASDvsCTRL_mat)

# Initial PROGENy calculation
PathwayActivity_counts <- progeny(normcounts_mat, scale = TRUE,
                                  organism = "Human", top = 100)
Activity_counts <- as.vector(PathwayActivity_counts)

paletteLength <- 100
myColor <- colorRampPalette(c("darkblue", "whitesmoke","indianred"))(paletteLength)

progenyBreaks <- c(seq(min(Activity_counts), 0, 
                       length.out=ceiling(paletteLength/2) + 1),
                   seq(max(Activity_counts)/paletteLength, 
                       max(Activity_counts), 
                       length.out=floor(paletteLength/2)))


# A. Extract sample IDs (row names) and determine groups (CTRL/AD)
sample_ids <- rownames(PathwayActivity_counts)

# Assuming 'CD28' samples are CTRL and 'AD28' samples are AD based on your data.
sample_groups <- ifelse(grepl("^CD28", sample_ids), "CTRL", "AD")
names(sample_groups) <- sample_ids

# B. Create annotation data frame and define the ordering
# Set 'CTRL' as the first level to ensure they appear first on the heatmap.
sample_info <- data.frame(Group = factor(sample_groups, levels = c("CTRL", "AD")))
rownames(sample_info) <- sample_ids

# C. Determine the manual order (CTRL first, then AD)
new_order <- rownames(sample_info)[order(sample_info$Group)]

# D. Reorder the matrix rows
PathwayActivity_counts_ordered <- PathwayActivity_counts[new_order, ]

# E. Define colors for the annotation bar
group_colors <- list(
  Group = c(CTRL = "turquoise", AD = "purple") # Adjust colors as desired
)

progeny_hmap <- pheatmap(t(PathwayActivity_counts_ordered), # **Use the ordered matrix**
                         fontsize=10, 
                         fontsize_row = 10, fontsize_col = 10, 
                         color=myColor, breaks = progenyBreaks, 
                         main = "PROGENy (100)", angle_col = 45,
                         treeheight_col = 0, border_color = NA,
                         cluster_cols = FALSE, # **CRITICAL: Prevents samples from re-clustering**
                         annotation_col = sample_info, # Adds the group color bar
                         annotation_colors = group_colors)
ggsave(
  filename = "PROGENY_heatmap.png",
  plot = progeny_hmap,
  width = 8,
  height = 8,
  dpi = 300
)

head(DEGASDvsCTRL_mat)
#enrichment analysis using competitive permutation approach to assess
#significance of pathway activity to end with NES for each pathway
PathwayActivity_zscore <- progeny(DEGASDvsCTRL_mat, 
                                   scale=TRUE, organism="Human", top = 100, perm = 10000, z_scores = TRUE) %>%
  t()
colnames(PathwayActivity_zscore) <- "NES"

PathwayActivity_zscore_df <- as.data.frame(PathwayActivity_zscore) %>% 
  rownames_to_column(var = "Pathway") %>%
  dplyr::arrange(NES) %>%
  dplyr::mutate(Pathway = factor(Pathway))

NES_pathway <- ggplot(PathwayActivity_zscore_df,aes(x = reorder(Pathway, NES), y = NES)) + 
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
  filename = "PROGENY_NES_pathways.png",
  plot = NES_pathway,
  width = 8,
  height = 8,
  dpi = 300
)

#JAK.STAT pathway most active - visualise most responsive genes (progeny weights)
#along with their statistic values to interpret results
#scatterplot shows which genes contributing most to the pathway enrichment

prog_matrix <- getModel("Human", top=100) %>% 
  as.data.frame()  %>%
  tibble::rownames_to_column("GeneID")

DEGASDvsCTRL_df <- DEGASDvsCTRL_mat %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("GeneID")

scat_plots <- progeny::progenyScatter(df = DEGASDvsCTRL_df, 
                                       weight_matrix = prog_matrix, 
                                       statName = "stat", verbose = FALSE)
plot(scat_plots[[1]]$`JAK-STAT`) 

#progeny results as input for CARNIVAL
#CARNIVAL sets weights based on PROGENy scores in each pathway-related node
#to find more relevant solutions
#set progeny zscores false to return pathway activity values between 1 and -1

PathwayActivity_CARNIVALinput <- progeny(DEGASDvsCTRL_mat, 
                                          scale=TRUE, organism="Human", top = 100, perm = 10000, z_scores = FALSE) %>%
  t () %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "Pathway") 
colnames(PathwayActivity_CARNIVALinput)[2] <- "score"
write_csv(PathwayActivity_CARNIVALinput, 
          "~/hnrnpu-causal-multiomics/processeddata/PathwayActivity_CARNIVALinput.csv")

##### CARNIVAL #####
library(progeny)
library(dorothea)
library(CARNIVAL)
library(OmnipathR)
library(readr)
library(tibble)
library(tidyr)
library(dplyr)
library(visNetwork)
library(ggplot2)
library(pheatmap)
## For the volcano plot (related to support functions)
library(ggrepel)

source("assignPROGENyScores.r")
source("generateTFList.r")
source("carnival_visNetwork.r")

shiny::runGitHub(repo = "ShinyFUNKI", username = "saezlab", subdir = "FUNKI")
