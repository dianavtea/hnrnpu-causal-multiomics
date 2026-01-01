shiny::runGitHub(repo = "ShinyFUNKI", username = "saezlab", subdir = "FUNKI")

##### libraries #####
library(OmnipathR)
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
library(progeny)
library(dorothea)
library(tibble)
library(tidyr)
library(viper)
##### RNA prep #####
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

##### CollecTRI #####
net <- decoupleR::get_collectri(organism='human', split_complexes=FALSE) #43159
net

DESeq_res <- res.annot[, c("hgnc_symbol", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
head(DESeq_res)

# Remove rows with blank hgnc_symbol
DESeq_res <- DESeq_res[DESeq_res$hgnc_symbol != "", ]
# Remove duplicated hgnc_symbol rows (keeping the first one)
DESeq_res <- DESeq_res[!duplicated(DESeq_res$hgnc_symbol), ]
row.names(DESeq_res) <- DESeq_res$hgnc_symbol

names(DESeq_res)[1] <- "ID"
head(DESeq_res)
DESeq_res <- unique(DESeq_res)
write.csv(DESeq_res, "DESeq_results_formatted.csv")
#12878 transcripts

eset <- DESeq_res$stat
names(eset) <- DESeq_res$ID

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
head(viperRes) #217 TFs
write.csv(viperRes, "TF_scores.csv", row.names = TRUE) #217

n_tfs <- 25

## top 25 TF in general (abs)
# Filter for top 25 most active TFs (by absolute score)
top_tfs <- viperRes_decoupleR %>%
  arrange(desc(abs(score))) %>%
  slice_head(n = 25)

# Plot them directly 
ggplot(top_tfs, aes(x = reorder(source, score), y = score, fill = score > 0)) +
  geom_col() +
  coord_flip() + # Makes it horizontal
  theme_minimal() +
  labs(title = "Top 25 TF Activities (ASD vs CTRL)", 
       y = "Normalized Enrichment Score (NES)", 
       x = "Transcription Factor") +
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "blue"), 
                    labels = c("Activated", "Inhibited"),
                    name = "Status")

##top 15 for activated and top 15 for inhibited TFs
# Select Top 15 Activated (Highest Positive Scores)
top_up <- viperRes_decoupleR %>%
  arrange(desc(score)) %>%
  slice_head(n = 15)

# Select Top 15 Inhibited (Lowest Negative Scores)
top_down <- viperRes_decoupleR %>%
  arrange(score) %>%
  slice_head(n = 15)

# Combine them into one dataframe
top_30_balanced <- bind_rows(top_up, top_down)

# We use reorder() to ensure they appear sorted by score, not alphabetically
top30 <- ggplot(top_30_balanced, aes(x = reorder(source, score), y = score, fill = score > 0)) +
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
  filename = "top30TFs.png",      # use PDF or PNG
  plot = top30,           
  width = 8, height = 8,          # size in inches
  dpi = 300                       # resolution for PNG (ignored for PDF)
)

measObj <- viperRes$NES
names(measObj) <- viperRes$ID

is.numeric(measObj)
head(measObj)
head(names(measObj))

weightObj <- PROGENY$score
names(weightObj) <- PROGENY$Pathway

##### CARNIVAL #####
PROGENY <- read.csv("processeddata/PathwayActivity_CARNIVALinputf.csv")
sif_clean <- read.csv("clean_omnipath_PKN.csv")

?runCARNIVAL
carnival_res <- runCARNIVAL(
  inputObj  = NULL,             
  measObj   = measObj,             # TF activity measurements
  netObj    = sif_clean,                    # prior knowledge network
  weightObj = weightObj,         # pathway activity scores (or progeny_matrix)
  solverPath = "C:/Program Files/IBM/ILOG/CPLEX_Studio2211/cplex/bin/x64_win64/cplex.exe", 
  solver     = "cplex",
  timelimit  = 7200, 
  mipGAP     = 0.05,
  poolrelGAP = 0
)
saveRDS(carnival_res,"~/hnrnpu-causal-multiomics/processeddata/carnival_res.rds")

head(viperResdf)

#inputObj: data frame of list fir target of perturbation, OPTIONAL or default set to NULL to run invCARNIVAL when inputs are not known

##### metabolic prep #####
library(readr)
library(limma)


limmaRes <- runLimma(raw_metabolomic, targets_1, comparisons = list(c(1,-2)))
ttop_tumour_vs_healthy <- ttopFormatter(
  topTable(limmaRes[[1]], coef = 1, number = nrow(raw_metabolomic), adjust.method = "fdr")
)
write_csv(ttop_tumour_vs_healthy, "results/metab_ttop_tumour_vs_healthy.csv")
