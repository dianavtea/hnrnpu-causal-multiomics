#Another trial for data setup 01.12.25
#DIFFERENCE = filtered on merged gene counts making file
#from 63677 obs to 17636 genes, typical amount for cleaned up dataframe
#still transcriptutorial

library(readr)
library(ggplot2)
library(reshape)
library(pheatmap)
library(gridExtra)
library(grid)
library(cowplot)
library(ggrepel)
BiocManager::install("hexbin")
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
#63677 genes

#simple summary of experimental design
targets <- as.data.frame(matrix(NA, length(names(merged_gene_counts)),2))
names(targets) <- c("sample", "condition")
targets$sample <- names(merged_gene_counts)
targets$condition <- c("CTRL", "CTRL", "CTRL", "ASD", "ASD", "ASD")
head(targets)

#merged gene counts without NA
# Keep genes expressed (>20 counts) in at least 2/3 replicates per group (CTRL vs ASD)
gene_counts_CTRL <- rowSums(merged_gene_counts[,1:3] > 20) >= 2
gene_counts_ASD <- rowSums(merged_gene_counts[,4:6] > 20) >= 2
merged_gene_counts_filtered <- merged_gene_counts[gene_counts_CTRL | gene_counts_ASD, ] #17636 genes
any(is.na(merged_gene_counts_filtered)) #FALSE, no NA

ddsf <- DESeqDataSetFromMatrix(countData = merged_gene_counts_filtered,
                              colData = targets,
                              design= ~ condition)

ddsf$condition<- relevel( ddsf$condition, "CTRL" )
ddsf <- DESeq(ddsf)
resultsNames(ddsf)
resf <- results(ddsf)

resf.mat <- cbind(counts(ddsf, normalized=TRUE), resf$log2FoldChange, resf$padj)
mcols(resf, use.names = TRUE)
saveRDS(ddsf, file="dds_DEG_ASDvsCTRL_D28_filt.rds")
head(resf)

##### annotations #####
#Connect to GRCh37 Ensembl
# Use rownames of res directly
ensembl_ids <- rownames(resf)

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

resf.df <- as.data.frame(resf)
resf.df$ensembl <- ensembl_ids

resf.annot <- merge(
  resf.df,
  annot,
  by.x = "ensembl",
  by.y = "ensembl_gene_id",
  all.x = TRUE,
  sort = FALSE
)

resf.annot <- resf.annot[match(ensembl_ids, resf.annot$ensembl), ]

colnames(resf.annot)[colnames(resf.annot) == "entrezgene_id"] <- "entrez"
colnames(resf.annot)[colnames(resf.annot) == "hgnc_symbol"]   <- "hgnc_symbol"
head(resf.annot)
write.csv(resf.annot, "DESeq2_results.csv", row.names = FALSE)
##### VST #####
vst_ddsf <- vst(ddsf, blind = FALSE)
vst_matf <- assay(vst_ddsf)
#105816 elements
dim(vst_matf) #17636 genes x 6 samples

# Ensure rows of res.annot match vst_mat
resf.annot <- resf.annot[match(rownames(vst_matf), resf.annot$ensembl), ]

vst_annotf <- cbind(resf.annot[, c("ensembl", "entrez", "hgnc_symbol")], vst_matf)

write.csv(vst_annotf, "VST_counts_with_annotation_filt.csv", row.names = FALSE)

#vst for progeny
vst_progenyf <- vst_matf
rownames(vst_progenyf) <- resf.annot$hgnc_symbol
head(vst_progenyf)
any(is.na(vst_progenyf)) #false

rownames(vst_progenyf)[duplicated(rownames(vst_progenyf))] 
vst_progenyf <- vst_progenyf[rownames(vst_progenyf) != "", ] #remove blanks
vst_progenyf <- vst_progenyf[!duplicated(rownames(vst_progenyf)), ] #remove dups

any(rownames(vst_progenyf) == "")       # should be FALSE
any(duplicated(rownames(vst_progenyf))) # should be FALSE

write.csv(vst_progenyf, "VST_progeny_filt.csv", row.names = TRUE)
write.csv(targets, "targets.csv", row.names = FALSE)

plotPCA(vst_ddsf, intgroup = "condition")


# PCA plot
pca_plot <- plotPCA(vst_ddsf, intgroup = "condition") +
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
  filename = "PCA_plot.pdf",      # use PDF or PNG
  plot = pca_plot,           
  width = 6, height = 5,          # size in inches
  dpi = 300                       # resolution for PNG (ignored for PDF)
)

##### prep for infer TF activities #####
#open file which is currently excel format
DEG_ASD_vs_CTRL <- read.csv("~/hnrnpu-causal-multiomics/DESeq2_results.csv")
head(DEG_ASD_vs_CTRL)
#clean file w symbol, stat, logFC, pvalue and padj 
DEG_allstatsf <- DEG_ASD_vs_CTRL %>%
  # Keep only the columns you want
  dplyr::select(hgnc_symbol, stat, log2FoldChange, pvalue, padj) %>%
  # Make into numeric and not chr
  mutate(
    stat = as.numeric(stat), 
    log2FoldChange = as.numeric(log2FoldChange),
    pvalue = as.numeric(pvalue),
    padj = as.numeric(padj)
  ) %>%
  rename(
    ID = hgnc_symbol,
    statistic = stat
  ) %>%
  # Remove rows with missing values
  filter(!is.na(ID), !is.na(statistic), !is.na(log2FoldChange), !is.na(pvalue),!is.na(padj))
head(DEG_allstatsf)
any(is.na(DEG_allstatsf)) #check if any NA in DEG_allstats
write.csv(DEG_allstatsf, "DEG_allstatsf.csv", row.names = FALSE)
#17627 genes now


##### important outputs #####
ddsf
resf
vst_annotf
vst_progenyf
targets
DEG_allstatsf