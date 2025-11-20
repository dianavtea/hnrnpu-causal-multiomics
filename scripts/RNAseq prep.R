##### 20.11.2025 RNAseq DESeq2 and vst #####
##### packages #####
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("ggpubr")
BiocManager::install("org.Hs.eg.db")
library(DESeq2)
library(ggplot2)
library(ggpubr)
library(biomaRt)
library(RSQLite)
library(org.Hs.eg.db)

#clean environment
rm(list=ls())

##### DESeq2 #####
# read data
cr <- read.delim("~/Masters/RP2 Tammimies/hnrnpu-causal-multiomics/rawdata/merged_gene_counts 1.txt", header=TRUE, stringsAsFactors = FALSE,row.names=1)
md <- read.delim("~/Masters/RP2 Tammimies/hnrnpu-causal-multiomics/rawdata/columndata 1.txt", header=TRUE, sep="\t", row.names = 1,stringsAsFactors = FALSE)
#check that col names of cr are the same as raw names in md
all.equal(colnames(cr),rownames(md))

#check if data is correctly imported
head(colnames(md))
all(rownames(md) %in% colnames(cr))
all(rownames(md) == colnames(cr))

#create group for multiple factor analysis
md$group <- factor(paste0(md$StatusCells, md$Day))
md

#construct a deseq dataset
dds <- DESeqDataSetFromMatrix(countData = cr,
                              colData = md,
                              design= ~ group)

dds$group<- relevel( dds$group, "CTRLD28" )

dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds)

res.mat <- cbind(counts(dds, normalized=TRUE), res$log2FoldChange, res$padj)

mcols(res, use.names = TRUE)
saveRDS(dds, file="dds_DEG_ASDvsCTRL_D28.rds")
head(res)

##### annotations #####
#Connect to GRCh37 Ensembl
library(biomaRt)
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

##### VST #####
vst_dds <- vst(dds, blind = FALSE)
vst_mat <- assay(vst_dds)
# Ensure rows of res.annot match vst_mat
res.annot <- res.annot[match(rownames(vst_mat), res.annot$ensembl), ]

vst_annot <- cbind(res.annot[, c("ensembl", "entrez", "hgnc_symbol")], vst_mat)

write.csv(vst_annot, "VST_counts_with_annotation.csv", row.names = FALSE)

##### prep for infer TF activities (ID and stat) #####

