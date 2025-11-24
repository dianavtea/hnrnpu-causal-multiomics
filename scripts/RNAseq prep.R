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
library(readxl)

#clean environment
rm(list=ls())

##### DESeq2 #####
# read data
cr <- read.delim("~/Masters/RP2/hnrnpu-causal-multiomics/rawdata/merged_gene_counts 1.txt", header=TRUE, stringsAsFactors = FALSE,row.names=1)
md <- read.delim("~/Masters/RP2/hnrnpu-causal-multiomics/rawdata/columndata 1.txt", header=TRUE, sep="\t", row.names = 1,stringsAsFactors = FALSE)
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

#OPTIONAL - not very valuable for only 6 samples 
#dists <- dist(t(vst_mat))  # transpose so samples are columns
#hc <- hclust(dists)
#plot(hc, main = "Hierarchical clustering of samples (VST)")
#pheatmap(as.matrix(dists),
#         clustering_distance_rows = "euclidean",
#         clustering_distance_cols = "euclidean",
#         main = "Sample-to-sample distances (VST)")

# Ensure rows of res.annot match vst_mat
res.annot <- res.annot[match(rownames(vst_mat), res.annot$ensembl), ]

vst_annot <- cbind(res.annot[, c("ensembl", "entrez", "hgnc_symbol")], vst_mat)

write.csv(vst_annot, "VST_counts_with_annotation.csv", row.names = FALSE)

head(vst_mat)
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

write.csv(vst_progeny, "VST_progeny.csv", row.names = TRUE)

##### prep for infer TF activities (ID and stat) #####
#open file which is currently excel format
DEG_ASD12D28_vs_CTRL9D28 <- read_excel("~/Masters/RP2/hnrnpu-causal-multiomics/rawdata/bulkRNAseq/DEG_ASD12D28_vs_CTRL9D28 .xlsx")
View(DEG_ASD12D28_vs_CTRL9D28)

#clean file w symbol, stat, logFC and padj (filtered later)
DEG_allstats <- DEG_ASD12D28_vs_CTRL9D28 %>%
  # Keep only the columns you want
  dplyr::select(hgnc_symbol, stat, log2FoldChange, padj) %>%
  # Make into numeric and not chr
  mutate(
    stat = as.numeric(stat), 
    log2FoldChange = as.numeric(log2FoldChange),
    padj = as.numeric(padj)
  ) %>%
  # Rename columns to how FUNKI wants it
  rename(
    ID = hgnc_symbol,
    statistic = stat
  ) %>%
  # Remove rows with missing values
  filter(!is.na(ID), !is.na(statistic), !is.na(log2FoldChange), !is.na(padj))
head(DEG_allstats)
any(is.na(DEG_allstats)) #check if any NA in DEG_allstats
write.csv(DEG_allstats, "DEG_allstats.csv", row.names = FALSE)
#can put this csv put into shinyFUNKI (optional)

#OPTIONAL: run shinyFUNKI locally (setup available in other R file) - optional
#shiny::runGitHub(repo = "ShinyFUNKI", username = "saezlab", subdir = "FUNKI")

#remove any random or dup ids with numbers only
deg_clean <- DEG_allstats %>%
  dplyr::filter(!grepl("^[0-9]+$", ID))

deg_mat <- deg_clean %>%
  tibble::column_to_rownames("ID") %>%
  as.matrix()

#need to be in matrix format
head(deg_mat)

#filter for only one statistic - i choose wald statistic
degstat <- DEG_allstats %>%
  dplyr::select(ID, statistic)
head(degstat)

deg_clean1 <- degstat %>%
  dplyr::filter(!grepl("^[0-9]+$", ID)) %>%
  dplyr::distinct(ID, .keep_all = TRUE)

degstat_mat <- deg_clean1 %>%
  tibble::column_to_rownames("ID") %>%
  as.matrix()
head(degstat_mat)
any(is.na(degstat_mat))
write.csv(degstat_mat, "DEG_statistic.csv", row.names = TRUE)

##### most important outputs here #####
dds #DESeq2
res #results of DESeq2
vst_annot
vst_progeny
DEG_allstats
degstat_mat #wald stat only and ID 'ID', 'statistic'