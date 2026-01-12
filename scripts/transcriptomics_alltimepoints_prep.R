##### libraries #####
library(dplyr)
library(DESeq2)
library(biomaRt)

# repeat code for day 0, 5 and 28

##### RNAseq data prep #####
merged_gene_counts_all <- as.data.frame(
  read.delim("~/hnrnpu-causal-multiomics/rawdata/bulkRNAseq/Complete_merged_gene_counts_more20count.txt", header=TRUE, sep = ',', stringsAsFactors = FALSE,row.names=1))
#21262 transcripts

meta <- as.data.frame(
  read.delim("~/hnrnpu-causal-multiomics/rawdata/bulkRNAseq/Complete_columndata.txt", header=TRUE, stringsAsFactors = FALSE, row.names=1))

#extract relevant metadata
meta_filtered <- meta %>%
  filter(
    !grepl("^si", Sample_Name),               # removes siNTC & siHNRNPU
    grepl("^CTRL|^HNRNPU", Sample_Name),       # keeps CTRL & HNRNPUdel
    Day %in% c("D0", "D5", "D28")
  )

#counts = metadata
counts_filtered <- merged_gene_counts_all[, colnames(merged_gene_counts_all) %in% meta_filtered$Sample_ID]
all(colnames(counts_filtered) == meta_filtered$Sample_ID) #TRUE

meta_filtered$Day <- factor(
  meta_filtered$Day,
  levels = c("D0", "D5", "D28")
)

meta_filtered$condition <- factor(
  meta_filtered$StatusHNRNPU,
  levels = c("HNRNPU_WT", "HNRNPU_KD"),
  labels = c("CTRL", "HNRNPUdel")
)

head(counts_filtered)
head(meta_filtered)
any(is.na(counts_filtered)) #FALSE

#DESeq2
dds <- DESeqDataSetFromMatrix(countData = counts_filtered,
                                colData = meta_filtered,
                                design= ~ Day + condition + Day:condition)

dds <- DESeq(dds)
resultsNames(dds)
saveRDS(dds, file="DEG_ASDvsCTRL.rds")
dds <- readRDS("DEG_ASDvsCTRL.rds")

#ASD vs CTRL at D0
resD0 <- results(dds, name = "condition_HNRNPUdel_vs_CTRL")
head(resD0)

# ASD vs CTRL specifically at Day 5
res_D5_only <- results(dds, 
                           contrast = list( c("condition_HNRNPUdel_vs_CTRL", "DayD5.conditionHNRNPUdel") ))
summary(res_D5_only)
head(res_D5_only)

# ASD vs CTRL specifically at Day 28
res_D28_only <- results(dds, 
                       contrast = list( c("condition_HNRNPUdel_vs_CTRL", "DayD28.conditionHNRNPUdel") ))
summary(res_D28_only)
head(res_D28_only)

##### D28 only #####

#Connect to GRCh37 Ensembl - annotations
# formatting!

# Use rownames of res directly
ensembl_ids <- rownames(res_D28_only)

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

res_D28_only.df <- as.data.frame(res_D28_only)
res_D28_only.df$ensembl <- ensembl_ids

res_D28_only.annot <- merge(
  res_D28_only.df,
  annot,
  by.x = "ensembl",
  by.y = "ensembl_gene_id",
  all.x = TRUE,
  sort = FALSE
)

res_D28_only.annot <- res_D28_only.annot[match(ensembl_ids, res_D28_only.annot$ensembl), ]

colnames(res_D28_only.annot)[colnames(res_D28_only.annot) == "entrezgene_id"] <- "entrez"
colnames(res_D28_only.annot)[colnames(res_D28_only.annot) == "hgnc_symbol"]   <- "hgnc_symbol"
head(res_D28_only.annot)
#write.csv(res_D28_only.annot, "DESeq2_res_D28_only.csv", row.names = FALSE)


res_D28only_clean <- res_D28_only.annot %>%
  dplyr::filter(!is.na(hgnc_symbol), hgnc_symbol != "") %>%
  dplyr::arrange(desc(abs(stat))) %>%   
  dplyr::distinct(hgnc_symbol, .keep_all = TRUE)
#17017 

write.csv(res_D28only_clean, "DESeq2_res_D28only_unique.csv", row.names = FALSE) #for CARNIVAL
res_D28only_clean <- read.csv("DESeq2_res_D28only_unique.csv")

##### D5 only #####
ensembl_ids <- rownames(res_D5_only)

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

res_D5_only.df <- as.data.frame(res_D5_only)
res_D5_only.df$ensembl <- ensembl_ids

res_D5_only.annot <- merge(
  res_D5_only.df,
  annot,
  by.x = "ensembl",
  by.y = "ensembl_gene_id",
  all.x = TRUE,
  sort = FALSE
)

res_D5_only.annot <- res_D5_only.annot[match(ensembl_ids, res_D5_only.annot$ensembl), ]

colnames(res_D5_only.annot)[colnames(res_D5_only.annot) == "entrezgene_id"] <- "entrez"
colnames(res_D5_only.annot)[colnames(res_D5_only.annot) == "hgnc_symbol"]   <- "hgnc_symbol"
head(res_D5_only.annot)
#write.csv(res_D5_only.annot, "DESeq2_res_D5_only.csv", row.names = FALSE)

#addition 01.01.26
res_D5only_clean <- res_D5_only.annot %>%
  dplyr::filter(!is.na(hgnc_symbol), hgnc_symbol != "") %>%
  dplyr::arrange(desc(abs(stat))) %>%   
  dplyr::distinct(hgnc_symbol, .keep_all = TRUE)

write.csv(res_D5only_clean, "DESeq2_res_D5only_unique.csv", row.names = FALSE) #for CARNIVAL
res_D5only_clean <- read.csv("DESeq2_res_D5only_unique.csv")

##### D0 #####
ensembl_ids <- rownames(resD0)

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

resD0.df <- as.data.frame(resD0)
resD0.df$ensembl <- ensembl_ids

resD0.annot <- merge(
  resD0.df,
  annot,
  by.x = "ensembl",
  by.y = "ensembl_gene_id",
  all.x = TRUE,
  sort = FALSE
)

resD0.annot <- resD0.annot[match(ensembl_ids, resD0.annot$ensembl), ]

colnames(resD0.annot)[colnames(resD0.annot) == "entrezgene_id"] <- "entrez"
colnames(resD0.annot)[colnames(resD0.annot) == "hgnc_symbol"]   <- "hgnc_symbol"
head(resD0.annot)
#write.csv(res_D0_only.annot, "DESeq2_res_D0_only.csv", row.names = FALSE)

resD0_clean <- resD0.annot %>%
  dplyr::filter(!is.na(hgnc_symbol), hgnc_symbol != "") %>%
  dplyr::arrange(desc(abs(stat))) %>%   
  dplyr::distinct(hgnc_symbol, .keep_all = TRUE)

write.csv(resD0_clean, "DESeq2_res_D0only_unique.csv", row.names = FALSE) #for CARNIVAL
resD0_clean <- read.csv("DESeq2_res_D0only_unique.csv")
