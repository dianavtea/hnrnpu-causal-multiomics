library(dplyr)
library(DESeq2)
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
#ASD vs CTRL at D0
resD0 <- results(dds, name = "condition_HNRNPUdel_vs_CTRL")
head(resD0)

#ASD vs CTRL at D5 relative to D0
resD5 <- results(dds, name = "DayD5.conditionHNRNPUdel")
summary(resD5)
head(resD5)

# ASD vs CTRL specifically at Day 5
res_D5_only <- results(dds, 
                           contrast = list( c("condition_HNRNPUdel_vs_CTRL", "DayD5.conditionHNRNPUdel") ))
summary(res_D5_only)
head(res_D5_only)

#ASD vs CTRL at D28 relative to D0
resD28 <- results(dds, name = "DayD28.conditionHNRNPUdel")
summary(resD28)
head(resD28)

# ASD vs CTRL specifically at Day 28
res_D28_only <- results(dds, 
                       contrast = list( c("condition_HNRNPUdel_vs_CTRL", "DayD28.conditionHNRNPUdel") ))
summary(res_D28_only)
head(res_D28_only)

##### D28 only #####
library(biomaRt)
#Connect to GRCh37 Ensembl - annotations
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

#ADDITION 28.12.25
# Remove rows where symbol is NA or blank
#res_D28only_clean <- res_D28_only.annot[ !is.na(res_D28_only.annot$hgnc_symbol) & res_D28_only.annot$hgnc_symbol != "", ]
#17017
# Sort by baseMean (highest expression first)
#res_D28only_clean <- res_D28only_clean[order(res_D28only_clean$baseMean, decreasing = TRUE), ]
# Remove duplicates
# !duplicated() keeps the FIRST instance (which is now the highest expressed one)
#res_D28only_clean <- res_D28only_clean[!duplicated(res_D28only_clean$hgnc_symbol), ]

write.csv(res_D28only_clean, "DESeq2_res_D28only_unique.csv", row.names = FALSE) #for CARNIVAL

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

vstD28only_annot <- cbind(res_D28only_clean[, c("ensembl", "entrez", "hgnc_symbol")], vst_D28only_clean)
head(vstD28only_annot)
write.csv(vstD28only_annot, "VSTcounts_annot_D28only.csv", row.names = FALSE)

##### D28 vs D0  #####
##### D5 only #####
#Connect to GRCh37 Ensembl - annotations
# Use rownames of res directly
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

#ADDITION 28.12.25
# Remove rows where symbol is NA or blank
#res_D5only_clean <- res_D5_only.annot[ !is.na(res_D5_only.annot$hgnc_symbol) & res_D5_only.annot$hgnc_symbol != "", ]
#17017
# Sort by baseMean (highest expression first)
#res_D5only_clean <- res_D5only_clean[order(res_D5only_clean$baseMean, decreasing = TRUE), ]
# Remove duplicates
# !duplicated() keeps the FIRST instance (which is now the highest expressed one)
#res_D5only_clean <- res_D5only_clean[!duplicated(res_D5only_clean$hgnc_symbol), ]

write.csv(res_D5only_clean, "DESeq2_res_D5only_unique.csv", row.names = FALSE) #for CARNIVAL

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

vstD5only_annot <- cbind(res_D5only_clean[, c("ensembl", "entrez", "hgnc_symbol")], vst_D5only_clean)
head(vstD5only_annot)
write.csv(vstD5only_annot, "VSTcounts_annot_D5only.csv", row.names = FALSE)

##### D5 vs D0 #####
##### D0 #####
#Connect to GRCh37 Ensembl - annotations
# Use rownames of res directly
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

#ADDITION 28.12.25
# Remove rows where symbol is NA or blank
#resD0_clean <- resD0.annot[ !is.na(resD0.annot$hgnc_symbol) & resD0.annot$hgnc_symbol != "", ]
#17017
# Sort by baseMean (highest expression first)
#resD0_clean <- resD0_clean[order(resD0_clean$baseMean, decreasing = TRUE), ]
# Remove duplicates
# !duplicated() keeps the FIRST instance (which is now the highest expressed one)
#resD0_clean <- resD0_clean[!duplicated(resD0_clean$hgnc_symbol), ]

write.csv(resD0_clean, "DESeq2_res_D0only_unique.csv", row.names = FALSE) #for CARNIVAL

###### VST ######
vst_all <- vst(dds, blind = FALSE)
vstD0 <- vst_all[, colData(dds)$Day == "D0"]
vstmatD0 <- assay(vstD0) #127572 elements
dim(vstmatD0) #21262 transcripts x 6 samples

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

vstD0_annot <- cbind(resD0_clean[, c("ensembl", "entrez", "hgnc_symbol")], vst_D0_clean)
head(vstD0_annot)
write.csv(vstD0_annot, "VSTcounts_annot_D0.csv", row.names = FALSE)

