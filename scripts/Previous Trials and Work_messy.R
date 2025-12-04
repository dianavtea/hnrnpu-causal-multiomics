##### 04.12.25 trials and rough work to clean up real files #####

##### 20.11.2025 RNAseq DESeq2 and vst #####
##### packages #####
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("ggpubr")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("biomaRt")
install.packages("dplyr")
library(DESeq2)
library(ggplot2)
library(ggpubr)
library(biomaRt)
library(RSQLite)
library(org.Hs.eg.db)
library(readxl)
library(dplyr)

#clean environment
rm(list=ls())

##### DESeq2 #####
# read data
cr <- read.delim("~/hnrnpu-causal-multiomics/rawdata/merged_gene_counts 1.txt", header=TRUE, stringsAsFactors = FALSE,row.names=1)
md <- read.delim("~/hnrnpu-causal-multiomics/rawdata/columndata 1.txt", header=TRUE, sep="\t", row.names = 1,stringsAsFactors = FALSE)
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
vst_progeny <- vst_mat #382062 elements
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
DEG_ASD12D28_vs_CTRL9D28 <- read_excel("~/hnrnpu-causal-multiomics/rawdata/bulkRNAseq/DEG_ASD12D28_vs_CTRL9D28 .xlsx")
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

##### PROGENy #####
##### 21.11.2025 #####
##### library installation #####
BiocManager::install("decoupleR")
BiocManager::install("progeny")
BiocManager::install("dorothea")
library(Biobase)
library(stats)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(decoupleR)
library(reshape2)
library(progeny)
library(dorothea)
library(tibble)
library(pheatmap)
library(readr)

##### trials and tests of PROGENy #####
#This function is designed for getting a model matrix with top significant genes for each pathway
progenymodel <- getModel(organism = "Human", top = 100, decoupleR = T)
head(progenymodel)
#top = top n of genes according to significance
#decoupleR true = compatible w decoupleR

model_human_full
#22479 genes, associated pathways, weight and the p-value
get("model_human_full", envir = .GlobalEnv)

##### Trial 1: Progeny with 'progeny' document as tutorial #####
#hgnc symbols needed in rows, vst counts in columns

pathwayscores <- progeny(vst_progeny, 
                         scale = TRUE, 
                         organism = "Human", 
                         top = 100, 
                         perm = 10000, 
                         verbose = TRUE, 
                         z_scores = TRUE, 
                         get_nulldist = TRUE, 
                         assay_name = "RNA", 
                         return_assay = FALSE)

pathwayscores[[1]] #individual pathway scores for each

#perm = 10000 to ensure stat significance accurately determined, recommended by saez
head(pathwayscores)

##### visualisation #####
# Compute average pathway score per group
ctrl_avg <- rowMeans(pathwayscores[[1]][, c("CD28_146", "CD28_147", "CD28_148")])
asd_avg  <- rowMeans(pathwayscores[[1]][, c("AD28_144", "AD28_145", "AD28_149")])

# Create a data frame for plotting
df <- data.frame(
  pathway = names(ctrl_avg),
  Ctrl = ctrl_avg,
  ASD  = asd_avg
)
head(df)

#scatterplot
scatterplot <- ggplot(df, aes(x = Ctrl, y = ASD, label = pathway)) +
  geom_point(color = "darkred", size = 3) +
  geom_text(vjust = -0.5, size = 3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(
    title = "PROGENY Pathway Activity: Ctrl vs ASD",
    x = "Ctrl (average z-score)",
    y = "ASD (average z-score)"
  )
ggsave(
  filename = "PROGENY_scatterplot.png",
  plot = scatterplot,
  width = 8,
  height = 8,
  dpi = 300
)

#barplot
df_long <- data.frame(
  pathway = rep(names(ctrl_avg), 2),
  value = c(ctrl_avg, asd_avg),
  group = rep(c("Ctrl", "ASD"), each = length(ctrl_avg))
)

activitybarplot <- ggplot(df_long, aes(x = pathway, y = value, fill = group)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(
    title = "PROGENY Pathway Activity",
    y = "Average z-score",
    x = "Pathway"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(
  filename = "PROGENY_pathway_activity.png",
  plot = activitybarplot,
  width = 8,
  height = 8,
  dpi = 300
)


#per sample analysis
# Extract the pathway scores matrix
scores <- as.matrix(pathwayscores[[1]])  # pathways x samples

# Define sample groups
ctrl_samples <- c("CD28_146", "CD28_147", "CD28_148")
asd_samples  <- c("AD28_144", "AD28_145", "AD28_149")

#Convert to long format for ggplot
df_long <- data.frame(
  pathway = rep(rownames(scores), times = ncol(scores)),
  sample = rep(colnames(scores), each = nrow(scores)),
  value  = as.vector(scores),
  group  = rep(c(rep("Ctrl", length(ctrl_samples)), rep("ASD", length(asd_samples))),
               each = nrow(scores))
)

#Plot boxplots + jitter for all pathways
boxplot <- ggplot(df_long, aes(x = pathway, y = value, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.4, position = position_dodge(width = 0.8)) +
  geom_jitter(aes(color = group), size = 2, position = position_dodge(width = 0.8)) +
  theme_minimal() +
  labs(title = "PROGENY Pathway Scores per Sample",
       x = "Pathway",
       y = "PROGENY z-score") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(
  filename = "PROGENY_boxplotpersample.png",
  plot = boxplot,
  width = 8,
  height = 8,
  dpi = 300
)

##### Trial 2: progeny with transcriptutorial #####
PathwayActivity_counts <- progeny(vst_progeny, scale=TRUE, 
                                  organism="Human", top = 100)
Activity_counts <- as.vector(PathwayActivity_counts)

#heatmap
paletteLength <- 100
myColor <- 
  colorRampPalette(c("darkblue", "whitesmoke","indianred"))(paletteLength)

progenyBreaks <- c(seq(min(Activity_counts), 0, 
                       length.out=ceiling(paletteLength/2) + 1),
                   seq(max(Activity_counts)/paletteLength, 
                       max(Activity_counts), 
                       length.out=floor(paletteLength/2)))

progeny_hmap <- pheatmap(t(PathwayActivity_counts),fontsize=14, 
                         fontsize_row = 10, fontsize_col = 10, 
                         color=myColor, breaks = progenyBreaks, 
                         main = "PROGENy (Top 100)", angle_col = 45,
                         treeheight_col = 0,  border_color = NA)
ggsave(
  filename = "PROGENY_heatmap_t2.png",
  plot = progeny_hmap,
  width = 8,
  height = 8,
  dpi = 300
)

#Now, we run an enrichment analysis using a competitive permutation approach to 
#assess the significance of the pathway activity. 
#We end up with Normalised Enrichment Scores (NES) for each pathway.

PathwayActivity_zscore <- progeny(degstat_mat, 
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

#TGFb most active pathway upon perturbation (HNRNPUdel vs CTRL)
#visualise MAPK most responsive genes (progeny_weights) with t vals to interpret results
#in scatterplot we can see genes contributing most to pathway enrichment
prog_matrix <- getModel("Human", top=100) %>% 
  as.data.frame()  %>%
  tibble::rownames_to_column("GeneID")

degstat_mat_df <- degstat_mat %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("ID")

scat_plots <- progeny::progenyScatter(df = degstat_mat_df, 
                                      weight_matrix = prog_matrix, 
                                      statName = "statistic", verbose = FALSE)
plot(scat_plots[[1]]$`TGFb`) 

##### progeny input for CARNIVAL #####
#CARNIVAL sets weights based on **PROGENy** scores in each pathway-related node in order
#to find more relevant solutions. We therefore run **PROGENy** again with 
#slightly different parameters, setting `z_scores = FALSE` so that **PROGENy** returns pathway activity 
#values between 1 and -1, rather than converting to Z-Scores.

PathwayActivity_CARNIVALinput <- progeny(degstat_mat, 
                                         scale=TRUE, organism="Human", top = 100, perm = 10000, z_scores = FALSE) %>%
  t () %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "Pathway") 
colnames(PathwayActivity_CARNIVALinput)[2] <- "score"

write.csv(PathwayActivity_CARNIVALinput, "PathwayActivity_CARNIVALinput.csv")

##### most important outputs here #####
#trial 1 'progeny' doc
pathwayscores[[1]]
pathwayscores
df #avg CTRL and ASD pathway activity scores

#visualisation
scatterplot
activitybarplot
boxplot #per sample boxplot
progeny_hmap
NES_pathway
plot(scat_plots[[1]]$`TGFb`) 

#trial 2 'transcriptutorial'
PathwayActivity_counts
PathwayActivity_zscore_df
PathwayActivity_CARNIVALinput 

##### infer TF activities #####
##### 21.11.2025 #####

##### packages #####
BiocManager::install("OmnipathR")
BiocManager::install("viper")
library(decoupleR)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(progeny)
library(dorothea)
library(readr)
library(OmnipathR)
library(viper)

##### transcriptutorial decoupleR, ULM, collecTRI #####

#CollecTRI is a comprehensive resource containing a curated collection of TFs and their transcriptional 
#targets compiled from 12 different resources. 
#This collection provides an increased coverage of transcription factors and a superior performance in 
#identifying perturbed TFs compared to our previous DoRothEA network and other literature based GRNs. 
#Similar to DoRothEA, interactions are weighted by their mode of regulation (activation or inhibition

net <- decoupleR::get_collectri(organism = 'human', 
                                split_complexes = FALSE)
net

#To infer TF enrichment scores we will run the Univariate Linear Model (ulm) method. 
#For each sample in our dataset (mat) and each TF in our network (net), it fits a linear model that predicts 
#the observed gene expression based solely on the TF’s TF-Gene interaction weights. 
#Once fitted, the obtained t-value of the slope is the score. If it is positive, we interpret that the TF 
#is active and if it is negative we interpret that it is inactive.

#To run decoupleR methods, we need an input matrix (mat), an input prior knowledge network/resource (net), 
#and the name of the columns of net that we want to use.

# Run ulm
sample_acts <- decoupleR::run_ulm(mat = vst_progeny, 
                                  net = net, 
                                  .source = 'source', 
                                  .target = 'target',
                                  .mor = 'mor', 
                                  minsize = 5)
sample_acts

#visualisation
#From the obtained results we will observe the most variable activities across samples in a heat-map
n_tfs <- 25

# Transform to wide matrix
sample_acts_mat <- sample_acts %>%
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
top25_heatmap <- pheatmap::pheatmap(mat = sample_acts_mat,
                                    color = colors.use,
                                    border_color = "white",
                                    breaks = my_breaks,
                                    cellwidth = 15,
                                    cellheight = 15,
                                    treeheight_row = 20,
                                    treeheight_col = 20)
ggsave(
  filename = "top25_heatmap_ulm.png",
  plot = top25_heatmap,
  width = 12,
  height = 8,
  dpi = 300
)

#We can also infer TF activities from the statistics of the DEGs between ASD and CTRL:
# Run ulm
contrast_acts <- decoupleR::run_ulm(mat = degstat_mat[, 'statistic', drop = FALSE], 
                                    net = net, 
                                    .source = 'source', 
                                    .target = 'target',
                                    .mor='mor', 
                                    minsize = 5)
contrast_acts

# Filter top TFs in both signs
f_contrast_acts <- contrast_acts %>%
  dplyr::mutate(rnk = NA)

msk <- f_contrast_acts$score > 0

f_contrast_acts[msk, 'rnk'] <- rank(-f_contrast_acts[msk, 'score'])
f_contrast_acts[!msk, 'rnk'] <- rank(-abs(f_contrast_acts[!msk, 'score']))

tfs <- f_contrast_acts %>%
  dplyr::arrange(rnk) %>%
  head(n_tfs) %>%
  dplyr::pull(source)

f_contrast_acts <- f_contrast_acts %>%
  filter(source %in% tfs)

colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")[c(2, 10)])

p <- ggplot2::ggplot(data = f_contrast_acts, 
                     mapping = ggplot2::aes(x = stats::reorder(source, score), 
                                            y = score)) + 
  ggplot2::geom_bar(mapping = ggplot2::aes(fill = score),
                    color = "black",
                    stat = "identity") +
  ggplot2::scale_fill_gradient2(low = colors[1], 
                                mid = "whitesmoke", 
                                high = colors[2], 
                                midpoint = 0) + 
  ggplot2::theme_minimal() +
  ggplot2::theme(axis.title = element_text(face = "bold", size = 12),
                 axis.text.x = ggplot2::element_text(angle = 45, 
                                                     hjust = 1, 
                                                     size = 10, 
                                                     face = "bold"),
                 axis.text.y = ggplot2::element_text(size = 10, 
                                                     face = "bold"),
                 panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank()) +
  ggplot2::xlab("TFs")

p
ggsave(
  filename = "barplot_ASDvsCTRL.png",
  plot = p,
  width = 12,
  height = 8,
  dpi = 300
)

#We can further visualize the most differential target genes in each TF along their p-values to 
#interpret the results. For example, let’s see the genes that are belong to SP1
tf <- 'SP1'

df1 <- net %>%
  dplyr::filter(source == tf) %>%
  dplyr::arrange(target) %>%
  dplyr::mutate(ID = target, color = "3") %>%
  tibble::column_to_rownames('target')

inter <- sort(dplyr::intersect(rownames(deg_mat), rownames(df1)))

df1 <- df1[inter, ]

df1[,c('logfc', 't_value', 'p_value')] <- deg_mat[inter, ]

df1 <- df1 %>%
  dplyr::mutate(color = dplyr::if_else(mor > 0 & t_value > 0, '1', color)) %>%
  dplyr::mutate(color = dplyr::if_else(mor > 0 & t_value < 0, '2', color)) %>%
  dplyr::mutate(color = dplyr::if_else(mor < 0 & t_value > 0, '2', color)) %>%
  dplyr::mutate(color = dplyr::if_else(mor < 0 & t_value < 0, '1', color))

colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")[c(2, 10)])

p1 <- ggplot2::ggplot(data = df1, 
                      mapping = ggplot2::aes(x = logfc, 
                                             y = -log10(p_value), 
                                             color = color,
                                             size = abs(mor))) + 
  ggplot2::geom_point(size = 2.5, 
                      color = "black") + 
  ggplot2::geom_point(size = 1.5) +
  ggplot2::scale_colour_manual(values = c(colors[2], colors[1], "grey")) +
  ggrepel::geom_label_repel(mapping = ggplot2::aes(label = ID,
                                                   size = 1)) + 
  ggplot2::theme_minimal() +
  ggplot2::theme(legend.position = "none") +
  ggplot2::geom_vline(xintercept = 0, linetype = 'dotted') +
  ggplot2::geom_hline(yintercept = 0, linetype = 'dotted') +
  ggplot2::ggtitle(tf)

p1
ggsave(
  filename = "SP1_volcanoplot.png",
  plot = p1,
  width = 8,
  height = 8,
  dpi = 300
)

##### transcriptutorial with DOROTHEA 24.11.2025 #####
## We load Dorothea Regulons
data(dorothea_hs, package = "dorothea")
regulons <- dorothea_hs %>%
  dplyr::filter(confidence %in% c("A", "B","C"))

tf_activities_stat <- dorothea::run_viper(degstat_mat, regulons,
                                          options =  list(minsize = 5, eset.filter = FALSE, 
                                                          cores = 1, verbose = FALSE, nes = TRUE))
tf_activities_stat_top25 <- tf_activities_stat %>%
  as.data.frame() %>% 
  rownames_to_column(var = "GeneID") %>%
  dplyr::rename(NES = "statistic") %>%
  dplyr::top_n(25, wt = abs(NES)) %>%
  dplyr::arrange(NES) %>% 
  dplyr::mutate(GeneID = factor(GeneID))

tftop25 <- ggplot(tf_activities_stat_top25,aes(x = reorder(GeneID, NES), y = NES)) + 
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
  xlab("Transcription Factors")
ggsave(
  filename = "TFtop25_NES_dorothea.png",
  plot = tftop25,
  width = 10,
  height = 8,
  dpi = 300
)

targets_SP1 <- regulons$target[regulons$tf == "SP1"]
SP1 <- volcano_nice(as.data.frame(DEG_allstats[DEG_allstats$ID %in% targets_SP1,]), 
                    FCIndex = 3, pValIndex = 4, IDIndex = 1,nlabels = 20, label = TRUE, 
                    straight = FALSE) 
ggsave(
  filename = "SP1_dorothea.png",
  plot = SP1,
  width = 8,
  height = 12,
  dpi = 300
)

targets_NANOG <- regulons$target[regulons$tf == "NANOG"]
NANOG <- volcano_nice(as.data.frame(DEG_allstats[DEG_allstats$ID %in% targets_NANOG,]), 
                      FCIndex = 3, pValIndex = 4, IDIndex = 1,nlabels = 20, label = TRUE, 
                      straight = FALSE) 

tf_activities_CARNIVALinput<- tf_activities_stat %>%
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "TF") 
write.csv(tf_activities_CARNIVALinput, "tf_activities_CARNIVALinput.csv")

#tf activities per sample
tf_activities_counts <- 
  dorothea::run_viper(vst_progeny, regulons,
                      options =  list(minsize = 5, eset.filter = FALSE, 
                                      cores = 1, verbose = FALSE, method = c("scale")))

tf_activities_counts_filter <- tf_activities_counts %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "GeneID") %>%
  dplyr::filter(GeneID %in% tf_activities_stat_top25$GeneID) %>%
  column_to_rownames(var = "GeneID") %>%
  as.matrix()
tf_activities_vector <- as.vector(tf_activities_counts_filter)

paletteLength <- 100
myColor <- 
  colorRampPalette(c("darkblue", "whitesmoke","indianred"))(paletteLength)

dorotheaBreaks <- c(seq(min(tf_activities_vector), 0, 
                        length.out=ceiling(paletteLength/2) + 1),
                    seq(max(tf_activities_vector)/paletteLength, 
                        max(tf_activities_vector), 
                        length.out=floor(paletteLength/2)))
dorothea_hmap <- pheatmap(tf_activities_counts_filter,
                          fontsize=14, fontsize_row = 8, fontsize_col = 8, 
                          color=myColor, breaks = dorotheaBreaks,
                          main = "Dorothea ABC", angle_col = 45,
                          treeheight_col = 0,  border_color = NA)
ggsave(
  filename = "hmap_dorothea.png",
  plot = dorothea_hmap,
  width = 8,
  height = 10,
  dpi = 300
)
##### most important outputs here #####
net
sample_acts
contrast_acts
tf_activities_CARNIVALinput

#visualisation
top25_heatmap
p #barplot DEG ASD vs CTRL TFs
p1 #volcano plot SP1
tftop25
SP1
dorothea_hmap


##### support code #####
volcano_nice <- function (df, hAss = 0.05, FCIndex, pValIndex, IDIndex, vAss = NULL,
                          label = FALSE, straight = FALSE, nlabels, manual_labels = NA)
{
  df <- df[complete.cases(df), ]
  names(df)[1] <- "X1"
  hAssOri <- hAss
  hAss <- -log(hAss)
  names(df) <- gsub("adj.P.Val", "FDR", names(df))
  names(df)[FCIndex] <- "logFC"
  names(df)[pValIndex] <- "adj.P.Val"
  if (max(abs(df[, FCIndex])) >= 1) {
    xlimAbs <- ceiling(max(abs(df[, FCIndex])))
    ylimAbs <- ceiling(max(abs(-log(df[, pValIndex]))))
  }
  else {
    xlimAbs <- max(abs(df[, FCIndex]))
    ylimAbs <- max(abs(-log(df[, pValIndex])))
  }
  if (is.null(vAss)) {
    vAss <- xlimAbs/10
  }
  xneg <- function(x) abs(hAss - 1 + x/(x + vAss))
  xpos <- function(x) abs(hAss - 1 + x/(x - vAss))
  test <- function(x, y, vAss) {
    if (x < -vAss) {
      if (xneg(x) < -log(y)) {
        return("1")
      }
      else {
        return("0")
      }
    }
    else {
      if (x > vAss) {
        if (xpos(x) < -log(y)) {
          return("1")
        }
        else {
          return("0")
        }
      }
      else {
        return("0")
      }
    }
  }
  if (straight) {
    df$couleur <- ifelse(abs(df$logFC) >= vAss & df$adj.P.Val <=
                           hAssOri, "1", "0")
  }
  else {
    df$couleur <- "0"
    df$couleur <- apply(df, 1, FUN = function(x) test(as.numeric(x[FCIndex]),
                                                      as.numeric(x[pValIndex]), vAss))
  }
  df <- df[order(df$adj.P.Val, decreasing = F), ]
  df$condLabel <- df[, IDIndex]
  df[df$couleur == "0", "condLabel"] <- NA
  labels_to_keep <- c(df[c(1:nlabels), "condLabel"],manual_labels)
  df[!(df$condLabel %in% labels_to_keep), "condLabel"] <- NA
  df$couleur <- ifelse(df$couleur == "1" & df$logFC < 0, "2",
                       df$couleur)
  if (label) {
    a <- ggplot(df, aes(x = logFC, y = -log(adj.P.Val), color = couleur)) +
      geom_point(alpha = 0.5) + geom_label_repel(aes(label = condLabel)) +
      stat_function(fun = xneg, xlim = c(-xlimAbs, -vAss),
                    color = "black", alpha = 0.7) + ylim(c(0, ylimAbs)) +
      xlim(c(-xlimAbs, xlimAbs)) + stat_function(fun = xpos,
                                                 xlim = c(vAss, xlimAbs), color = "black", alpha = 0.7) +
      scale_colour_manual(values = c("grey30", "red",
                                     "royalblue3")) + theme_minimal() + theme(legend.position = "none")
  }
  else {
    if (straight) {
      a <- ggplot(df, aes(x = logFC, y = -log(adj.P.Val),
                          color = couleur)) + geom_point(alpha = 0.5) +
        geom_vline(xintercept = -vAss, color = "blue") +
        geom_vline(xintercept = vAss, color = "blue") +
        ylim(c(0, ylimAbs)) + xlim(c(-xlimAbs, xlimAbs)) +
        geom_hline(yintercept = hAss, color = "red") +
        scale_colour_manual(values = c("grey30", "red",
                                       "royalblue3")) + theme_minimal() + theme(legend.position = "none")
    }
    else {
      a <- ggplot(df, aes(x = logFC, y = -log(adj.P.Val),
                          color = couleur)) + geom_point(alpha = 0.5) +
        stat_function(fun = xneg, xlim = c(-xlimAbs,
                                           -vAss), color = "black", alpha = 0.7) + ylim(c(0,
                                                                                          ylimAbs)) + xlim(c(-xlimAbs, xlimAbs)) + stat_function(fun = xpos,
                                                                                                                                                 xlim = c(vAss, xlimAbs), color = "black", alpha = 0.7) +
        scale_colour_manual(values = c("grey30", "red",
                                       "royalblue3")) + theme_minimal() + theme(legend.position = "none")
    }
  }
  return(a)
}

##### CARNIVAL #####
##### 24.11.2025 #####
##### CARNIVAL transcriptutorial non-adapted or modified TEST! #####
#network reconstruction
source("support_functions.R")
BiocManager::install("CARNIVAL")
if (!require(remotes))
  install.packages("remotes")
remotes::install_github("dirkschumacher/rcbc")

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
library(ggrepel)
library(rcbc)
Sys.setenv(CBC_HOME = "C:/Users/diana/Downloads/Cbc-releases.2.10.12-windows-2022-msvs-v17-Release-x64/bin/cbc.exe") 
Sys.setenv(PATH = paste(Sys.getenv("PATH"), "C:/Users/diana/Downloads/Cbc-releases.2.10.12-windows-2022-msvs-v17-Release-x64/bin/cbc.exe", sep=";"))
Sys.getenv("CBC_HOME")
ls("package:rcbc")


## We also load the support functions
source("assignPROGENyScores.r")
source("generateTFList.r")
source("carnival_visNetwork.r")

## We read the normalised counts and the experimental design 
tf_activities <- read_csv("~/hnrnpu-causal-multiomics/processeddata/tf_activities_CARNIVALinput.csv")
PathwayActivity <- read_csv("~/hnrnpu-causal-multiomics/processeddata/PathwayActivity_CARNIVALinput.csv")

#create or upload scaffold network w omnipath
# need sif table format (node1, interaction, node2)
#therefore we use the consensus columns of direction (consensus_direction) and sign 
#(consensus_stimulation and consensus_inhibition) to extract it.

omniR <- import_omnipath_interactions()
#85217 interactions

# signed and directed
omnipath_sd <- omniR %>% dplyr::filter(consensus_direction == 1 &
                                         (consensus_stimulation == 1 | 
                                            consensus_inhibition == 1
                                         ))

# changing 0/1 criteria in consensus_stimulation/inhibition to -1/1
omnipath_sd$consensus_stimulation[which( omnipath_sd$consensus_stimulation == 0)] = -1
omnipath_sd$consensus_inhibition[which( omnipath_sd$consensus_inhibition == 1)] = -1
omnipath_sd$consensus_inhibition[which( omnipath_sd$consensus_inhibition == 0)] = 1

# check consistency on consensus sign and select only those in a SIF format
sif <- omnipath_sd[,c('source_genesymbol', 'consensus_stimulation', 'consensus_inhibition', 'target_genesymbol')] %>%
  dplyr::filter(consensus_stimulation==consensus_inhibition) %>%
  unique.data.frame()

sif$consensus_stimulation <- NULL
colnames(sif) <- c('source', 'interaction', 'target')

# remove complexes
sif$source <- gsub(":", "_", sif$source)
sif$target <- gsub(":", "_", sif$target)

#save SIF
write_tsv(sif, "~/Masters/RP2/hnrnpu-causal-multiomics/processeddata/omnipath_carnival.tsv")

#tf and pathway activities CARNIVAL
#We use the supplementary functions generateTFList.r and assignPROGENyScores.r to 
#shift the formats of tf_activities and PathwayActivity to the one required by CARNIVAL.

# dorothea for CARNIVAL
tf_activities_carnival <- data.frame(tf_activities, stringsAsFactors = F)
rownames(tf_activities_carnival) <- tf_activities$TF
tf_activities_carnival$TF <- NULL
tfList = generateTFList(tf_activities_carnival, top=50, access_idx = 1)
head(tfList)

# progeny for CARNIVAL
load(file = system.file("progenyMembers.RData",package="CARNIVAL"))

PathwayActivity_carnival <- data.frame(PathwayActivity, stringsAsFactors = F)
rownames(PathwayActivity_carnival) <- PathwayActivity_carnival$Pathway
PathwayActivity_carnival$Pathway <- NULL
progenylist = assignPROGENyScores(progeny = t(PathwayActivity_carnival), 
                                  progenyMembers = progenyMembers, 
                                  id = "gene", 
                                  access_idx = 1)

#CARNIVAL has been developed to find the causal link between the activities of the transcription factors (TFs) 
#and the ‘perturbed’ nodes. 
#3 inputs
#measObj: The TFs’ activities (like the ones we have obtained from DoRothEA)
#inputObj: The ‘perturbed’ nodes we want that CARNIVAL connects with the activity of TFs. 
#There are 3 ways of using it:
# 1. Give the name and sign of the selected nodes;
# 2. Give the name only, so the algorithm will select the sign that best fit the models,
# 3. Give NULL as value will create a “Perturbation” node that will try both signs for all ‘initial’ nodes of the given network ( netObj ).

#netObj: The network that will serve as map to connect the TFs’ activities ( measObj ) and the perturbed nodes ( inputObj )
#Although it is not required, a fourth object called weightObj can be also given. This object gives values ranged from -1 to 1 
#for a set of nodes of the network. The aim of weightObj is helping the solver to find optimal solutions faster.

#In the present example, we use assign as perturbation nodes all the “initial” nodes (option 2), 
#and as weightObj the PROGENy scores assigned to the most representative genes of the calculated pathways,

# get initial nodes
iniMTX = base::setdiff(sif$source, sif$target)
iniciators = base::data.frame(base::matrix(data = NaN, nrow = 1, ncol = length(iniMTX)), stringsAsFactors = F)
colnames(iniciators) = iniMTX

# If you used top=50 and access_idx = 1
tf_df <- tfList[[1]]              # first element of the list
measVec <- as.numeric(tf_df[1, ]) # extract numeric values
names(measVec) <- colnames(tf_df) # TF names

carnival_result = runCARNIVAL(inputObj= iniciators,
                              measObj = measVec, 
                              netObj = sif, 
                              weightObj = progenylist$score, 
                              solverPath = "C:/Users/diana/Downloads/Cbc-releases.2.10.12-windows-2022-msvs-v17-Release-x64/bin/cbc.exe", 
                              solver = "cbc",
                              timelimit=7200,
                              mipGAP=0,
                              poolrelGAP=0)




#transform to data.frame
carnival_result$weightedSIF <- data.frame(carnival_result$weightedSIF, stringsAsFactors = F)
carnival_result$weightedSIF$Sign <- as.numeric(carnival_result$weightedSIF$Sign)
carnival_result$weightedSIF$Weight <- as.numeric(carnival_result$weightedSIF$Weight)

carnival_result$nodesAttributes <- data.frame(carnival_result$nodesAttributes, stringsAsFactors = F)
carnival_result$nodesAttributes$ZeroAct <- as.numeric(carnival_result$nodesAttributes$ZeroAct)
carnival_result$nodesAttributes$UpAct <- as.numeric(carnival_result$nodesAttributes$UpAct)
carnival_result$nodesAttributes$DownAct <- as.numeric(carnival_result$nodesAttributes$DownAct)
carnival_result$nodesAttributes$AvgAct <- as.numeric(carnival_result$nodesAttributes$AvgAct)

saveRDS(carnival_result,"~/hnrnpu-causal-multiomics/processeddata/carnival_result.rds")

# visualization
visNet = carnival_visNet(evis = carnival_result$weightedSIF,
                         nvis = carnival_result$nodesAttributes)
#visNet
visSave(visNet, file = paste0('carnival_visualization_visNetwork.html'), selfcontained = TRUE)

##### analysis of CARNIVAL results #####
BiocManager::install("piano")
BiocManager::install("GSEABase")
library(readr)
library(piano)
library(dplyr)
library(ggplot2)
library(tibble)
library(tidyr)
library(dplyr)
library(scales)
library(plyr)
library(GSEABase)
library(network)
library(reshape2)
library(cowplot)
library(pheatmap)
library(ggraph)
library(tidygraph)


## We also load the support functions
source("support_enrichment.r")
source("support_networks.r")

## and the data

#read CARNIVAL results
carnival_result = readRDS("~/hnrnpu-causal-multiomics/processeddata/carnival_result.rds")
pkn = read_tsv("~/hnrnpu-causal-multiomics/processeddata/omnipath_carnival.tsv")
carnival_sample_resolution = readRDS("~/hnrnpu-causal-multiomics/carnival_sample_resolution.rds")

# Load pathways
pathways = gmt_to_csv("~/hnrnpu-causal-multiomics/c2.cp.v2025.1.Hs.symbols.gmt")

# Extract nodes and background
nodes_carnival = extractCARNIVALnodes(carnival_result)

# Run GSA hyper Geometric test
sig_pathways <- runGSAhyper(genes = nodes_carnival$sucesses, 
                            universe = nodes_carnival$bg, gsc = loadGSC(pathways))
sig_pathways_df <- as.data.frame(sig_pathways$resTab)  %>% 
  tibble::rownames_to_column(var = "pathway") 

#data for plotting
PathwaysSelect <- sig_pathways_df %>%
  dplyr::select(pathway, `p-value`, `Adjusted p-value`) %>%
  dplyr::filter(`Adjusted p-value` <= 0.05) %>%
  dplyr::rename(pvalue = `p-value`, AdjPvalu = `Adjusted p-value`) %>% 
  dplyr::mutate(pathway = as.factor(pathway))

PathwaysSelect <- data.frame(t(apply(PathwaysSelect, 1, function(r){
  aux = unlist(strsplit( sub("_",";", r["pathway"]), ";" ))
  r["pathway"] = gsub("_", " ", aux[2])
  return(c(r, "source" = aux[1]))
})))

colnames(PathwaysSelect) = c("pathway", "pvalue", "AdjPvalu", "source")
PathwaysSelect$AdjPvalu = as.numeric(PathwaysSelect$AdjPvalu)

ggdata = PathwaysSelect %>% 
  dplyr::filter(AdjPvalu <= 0.05) %>% 
  dplyr::group_by(source) %>% 
  dplyr::arrange(AdjPvalu) %>%
  dplyr::slice(1:5)


# Visualize top results
ggplot(ggdata, aes(y = reorder(pathway, AdjPvalu), x = -log10(AdjPvalu)), color = source) + 
  geom_bar(stat = "identity") +
  facet_grid(source ~ ., scales="free_y") +
  scale_x_continuous(
    expand = c(0.01, 0.01),
    limits = c(0, ceiling(max(-log10(PathwaysSelect$AdjPvalu)))),
    breaks = seq(floor(min(-log10(PathwaysSelect$AdjPvalu))), ceiling(max(-log10(PathwaysSelect$AdjPvalu))), 1),
    labels = math_format(10^-.x)
  ) +
  annotation_logticks(sides = "bt") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.y = element_text(size = 6)) +
  ylab("")

#trying to fix code
length(nodes_carnival$sucesses)
head(nodes_carnival$sucesses)

length(nodes_carnival$bg)
head(nodes_carnival$bg)