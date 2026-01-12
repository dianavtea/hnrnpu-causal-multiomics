##### D0 analysis #####

metafiltD0 <- meta %>%
  filter(
    !grepl("^si", Sample_Name),               # removes siNTC & siHNRNPU
    grepl("^CTRL|^HNRNPU", Sample_Name),       # keeps CTRL & HNRNPUdel
    Day %in% "D0"
  )
metafiltD0$condition <- factor(
  metafiltD0$StatusHNRNPU,
  levels = c("HNRNPU_WT", "HNRNPU_KD"),
  labels = c("CTRL", "HNRNPUdel")
)
countsfiltD0 <- merged_gene_counts_all[, colnames(merged_gene_counts_all) %in% metafiltD0$Sample_ID]

any(is.na(countsfiltD0)) #FALSE

library(DESeq2)
ddsD0 <- DESeqDataSetFromMatrix(countData = countsfiltD0,
                                colData = metafiltD0,
                                design= ~ condition)

ddsD0$condition<- relevel(ddsD0$condition, "CTRL")
ddsD0 <- DESeq(ddsD0)
resultsNames(ddsD0)
resD0 <- results(ddsD0)

resD0.mat <- cbind(counts(ddsD0, normalized=TRUE), resD0$log2FoldChange, resD0$padj)
mcols(resD0, use.names = TRUE)
saveRDS(ddsD0, file="ddsD0_ASDvsCTRL_D0.rds")
head(resD0)


resD0.mat <- cbind(
  norm_counts,                        # all normalized counts
  log2FC = resD0$log2FoldChange,
  padj  = resD0$padj
)
mcols(resD0, use.names = TRUE)
head(resD0.mat)

##### IGNORE ABOVE #####
norm_counts <- counts(dds, normalized=TRUE)

saveRDS(ddsD0, file="ddsD0_ASDvsCTRL_D0.rds")

library(biomaRt)
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
write.csv(resD0.annot, "DESeq2_resD0.csv", row.names = FALSE)


###### VST ######
vst_all <- vst(dds, blind = FALSE)
vstD0 <- vst_all[, colData(dds)$Day == "D0"]
vstmatD0 <- assay(vstD0) #212620 elements
dim(vstmatD0) #21262 transcripts x 10 samples

# Ensure rows of res.annot match vstmat
resD0.annot <- resD0.annot[match(rownames(vstmatD0), resD0.annot$ensembl), ]

vstD0_annot <- cbind(resD0.annot[, c("ensembl", "entrez", "hgnc_symbol")], vstmatD0)
head(vstD0_annot)
write.csv(vstD0_annot, "VSTcounts_annot_D0.csv", row.names = FALSE)

head(vstmatD0)
#vst for progeny
vst_progD0 <- vstmatD0
rownames(vst_progD0) <- resD0.annot$hgnc_symbol
head(vst_progD0)
any(is.na(vst_progD0)) #false

rownames(vst_progD0)[duplicated(rownames(vst_progD0))] 
vst_progD0 <- vst_progD0[rownames(vst_progD0) != "", ] #remove blanks
vst_progD0 <- vst_progD0[!duplicated(rownames(vst_progD0)), ] #remove dups

any(rownames(vst_progD0) == "")       # should be FALSE
any(duplicated(rownames(vst_progD0))) # should be FALSE


#num_samples <- ncol(vst_progdf_D0) - 1 

# Select the last column (the ID column) first, followed by the original sample columns (1 to num_samples)
#vst_progdf_D0 <- vst_progdf_D0[, c(ncol(vst_progdf_D0), 1:num_samples)]
head(vst_progD0)
write.csv(vst_progD0, "~/hnrnpu-causal-multiomics/processeddata/VST_progD0.csv", row.names = TRUE)
write.csv(metafiltD0, "~/hnrnpu-causal-multiomics/processeddata/metadataD0.csv", row.names = FALSE)

plotPCA(vstD0, intgroup = "condition")
#73% variance explained in PC1, 18% in PC2, distinction between ASD and CTRL groups

library(ggplot2)
# PCA plot
pca_D0 <- plotPCA(vstD0, intgroup = "condition") +
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
  filename = "PCA_d0.png",      # use PDF or PNG
  plot = pca_D0,           
  width = 6, height = 5,          # size in inches
  dpi = 300                       # resolution for PNG (ignored for PDF)
)

##### prep for inferring TF activities #####

##### collectri #####
net <- decoupleR::get_collectri(organism='human', split_complexes=FALSE)
net #43159 - source, target, mor

resD0.annotmat <- resD0.annot[, c("hgnc_symbol", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
head(resD0.annotmat)

#Remove rows with blank hgnc_symbol
resD0.annotmat <- resD0.annotmat[resD0.annotmat$hgnc_symbol != "", ]
# Remove duplicated hgnc_symbol rows (keeping the first one)
resD0.annotmat <- resD0.annotmat[!duplicated(resD0.annotmat$hgnc_symbol), ]
row.names(resD0.annotmat) <- resD0.annotmat$hgnc_symbol
#17017 transcripts

names(resD0.annotmat)[1] <- "ID"
head(resD0.annotmat)

resD0.annotmat <- unique(resD0.annotmat)

#make into numeric factor
resD0.filt <- resD0.annotmat$stat
names(resD0.filt) <- resD0.annotmat$ID

resD0.filt <- resD0.filt[!is.na(resD0.filt)]
any(is.na(resD0.filt))
head(resD0.filt)

TFscores_D0 <- decoupleR::run_viper(
  mat = resD0.filt, #statistic matrix 
  net = net,                 
  .source = 'source', 
  .target = 'target', 
  .mor = 'mor', 
  minsize = 25
)
TFscores_D0$condition <- "HNRNPUdel_vs_CTRL"
#278 TF scores
head(TFscores_D0)

TFscoresD0 <- TFscores_D0[, c('source', 'score')]
names(TFscoresD0) <- c("TF", "score")
head(TFscoresD0)

TFscoresD0 <- as.data.frame(TFscoresD0)
row.names(TFscoresD0) <- TFscoresD0$TF
TFscoresD0 <- TFscoresD0 %>%
  tibble::column_to_rownames("TF") %>%
  as.matrix()
head(TFscoresD0)
write.csv(TFscoresD0, "~/hnrnpu-causal-multiomics/processeddata/TFscoresD0.csv", row.names = TRUE) #278

n_tfs <- 25

library(dplyr)
## top 25 TF in general (abs)
# Filter for top 25 most active TFs (by absolute score)
top_tfsD0 <- TFscores_D0 %>%
  arrange(desc(abs(score))) %>%
  slice_head(n = 25)

# Plot them directly 
top25allD0 <- ggplot(top_tfsD0, aes(x = reorder(source, score), y = score, fill = score > 0)) +
  geom_col() +
  coord_flip() + # Makes it horizontal
  theme_minimal() +
  labs(title = "Top 25 TF Activities (ASD vs CTRL)", 
       y = "Normalized Enrichment Score (NES)", 
       x = "Transcription Factor") +
  scale_fill_manual(values = c("TRUE" = "#E41A1C", "FALSE" = "#377EB8"), 
                    labels = c("Activated", "Inhibited"),
                    name = "Status")

ggsave(
  filename = "top25abs_D0.png",      # use PDF or PNG
  plot = top25allD0,           
  width = 8, height = 10,          # size in inches
  dpi = 300                       # resolution for PNG (ignored for PDF)
)

##top 15 for activated and top 15 for inhibited TFs
# Select Top 15 Activated (Highest Positive Scores)
top_up <- TFscores_D0 %>%
  arrange(desc(score)) %>%
  slice_head(n = 15)

# Select Top 15 Inhibited (Lowest Negative Scores)
top_down <- TFscores_D0 %>%
  arrange(score) %>%
  slice_head(n = 15)

# Combine them into one dataframe
top_30_balanced <- bind_rows(top_up, top_down)

# We use reorder() to ensure they appear sorted by score, not alphabetically
top30D0 <- ggplot(top_30_balanced, aes(x = reorder(source, score), y = score, fill = score > 0)) +
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
  filename = "top30TFs_D0.png",      # use PDF or PNG
  plot = top30D0,           
  width = 10, height = 8,          # size in inches
  dpi = 300                       # resolution for PNG (ignored for PDF)
)

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
Normalised_counts <- read_csv("~/hnrnpu-causal-multiomics/processeddata/VST_progD0.csv") 
Experimental_design <- read_csv("~/hnrnpu-causal-multiomics/processeddata/metadataD0.csv")

## We read the results from the differential analysis. 
DEAD0 <- read_csv("~/hnrnpu-causal-multiomics/processeddata/DESeq2_resD0.csv") #21262 genes

Normalised_counts_matrix <- as.matrix(
  column_to_rownames(
    as.data.frame(Normalised_counts[, -2]),
    var = "...1"
  )
)

head(Normalised_counts_matrix)

DEAD0 <- DEAD0 %>%
  mutate(hgnc_symbol = as.character(hgnc_symbol)) %>%
  mutate(hgnc_symbol = trimws(hgnc_symbol))

DEAD0_mat <- as.matrix(DEAD0$stat)
rownames(DEAD0_mat) <- DEAD0$hgnc_symbol
colnames(DEAD0_mat) <- "statistic"
head(DEAD0_mat)

PathwayActivity_countsD0 <- progeny(Normalised_counts_matrix, scale = TRUE,
                                   organism = "Human", top = 100)
Activity_countsD0 <- as.vector(PathwayActivity_countsD0)

paletteLength <- 100
myColor <- 
  colorRampPalette(c("darkblue", "whitesmoke","indianred"))(paletteLength)

progenyBreaksD0 <- c(seq(min(Activity_countsD0), 0, 
                        length.out=ceiling(paletteLength/2) + 1),
                    seq(max(Activity_countsD0)/paletteLength, 
                        max(Activity_countsD0), 
                        length.out=floor(paletteLength/2)))

annotation_col <- Experimental_design %>%
  dplyr::select(Sample_ID, condition) %>%
  as.data.frame()

rownames(annotation_col) <- annotation_col$Sample_ID
annotation_col$Sample_ID <- NULL

progeny_hmapD0 <- pheatmap(
  t(PathwayActivity_countsD0),
  fontsize = 14,
  fontsize_row = 10,
  fontsize_col = 10,
  color = myColor,
  breaks = progenyBreaksD0,
  main = "PROGENy (100)",
  angle_col = 45,
  treeheight_col = 0,
  border_color = NA,
  annotation_col = annotation_col
)
ggsave(
  filename = "PROGENY_heatmapD0.png",
  plot = progeny_hmapD0,
  width = 8,
  height = 12,
  dpi = 300
)

#enrichment analysis using competitive permutation approach to assess
#significance of pathway activity to end with NES for each pathway
PathwayActivity_zscoreD0 <- progeny(DEAD0_mat, 
                                   scale=TRUE, organism="Human", top = 100, perm = 10000, z_scores = TRUE) %>%
  t()
colnames(PathwayActivity_zscoreD0) <- "NES"

PathwayActivity_zscoreD0_df <- as.data.frame(PathwayActivity_zscoreD0) %>% 
  rownames_to_column(var = "Pathway") %>%
  dplyr::arrange(NES) %>%
  dplyr::mutate(Pathway = factor(Pathway))

NES_pathwayD0 <- ggplot(PathwayActivity_zscoreD0_df,aes(x = reorder(Pathway, NES), y = NES)) + 
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
  filename = "PROGENY_NES_pathwaysD0.png",
  plot = NES_pathwayD0,
  width = 10,
  height = 8,
  dpi = 300
)

# pathway most active - visualise most responsive genes (progeny weights)
#along with their statistic values to interpret results
#scatterplot shows which genes contributing most to the pathway enrichment

prog_matrixD0 <- getModel("Human", top=100) %>% 
  as.data.frame()  %>%
  tibble::rownames_to_column("GeneID")

DEAD0_df <- DEAD0_mat %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("GeneID")

sum(is.na(DEAD0_df$statistic))
DEAD0_df <- DEAD0_df %>%
  dplyr::filter(!is.na(statistic))

scat_plotsD0 <- progeny::progenyScatter(df = DEAD0_df, 
                                       weight_matrix = prog_matrixD0, 
                                       statName = "statistic", verbose = FALSE)
plot(scat_plotsD0[[1]]$`TGFb`) 

#progeny results as input for CARNIVAL
#CARNIVAL sets weights based on PROGENy scores in each pathway-related node
#to find more relevant solutions
#set progeny zscores false to return pathway activity values between 1 and -1

PathwayActivity_CARNIVALinputD0 <- progeny(DEAD0_mat, 
                                          scale=TRUE, organism="Human", top = 100, perm = 10000, z_scores = FALSE) %>%
  t () %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "Pathway") 
colnames(PathwayActivity_CARNIVALinputD0)[2] <- "score"
head(PathwayActivity_CARNIVALinputD0)
write_csv(PathwayActivity_CARNIVALinputD0, 
          "~/hnrnpu-causal-multiomics/processeddata/PathwayActivity_CARNIVALinputD0.csv")


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

measObj <- viperRes$NES
names(measObj) <- viperRes$ID
is.numeric(measObj)
head(measObj)
head(names(measObj))

weightObj <- PROGENY$score
names(weightObj) <- PROGENY$Pathway

###### omnipath ######

library(OmnipathR)

full_pkn <- as.data.frame(omnipath_interactions()) #85217 interactions
full_pkn <- full_pkn[!is.na(full_pkn$references),]
head(full_pkn)
clean_PKN <- full_pkn[full_pkn$consensus_stimulation == 1 | full_pkn$consensus_inhibition == 1,]

clean_PKN$sign <- clean_PKN$consensus_stimulation - clean_PKN$consensus_inhibition

clean_PKN <- clean_PKN[,c(3,4,16)]
head(clean_PKN)

clean_PKN_supp <- clean_PKN[clean_PKN$sign == 0,]
clean_PKN_supp$sign <- -1
clean_PKN[clean_PKN$sign == 0,"sign"] <- 1

clean_PKN <- as.data.frame(rbind(clean_PKN, clean_PKN_supp))
#71273 interactions

names(clean_PKN) <- c("source","target","interaction")
write_csv(x = clean_PKN, file = "clean_omnipath_PKN.csv")


head(TFscoresD0)
TFscoresD0_CARNIVAL <- TFscoresD0[, "score"]
class(TFscoresD0_CARNIVAL)
# numeric
names(TFscoresD0_CARNIVAL)
# "AHR" "AP1" "APEX1" ...


head(PathwayActivity_CARNIVALinputD0)
Pathway_CARNIVALD0 <- PathwayActivity_CARNIVALinputD0$score
names(Pathway_CARNIVALD0) <- PathwayActivity_CARNIVALinputD0$Pathway
head(Pathway_CARNIVALD0)
class(Pathway_CARNIVALD0)
# numeric

names(Pathway_CARNIVALD0)
# "Androgen" "EGFR" "Estrogen" "Hypoxia" "JAK.STAT" "MAPK"

##### CARNIVAL #####

sif_clean <- read.csv("clean_omnipath_PKN.csv")
clean_PKN

?runCARNIVAL
carnival_D0 <- runCARNIVAL(
  inputObj  = NULL,             
  measObj   = TFscoresD0_CARNIVAL,             # TF activity measurements
  netObj    = clean_PKN,                    # prior knowledge network
  weightObj = Pathway_CARNIVALD0,         # pathway activity scores (or progeny_matrix)
  solverPath = "C:/Program Files/IBM/ILOG/CPLEX_Studio2211/cplex/bin/x64_win64/cplex.exe", 
  solver     = "cplex",
  timelimit  = 7200, 
  mipGAP     = 0.05,
  poolrelGAP = 0
)
saveRDS(carnival_D0,"~/hnrnpu-causal-multiomics/processeddata/carnival_D0.rds")

head(carnival_D0)

#inputObj: data frame of list fir target of perturbation, OPTIONAL or default set to NULL to run invCARNIVAL when inputs are not known

## trial
# Load libraries
library(igraph)
library(ggraph)
library(tidyverse)
library(dplyr)

# 1. Prepare edges
edges <- carnival_D0$weightedSIF %>%
  filter(!is.na(Node1) & !is.na(Node2)) %>% # remove incomplete rows
  dplyr::select(Node1, Node2)

# 2. Prepare node attributes
nodes <- carnival_D0$nodesAttributes %>%
  mutate(Node = Node) %>%
  dplyr::select(Node, AvgAct)

# 3. Create igraph object
g <- graph_from_data_frame(d = edges, vertices = nodes, directed = TRUE)

# 4. Compute node degree (for sizing)
V(g)$degree <- degree(g, mode = "all")

# 5. Map activity to colors
# AvgAct: 0 = inactive (blue), 1 = active (red), interpolate if needed
V(g)$color <- scales::col_numeric(
  palette = c("blue", "white", "red"),
  domain = c(0, 1)
)(V(g)$AvgAct)

# 6. Map degree to node size
V(g)$size <- scales::rescale(V(g)$degree, to = c(3, 10))  # smaller to larger

# 7. Plot network using ggraph
ggraph(g, layout = "fr") +  # Fruchterman-Reingold layout
  geom_edge_link(alpha = 0.4, arrow = arrow(length = unit(4, 'mm')), end_cap = circle(3, 'mm')) +
  geom_node_point(aes(size = size, color = color)) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  scale_color_identity() +
  theme_void() +
  ggtitle("CARNIVAL Network: Node Activity Visualization")

### trial 2
top_nodes <- carnival_D0$nodesAttributes %>%
  arrange(desc(AvgAct)) %>%
  dplyr::slice(1:50) %>%
  pull(Node)

sub_edges <- carnival_D0$weightedSIF %>%
  filter(Node1 %in% top_nodes | Node2 %in% top_nodes) %>%
  dplyr::select(Node1, Node2)

g_sub <- graph_from_data_frame(sub_edges, directed = TRUE)

nodes_sub <- carnival_D0$nodesAttributes %>%
  filter(Node %in% V(g_sub)$name) %>%
  dplyr::select(Node, AvgAct)

# 5️⃣ Add attributes to igraph
V(g_sub)$AvgAct <- nodes_sub$AvgAct[match(V(g_sub)$name, nodes_sub$Node)]
V(g_sub)$degree <- degree(g_sub, mode = "all")
V(g_sub)$color <- scales::col_numeric(c("blue", "white", "red"), domain = c(0,1))(V(g_sub)$AvgAct)
V(g_sub)$size <- scales::rescale(V(g_sub)$degree, to = c(4,10))

# 6️⃣ Plot subnetwork
ggraph(g_sub, layout = "fr") +  # Fruchterman-Reingold layout
  geom_edge_link(alpha = 0.4, arrow = arrow(length = unit(4, 'mm')), end_cap = circle(3, 'mm')) +
  geom_node_point(aes(size = size, color = color)) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  scale_color_identity() +
  theme_void() +
  ggtitle("Top 50 Active Nodes Subnetwork")

#### cytoscape prep
library(dplyr)

library(dplyr)
library(stringr)

# Clean edges
edges <- carnival_D0$weightedSIF %>%
  dplyr::select(Node1, Node2, Weight) %>%
  mutate(
    Node1 = str_trim(Node1),
    Node2 = str_trim(Node2)
  )

write.csv(edges, "carnival_edgesD0.csv", row.names = FALSE)

# Clean nodes and rename AvgAct
nodes <- carnival_D0$nodesAttributes %>%
  dplyr::select(Node, AvgAct) %>%
  mutate(Node = str_trim(Node)) %>%  # remove extra spaces
  dplyr::rename(id = Node,
         activity_D0 = AvgAct)

write.csv(nodes, "carnival_nodesD0.csv", row.names = FALSE)



##### MOON #####
library(cosmosR)

#TF scores for MOON
TFscoresMOON <- TFscores_D0[, c('source', 'score')]
names(TFscoresMOON) <- c("TF", "score")
head(TFscoresMOON)
TFscoresMOON <- as.data.frame(TFscoresMOON)
row.names(TFscoresMOON) <- TFscoresMOON$TF
head(TFscoresMOON)

clean_PKN

#DEA for MOON
DEA_MOON <- resD0.annot[, c("hgnc_symbol", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
head(DEA_MOON)
#data frame, 13821 transcripts

# Remove rows with blank hgnc_symbol
DEA_MOON <- DEA_MOON[DEA_MOON$hgnc_symbol != "", ]
# Remove duplicated hgnc_symbol rows (keeping the first one)
DEA_MOON <- DEA_MOON[!duplicated(DEA_MOON$hgnc_symbol), ]
row.names(DEA_MOON) <- DEA_MOON$hgnc_symbol
#12878 transcripts

names(DEA_MOON)[1] <- "gene"
row.names(DEA_MOON) <- DEA_MOON$gene
head(DEA_MOON)
DEA_clean_MOON <- DEA_MOON[, "stat", drop = FALSE] 

# Rename the column to represent comparison
names(DEA_clean_MOON) <- "ASD_vs_CTRL"
DEA_moon <- DEA_clean_MOON
head(DEA_moon)
conditions_of_interest <- "ASD_vs_CTRL"

head(DEA_moon)

PKND0 <- cosmosR:::filter_pkn_expressed_genes(row.names(DEA_moon), meta_pkn = clean_PKN)

#####
DEA_moon 
obsgenes <- TFscoresMOON$TF
PKN

head(names(TFscoresMOON))
head(obsgenes)

if(!"source" %in% names(PKND0)) { names(PKND0)[1] <- "source" }
if(!"target" %in% names(PKND0)) { names(PKND0)[2] <- "target" }

network_nodesD0 <- unique(c(PKND0$source, PKND0$target))
overlap_countD0 <- sum(obsgenes %in% network_nodesD0)
print(paste("Original genes to observe:", length(obsgenes)))
print(paste("Genes found in the network:", overlap_countD0))

obsgenes_D0 <- obsgenes[obsgenes %in% network_nodesD0]
# Prune the Network (Keep only nodes upstream of your TFs)
n_steps <- 6
PKNfiltD0 <- cosmosR:::keep_observable_neighbours(PKND0, n_steps, obsgenes_D0)

# Compress the network (Simplify redundant paths)
meta_network_compressed_listD0 <- compress_same_children(PKNfiltD0, sig_input = NULL, metab_input = obsgenes_D0)
meta_network_compressed_D0 <- meta_network_compressed_listD0$compressed_network

node_signaturesD0 <- meta_network_compressed_listD0$node_signaturesD0
duplicated_parentsD0 <- meta_network_compressed_listD0$duplicated_signatures
meta_network_compressed_D0 <- meta_network_cleanup(meta_network_compressed_D0)


# 7. The Main Loop (The "MOON" Algorithm)
# This iterates to find the best signaling path explaining your data
meta_network_TFD0 <- meta_network_compressed_D0
i <- 1
before <- 1
after <- 0

load("collectri_regulon_R.RData")

TFscoresMOON_vec <- TFscoresMOON$score
names(TFscoresMOON_vec) <- TFscoresMOON$TF

print("Starting COSMOS optimization...")
while (before != after & i < 10) {
  before <- length(meta_network_TFD0[,1])
  
  # Run MOON: Upstream=NULL because we don't know the drug target
  moon_resD0 <- moon(upstream_input = NULL, 
                   downstream_input = TFscoresMOON_vec, 
                   meta_network = meta_network_TFD0, 
                   n_layers = n_steps, 
                   statistic = "ulm") 
  
  # Prune incoherent links (contradictions between Kinase -> TF -> Gene)
  meta_network_TFD0 <- filter_incohrent_TF_target(moon_resD0, net, meta_network_TFD0, DEA_moon)
  
  after <- length(meta_network_TFD0[,1])
  i <- i + 1
  print(paste("Iteration:", i, "- Network size:", after))
}

getAnywhere("filter_incohrent_TF_target")

nrow(meta_network_TFD0)
head(meta_network_TFD0)


# Look at the top regulators (Kinases/Signaling Proteins)
top_regulators <- moon_resD0 %>%
  dplyr::arrange(desc(abs(score))) # Sort by magnitude of score

# Display the top 10 potential drivers of your disease
print(head(top_regulators, 10))

# Visualize it nicely
ggplot(head(top_regulators, 20), aes(x = reorder(source, score), y = score, fill = score > 0)) +
  geom_col() +
  coord_flip() +
  theme_minimal() +
  labs(title = "Top Upstream Regulators (Disease vs Control)",
       x = "Signaling Protein",
       y = "Predicted Activity Score") +
  scale_fill_manual(values = c("TRUE" = "firebrick", "FALSE" = "steelblue"), 
                    labels = c("Activated", "Inhibited"))


moonMAP2K2 <- get_moon_scoring_network("MAP2K2", meta_network_TFD0, moon_resD0)
library(igraph)


# Create igraph object
g_moon <- graph_from_data_frame(
  d = moonMAP2K2$SIF,
  directed = TRUE,
  vertices = moonMAP2K2$ATT  # optional, adds attributes
)

plot(
  g_moon,
  vertex.size = 6,
  vertex.label.cex = 0.7,
  vertex.color = "steelblue",
  edge.arrow.size = 0.3
)

library(ggraph)
library(ggplot2)

# Optionally, color nodes by MOON score
V(g_moon)$score <- moonTUB$ATT$score[match(V(g_moon)$name, moonTUB$ATT$source)]

ggraph(g_moon, layout = "fr") +
  geom_edge_link(aes(edge_alpha = 0.5, edge_color = factor(E(g_moon)$interaction))) +
  geom_node_point(aes(color = score), size = 5) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  scale_color_gradient2(low = "steelblue", mid = "white", high = "firebrick", midpoint = 0) +
  scale_edge_color_manual(values = c("1" = "darkgreen", "-1" = "red")) +
  theme_void() +
  labs(title = "MOON Network for MAP2K2")


