##### D28 analysis #####

library(biomaRt)
#Connect to GRCh37 Ensembl - annotations
# Use rownames of res directly
ensembl_ids <- rownames(resD28)

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

resD28.df <- as.data.frame(resD28)
resD28.df$ensembl <- ensembl_ids

resD28.annot <- merge(
  resD28.df,
  annot,
  by.x = "ensembl",
  by.y = "ensembl_gene_id",
  all.x = TRUE,
  sort = FALSE
)

resD28.annot <- resD28.annot[match(ensembl_ids, resD28.annot$ensembl), ]

colnames(resD28.annot)[colnames(resD28.annot) == "entrezgene_id"] <- "entrez"
colnames(resD28.annot)[colnames(resD28.annot) == "hgnc_symbol"]   <- "hgnc_symbol"
head(resD28.annot)
write.csv(resD28.annot, "DESeq2_resD28.csv", row.names = FALSE)


###### VST ######
vst_all <- vst(dds, blind = FALSE)
vstD28 <- vst_all[, colData(dds)$Day == "D28"]
vstmatD28 <- assay(vstD28) #212620 elements
dim(vstmatD28) #21262 transcripts x 6 samples

# Ensure rows of res.annot match vstmat
resD28.annot <- resD28.annot[match(rownames(vstmatD28), resD28.annot$ensembl), ]

vstD28_annot <- cbind(resD28.annot[, c("ensembl", "entrez", "hgnc_symbol")], vstmatD28)
head(vstD28_annot)
write.csv(vstD28_annot, "VSTcounts_annot_D28.csv", row.names = FALSE)

head(vstmatD28)
#vst for progeny
vst_progD28 <- vstmatD28
rownames(vst_progD28) <- resD28.annot$hgnc_symbol
head(vst_progD28)
any(is.na(vst_progD28)) #false

rownames(vst_progD28)[duplicated(rownames(vst_progD28))] 
vst_progD28 <- vst_progD28[rownames(vst_progD28) != "", ] #remove blanks
vst_progD28 <- vst_progD28[!duplicated(rownames(vst_progD28)), ] #remove dups

any(rownames(vst_progD28) == "")       # should be FALSE
any(duplicated(rownames(vst_progD28))) # should be FALSE


#num_samples <- ncol(vst_progdf_D28) - 1 

# Select the last column (the ID column) first, followed by the original sample columns (1 to num_samples)
#vst_progdf_D28 <- vst_progdf_D28[, c(ncol(vst_progdf_D28), 1:num_samples)]
head(vst_progD28)
write.csv(vst_progD28, "~/hnrnpu-causal-multiomics/processeddata/VST_progD28.csv", row.names = TRUE)


plotPCA(vstD28, intgroup = "condition")
#82% variance explained in PC1, 10% in PC2, distinction between ASD and CTRL groups

library(ggplot2)
# PCA plot
pca_D28 <- plotPCA(vstD28, intgroup = "condition") +
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
  filename = "PCA_D28.png",      # use PDF or PNG
  plot = pca_D28,           
  width = 6, height = 5,          # size in inches
  dpi = 300                       # resolution for PNG (ignored for PDF)
)

##### prep for inferring TF activities #####

##### collectri #####
net <- decoupleR::get_collectri(organism='human', split_complexes=FALSE)
net #43159 - source, target, mor

resD28.annotmat <- resD28.annot[, c("hgnc_symbol", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
head(resD28.annotmat)

#Remove rows with blank hgnc_symbol
resD28.annotmat <- resD28.annotmat[resD28.annotmat$hgnc_symbol != "", ]
# Remove duplicated hgnc_symbol rows (keeping the first one)
resD28.annotmat <- resD28.annotmat[!duplicated(resD28.annotmat$hgnc_symbol), ]
row.names(resD28.annotmat) <- resD28.annotmat$hgnc_symbol
#17017 transcripts

names(resD28.annotmat)[1] <- "ID"
head(resD28.annotmat)

resD28.annotmat <- unique(resD28.annotmat)

#make into numeric factor
resD28.filt <- resD28.annotmat$stat
names(resD28.filt) <- resD28.annotmat$ID

resD28.filt <- resD28.filt[!is.na(resD28.filt)]
any(is.na(resD28.filt))
head(resD28.filt)

TFscores_D28 <- decoupleR::run_viper(
  mat = resD28.filt, #statistic matrix 
  net = net,                 
  .source = 'source', 
  .target = 'target', 
  .mor = 'mor', 
  minsize = 25
)
TFscores_D28$condition <- "HNRNPUdel_vs_CTRL"
#280 TF scores
head(TFscores_D28)

TFscoresD28 <- TFscores_D28[, c('source', 'score')]
names(TFscoresD28) <- c("TF", "score")
head(TFscoresD28)

TFscoresD28 <- as.data.frame(TFscoresD28)
TFscoresD28 <- TFscoresD28 %>%
  tibble::column_to_rownames("TF") %>%
  as.matrix()
head(TFscoresD28)
write.csv(TFscoresD28, "~/hnrnpu-causal-multiomics/processeddata/TFscoresD28.csv", row.names = TRUE) #278

n_tfs <- 25

library(dplyr)
## top 25 TF in general (abs)
# Filter for top 25 most active TFs (by absolute score)
top_tfsD28 <- TFscores_D28 %>%
  arrange(desc(abs(score))) %>%
  slice_head(n = 25)

# Plot them directly 
top25allD28 <- ggplot(top_tfsD28, aes(x = reorder(source, score), y = score, fill = score > 0)) +
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
  filename = "top25abs_D28.png",      # use PDF or PNG
  plot = top25allD28,           
  width = 8, height = 10,          # size in inches
  dpi = 300                       # resolution for PNG (ignored for PDF)
)

##top 15 for activated and top 15 for inhibited TFs
# Select Top 15 Activated (Highest Positive Scores)
top_up <- TFscores_D28 %>%
  arrange(desc(score)) %>%
  slice_head(n = 15)

# Select Top 15 Inhibited (Lowest Negative Scores)
top_down <- TFscores_D28 %>%
  arrange(score) %>%
  slice_head(n = 15)

# Combine them into one dataframe
top_30_balanced <- bind_rows(top_up, top_down)

# We use reorder() to ensure they appear sorted by score, not alphabetically
top30D28 <- ggplot(top_30_balanced, aes(x = reorder(source, score), y = score, fill = score > 0)) +
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
  filename = "top30TFs_D28.png",      # use PDF or PNG
  plot = top30D28,           
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
Normalised_counts <- read_csv("~/hnrnpu-causal-multiomics/processeddata/VST_progD28.csv") 

metafiltD28 <- meta_filtered[21:26,]
Experimental_design <- metafiltD28

## We read the results from the differential analysis. 
DEAD28 <- read_csv("~/hnrnpu-causal-multiomics/DESeq2_resD28.csv") #21262 genes

Normalised_counts_matrix <- as.matrix(
  data.frame(Normalised_counts, row.names = Normalised_counts$...1)[, -1]
)


head(Normalised_counts_matrix)

DEAD28 <- DEAD28 %>%
  mutate(hgnc_symbol = as.character(hgnc_symbol)) %>%
  mutate(hgnc_symbol = trimws(hgnc_symbol))

DEAD28_mat <- as.matrix(DEAD28$stat)
rownames(DEAD28_mat) <- DEAD28$hgnc_symbol
colnames(DEAD28_mat) <- "statistic"
head(DEAD28_mat)

PathwayActivity_countsD28 <- progeny(Normalised_counts_matrix, scale = TRUE,
                                    organism = "Human", top = 100)
Activity_countsD28 <- as.vector(PathwayActivity_countsD28)

paletteLength <- 100
myColor <- 
  colorRampPalette(c("darkblue", "whitesmoke","indianred"))(paletteLength)

progenyBreaksD28 <- c(seq(min(Activity_countsD28), 0, 
                         length.out=ceiling(paletteLength/2) + 1),
                     seq(max(Activity_countsD28)/paletteLength, 
                         max(Activity_countsD28), 
                         length.out=floor(paletteLength/2)))

annotation_col <- Experimental_design %>%
  dplyr::select(Sample_ID, condition) %>%
  as.data.frame()

rownames(annotation_col) <- annotation_col$Sample_ID
annotation_col$Sample_ID <- NULL

# force condition order
annotation_col$condition <- factor(
  annotation_col$condition,
  levels = c("CTRL", "HNRNPUdel")
)

# order samples
sample_order <- rownames(annotation_col)[order(annotation_col$condition)]

# reorder matrix
mat <- t(PathwayActivity_countsD28)[, sample_order, drop = FALSE]

# reorder annotation (IMPORTANT: keep as data.frame)
annotation_col2 <- annotation_col[sample_order, , drop = FALSE]

# plot
progeny_hmapD28 <- pheatmap(
  mat,
  fontsize = 14,
  fontsize_row = 10,
  fontsize_col = 10,
  color = myColor,
  breaks = progenyBreaksD28,
  main = "PROGENy (100)",
  angle_col = 45,
  treeheight_col = 0,
  border_color = NA,
  annotation_col = annotation_col2,
  cluster_cols = FALSE,
  gaps_col = sum(annotation_col2$condition == "CTRL")
)


ggsave(
  filename = "PROGENY_heatmapD28.png",
  plot = progeny_hmapD28,
  width = 8,
  height = 12,
  dpi = 300
)

#enrichment analysis using competitive permutation approach to assess
#significance of pathway activity to end with NES for each pathway
PathwayActivity_zscoreD28 <- progeny(DEAD28_mat, 
                                    scale=TRUE, organism="Human", top = 100, perm = 10000, z_scores = TRUE) %>%
  t()
colnames(PathwayActivity_zscoreD28) <- "NES"

PathwayActivity_zscoreD28_df <- as.data.frame(PathwayActivity_zscoreD28) %>% 
  rownames_to_column(var = "Pathway") %>%
  dplyr::arrange(NES) %>%
  dplyr::mutate(Pathway = factor(Pathway))

NES_pathwayD28 <- ggplot(PathwayActivity_zscoreD28_df,aes(x = reorder(Pathway, NES), y = NES)) + 
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
  filename = "PROGENY_NES_pathwaysD28.png",
  plot = NES_pathwayD28,
  width = 10,
  height = 8,
  dpi = 300
)

# pathway most active - visualise most responsive genes (progeny weights)
#along with their statistic values to interpret results
#scatterplot shows which genes contributing most to the pathway enrichment

prog_matrixD28 <- getModel("Human", top=100) %>% 
  as.data.frame()  %>%
  tibble::rownames_to_column("GeneID")

DEAD28_df <- DEAD28_mat %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("GeneID")

sum(is.na(DEAD28_df$statistic))
DEAD28_df <- DEAD28_df %>%
  dplyr::filter(!is.na(statistic))

scat_plotsD28 <- progeny::progenyScatter(df = DEAD28_df, 
                                        weight_matrix = prog_matrixD28, 
                                        statName = "statistic", verbose = FALSE)
plot(scat_plotsD28[[1]]$`p53`) 

#progeny results as input for CARNIVAL
#CARNIVAL sets weights based on PROGENy scores in each pathway-related node
#to find more relevant solutions
#set progeny zscores false to return pathway activity values between 1 and -1

PathwayActivity_CARNIVALinputD28 <- progeny(DEAD28_mat, 
                                           scale=TRUE, organism="Human", top = 100, perm = 10000, z_scores = FALSE) %>%
  t () %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "Pathway") 
colnames(PathwayActivity_CARNIVALinputD28)[2] <- "score"
head(PathwayActivity_CARNIVALinputD28)
write_csv(PathwayActivity_CARNIVALinputD28, 
          "~/hnrnpu-causal-multiomics/processeddata/PathwayActivity_CARNIVALinputD28.csv")


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


head(TFscoresD28)
TFscoresD28_CARNIVAL <- TFscoresD28[, "score"]
class(TFscoresD28_CARNIVAL)
# numeric
names(TFscoresD28_CARNIVAL)
# "AHR" "AP1" "APEX1" ...


head(PathwayActivity_CARNIVALinputD28)
Pathway_CARNIVALD28 <- PathwayActivity_CARNIVALinputD28$score
names(Pathway_CARNIVALD28) <- PathwayActivity_CARNIVALinputD28$Pathway
head(Pathway_CARNIVALD28)
class(Pathway_CARNIVALD28)
# numeric

names(Pathway_CARNIVALD28)
# "Androgen" "EGFR" "Estrogen" "Hypoxia" "JAK.STAT" "MAPK"

##### CARNIVAL #####

sif_clean <- read.csv("clean_omnipath_PKN.csv")
clean_PKN

?runCARNIVAL
carnival_D28 <- runCARNIVAL(
  inputObj  = NULL,             
  measObj   = TFscoresD28_CARNIVAL,             # TF activity measurements
  netObj    = clean_PKN,                    # prior knowledge network
  weightObj = Pathway_CARNIVALD28,         # pathway activity scores (or progeny_matrix)
  solverPath = "C:/Program Files/IBM/ILOG/CPLEX_Studio2211/cplex/bin/x64_win64/cplex.exe", 
  solver     = "cplex",
  timelimit  = 7200, 
  mipGAP     = 0.05,
  poolrelGAP = 0
)
saveRDS(carnival_D28,"~/hnrnpu-causal-multiomics/processeddata/carnival_D28.rds")

head(carnival_D28)

#inputObj: data frame of list fir target of perturbation, OPTIONAL or default set to NULL to run invCARNIVAL when inputs are not known

## trial
# Load libraries
library(igraph)
library(ggraph)
library(tidyverse)
library(dplyr)

# 1. Prepare edges
edges <- carnival_D28$weightedSIF %>%
  filter(!is.na(Node1) & !is.na(Node2)) %>% # remove incomplete rows
  dplyr::select(Node1, Node2)

# 2. Prepare node attributes
nodes <- carnival_D28$nodesAttributes %>%
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
top_nodes <- carnival_D28$nodesAttributes %>%
  arrange(desc(AvgAct)) %>%
  dplyr::slice(1:50) %>%
  pull(Node)

sub_edges <- carnival_D28$weightedSIF %>%
  filter(Node1 %in% top_nodes | Node2 %in% top_nodes) %>%
  dplyr::select(Node1, Node2)

g_sub <- graph_from_data_frame(sub_edges, directed = TRUE)

nodes_sub <- carnival_D28$nodesAttributes %>%
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
library(stringr)

# Clean edges
edges <- carnival_D28$weightedSIF %>%
  dplyr::select(Node1, Node2, Weight) %>%
  mutate(
    Node1 = str_trim(Node1),
    Node2 = str_trim(Node2)
  )

write.csv(edges, "carnival_edgesD28.csv", row.names = FALSE)

# Clean nodes and rename AvgAct
nodes <- carnival_D28$nodesAttributes %>%
  dplyr::select(Node, AvgAct) %>%
  mutate(Node = str_trim(Node)) %>%  # remove extra spaces
  dplyr::rename(id = Node,
                activity_D28 = AvgAct)

write.csv(nodes, "carnival_nodesD28.csv", row.names = FALSE)



##### MOON #####
library(cosmosR)

#TF scores for MOON
TFscoresMOON <- TFscores_D28[, c('source', 'score')]
names(TFscoresMOON) <- c("TF", "score")
head(TFscoresMOON)
TFscoresMOON <- as.data.frame(TFscoresMOON)
row.names(TFscoresMOON) <- TFscoresMOON$TF
head(TFscoresMOON)

clean_PKN

#DEA for MOON
DEA_MOON <- resD28.annot[, c("hgnc_symbol", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
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

PKND28 <- cosmosR:::filter_pkn_expressed_genes(row.names(DEA_moon), meta_pkn = clean_PKN)

#####
DEA_moon 
obsgenes <- TFscoresMOON$TF
PKND28

head(names(TFscoresMOON))
head(obsgenes)

if(!"source" %in% names(PKND28)) { names(PKND28)[1] <- "source" }
if(!"target" %in% names(PKND28)) { names(PKND28)[2] <- "target" }

network_nodesD28 <- unique(c(PKND28$source, PKND28$target))
overlap_countD28 <- sum(obsgenes %in% network_nodesD28)
print(paste("Original genes to observe:", length(obsgenes)))
print(paste("Genes found in the network:", overlap_countD28))

obsgenes_D28 <- obsgenes[obsgenes %in% network_nodesD28]
# Prune the Network (Keep only nodes upstream of your TFs)
n_steps <- 6
PKNfiltD28 <- cosmosR:::keep_observable_neighbours(PKND28, n_steps, obsgenes_D28)

# Compress the network (Simplify redundant paths)
meta_network_compressed_listD28 <- compress_same_children(PKNfiltD28, sig_input = NULL, metab_input = obsgenes_D28)
meta_network_compressed_D28 <- meta_network_compressed_listD28$compressed_network

node_signaturesD28 <- meta_network_compressed_listD28$node_signaturesD28
duplicated_parentsD28 <- meta_network_compressed_listD28$duplicated_signatures
meta_network_compressed_D28 <- meta_network_cleanup(meta_network_compressed_D28)


# 7. The Main Loop (The "MOON" Algorithm)
# This iterates to find the best signaling path explaining your data
meta_network_TFD28 <- meta_network_compressed_D28
i <- 1
before <- 1
after <- 0

load("collectri_regulon_R.RData")

TFscoresMOON_vec <- TFscoresMOON$score
names(TFscoresMOON_vec) <- TFscoresMOON$TF

print("Starting COSMOS optimization...")
while (before != after & i < 10) {
  before <- length(meta_network_TFD28[,1])
  
  # Run MOON: Upstream=NULL because we don't know the drug target
  moon_resD28 <- moon(upstream_input = NULL, 
                     downstream_input = TFscoresMOON_vec, 
                     meta_network = meta_network_TFD28, 
                     n_layers = n_steps, 
                     statistic = "ulm") 
  
  # Prune incoherent links (contradictions between Kinase -> TF -> Gene)
  meta_network_TFD28 <- filter_incohrent_TF_target(moon_resD28, net, meta_network_TFD28, DEA_moon)
  
  after <- length(meta_network_TFD28[,1])
  i <- i + 1
  print(paste("Iteration:", i, "- Network size:", after))
}

getAnywhere("filter_incohrent_TF_target")

nrow(meta_network_TFD28)
head(meta_network_TFD28)


# Look at the top regulators (Kinases/Signaling Proteins)
top_regulators <- moon_resD28 %>%
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


moonGLI3 <- get_moon_scoring_network("RBL2", meta_network_TFD28, moon_resD28)
library(igraph)


# Create igraph object
g_moon <- graph_from_data_frame(
  d = moonE2F4$SIF,
  directed = TRUE,
  vertices = moonE2F4$ATT  # optional, adds attributes
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
V(g_moon)$score <- moonE2F4$ATT$score[match(V(g_moon)$name, moonE2F4$ATT$source)]

ggraph(g_moon, layout = "fr") +
  geom_edge_link(aes(edge_alpha = 0.5, edge_color = factor(E(g_moon)$interaction))) +
  geom_node_point(aes(color = score), size = 5) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  scale_color_gradient2(low = "steelblue", mid = "white", high = "firebrick", midpoint = 0) +
  scale_edge_color_manual(values = c("1" = "darkgreen", "-1" = "red")) +
  theme_void() +
  labs(title = "MOON Network for MAP2K2")

###ignore
nrow(carnival_D28$tfActivities)

