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
PKND0

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
topregD0 <- ggplot(head(top_regulators, 20), aes(x = reorder(source, score), y = score, fill = score > 0)) +
  geom_col() +
  coord_flip() +
  theme_minimal() +
  labs(title = "Top Upstream Regulators (Disease vs Control)",
       x = "Signaling Protein",
       y = "Predicted Activity Score") +
  scale_fill_manual(values = c("TRUE" = "firebrick", "FALSE" = "steelblue"), 
                    labels = c("Activated", "Inhibited"))

ggsave(
  filename = "top_upstream_regulatorsD0.png",
  plot = topregD0,
  width = 8,
  height = 6,
  dpi = 300
)

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
V(g_moon)$score <- moonMAP2K2$ATT$score[match(V(g_moon)$name, moonMAP2K2$ATT$source)]

ggraph(g_moon, layout = "fr") +
  geom_edge_link(aes(edge_alpha = 0.5, edge_color = factor(E(g_moon)$interaction))) +
  geom_node_point(aes(color = score), size = 5) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  scale_color_gradient2(low = "steelblue", mid = "white", high = "firebrick", midpoint = 0) +
  scale_edge_color_manual(values = c("1" = "darkgreen", "-1" = "red")) +
  theme_void() +
  labs(title = "MOON Network for MAP2K2")


##### D5 #####
library(cosmosR)

#TF scores for MOON
TFscoresMOOND5 <- TFscores_D5[, c('source', 'score')]
names(TFscoresMOOND5) <- c("TF", "score")
head(TFscoresMOOND5)
TFscoresMOOND5 <- as.data.frame(TFscoresMOOND5)
row.names(TFscoresMOOND5) <- TFscoresMOOND5$TF
head(TFscoresMOOND5)

clean_PKN

#DEA for MOON
DEA_MOOND5 <- res_D5_only.annot[, c("hgnc_symbol", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
head(DEA_MOOND5)


# Remove rows with blank hgnc_symbol
DEA_MOOND5 <- DEA_MOOND5[DEA_MOOND5$hgnc_symbol != "", ]
# Remove duplicated hgnc_symbol rows (keeping the first one)
DEA_MOOND5 <- DEA_MOOND5[!duplicated(DEA_MOOND5$hgnc_symbol), ]
row.names(DEA_MOOND5) <- DEA_MOOND5$hgnc_symbol
#12878 transcripts

names(DEA_MOOND5)[1] <- "gene"
row.names(DEA_MOOND5) <- DEA_MOOND5$gene
head(DEA_MOOND5)
DEA_clean_MOOND5 <- DEA_MOOND5[, "stat", drop = FALSE] 

# Rename the column to represent comparison
names(DEA_clean_MOOND5) <- "ASD_vs_CTRL"
DEA_moonD5 <- DEA_clean_MOOND5
head(DEA_moonD5)
conditions_of_interest <- "ASD_vs_CTRL"

head(DEA_moonD5)

PKND5 <- cosmosR:::filter_pkn_expressed_genes(row.names(DEA_moonD5), meta_pkn = clean_PKN)

#####
DEA_moonD5
obsgenesD5 <- TFscoresMOOND5$TF
PKND5

head(names(TFscoresMOOND5))
head(obsgenesD5)

if(!"source" %in% names(PKND5)) { names(PKND5)[1] <- "source" }
if(!"target" %in% names(PKND5)) { names(PKND5)[2] <- "target" }

network_nodesD5 <- unique(c(PKND5$source, PKND5$target))
overlap_countD5 <- sum(obsgenesD5 %in% network_nodesD5)
print(paste("Original genes to observe:", length(obsgenesD5)))
print(paste("Genes found in the network:", overlap_countD5))

obsgenes_D5 <- obsgenesD5[obsgenesD5 %in% network_nodesD5]
# Prune the Network (Keep only nodes upstream of your TFs)
n_steps <- 6
PKNfiltD5 <- cosmosR:::keep_observable_neighbours(PKND5, n_steps, obsgenes_D5)

# Compress the network (Simplify redundant paths)
meta_network_compressed_listD5 <- compress_same_children(PKNfiltD5, sig_input = NULL, metab_input = obsgenes_D5)
meta_network_compressed_D5 <- meta_network_compressed_listD5$compressed_network

node_signaturesD5 <- meta_network_compressed_listD5$node_signaturesD5
duplicated_parentsD5 <- meta_network_compressed_listD5$duplicated_signatures
meta_network_compressed_D5 <- meta_network_cleanup(meta_network_compressed_D5)


# 7. The Main Loop (The "MOON" Algorithm)
# This iterates to find the best signaling path explaining your data
meta_network_TFD5 <- meta_network_compressed_D5
i <- 1
before <- 1
after <- 0

load("collectri_regulon_R.RData")

TFscoresMOON_vecD5 <- TFscoresMOOND5$score
names(TFscoresMOON_vecD5) <- TFscoresMOOND5$TF

print("Starting COSMOS optimization...")
while (before != after & i < 10) {
  before <- length(meta_network_TFD5[,1])
  
  # Run MOON: Upstream=NULL because we don't know the drug target
  moon_resD5 <- moon(upstream_input = NULL, 
                     downstream_input = TFscoresMOON_vecD5, 
                     meta_network = meta_network_TFD5, 
                     n_layers = n_steps, 
                     statistic = "ulm") 
  
  # Prune incoherent links (contradictions between Kinase -> TF -> Gene)
  meta_network_TFD5 <- filter_incohrent_TF_target(moon_resD5, net, meta_network_TFD5, DEA_moonD5)
  
  after <- length(meta_network_TFD5[,1])
  i <- i + 1
  print(paste("Iteration:", i, "- Network size:", after))
}

getAnywhere("filter_incohrent_TF_target")

nrow(meta_network_TFD5)
head(meta_network_TFD5)


# Look at the top regulators (Kinases/Signaling Proteins)
top_regulatorsD5 <- moon_resD5 %>%
  dplyr::arrange(desc(abs(score))) # Sort by magnitude of score

# Display the top 10 potential drivers of your disease
print(head(top_regulatorsD5, 10))

# Visualize it nicely
topregD5 <- ggplot(head(top_regulatorsD5, 20), aes(x = reorder(source, score), y = score, fill = score > 0)) +
  geom_col() +
  coord_flip() +
  theme_minimal() +
  labs(title = "Top Upstream Regulators (Disease vs Control)",
       x = "Signaling Protein",
       y = "Predicted Activity Score") +
  scale_fill_manual(values = c("TRUE" = "firebrick", "FALSE" = "steelblue"), 
                    labels = c("Activated", "Inhibited"))

ggsave(
  filename = "top_upstream_regulatorsD5.png",
  plot = topregD5,
  width = 8,
  height = 6,
  dpi = 300
)
moonMAP2K2 <- get_moon_scoring_network("MAP2K2", meta_network_TFD5, moon_resD5)
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

##### D28 #####
library(cosmosR)

#TF scores for MOON
TFscoresMOOND28 <- TFscores_D28[, c('source', 'score')]
names(TFscoresMOOND28) <- c("TF", "score")
head(TFscoresMOOND28)
TFscoresMOOND28 <- as.data.frame(TFscoresMOOND28)
row.names(TFscoresMOOND28) <- TFscoresMOOND28$TF
head(TFscoresMOOND28)

clean_PKN

#DEA for MOON
DEA_MOOND28 <- res_D28_only.annot[, c("hgnc_symbol", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
head(DEA_MOOND28)


# Remove rows with blank hgnc_symbol
DEA_MOOND28 <- DEA_MOOND28[DEA_MOOND28$hgnc_symbol != "", ]
# Remove duplicated hgnc_symbol rows (keeping the first one)
DEA_MOOND28 <- DEA_MOOND28[!duplicated(DEA_MOOND28$hgnc_symbol), ]
row.names(DEA_MOOND28) <- DEA_MOOND28$hgnc_symbol
#12878 transcripts

names(DEA_MOOND28)[1] <- "gene"
row.names(DEA_MOOND28) <- DEA_MOOND28$gene
head(DEA_MOOND28)
DEA_clean_MOOND28 <- DEA_MOOND28[, "stat", drop = FALSE] 

# Rename the column to represent comparison
names(DEA_clean_MOOND28) <- "ASD_vs_CTRL"
DEA_moonD28 <- DEA_clean_MOOND28
head(DEA_moonD28)
conditions_of_interest <- "ASD_vs_CTRL"

head(DEA_moonD28)

PKND28 <- cosmosR:::filter_pkn_expressed_genes(row.names(DEA_moonD28), meta_pkn = clean_PKN)

#####
DEA_moonD28
obsgenesD28 <- TFscoresMOOND28$TF
PKND28

head(names(TFscoresMOOND28))
head(obsgenesD28)

if(!"source" %in% names(PKND28)) { names(PKND28)[1] <- "source" }
if(!"target" %in% names(PKND28)) { names(PKND28)[2] <- "target" }

network_nodesD28 <- unique(c(PKND28$source, PKND28$target))
overlap_countD28 <- sum(obsgenesD28 %in% network_nodesD28)
print(paste("Original genes to observe:", length(obsgenesD28)))
print(paste("Genes found in the network:", overlap_countD28))

obsgenes_D28 <- obsgenesD28[obsgenesD28 %in% network_nodesD28]
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

TFscoresMOON_vecD28 <- TFscoresMOOND28$score
names(TFscoresMOON_vecD28) <- TFscoresMOOND28$TF

print("Starting COSMOS optimization...")
while (before != after & i < 10) {
  before <- length(meta_network_TFD28[,1])
  
  # Run MOON: Upstream=NULL because we don't know the drug target
  moon_resD28 <- moon(upstream_input = NULL, 
                     downstream_input = TFscoresMOON_vecD28, 
                     meta_network = meta_network_TFD28, 
                     n_layers = n_steps, 
                     statistic = "ulm") 
  
  # Prune incoherent links (contradictions between Kinase -> TF -> Gene)
  meta_network_TFD28 <- filter_incohrent_TF_target(moon_resD28, net, meta_network_TFD28, DEA_moonD28)
  
  after <- length(meta_network_TFD28[,1])
  i <- i + 1
  print(paste("Iteration:", i, "- Network size:", after))
}

getAnywhere("filter_incohrent_TF_target")

nrow(meta_network_TFD28)
head(meta_network_TFD28)


# Look at the top regulators (Kinases/Signaling Proteins)
top_regulatorsD28 <- moon_resD28 %>%
  dplyr::arrange(desc(abs(score))) # Sort by magnitude of score

# Display the top 10 potential drivers of your disease
print(head(top_regulatorsD28, 10))

# Visualize it nicely
topregD28 <- ggplot(head(top_regulatorsD28, 20), aes(x = reorder(source, score), y = score, fill = score > 0)) +
  geom_col() +
  coord_flip() +
  theme_minimal() +
  labs(title = "Top Upstream Regulators (Disease vs Control)",
       x = "Signaling Protein",
       y = "Predicted Activity Score") +
  scale_fill_manual(values = c("TRUE" = "firebrick", "FALSE" = "steelblue"), 
                    labels = c("Activated", "Inhibited"))
ggsave(
  filename = "top_upstream_regulatorsD28.png",
  plot = topregD28,
  width = 8,
  height = 6,
  dpi = 300
)

moonMAP2K2 <- get_moon_scoring_network("MAP2K2", meta_network_TFD28, moon_resD28)
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