##### COSMOS_GDPx2 adaptation #####
##### collecTRI TF scores 10.12.25 #####
# We advise to instal from github to get the latest version of the tool.
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")

devtools::install_github("saezlab/cosmosR")
library(readr)
library(decoupleR)
library(reshape2)
library(progeny)
library(dorothea)
library(tibble)
library(tidyr)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(viper)
library(ggrepel)
source("support_functions.R")

#res.annot in "herewegoagain.R"
DEA <- res.annot[, c("hgnc_symbol", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
head(DEA)
#data frame, 13821 transcripts

net <- decoupleR::get_collectri(organism='human', split_complexes=FALSE)
#43159, source, target, mor

# Remove rows with blank hgnc_symbol
DEA <- DEA[DEA$hgnc_symbol != "", ]
# Remove duplicated hgnc_symbol rows (keeping the first one)
DEA <- DEA[!duplicated(DEA$hgnc_symbol), ]
row.names(DEA) <- DEA$hgnc_symbol
#12878 transcripts

names(DEA)[1] <- "gene"

DEA_vec <- DEA$stat
names(DEA_vec) <- DEA$gene
head(DEA_vec)

TFscores_res <- decoupleR::run_viper(
  mat = DEA_vec, #statistic matrix 
  net = net,                 
  .source = 'source', 
  .target = 'target', 
  .mor = 'mor', 
  minsize = 25
)
#217 TF scores
TFscores <- TFscores_res[, c('source', 'score')]
names(TFscores) <- c("TF", "score")
head(TFscores)

TFscores <- as.data.frame(TFscores)
row.names(TFscores) <- TFscores$TF
head(TFscores)
write.csv(TFscores, "TFscores.csv", row.names = TRUE) #217

## top 25 TF in general (abs)
# Filter for top 25 most active TFs (by absolute score)
top_tfs <- TFscores_res %>%
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
top_up <- TFscores_res %>%
  arrange(desc(score)) %>%
  slice_head(n = 15)

# Select Top 15 Inhibited (Lowest Negative Scores)
top_down <- TFscores_res %>%
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

##### omnipath PKN #####
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

##### run MOON #####
library(cosmosR)
TFscores_moon <- TFscores
clean_PKN
DEA_moon <- DEA
row.names(DEA_moon) <- DEA_moon$gene
head(DEA_moon)

PKN <- cosmosR:::filter_pkn_expressed_genes(row.names(DEA_moon), meta_pkn = clean_PKN)

# Where the function is defined (works even if it's inside a namespace)
getAnywhere("moon")

# Select ONLY the 'stat' column 
# We drop=FALSE to keep it as a dataframe, not a vector
DEA_clean <- DEA_moon[, "stat", drop = FALSE] 

# Rename the column to represent comparison
names(DEA_clean) <- "ASD_vs_CTRL"
DEA_moon <- DEA_clean
head(DEA_moon)
conditions_of_interest <- "ASD_vs_CTRL"

head(DEA_moon)

#####
DEA_moon 
obsgenes <- TFscores$TF
PKN

head(names(TFscores))
head(obsgenes)

if(!"source" %in% names(PKN)) { names(PKN)[1] <- "source" }
if(!"target" %in% names(PKN)) { names(PKN)[2] <- "target" }

network_nodes <- unique(c(PKN$source, PKN$target))
overlap_count <- sum(obsgenes %in% network_nodes)
print(paste("Original genes to observe:", length(obsgenes)))
print(paste("Genes found in the network:", overlap_count))

obsgenes_final <- obsgenes[obsgenes %in% network_nodes]
# Prune the Network (Keep only nodes upstream of your TFs)
n_steps <- 6
PKNfilt <- cosmosR:::keep_observable_neighbours(PKN, n_steps, obsgenes_final)

# Compress the network (Simplify redundant paths)
meta_network_compressed_list <- compress_same_children(PKNfilt, sig_input = NULL, metab_input = obsgenes_final)
meta_network_compressed <- meta_network_compressed_list$compressed_network

node_signatures <- meta_network_compressed_list$node_signatures
duplicated_parents <- meta_network_compressed_list$duplicated_signatures
meta_network_compressed <- meta_network_cleanup(meta_network_compressed)


# 7. The Main Loop (The "MOON" Algorithm)
# This iterates to find the best signaling path explaining your data
meta_network_TF <- meta_network_compressed
i <- 1
before <- 1
after <- 0

load("collectri_regulon_R.RData")

TFscores_vec <- TFscores$score
names(TFscores_vec) <- TFscores$TF

print("Starting COSMOS optimization...")
while (before != after & i < 10) {
  before <- length(meta_network_TF[,1])
  
  # Run MOON: Upstream=NULL because we don't know the drug target
  moon_res <- moon(upstream_input = NULL, 
                   downstream_input = TFscores_vec, 
                   meta_network = meta_network_TF, 
                   n_layers = n_steps, 
                   statistic = "ulm") 
  
  # Prune incoherent links (contradictions between Kinase -> TF -> Gene)
  meta_network_TF <- filter_incohrent_TF_target(moon_res, net, meta_network_TF, DEA_moon)
  
  after <- length(meta_network_TF[,1])
  i <- i + 1
  print(paste("Iteration:", i, "- Network size:", after))
}

getAnywhere("filter_incohrent_TF_target")


# Look at the top regulators (Kinases/Signaling Proteins)
top_regulators <- moon_res %>%
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

# Visualize top upstream regulators and paths
plot_moon_network(top_regulators, top_n = 20)
