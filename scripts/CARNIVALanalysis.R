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

##### CARNIVAL #####
##### D0 #####

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

unique(c(sif$source, sif$target))
#save SIF
write_tsv(sif, "~/hnrnpu-causal-multiomics/processeddata/omnipath_carnival.tsv")
#70565

#top 100
top100tfs_D0 <- TFscores_D0 %>%
  arrange(desc(abs(score))) %>%
  slice_head(n = 100)
head(top100tfs_D0)

top100tfs_D0 <- top100tfs_D0 %>%
  dplyr::select(source, score) %>%
  deframe()

head(PathwayActivity_CARNIVALinputD0only)
Pathway_CARNIVALD0 <- PathwayActivity_CARNIVALinputD0only$score
names(Pathway_CARNIVALD0) <- PathwayActivity_CARNIVALinputD0only$Pathway
head(Pathway_CARNIVALD0)
class(Pathway_CARNIVALD0)
# numeric
names(Pathway_CARNIVALD0)
# "Androgen" "EGFR" "Estrogen" "Hypoxia" "JAK.STAT" "MAPK"

iniMTX = base::setdiff(sif_clean$source, sif_clean$target)
iniciators = base::data.frame(base::matrix(data = NaN, nrow = 1, ncol = length(iniMTX)), stringsAsFactors = F)
colnames(iniciators) = iniMTX

carnival_D0top100 <- runCARNIVAL(
  inputObj  = iniciators,             
  measObj   = top100tfs_D0,             # TF activity measurements
  netObj    = sif,                    # prior knowledge network
  weightObj = Pathway_CARNIVALD0,         # pathway activity scores 
  solverPath = "C:/Program Files/IBM/ILOG/CPLEX_Studio2211/cplex/bin/x64_win64/cplex.exe", #cplex link
  solver     = "cplex",
  timelimit  = 900, 
  mipGAP     = 0.05,
  poolrelGAP = 0
)
saveRDS(carnival_D0top100,"~/hnrnpu-causal-multiomics/processeddata/carnival_D0top100.rds")
carnival_D0top100 <- readRDS("~/hnrnpu-causal-multiomics/processeddata/carnival_D0top100.rds")

clean_edgesD0_100 <- carnival_D0top100$weightedSIF %>%
  filter(Weight > 0)

active_node_namesD0_100 <- unique(c(clean_edgesD0_100$Node1, clean_edgesD0_100$Node2))

clean_nodesD0_100 <- carnival_D0top100$nodesAttributes %>%
  filter(Node %in% active_node_namesD0_100)

write.csv(clean_edgesD0_100, "Clean_CARNIVAL_EdgesD0_100.csv", row.names = FALSE)
write.csv(clean_nodesD0_100, "Clean_CARNIVAL_NodesD0_100.csv", row.names = FALSE)

cytoscapePing() 
edges_for_cyD0_100 <- clean_edgesD0_100
colnames(edges_for_cyD0_100) <- c("source", "sign", "target", "weight") 

nodes_for_cyD0_100 <- clean_nodesD0_100
colnames(nodes_for_cyD0_100)[1] <- "id" # Ensure first column is 'id'

# Send to Cytoscape
createNetworkFromDataFrames(nodes = nodes_for_cyD0_100, 
                            edges = edges_for_cyD0_100, 
                            title = "CARNIVAL_D0_top100", 
                            collection = "CARNIVAL_Analysis_D0_top100")
#142 nodes 154 edges

##### D5 #####
head(PathwayActivity_CARNIVALinputD5only)
Pathway_CARNIVALD5 <- PathwayActivity_CARNIVALinputD5only$score
names(Pathway_CARNIVALD5) <- PathwayActivity_CARNIVALinputD5only$Pathway
head(Pathway_CARNIVALD5)
class(Pathway_CARNIVALD5)
# numeric
names(Pathway_CARNIVALD5)

#top 100
top100tfs_D5 <- TFscores_D5 %>%
  arrange(desc(abs(score))) %>%
  slice_head(n = 100)
head(top100tfs_D5)

top100tfs_D5 <- top100tfs_D5 %>%
  dplyr::select(source, score) %>%
  deframe()

carnival_D5top100 <- runCARNIVAL(
  inputObj  = iniciators,             
  measObj   = top100tfs_D5,             # TF activity measurements
  netObj    = sif,                    # prior knowledge network
  weightObj = Pathway_CARNIVALD5,         # pathway activity scores 
  solverPath = "C:/Program Files/IBM/ILOG/CPLEX_Studio2211/cplex/bin/x64_win64/cplex.exe", 
  solver     = "cplex",
  timelimit  = 900, 
  mipGAP     = 0.05,
  poolrelGAP = 0
)
saveRDS(carnival_D5top100,"~/hnrnpu-causal-multiomics/processeddata/carnival_D5top100.rds")
carnival_D5top100 <- readRDS("~/hnrnpu-causal-multiomics/processeddata/carnival_D5top100.rds")

clean_edgesD5_100 <- carnival_D5top100$weightedSIF %>%
  filter(Weight > 0)

active_node_namesD5_100 <- unique(c(clean_edgesD5_100$Node1, clean_edgesD5_100$Node2))

clean_nodesD5_100 <- carnival_D5top100$nodesAttributes %>%
  filter(Node %in% active_node_namesD5_100)

write.csv(clean_edgesD5_100, "Clean_CARNIVAL_EdgesD5_100.csv", row.names = FALSE)
write.csv(clean_nodesD5_100, "Clean_CARNIVAL_NodesD5_100.csv", row.names = FALSE)

cytoscapePing() 

edges_for_cyD5_100 <- clean_edgesD5_100
colnames(edges_for_cyD5_100) <- c("source", "sign", "target", "weight") # Adjust based on your actual column order!

nodes_for_cyD5_100 <- clean_nodesD5_100
colnames(nodes_for_cyD5_100)[1] <- "id" # Ensure first column is 'id'

# Send to Cytoscape
createNetworkFromDataFrames(nodes = nodes_for_cyD5_100, 
                            edges = edges_for_cyD5_100, 
                            title = "CARNIVAL_D5_top100", 
                            collection = "CARNIVAL_Analysis_D5_top100")
#137 nodes 184 edges

##### D28 #####

head(PathwayActivity_CARNIVALinputD28only)
Pathway_CARNIVALD28 <- PathwayActivity_CARNIVALinputD28only$score
names(Pathway_CARNIVALD28) <- PathwayActivity_CARNIVALinputD28only$Pathway
head(Pathway_CARNIVALD28)
class(Pathway_CARNIVALD28)
# numeric
names(Pathway_CARNIVALD28)

#top 100
top100tfs_D28 <- TFscores_D28 %>%
  arrange(desc(abs(score))) %>%
  slice_head(n = 100)
head(top100tfs_D28)

top100tfs_D28 <- top100tfs_D28 %>%
  dplyr::select(source, score) %>%
  deframe()

carnival_D28top100 <- runCARNIVAL(
  inputObj  = iniciators,             
  measObj   = top100tfs_D28,             # TF activity measurements
  netObj    = sif,                    # prior knowledge network
  weightObj = Pathway_CARNIVALD28,         # pathway activity scores 
  solverPath = "C:/Program Files/IBM/ILOG/CPLEX_Studio2211/cplex/bin/x64_win64/cplex.exe", 
  solver     = "cplex",
  timelimit  = 900, 
  mipGAP     = 0.05,
  poolrelGAP = 0
)
saveRDS(carnival_D28top100,"~/hnrnpu-causal-multiomics/processeddata/carnival_D28top100.rds")
carnival_D28top100 <- readRDS("~/hnrnpu-causal-multiomics/processeddata/carnival_D28top100.rds")

clean_edgesD28_100 <- carnival_D28top100$weightedSIF %>%
  filter(Weight > 0)

active_node_namesD28_100 <- unique(c(clean_edgesD28_100$Node1, clean_edgesD28_100$Node2))

clean_nodesD28_100 <- carnival_D28top100$nodesAttributes %>%
  filter(Node %in% active_node_namesD28_100)

write.csv(clean_edgesD28_100, "Clean_CARNIVAL_EdgesD28_100.csv", row.names = FALSE)
write.csv(clean_nodesD28_100, "Clean_CARNIVAL_NodesD28_100.csv", row.names = FALSE)

cytoscapePing() 
edges_for_cyD28_100 <- clean_edgesD28_100
colnames(edges_for_cyD28_100) <- c("source", "sign", "target", "weight") # Adjust based on your actual column order!

nodes_for_cyD28_100 <- clean_nodesD28_100
colnames(nodes_for_cyD28_100)[1] <- "id" # Ensure first column is 'id'

# Send to Cytoscape
createNetworkFromDataFrames(nodes = nodes_for_cyD28_100, 
                            edges = edges_for_cyD28_100, 
                            title = "CARNIVAL_D28_top100", 
                            collection = "CARNIVAL_Analysis_D28_top100")
#211 nodes 420 edges

##### node venn diagram #####
list_d0  <- unique(clean_nodesD0_100$Node)
list_d5  <- unique(clean_nodesD5_100$Node)
list_d28 <- unique(clean_nodesD28_100$Node)

# Create a named list for the plotting function
venn_data <- list(
  "Day 0" = list_d0,
  "Day 5" = list_d5,
  "Day 28" = list_d28
)

#Present in ALL 3 days
common_all <- Reduce(intersect, venn_data)

#Unique to Each Day 
unique_d0  <- setdiff(list_d0,  union(list_d5, list_d28))
unique_d5  <- setdiff(list_d5,  union(list_d0, list_d28))
unique_d28 <- setdiff(list_d28, union(list_d0, list_d5))

#present in D5 and D28 but NOT D0 (
shared_d5_d28 <- setdiff(intersect(list_d5, list_d28), list_d0)

cat("\n=== NETWORK OVERLAP SUMMARY ===\n")
cat("Total Nodes in Day 0: ", length(list_d0), "\n")
cat("Total Nodes in Day 5: ", length(list_d5), "\n")
cat("Total Nodes in Day 28:", length(list_d28), "\n\n")

cat("--- CORE NODES (Preserved across all stages) ---\n")
print(common_all)

cat("\n--- UNIQUE NODES (Specific to Day 0) ---\n")
print(unique_d0)

cat("\n--- UNIQUE NODES (Specific to Day 5) ---\n")
print(unique_d5)

cat("\n--- UNIQUE NODES (Specific to Day 28) ---\n")
print(unique_d28)

ggvenn(
  venn_data, 
  fill_color = c("#E41A1C", "#377EB8", "#4DAF4A"), # Red, Blue, Green
  stroke_size = 0.5, 
  set_name_size = 5
) +
  ggtitle("Overlap of Signaling Nodes Across Differentiation") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
ggsave("Node_Overlap_Venn.png", width = 8, height = 8, dpi = 300)

max_len <- max(length(common_all), length(unique_d0), length(unique_d5), length(unique_d28))
pad_na <- function(x, n) { c(x, rep(NA, n - length(x))) }
export_df <- data.frame(
  Common_Nodes = pad_na(common_all, max_len),
  Unique_D0    = pad_na(unique_d0, max_len),
  Unique_D5    = pad_na(unique_d5, max_len),
  Unique_D28   = pad_na(unique_d28, max_len),
  Shared_d5_d28 = pad_na(shared_d5_d28, max_len)
)
write.csv(export_df, "Node_Comparison_Lists.csv", row.names = FALSE)