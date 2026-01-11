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

head(PathwayActivity_CARNIVALinputD0only)
Pathway_CARNIVALD0 <- PathwayActivity_CARNIVALinputD0only$score
names(Pathway_CARNIVALD0) <- PathwayActivity_CARNIVALinputD0only$Pathway
head(Pathway_CARNIVALD0)
class(Pathway_CARNIVALD0)
# numeric
names(Pathway_CARNIVALD0)
# "Androgen" "EGFR" "Estrogen" "Hypoxia" "JAK.STAT" "MAPK"


##### CARNIVAL #####
sif_clean <- read.csv("clean_omnipath_PKN.csv")

iniMTX = base::setdiff(sif_clean$source, sif_clean$target)
iniciators = base::data.frame(base::matrix(data = NaN, nrow = 1, ncol = length(iniMTX)), stringsAsFactors = F)
colnames(iniciators) = iniMTX

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
carnival_D0 <- readRDS("~/hnrnpu-causal-multiomics/processeddata/carnival_D0.rds")
head(carnival_D0)

#inputObj: data frame of list fir target of perturbation, OPTIONAL or default set to NULL to run invCARNIVAL when inputs are not known

##### approach 2 #####
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
carnival_D0 <- runCARNIVAL(
  inputObj  = iniciators,             
  measObj   = TFscoresD0_CARNIVAL,             # TF activity measurements
  netObj    = sif,                    # prior knowledge network
  weightObj = Pathway_CARNIVALD0,         # pathway activity scores (or progeny_matrix)
  solverPath = "C:/Program Files/IBM/ILOG/CPLEX_Studio2211/cplex/bin/x64_win64/cplex.exe", 
  solver     = "cplex",
  timelimit  = 7200, 
  mipGAP     = 0.05,
  poolrelGAP = 0
)
saveRDS(carnival_D0,"~/hnrnpu-causal-multiomics/processeddata/carnival_D0.rds")

clean_edgesD0 <- carnival_D0$weightedSIF %>%
  filter(Weight > 0)

active_node_namesD0 <- unique(c(clean_edgesD0$Node1, clean_edgesD0$Node2))

clean_nodesD0 <- carnival_D0$nodesAttributes %>%
  filter(Node %in% active_node_namesD0)

###### visual ######
library(RCy3)
library(igraph)
library(dplyr)
write.csv(clean_edgesD0, "Clean_CARNIVAL_EdgesD0.csv", row.names = FALSE)
write.csv(clean_nodesD0, "Clean_CARNIVAL_NodesD0.csv", row.names = FALSE)

cytoscapePing() 

# Rename columns to standard RCy3 formats to be safe
# (RCy3 prefers 'source', 'target', 'interaction' but handles others too)
edges_for_cyD0 <- clean_edgesD0
colnames(edges_for_cyD0) <- c("source", "sign", "target", "weight") # Adjust based on your actual column order!

nodes_for_cyD0 <- clean_nodesD0
colnames(nodes_for_cyD0)[1] <- "id" # Ensure first column is 'id'

# Send to Cytoscape
createNetworkFromDataFrames(nodes = nodes_for_cyD0, 
                            edges = edges_for_cyD0, 
                            title = "CARNIVAL_Automated", 
                            collection = "CARNIVAL_Analysis")

setNodeColorMapping(table.column = "AvgAct", 
                    table.column.values = c(-100, 0, 100), 
                    colors = c("#0000FF", "#FFFFFF", "#FF0000"), 
                    mapping.type = "continuous")

## top 30
top30tfs_D0 <- TFscores_D0 %>%
  arrange(desc(abs(score))) %>%
  slice_head(n = 30)
head(top30tfs_D0)

top30tfs_D0 <- top30tfs_D0 %>%
  dplyr::select(source, score) %>%
  deframe()

carnival_D0top30 <- runCARNIVAL(
  inputObj  = iniciators,             
  measObj   = top30tfs_D0,             # TF activity measurements
  netObj    = sif,                    # prior knowledge network
  weightObj = Pathway_CARNIVALD0,         # pathway activity scores (or progeny_matrix)
  solverPath = "C:/Program Files/IBM/ILOG/CPLEX_Studio2211/cplex/bin/x64_win64/cplex.exe", 
  solver     = "cplex",
  timelimit  = 900, 
  mipGAP     = 0.05,
  poolrelGAP = 0
)
saveRDS(carnival_D0top30,"~/hnrnpu-causal-multiomics/processeddata/carnival_D0top30.rds")
carnival_D0top30 <- readRDS("~/hnrnpu-causal-multiomics/processeddata/carnival_D0top30.rds")

clean_edgesD0_30 <- carnival_D0top30$weightedSIF %>%
  filter(Weight > 0)

active_node_namesD0_30 <- unique(c(clean_edgesD0_30$Node1, clean_edgesD0_30$Node2))

clean_nodesD0_30 <- carnival_D0top30$nodesAttributes %>%
  filter(Node %in% active_node_namesD0_30)

write.csv(clean_edgesD0_30, "Clean_CARNIVAL_EdgesD28_30.csv", row.names = FALSE)
write.csv(clean_nodesD0_30, "Clean_CARNIVAL_NodesD28_30.csv", row.names = FALSE)

cytoscapePing() 

# Rename columns to standard RCy3 formats to be safe
# (RCy3 prefers 'source', 'target', 'interaction' but handles others too)
edges_for_cyD0_30 <- clean_edgesD0_30
colnames(edges_for_cyD0_30) <- c("source", "sign", "target", "weight") # Adjust based on your actual column order!

nodes_for_cyD0_30 <- clean_nodesD0_30
colnames(nodes_for_cyD0_30)[1] <- "id" # Ensure first column is 'id'

# Send to Cytoscape
createNetworkFromDataFrames(nodes = nodes_for_cyD0_30, 
                            edges = edges_for_cyD0_30, 
                            title = "CARNIVAL_D0_top30", 
                            collection = "CARNIVAL_Analysis_D0_top30")

#top 100
top100tfs_D0 <- TFscores_D0 %>%
  arrange(desc(abs(score))) %>%
  slice_head(n = 100)
head(top100tfs_D0)

top100tfs_D0 <- top100tfs_D0 %>%
  dplyr::select(source, score) %>%
  deframe()

carnival_D0top100 <- runCARNIVAL(
  inputObj  = iniciators,             
  measObj   = top100tfs_D0,             # TF activity measurements
  netObj    = sif,                    # prior knowledge network
  weightObj = Pathway_CARNIVALD0,         # pathway activity scores (or progeny_matrix)
  solverPath = "C:/Program Files/IBM/ILOG/CPLEX_Studio2211/cplex/bin/x64_win64/cplex.exe", 
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

# Rename columns to standard RCy3 formats to be safe
# (RCy3 prefers 'source', 'target', 'interaction' but handles others too)
edges_for_cyD0_100 <- clean_edgesD0_100
colnames(edges_for_cyD0_100) <- c("source", "sign", "target", "weight") # Adjust based on your actual column order!

nodes_for_cyD0_100 <- clean_nodesD0_100
colnames(nodes_for_cyD0_100)[1] <- "id" # Ensure first column is 'id'

# Send to Cytoscape
createNetworkFromDataFrames(nodes = nodes_for_cyD0_100, 
                            edges = edges_for_cyD0_100, 
                            title = "CARNIVAL_D0_top100", 
                            collection = "CARNIVAL_Analysis_D0_top100")
#142 nodes 154 edges
##### D5 #####
head(TFscoresD5)
TFscoresD5_CARNIVAL <- TFscoresD5[, "score"]
class(TFscoresD5_CARNIVAL)
# numeric
names(TFscoresD5_CARNIVAL)
# "AHR" "AP1" "APEX1" ...

head(PathwayActivity_CARNIVALinputD5only)
Pathway_CARNIVALD5 <- PathwayActivity_CARNIVALinputD5only$score
names(Pathway_CARNIVALD5) <- PathwayActivity_CARNIVALinputD5only$Pathway
head(Pathway_CARNIVALD5)
class(Pathway_CARNIVALD5)
# numeric
names(Pathway_CARNIVALD5)

carnival_D5 <- runCARNIVAL(
  inputObj  = iniciators,             
  measObj   = TFscoresD5_CARNIVAL,             # TF activity measurements
  netObj    = sif,                    # prior knowledge network
  weightObj = Pathway_CARNIVALD5,         # pathway activity scores (or progeny_matrix)
  solverPath = "C:/Program Files/IBM/ILOG/CPLEX_Studio2211/cplex/bin/x64_win64/cplex.exe", 
  solver     = "cplex",
  timelimit  = 7200, 
  mipGAP     = 0.05,
  poolrelGAP = 0
)
saveRDS(carnival_D5, "~/hnrnpu-causal-multiomics/processeddata/carnival_D5.rds")
carnival_D5 <- readRDS("~/hnrnpu-causal-multiomics/processeddata/carnival_D5.rds")
clean_edgesD5 <- carnival_D5$weightedSIF %>%
  filter(Weight > 0)

active_node_namesD5 <- unique(c(clean_edgesD5$Node1, clean_edgesD5$Node2))

clean_nodesD5 <- carnival_D5$nodesAttributes %>%
  filter(Node %in% active_node_namesD5)

###### visual ######
library(RCy3)
library(igraph)
library(dplyr)
write.csv(clean_edgesD5, "Clean_CARNIVAL_EdgesD5.csv", row.names = FALSE)
write.csv(clean_nodesD5, "Clean_CARNIVAL_NodesD5.csv", row.names = FALSE)

cytoscapePing() 

# Rename columns to standard RCy3 formats to be safe
# (RCy3 prefers 'source', 'target', 'interaction' but handles others too)
edges_for_cyD5 <- clean_edgesD5
colnames(edges_for_cyD5) <- c("source", "sign", "target", "weight") # Adjust based on your actual column order!

nodes_for_cyD5 <- clean_nodesD5
colnames(nodes_for_cyD5)[1] <- "id" # Ensure first column is 'id'

# Send to Cytoscape
createNetworkFromDataFrames(nodes = nodes_for_cyD5, 
                            edges = edges_for_cyD5, 
                            title = "CARNIVAL_D5", 
                            collection = "CARNIVAL_AnalysisD5")

setNodeColorMapping(table.column = "AvgAct", 
                    table.column.values = c(-100, 0, 100), 
                    colors = c("#0000FF", "#FFFFFF", "#FF0000"), 
                    mapping.type = "continuous")

## top 30
top30tfs_D5 <- TFscores_D5 %>%
  arrange(desc(abs(score))) %>%
  slice_head(n = 30)
head(top30tfs_D5)

top30tfs_D5 <- top30tfs_D5 %>%
  dplyr::select(source, score) %>%
  deframe()

carnival_D5top30 <- runCARNIVAL(
  inputObj  = iniciators,             
  measObj   = top30tfs_D5,             # TF activity measurements
  netObj    = sif,                    # prior knowledge network
  weightObj = Pathway_CARNIVALD5,         # pathway activity scores (or progeny_matrix)
  solverPath = "C:/Program Files/IBM/ILOG/CPLEX_Studio2211/cplex/bin/x64_win64/cplex.exe", 
  solver     = "cplex",
  timelimit  = 900, 
  mipGAP     = 0.05,
  poolrelGAP = 0
)
saveRDS(carnival_D5top30,"~/hnrnpu-causal-multiomics/processeddata/carnival_D5top30.rds")
carnival_D5top30 <- readRDS("~/hnrnpu-causal-multiomics/processeddata/carnival_D5top30.rds")

clean_edgesD5_30 <- carnival_D5top30$weightedSIF %>%
  filter(Weight > 0)

active_node_namesD5_30 <- unique(c(clean_edgesD5_30$Node1, clean_edgesD5_30$Node2))

clean_nodesD5_30 <- carnival_D5top30$nodesAttributes %>%
  filter(Node %in% active_node_namesD5_30)

write.csv(clean_edgesD5_30, "Clean_CARNIVAL_EdgesD28_30.csv", row.names = FALSE)
write.csv(clean_nodesD5_30, "Clean_CARNIVAL_NodesD28_30.csv", row.names = FALSE)

library(visNetwork)
source("~/hnrnpu-causal-multiomics/carnival_visNet.R")
carnival_visNet(evis = carnival_D5top30$weightedSIF,
                nvis = carnival_D5top30$nodesAttributes)

cytoscapePing() 

# Rename columns to standard RCy3 formats to be safe
# (RCy3 prefers 'source', 'target', 'interaction' but handles others too)
edges_for_cyD5_30 <- clean_edgesD5_30
colnames(edges_for_cyD5_30) <- c("source", "sign", "target", "weight") # Adjust based on your actual column order!

nodes_for_cyD5_30 <- clean_nodesD5_30
colnames(nodes_for_cyD5_30)[1] <- "id" # Ensure first column is 'id'

# Send to Cytoscape
createNetworkFromDataFrames(nodes = nodes_for_cyD5_30, 
                            edges = edges_for_cyD5_30, 
                            title = "CARNIVAL_D5_top30", 
                            collection = "CARNIVAL_Analysis_D5_top30")

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
  weightObj = Pathway_CARNIVALD5,         # pathway activity scores (or progeny_matrix)
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

# Rename columns to standard RCy3 formats to be safe
# (RCy3 prefers 'source', 'target', 'interaction' but handles others too)
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
head(TFscoresD28)
TFscoresD28_CARNIVAL <- TFscoresD28[, "score"]
class(TFscoresD28_CARNIVAL)
# numeric
names(TFscoresD28_CARNIVAL)
# "AHR" "AP1" "APEX1" ...

head(PathwayActivity_CARNIVALinputD28only)
Pathway_CARNIVALD28 <- PathwayActivity_CARNIVALinputD28only$score
names(Pathway_CARNIVALD28) <- PathwayActivity_CARNIVALinputD28only$Pathway
head(Pathway_CARNIVALD28)
class(Pathway_CARNIVALD28)
# numeric
names(Pathway_CARNIVALD28)

carnival_D28 <- runCARNIVAL(
  inputObj  = iniciators,             
  measObj   = TFscoresD28_CARNIVAL,             # TF activity measurements
  netObj    = sif,                    # prior knowledge network
  weightObj = Pathway_CARNIVALD28,         # pathway activity scores (or progeny_matrix)
  solverPath = "C:/Program Files/IBM/ILOG/CPLEX_Studio2211/cplex/bin/x64_win64/cplex.exe", 
  solver     = "cplex",
  timelimit  = 7200, 
  mipGAP     = 0.05,
  poolrelGAP = 0
)

clean_edgesD28 <- carnival_D28$weightedSIF %>%
  filter(Weight > 0)

active_node_namesD28 <- unique(c(clean_edgesD28$Node1, clean_edgesD28$Node2))

clean_nodesD28 <- carnival_D28$nodesAttributes %>%
  filter(Node %in% active_node_namesD28)

###### visual ######
library(RCy3)
library(igraph)
library(dplyr)
write.csv(clean_edgesD28, "Clean_CARNIVAL_EdgesD28.csv", row.names = FALSE)
write.csv(clean_nodesD28, "Clean_CARNIVAL_NodesD28.csv", row.names = FALSE)

cytoscapePing() 

# Rename columns to standard RCy3 formats to be safe
# (RCy3 prefers 'source', 'target', 'interaction' but handles others too)
edges_for_cyD28 <- clean_edgesD28
colnames(edges_for_cyD28) <- c("source", "sign", "target", "weight") # Adjust based on your actual column order!

nodes_for_cyD28 <- clean_nodesD28
colnames(nodes_for_cyD28)[1] <- "id" # Ensure first column is 'id'

# Send to Cytoscape
createNetworkFromDataFrames(nodes = nodes_for_cyD28, 
                            edges = edges_for_cyD28, 
                            title = "CARNIVAL_Automated", 
                            collection = "CARNIVAL_Analysis")

setNodeColorMapping(table.column = "AvgAct", 
                    table.column.values = c(-100, 0, 100), 
                    colors = c("#0000FF", "#FFFFFF", "#FF0000"), 
                    mapping.type = "continuous")

##### top 30 TFs CARNIVAL #####
top30tfs_D28 <- TFscores_D28 %>%
  arrange(desc(abs(score))) %>%
  slice_head(n = 30)
head(top30tfs_D28)

top30tfs_D28 <- top30tfs_D28 %>%
  dplyr::select(source, score) %>%
  deframe()

carnival_D28top30 <- runCARNIVAL(
  inputObj  = iniciators,             
  measObj   = top30tfs_D28,             # TF activity measurements
  netObj    = sif,                    # prior knowledge network
  weightObj = Pathway_CARNIVALD28,         # pathway activity scores (or progeny_matrix)
  solverPath = "C:/Program Files/IBM/ILOG/CPLEX_Studio2211/cplex/bin/x64_win64/cplex.exe", 
  solver     = "cplex",
  timelimit  = 900, 
  mipGAP     = 0.05,
  poolrelGAP = 0
)

#accidentally saved as carnival_D28 
saveRDS(carnival_D28, "~/hnrnpu-causal-multiomics/processeddata/carnival_D28top30.rds")
carnival_D28top30 <- readRDS("~/hnrnpu-causal-multiomics/processeddata/carnival_D28top30.rds")
clean_edgesD28_30 <- carnival_D28top30$weightedSIF %>%
  filter(Weight > 0)

active_node_namesD28_30 <- unique(c(clean_edgesD28_30$Node1, clean_edgesD28_30$Node2))

clean_nodesD28_30 <- carnival_D28top30$nodesAttributes %>%
  filter(Node %in% active_node_namesD28_30)

write.csv(clean_edgesD28_30, "Clean_CARNIVAL_EdgesD28_30.csv", row.names = FALSE)
write.csv(clean_nodesD28_30, "Clean_CARNIVAL_NodesD28_30.csv", row.names = FALSE)

cytoscapePing() 

# Rename columns to standard RCy3 formats to be safe
# (RCy3 prefers 'source', 'target', 'interaction' but handles others too)
edges_for_cyD28_30 <- clean_edgesD28_30
colnames(edges_for_cyD28_30) <- c("source", "sign", "target", "weight") # Adjust based on your actual column order!

nodes_for_cyD28_30 <- clean_nodesD28_30
colnames(nodes_for_cyD28_30)[1] <- "id" # Ensure first column is 'id'

# Send to Cytoscape
createNetworkFromDataFrames(nodes = nodes_for_cyD28_30, 
                            edges = edges_for_cyD28_30, 
                            title = "CARNIVAL_D28_top30", 
                            collection = "CARNIVAL_Analysis_D28_top30")

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
  weightObj = Pathway_CARNIVALD28,         # pathway activity scores (or progeny_matrix)
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

# Rename columns to standard RCy3 formats to be safe
# (RCy3 prefers 'source', 'target', 'interaction' but handles others too)
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