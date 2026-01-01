#CARNIVAL visualisation trial

#gemini approach
BiocManager::install("RCy3")

library(RCy3)
library(igraph)

library(dplyr)

clean_edges <- carnival_D0$weightedSIF %>%
  filter(Weight > 0)
head(clean_edges)
# 3. Filter Nodes to match
# We grab 'nodesAttributes' and keep only nodes present in our active edges
active_node_names <- unique(c(clean_edges$Node1, clean_edges$Node2))

clean_nodes <- carnival_D0$nodesAttributes %>%
  filter(Node %in% active_node_names)

# 4. Verify columns exist (You should see "Sign" here)
print(colnames(clean_edges))
print(paste("Final Edge Count:", nrow(clean_edges)))
write.csv(clean_edges, "Clean_CARNIVAL_Edges.csv", row.names = FALSE)
write.csv(clean_nodes, "Clean_CARNIVAL_Nodes.csv", row.names = FALSE)

#rcy3
library(RCy3)

# ---------------------------------------------------------
# STEP 1: Connect and Create Network
# ---------------------------------------------------------
# Make sure Cytoscape is open!
cytoscapePing() 

# Rename columns to standard RCy3 formats to be safe
# (RCy3 prefers 'source', 'target', 'interaction' but handles others too)
edges_for_cy <- clean_edges
colnames(edges_for_cy) <- c("source", "sign", "target", "weight") # Adjust based on your actual column order!

nodes_for_cy <- clean_nodes
colnames(nodes_for_cy)[1] <- "id" # Ensure first column is 'id'

# Send to Cytoscape
createNetworkFromDataFrames(nodes = nodes_for_cy, 
                            edges = edges_for_cy, 
                            title = "CARNIVAL_Automated", 
                            collection = "CARNIVAL_Analysis")

# ---------------------------------------------------------
# STEP 2: Style the Nodes (Activity Color)
# ---------------------------------------------------------
# Map 'activity_D0' to Fill Color: Blue (-100) -> White (0) -> Red (100)
setNodeColorMapping(table.column = "AvgAct", 
                    table.column.values = c(-100, 0, 100), 
                    colors = c("#0000FF", "#FFFFFF", "#FF0000"), 
                    mapping.type = "continuous")

# Lock node size (optional, makes it cleaner)
lockNodeDimensions(TRUE)
setNodeSizeDefault(40)

# ---------------------------------------------------------
# STEP 3: Style the Edges (Signs & Weights)
# ---------------------------------------------------------
# 1. Map 'sign' to Arrow Shape
# 1 = Delta (Arrow), -1 = T (Inhibitor)
setEdgeTargetArrowShapeMapping(table.column = "sign", 
                               table.column.values = c(1, -1), 
                               shapes = c("DELTA", "T"), 
                               mapping.type = "discrete")

# 2. Map 'sign' to Color (Blue=Activate, Red=Inhibit)
setEdgeColorMapping(table.column = "sign", 
                    table.column.values = c(1, -1), 
                    colors = c("#3399FF", "#FF3333"), 
                    mapping.type = "discrete")

# 3. Map 'weight' to Thickness
setEdgeLineWidthMapping(table.column = "weight", 
                        table.column.values = c(1, 100), 
                        widths = c(1, 5), 
                        mapping.type = "continuous")

# ---------------------------------------------------------
# STEP 4: Apply Layout
# ---------------------------------------------------------
# Force-directed layout to untangle the hairball
layoutNetwork(layout.name = "force-directed")

# Optional: If you have yFiles installed
layoutNetwork(layout.name = "yfiles-organic")

##### trial again 

library(dplyr)
library(stringr)

edges_D0 <- carnival_D0$weightedSIF %>%
  mutate(
    Node1 = str_trim(Node1),
    Node2 = str_trim(Node2)
  ) %>%
  filter(abs(Weight) >= 0.3)   # try 0.3–0.7 if needed

edges_D0 <- edges_D0 %>%
  filter(Node1 != "Perturbation" & Node2 != "Perturbation")
#### d28 

carnival_D28 <- readRDS("~/hnrnpu-causal-multiomics/processeddata/carnival_D28.rds")
carnival_res1 <- readRDS("~/hnrnpu-causal-multiomics/processeddata/Previous Trials/carnival_result1.rds")
edges <- carnival_res1$weightedSIF %>%
  mutate(
    Node1 = str_trim(Node1),
    Node2 = str_trim(Node2)
  ) %>%
  filter(abs(Weight) >= 0.3) 

#CARNIVAL visualisation trial

#gemini approach
BiocManager::install("RCy3")

library(RCy3)
library(igraph)

library(dplyr)

clean_edgesD28 <- carnival_D28$weightedSIF %>%
  filter(Weight > 0)
head(clean_edgesD28)
# 3. Filter Nodes to match
# We grab 'nodesAttributes' and keep only nodes present in our active edges
active_node_names <- unique(c(clean_edges$Node1, clean_edges$Node2))

clean_nodes <- carnival_D0$nodesAttributes %>%
  filter(Node %in% active_node_names)

# 4. Verify columns exist (You should see "Sign" here)
print(colnames(clean_edges))
print(paste("Final Edge Count:", nrow(clean_edges)))
write.csv(clean_edges, "Clean_CARNIVAL_Edges.csv", row.names = FALSE)
write.csv(clean_nodes, "Clean_CARNIVAL_Nodes.csv", row.names = FALSE)

#rcy3
library(RCy3)

# ---------------------------------------------------------
# STEP 1: Connect and Create Network
# ---------------------------------------------------------
# Make sure Cytoscape is open!
cytoscapePing() 

# Rename columns to standard RCy3 formats to be safe
# (RCy3 prefers 'source', 'target', 'interaction' but handles others too)
edges_for_cy <- clean_edges
colnames(edges_for_cy) <- c("source", "sign", "target", "weight") # Adjust based on your actual column order!

nodes_for_cy <- clean_nodes
colnames(nodes_for_cy)[1] <- "id" # Ensure first column is 'id'

# Send to Cytoscape
createNetworkFromDataFrames(nodes = nodes_for_cy, 
                            edges = edges_for_cy, 
                            title = "CARNIVAL_Automated", 
                            collection = "CARNIVAL_Analysis")

# ---------------------------------------------------------
# STEP 2: Style the Nodes (Activity Color)
# ---------------------------------------------------------
# Map 'activity_D0' to Fill Color: Blue (-100) -> White (0) -> Red (100)
setNodeColorMapping(table.column = "AvgAct", 
                    table.column.values = c(-100, 0, 100), 
                    colors = c("#0000FF", "#FFFFFF", "#FF0000"), 
                    mapping.type = "continuous")

# Lock node size (optional, makes it cleaner)
lockNodeDimensions(TRUE)
setNodeSizeDefault(40)

# ---------------------------------------------------------
# STEP 3: Style the Edges (Signs & Weights)
# ---------------------------------------------------------
# 1. Map 'sign' to Arrow Shape
# 1 = Delta (Arrow), -1 = T (Inhibitor)
setEdgeTargetArrowShapeMapping(table.column = "sign", 
                               table.column.values = c(1, -1), 
                               shapes = c("DELTA", "T"), 
                               mapping.type = "discrete")

# 2. Map 'sign' to Color (Blue=Activate, Red=Inhibit)
setEdgeColorMapping(table.column = "sign", 
                    table.column.values = c(1, -1), 
                    colors = c("#3399FF", "#FF3333"), 
                    mapping.type = "discrete")

# 3. Map 'weight' to Thickness
setEdgeLineWidthMapping(table.column = "weight", 
                        table.column.values = c(1, 100), 
                        widths = c(1, 5), 
                        mapping.type = "continuous")

# ---------------------------------------------------------
# STEP 4: Apply Layout
# ---------------------------------------------------------
# Force-directed layout to untangle the hairball
layoutNetwork(layout.name = "force-directed")

# Optional: If you have yFiles installed
layoutNetwork(layout.name = "yfiles-organic")
