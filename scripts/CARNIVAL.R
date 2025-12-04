

##### CARNIVAL transcriptutorial filtered set #####
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

source("assignPROGENyScores.r")
source("generateTFList.r")
source("carnival_visNetwork.r")

tf_activitiesf <- read_csv("~/hnrnpu-causal-multiomics/processeddata/TFActivity_CARNIVALinputf.csv")
#656 rows/TFs, 2 col (TF and score)
head(tf_activitiesf)
PathwayActivity <- read_csv("~/hnrnpu-causal-multiomics/processeddata/PathwayActivity_CARNIVALinputf.csv")
#14 rows and 2 columns

# creating OmniPath PKN scaffold

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
write_tsv(sif, "~/hnrnpu-causal-multiomics/processeddata/omnipath_carnival.tsv")

tf_vector <- tf_activitiesf %>%
  dplyr::slice_max(abs(score), n = 50) %>%   # same as top = 50
  tibble::column_to_rownames("TF") %>%
  pull(score)
head(tf_vector)

progeny_vector <- PathwayActivity %>%
  tibble::column_to_rownames("Pathway") %>%
  pull(score)

tf_activities_carnival <- data.frame(tf_activitiesf, stringsAsFactors = F)
rownames(tf_activities_carnival) <- tf_activitiesf$TF #656 TFs
tf_activities_carnival$TF <- NULL
tfList = generateTFList(tf_activities_carnival, top=50, access_idx = 1)

load(file = system.file("progenyMembers.RData",package="CARNIVAL"))

PathwayActivity_carnival <- data.frame(PathwayActivity, stringsAsFactors = F)
rownames(PathwayActivity_carnival) <- PathwayActivity_carnival$Pathway
PathwayActivity_carnival$Pathway <- NULL
progenylist = assignPROGENyScores(progeny = t(PathwayActivity_carnival), 
                                  progenyMembers = progenyMembers, 
                                  id = "gene", 
                                  access_idx = 1)

# get initial nodes
iniMTX = base::setdiff(sif$source, sif$target)
iniciators = base::data.frame(base::matrix(data = NaN, nrow = 1, ncol = length(iniMTX)), stringsAsFactors = F)
colnames(iniciators) = iniMTX

sum(names(tf_vector) %in% unique(c(sif$source, sif$target)))

# run carnival
carnival_result = runCARNIVAL( inputObj= iniciators,
                               measObj = tfList$score, 
                               netObj = sif, 
                               weightObj = progenylist$score, 
                               solverPath = "C:/Program Files/IBM/ILOG/CPLEX_Studio2211/cplex/bin/x64_win64/cplex.exe", 
                               solver = "cplex",
                               timelimit=7200,
                               mipGAP=0,
                               poolrelGAP=0 )

#transoform to data.frame
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


# ANALYSIS CARNIVAL results
packageDescription("CARNIVAL")
vignette("CARNIVAL")


##### trial again adapted to carnival 2.20.0

tf_measObj <- t(as.matrix(tf_vector))

tf_measObj <- as.data.frame(tf_measObj)
str(tf_measObj)
head(tf_measObj)
progeny_matrix <- t(as.matrix(progeny_vector))

iniMTX = base::setdiff(sif$source, sif$target)
iniciators = base::data.frame(base::matrix(data = NaN, nrow = 1, ncol = length(iniMTX)), stringsAsFactors = F)
colnames(iniciators) = iniMTX


sum(colnames(tf_measObj) %in% unique(c(sif$source, sif$target)))

carnival_result <- runCARNIVAL(
  inputObj  = iniciators,             # perturbation object
  measObj   = tf_measObj,             # TF activity measurements
  netObj    = sif,                    # prior knowledge network
  weightObj = progeny_vector,         # pathway activity scores (or progeny_matrix)
  solverPath = "C:/Program Files/IBM/ILOG/CPLEX_Studio2211/cplex/bin/x64_win64/cplex.exe", 
  solver     = "cplex",
  timelimit  = 7200, 
  mipGAP     = 0,
  poolrelGAP = 0
)
str(carnival_result)

# ANALYSIS
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

install.packages("msigdbr")
library(msigdbr)
msigdbr_collections()  
go_terms <- msigdbr(species = "Homo sapiens", collection = "C5") 
head(go_terms)

brain_go <- go_terms[grep("brain|neuro|synapse|axon|neuron|autism|neurodev", go_terms$gs_name, ignore.case = TRUE), ]
unique(brain_go$gs_name)
brain_gene_sets <- split(brain_go$gene_symbol, brain_go$gs_name)
head(brain_gene_sets)

nodes_carnival # your 624 CARNIVAL nodes
gsc_brain <- list(
  genesets = brain_gene_sets,
  gs_genes = unique(brain_go$gene_symbol),
  gs_descr = names(brain_gene_sets)
)
head(gsc_brain)

##
str(carnival_result$sifAll[[1]])
colnames(carnival_result$sifAll[[1]])
# Combine all edges from all solutions
all_edges <- do.call(rbind, lapply(carnival_result$sifAll, function(sol) {
  as.data.frame(sol)[, c("Node1", "Node2")]
}))

# Extract unique nodes for enrichment
nodes_carnival <- unique(c(all_edges$Node1, all_edges$Node2))

#THIS ONE HERE
nodes_carnival <- extractCARNIVALnodes(carnival_result)
# Check
length(nodes_carnival)
head(nodes_carnival)
str(nodes_carnival)



#ignore
nodes_background <- unique(c(carnival_result$nodesAttributes$Node))  # or whatever column has all nodes
length(nodes_background)

#sfari
sfari_df <- read.csv("SFARI-Gene_genes_10-23-2025release_12-03-2025export.csv", stringsAsFactors = FALSE)
#convert to GSC here
gsc_sfari <- list(
  genesets = list(SFARI = sfari_genes),
  gs_id <- "SFARI",
  gs_genes = sfari_genes,
  gs_descr = "SFARI curated autism gene set"
)
head(gsc_sfari)
class(gsc_sfari) <- "GSC"

colnames(sfari_df)
sfari_genes <- unique(sfari_df$gene.symbol)  # change column name accordingly
length(sfari_genes)
length(intersect(nodes_carnival$bg, sfari_genes))

length(gsc_sfari$genesets$SFARI)
gsc_sfari$gs_id
class(gsc_sfari)

## new approach
#create SFARI GMT file
sfari_genes <- unique(sfari_df$gene.symbol) #genes that exist in the universe
sfari_genes <- intersect(sfari_genes, nodes_carnival$bg)
#815
gsc_sfari <- list(
  genesets = list(SFARI = sfari_genes),
  gs_id    = "SFARI",
  gs_descr = "SFARI autism gene set",
  gs_genes = sfari_genes
)

length(gsc_sfari$genesets$SFARI)  # 815
length(intersect(nodes_carnival$bg, gsc_sfari$genesets$SFARI))  #815

sfari_gsc_df <- data.frame(
  gene_symbol = sfari_genes,
  gene_set = rep("SFARI", length(sfari_genes)),  # repeat "SFARI" for each gene
  stringsAsFactors = FALSE
)
head(sfari_gsc_df)
loadGSC(sfari_gsc_df)
# check
head(sfari_gsc_df)
sfari_genes <- intersect(sfari_df$gene.symbol, nodes_carnival$bg)

sig_sfari <- runGSAhyper(
  genes = nodes_carnival$sucesses,
  universe = nodes_carnival$bg,
  gsc = loadGSC(sfari_gsc_df)
)
head(sfari_df)
length(intersect(nodes_carnival$sucesses, nodes_carnival$bg))
## go from here after nodes_carnival
#sfari 
sig_pathways_sfari <- runGSAhyper(genes = nodes_carnival$sucesses, 
                            universe = nodes_carnival$bg, gsc = loadGSC(gsc_sfari))
sig_pathways_df <- as.data.frame(sig_pathways$resTab)  %>% 
  tibble::rownames_to_column(var = "pathway")

#gsc_brain