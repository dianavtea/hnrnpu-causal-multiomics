

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

#CollecTRI
tf_activitiesf <- read.csv("~/hnrnpu-causal-multiomics/processeddata/TFActivity_CARNIVALinputf.csv", row.names =1)
#656 rows/TFs, 2 col (TF and score)
head(tf_activitiesf)
measObj <- as.numeric(tf_activitiesf$score)
names(measObj) <- rownames(tf_activitiesf)

#PROGENy
PathwayActivity <- read.csv("~/hnrnpu-causal-multiomics/processeddata/PathwayActivity_CARNIVALinputf.csv", row.names=1)
#14 rows and 2 columns
weightObj <- as.numeric(PathwayActivity$score)
names(weightObj) <- rownames(PathwayActivity)

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


#get initial nodes
iniMTX = base::setdiff(sif$source, sif$target)
iniciators = base::data.frame(base::matrix(data = NaN, nrow = 1, ncol = length(iniMTX)), stringsAsFactors = F)
colnames(iniciators) = iniMTX

length(measObj) > 0
length(weightObj) > 0
unique(c(sif$source, sif$target))

carnival_result1 <- runCARNIVAL(
  inputObj  = iniciators,             
  measObj   = measObj,             # TF activity measurements
  netObj    = sif,                    # prior knowledge network
  weightObj = weightObj,         # pathway activity scores (or progeny_matrix)
  solverPath = "C:/Program Files/IBM/ILOG/CPLEX_Studio2211/cplex/bin/x64_win64/cplex.exe", 
  solver     = "cplex",
  timelimit  = 7200, 
  mipGAP     = 0,
  poolrelGAP = 0
)

carnival_result_test0_05 <- runCARNIVAL(
  inputObj  = iniciators,             
  measObj   = measObj,             # TF activity measurements
  netObj    = sif,                    # prior knowledge network
  weightObj = weightObj,         # pathway activity scores (or progeny_matrix)
  solverPath = "C:/Program Files/IBM/ILOG/CPLEX_Studio2211/cplex/bin/x64_win64/cplex.exe", 
  solver     = "cplex",
  timelimit  = 7200, 
  mipGAP     = 0.05,
)

carnival_result_test02 <- runCARNIVAL(
  inputObj  = iniciators,             
  measObj   = measObj,             # TF activity measurements
  netObj    = sif,                    # prior knowledge network
  weightObj = weightObj,         # pathway activity scores (or progeny_matrix)
  solverPath = "C:/Program Files/IBM/ILOG/CPLEX_Studio2211/cplex/bin/x64_win64/cplex.exe", 
  solver     = "cplex",
  timelimit  = 7200, 
  mipGAP     = 0.2,
)
str(carnival_result1)
#transform to data.frame
carnival_result1$weightedSIF <- data.frame(carnival_result1$weightedSIF, stringsAsFactors = F)
carnival_result1$weightedSIF$Sign <- as.numeric(carnival_result1$weightedSIF$Sign)
carnival_result1$weightedSIF$Weight <- as.numeric(carnival_result1$weightedSIF$Weight)

carnival_result1$nodesAttributes <- data.frame(carnival_result1$nodesAttributes, stringsAsFactors = F)
carnival_result1$nodesAttributes$ZeroAct <- as.numeric(carnival_result1$nodesAttributes$ZeroAct)
carnival_result1$nodesAttributes$UpAct <- as.numeric(carnival_result1$nodesAttributes$UpAct)
carnival_result1$nodesAttributes$DownAct <- as.numeric(carnival_result1$nodesAttributes$DownAct)
carnival_result1$nodesAttributes$AvgAct <- as.numeric(carnival_result1$nodesAttributes$AvgAct)

saveRDS(carnival_result1,"~/hnrnpu-causal-multiomics/processeddata/carnival_result1.rds")

# visualization
visNet1 = carnival_visNet(evis = carnival_result1$weightedSIF,
                         nvis = carnival_result1$nodesAttributes)
carnival_result1 <- readRDS("~/hnrnpu-causal-multiomics/processeddata/carnival_result1.rds")
#visNet
visSave(visNet1, file = paste0('carnival_visualization_visNetwork1.html'), selfcontained = TRUE)
?carnival_visNet
View(carnival_result1$weightedSIF)
View(carnival_result_test02$weightedSIF)


# 0.2 mipgap
str(carnival_result_test02)
#transform to data.frame
carnival_result_test02$weightedSIF <- data.frame(carnival_result_test02$weightedSIF, stringsAsFactors = F)
carnival_result_test02$weightedSIF$Sign <- as.numeric(carnival_result_test02$weightedSIF$Sign)
carnival_result_test02$weightedSIF$Weight <- as.numeric(carnival_result_test02$weightedSIF$Weight)

carnival_result_test02$nodesAttributes <- data.frame(carnival_result_test02$nodesAttributes, stringsAsFactors = F)
carnival_result_test02$nodesAttributes$ZeroAct <- as.numeric(carnival_result_test02$nodesAttributes$ZeroAct)
carnival_result_test02$nodesAttributes$UpAct <- as.numeric(carnival_result_test02$nodesAttributes$UpAct)
carnival_result_test02$nodesAttributes$DownAct <- as.numeric(carnival_result_test02$nodesAttributes$DownAct)
carnival_result_test02$nodesAttributes$AvgAct <- as.numeric(carnival_result_test02$nodesAttributes$AvgAct)

saveRDS(carnival_result_test02,"~/hnrnpu-causal-multiomics/processeddata/carnival_result_test02.rds")
readRDS("~/hnrnpu-causal-multiomics/processeddata/carnival_result_test02.rds")
carnival_result_test02 <- readRDS("~/hnrnpu-causal-multiomics/processeddata/carnival_result_test02.rds")

# visualization
visNet2 = carnival_visNet(evis = carnival_result_test02$weightedSIF,
                          nvis = carnival_result_test02$nodesAttributes)

#visNet
visSave(visNet2, file = paste0('carnival_visualization_visNetwork2.html'), selfcontained = TRUE)


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


##### tutorial
# Load pathways
pathways = gmt_to_csv("~/hnrnpu-causal-multiomics/c2.cp.v2025.1.Hs.symbols.gmt")

# Extract nodes and background
nodes_carnival = extractCARNIVALnodes(carnival_result_test02)

# Run GSA hyper Geometric test
sig_pathways <- runGSAhyper(genes = nodes_carnival$sucesses, 
                            universe = nodes_carnival$bg, gsc = loadGSC(pathways))
sig_pathways_df <- as.data.frame(sig_pathways$resTab)  %>% 
  tibble::rownames_to_column(var = "pathway") 

length(nodes_carnival$sucesses)
length(nodes_carnival$bg)
identical(sort(nodes_carnival$sucesses), sort(nodes_carnival$bg))

#data for plotting
library(dplyr)

PathwaysSelect <- sig_pathways_df %>%
  select(pathway, `p-value`, `Adjusted p-value`) %>%      # backticks for special names
  filter(`Adjusted p-value` <= 0.05) %>%                 # use backticks, no quotes
  rename(pvalue = `p-value`, AdjPvalue = `Adjusted p-value`) %>% 
  mutate(pathway = as.factor(pathway))                   # no quotes around pathway

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


##### trial 2
##### tutorial
# Load pathways
pathways = gmt_to_csv("~/hnrnpu-causal-multiomics/c2.cp.v2025.1.Hs.symbols.gmt")

# Extract nodes and background
nodes_carnival2 = extractCARNIVALnodes(carnival_result_test02)

# Run GSA hyper Geometric test
sig_pathways <- runGSAhyper(genes = nodes_carnival2$sucesses, 
                            universe = nodes_carnival2$bg, gsc = loadGSC(pathways))
sig_pathways_df <- as.data.frame(sig_pathways$resTab)  %>% 
  tibble::rownames_to_column(var = "pathway") 

length(nodes_carnival2$sucesses)
length(nodes_carnival2$bg)
identical(sort(nodes_carnival2$sucesses), sort(nodes_carnival2$bg))

#data for plotting
library(dplyr)

PathwaysSelect <- sig_pathways_df %>%
  select(pathway, `p-value`, `Adjusted p-value`) %>%      # backticks for special names
  filter(`Adjusted p-value` <= 0.05) %>%                 # use backticks, no quotes
  rename(pvalue = `p-value`, AdjPvalue = `Adjusted p-value`) %>% 
  mutate(pathway = as.factor(pathway))                   # no quotes around pathway

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