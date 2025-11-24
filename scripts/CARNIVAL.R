##### 24.11.2025 #####
##### CARNIVAL transcriptutorial non-adapted or modified TEST! #####
#network reconstruction
source("support_functions.R")

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
system("cbc")

## We also load the support functions
source("assignPROGENyScores.r")
source("generateTFList.r")
source("carnival_visNetwork.r")

## We read the normalised counts and the experimental design 
tf_activities <- read_csv("~/Masters/RP2/hnrnpu-causal-multiomics/processeddata/tf_activities_CARNIVALinput.csv")
PathwayActivity <- read_csv("~/Masters/RP2/hnrnpu-causal-multiomics/processeddata/PathwayActivity_CARNIVALinput.csv")

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

# run carnival
carnival_result = runCARNIVAL(inputObj= iniciators,
                               measObj = measVec, 
                               netObj = sif, 
                               weightObj = progenylist$score, 
                               solverPath = "~/IBM/ILOG/CPLEX_Studio_Community2212/cplex", 
                               solver = "cplex",
                               timelimit=7200,
                               mipGAP=0,
                               poolrelGAP=0 )


# If you used top=50 and access_idx = 1
tf_df <- tfList[[1]]              # first element of the list
measVec <- as.numeric(tf_df[1, ]) # extract numeric values
names(measVec) <- colnames(tf_df) # TF names

carnival_result = runCARNIVAL(inputObj= iniciators,
                              measObj = measVec, 
                              netObj = sif, 
                              weightObj = progenylist$score, 
                              solverPath = "C:/Cbc/bin/cbc.exe", 
                              solver = "cbc",
                              timelimit=7200,
                              mipGAP=0,
                              poolrelGAP=0)

#DID NOT WORK
carnival_result <- runCARNIVAL(
  inputObj = iniciators,
  measObj = measVec,
  netObj = sif,
  weightObj = progenylist$score,
  solver = "cbc",
  timelimit=7200,
  mipGAP=0,
  poolrelGAP=0 )

#if figure out above visualisation:
#transoform to data.frame
carnival_result$weightedSIF <- data.frame(carnival_result$weightedSIF, stringsAsFactors = F)
carnival_result$weightedSIF$Sign <- as.numeric(carnival_result$weightedSIF$Sign)
carnival_result$weightedSIF$Weight <- as.numeric(carnival_result$weightedSIF$Weight)

carnival_result$nodesAttributes <- data.frame(carnival_result$nodesAttributes, stringsAsFactors = F)
carnival_result$nodesAttributes$ZeroAct <- as.numeric(carnival_result$nodesAttributes$ZeroAct)
carnival_result$nodesAttributes$UpAct <- as.numeric(carnival_result$nodesAttributes$UpAct)
carnival_result$nodesAttributes$DownAct <- as.numeric(carnival_result$nodesAttributes$DownAct)
carnival_result$nodesAttributes$AvgAct <- as.numeric(carnival_result$nodesAttributes$AvgAct)

saveRDS(carnival_result,"~/Masters/RP2/hnrnpu-causal-multiomics/processeddata/carnival_result.rds")

# visualization
visNet = carnival_visNet(evis = carnival_result$weightedSIF,
                         nvis = carnival_result$nodesAttributes)
#visNet
visSave(visNet, file = paste0('carnival_visualization_visNetwork.html'), selfcontained = TRUE)

##### analysis of CARNIVAL results #####
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
carnival_result = readRDS("~/Masters/RP2/hnrnpu-causal-multiomics/processeddata/carnival_result.rds")



pkn = read_tsv("../results/omnipath_carnival.tsv")