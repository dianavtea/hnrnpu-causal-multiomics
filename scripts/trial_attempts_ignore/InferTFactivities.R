##### 02.12.25 filtered set trial #####
library(progeny)
library(dorothea)
library(tibble)
library(tidyr)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(readr)
library(viper)

## For the volcano plot (related to support functions)
library(ggrepel)

## We also load the support functions
source("support_functions.R")


# first compute a TF activity enrichment analysis using stats from DA
#select TFs whose activity varies with conditions under study

#Virtual Inference of Protein-activity by Enriched Regulon analysis
#CollecTRI is a comprehensive resource containing a curated collection of TFs and their transcriptional 
#targets compiled from 12 different resources. 
#This collection provides an increased coverage of transcription factors and a superior performance in 
#identifying perturbed TFs compared to our previous DoRothEA network and other literature based GRNs. 
#Similar to DoRothEA, interactions are weighted by their mode of regulation (activation or inhibition

net <- decoupleR::get_collectri(organism = 'human', 
                                split_complexes = FALSE)
net
#reprocessed large scale TF regulon consensus, more up to date and broader coverage
#We can also infer TF activities from the statistics of the DEGs between ASD and CTRL:

tf_activities_statf <- decoupleR::run_viper(mat = DEGASDvsCTRL_matrix[, 'statistic', drop = FALSE], 
                                      net = net, 
                                      .source = 'source', 
                                      .target = 'target',
                                      .mor='mor', 
                                      minsize = 5)
tf_activities_statf
head(rownames(DEGASDvsCTRL_matrix))
head(net$target)
length(intersect(rownames(DEGASDvsCTRL_matrix), net$target)) #4633

tf_activities_statf_top25 <- tf_activities_statf %>%
  as.data.frame() %>% 
  dplyr::filter(p_value < 0.05) %>%  # optional: only significant TFs
  dplyr::top_n(25, wt = abs(score)) %>%
  dplyr::arrange(score) %>% 
  dplyr::mutate(source = factor(source, levels = unique(source)))

contrast_tftop25f <- ggplot(tf_activities_statf_top25,aes(x = reorder(source, score), y = score)) + 
  geom_bar(aes(fill = score), stat = "identity") +
  scale_fill_gradient2(low = "darkblue", high = "indianred", 
                       mid = "whitesmoke", midpoint = 0) + 
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = 
          element_text(angle = 45, hjust = 1, size =10, face= "bold"),
        axis.text.y = element_text(size =10, face= "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("Transcription Factors")

ggsave(
  filename = "barplot_top25tfscores.png",
  plot = contrast_tftop25f,
  width = 12,
  height = 8,
  dpi = 300
  )

#To interpret the results, we can look at the expression of targets of one of the most deregulated TFs
targets_SP1 <- net$target[net$source == "SP1"]
DEG_SP1 <- DEGASDvsCTRL[DEGASDvsCTRL$ID %in% targets_SP1, ]
volcano_nice(as.data.frame(DEG_SP1), 
             FCIndex = 3,     # log2FoldChange column
             pValIndex = 5,   # padj column
             IDIndex = 1,     # gene ID column
             nlabels = 20, 
             label = TRUE, 
             straight = FALSE)
#input for CARNIVAL
#CARNIVAL tries to infer most likely upstream signaling events leading to 
#current TF activity results
tf_activities_CARNIVALinputf <- tf_activities_statf %>%
  as.data.frame() %>%
  dplyr::select(TF = source, score)  # keep only source and score, rename source to TF

write_csv(tf_activities_CARNIVALinputf, "~/hnrnpu-causal-multiomics/processeddata/TFActivity_CARNIVALinputf.csv")

#compute TF activities per sample (with the replicates) using normalised counts
tf_activities_countsf <- decoupleR::run_viper(
  mat = Normalised_counts_matrix,
  net = net,
  .source = "source",
  .target = "target",
  .mor = "mor",
  minsize = 5,
  method = "scale",
  verbose = FALSE
)

top_tfs <- tf_activities_statf_top25$source

# Filter and reshape
tf_activities_counts_filterf <- tf_activities_countsf %>%
  filter(source %in% top_tfs) %>%       # keep only top TFs
  select(source, condition, score) %>%  # keep relevant columns
  pivot_wider(names_from = condition, values_from = score) %>%  # TFs as rows, conditions as columns
  column_to_rownames(var = "source") %>% 
  as.matrix()

# Optional: flatten into a vector
tf_activities_vectorf <- as.vector(tf_activities_counts_filter)

paletteLength <- 100
myColor <- 
  colorRampPalette(c("darkblue", "whitesmoke","indianred"))(paletteLength)

decouplerBreaks <- c(seq(min(tf_activities_vectorf), 0, 
                        length.out=ceiling(paletteLength/2) + 1),
                    seq(max(tf_activities_vectorf)/paletteLength, 
                        max(tf_activities_vectorf), 
                        length.out=floor(paletteLength/2)))
decoupler_hmap <- pheatmap(tf_activities_counts_filterf,
                          fontsize=14, fontsize_row = 8, fontsize_col = 8, 
                          color=myColor, breaks = decouplerBreaks,
                          main = "CollecTRI", angle_col = 45,
                          treeheight_col = 0,  border_color = NA)
ggsave(
  filename = "collectripersampleheatmap.png",
  plot = decoupler_hmap,
  width = 12,
  height = 8,
  dpi = 300
)
