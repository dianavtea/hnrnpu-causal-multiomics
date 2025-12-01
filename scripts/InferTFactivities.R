##### 21.11.2025 #####

##### packages #####
BiocManager::install("OmnipathR")
BiocManager::install("viper")
library(decoupleR)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(progeny)
library(dorothea)
library(readr)
library(OmnipathR)
library(viper)

##### transcriptutorial decoupleR, ULM, collecTRI #####

#CollecTRI is a comprehensive resource containing a curated collection of TFs and their transcriptional 
#targets compiled from 12 different resources. 
#This collection provides an increased coverage of transcription factors and a superior performance in 
#identifying perturbed TFs compared to our previous DoRothEA network and other literature based GRNs. 
#Similar to DoRothEA, interactions are weighted by their mode of regulation (activation or inhibition

net <- decoupleR::get_collectri(organism = 'human', 
                                split_complexes = FALSE)
net

#To infer TF enrichment scores we will run the Univariate Linear Model (ulm) method. 
#For each sample in our dataset (mat) and each TF in our network (net), it fits a linear model that predicts 
#the observed gene expression based solely on the TF’s TF-Gene interaction weights. 
#Once fitted, the obtained t-value of the slope is the score. If it is positive, we interpret that the TF 
#is active and if it is negative we interpret that it is inactive.

#To run decoupleR methods, we need an input matrix (mat), an input prior knowledge network/resource (net), 
#and the name of the columns of net that we want to use.

# Run ulm
sample_acts <- decoupleR::run_ulm(mat = vst_progeny, 
                                  net = net, 
                                  .source = 'source', 
                                  .target = 'target',
                                  .mor = 'mor', 
                                  minsize = 5)
sample_acts

#visualisation
#From the obtained results we will observe the most variable activities across samples in a heat-map
n_tfs <- 25

# Transform to wide matrix
sample_acts_mat <- sample_acts %>%
  tidyr::pivot_wider(id_cols = 'condition', 
                     names_from = 'source',
                     values_from = 'score') %>%
  tibble::column_to_rownames('condition') %>%
  as.matrix()

# Get top tfs with more variable means across clusters
tfs <- sample_acts %>%
  dplyr::group_by(source) %>%
  dplyr::summarise(std = stats::sd(score)) %>%
  dplyr::arrange(-abs(std)) %>%
  head(n_tfs) %>%
  dplyr::pull(source)

sample_acts_mat <- sample_acts_mat[,tfs]

# Scale per sample
sample_acts_mat <- scale(sample_acts_mat)

# Choose color palette
colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))
colors.use <- grDevices::colorRampPalette(colors = colors)(100)

my_breaks <- c(seq(-2, 0, length.out = ceiling(100 / 2) + 1),
               seq(0.05, 2, length.out = floor(100 / 2)))

# Plot
top25_heatmap <- pheatmap::pheatmap(mat = sample_acts_mat,
                   color = colors.use,
                   border_color = "white",
                   breaks = my_breaks,
                   cellwidth = 15,
                   cellheight = 15,
                   treeheight_row = 20,
                   treeheight_col = 20)
ggsave(
  filename = "top25_heatmap_ulm.png",
  plot = top25_heatmap,
  width = 12,
  height = 8,
  dpi = 300
)

#We can also infer TF activities from the statistics of the DEGs between ASD and CTRL:
# Run ulm
contrast_acts <- decoupleR::run_ulm(mat = degstat_mat[, 'statistic', drop = FALSE], 
                                    net = net, 
                                    .source = 'source', 
                                    .target = 'target',
                                    .mor='mor', 
                                    minsize = 5)
contrast_acts

# Filter top TFs in both signs
f_contrast_acts <- contrast_acts %>%
  dplyr::mutate(rnk = NA)

msk <- f_contrast_acts$score > 0

f_contrast_acts[msk, 'rnk'] <- rank(-f_contrast_acts[msk, 'score'])
f_contrast_acts[!msk, 'rnk'] <- rank(-abs(f_contrast_acts[!msk, 'score']))

tfs <- f_contrast_acts %>%
  dplyr::arrange(rnk) %>%
  head(n_tfs) %>%
  dplyr::pull(source)

f_contrast_acts <- f_contrast_acts %>%
  filter(source %in% tfs)

colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")[c(2, 10)])

p <- ggplot2::ggplot(data = f_contrast_acts, 
                     mapping = ggplot2::aes(x = stats::reorder(source, score), 
                                            y = score)) + 
  ggplot2::geom_bar(mapping = ggplot2::aes(fill = score),
                    color = "black",
                    stat = "identity") +
  ggplot2::scale_fill_gradient2(low = colors[1], 
                                mid = "whitesmoke", 
                                high = colors[2], 
                                midpoint = 0) + 
  ggplot2::theme_minimal() +
  ggplot2::theme(axis.title = element_text(face = "bold", size = 12),
                 axis.text.x = ggplot2::element_text(angle = 45, 
                                                     hjust = 1, 
                                                     size = 10, 
                                                     face = "bold"),
                 axis.text.y = ggplot2::element_text(size = 10, 
                                                     face = "bold"),
                 panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank()) +
  ggplot2::xlab("TFs")

p
ggsave(
  filename = "barplot_ASDvsCTRL.png",
  plot = p,
  width = 12,
  height = 8,
  dpi = 300
)

#We can further visualize the most differential target genes in each TF along their p-values to 
#interpret the results. For example, let’s see the genes that are belong to SP1
tf <- 'SP1'

df1 <- net %>%
  dplyr::filter(source == tf) %>%
  dplyr::arrange(target) %>%
  dplyr::mutate(ID = target, color = "3") %>%
  tibble::column_to_rownames('target')

inter <- sort(dplyr::intersect(rownames(deg_mat), rownames(df1)))

df1 <- df1[inter, ]

df1[,c('logfc', 't_value', 'p_value')] <- deg_mat[inter, ]

df1 <- df1 %>%
  dplyr::mutate(color = dplyr::if_else(mor > 0 & t_value > 0, '1', color)) %>%
  dplyr::mutate(color = dplyr::if_else(mor > 0 & t_value < 0, '2', color)) %>%
  dplyr::mutate(color = dplyr::if_else(mor < 0 & t_value > 0, '2', color)) %>%
  dplyr::mutate(color = dplyr::if_else(mor < 0 & t_value < 0, '1', color))

colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")[c(2, 10)])

p1 <- ggplot2::ggplot(data = df1, 
                     mapping = ggplot2::aes(x = logfc, 
                                            y = -log10(p_value), 
                                            color = color,
                                            size = abs(mor))) + 
  ggplot2::geom_point(size = 2.5, 
                      color = "black") + 
  ggplot2::geom_point(size = 1.5) +
  ggplot2::scale_colour_manual(values = c(colors[2], colors[1], "grey")) +
  ggrepel::geom_label_repel(mapping = ggplot2::aes(label = ID,
                                                   size = 1)) + 
  ggplot2::theme_minimal() +
  ggplot2::theme(legend.position = "none") +
  ggplot2::geom_vline(xintercept = 0, linetype = 'dotted') +
  ggplot2::geom_hline(yintercept = 0, linetype = 'dotted') +
  ggplot2::ggtitle(tf)

p1
ggsave(
  filename = "SP1_volcanoplot.png",
  plot = p1,
  width = 8,
  height = 8,
  dpi = 300
)

##### transcriptutorial with DOROTHEA 24.11.2025 #####
## We load Dorothea Regulons
data(dorothea_hs, package = "dorothea")
regulons <- dorothea_hs %>%
  dplyr::filter(confidence %in% c("A", "B","C"))

tf_activities_stat <- dorothea::run_viper(degstat_mat, regulons,
                                          options =  list(minsize = 5, eset.filter = FALSE, 
                                                          cores = 1, verbose = FALSE, nes = TRUE))
tf_activities_stat_top25 <- tf_activities_stat %>%
  as.data.frame() %>% 
  rownames_to_column(var = "GeneID") %>%
  dplyr::rename(NES = "statistic") %>%
  dplyr::top_n(25, wt = abs(NES)) %>%
  dplyr::arrange(NES) %>% 
  dplyr::mutate(GeneID = factor(GeneID))

tftop25 <- ggplot(tf_activities_stat_top25,aes(x = reorder(GeneID, NES), y = NES)) + 
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
  xlab("Transcription Factors")
ggsave(
  filename = "TFtop25_NES_dorothea.png",
  plot = tftop25,
  width = 10,
  height = 8,
  dpi = 300
)

targets_SP1 <- regulons$target[regulons$tf == "SP1"]
SP1 <- volcano_nice(as.data.frame(DEG_allstats[DEG_allstats$ID %in% targets_SP1,]), 
             FCIndex = 3, pValIndex = 4, IDIndex = 1,nlabels = 20, label = TRUE, 
             straight = FALSE) 
ggsave(
  filename = "SP1_dorothea.png",
  plot = SP1,
  width = 8,
  height = 12,
  dpi = 300
)

targets_NANOG <- regulons$target[regulons$tf == "NANOG"]
NANOG <- volcano_nice(as.data.frame(DEG_allstats[DEG_allstats$ID %in% targets_NANOG,]), 
                    FCIndex = 3, pValIndex = 4, IDIndex = 1,nlabels = 20, label = TRUE, 
                    straight = FALSE) 

tf_activities_CARNIVALinput<- tf_activities_stat %>%
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "TF") 
write.csv(tf_activities_CARNIVALinput, "tf_activities_CARNIVALinput.csv")

#tf activities per sample
tf_activities_counts <- 
  dorothea::run_viper(vst_progeny, regulons,
                      options =  list(minsize = 5, eset.filter = FALSE, 
                                      cores = 1, verbose = FALSE, method = c("scale")))

tf_activities_counts_filter <- tf_activities_counts %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "GeneID") %>%
  dplyr::filter(GeneID %in% tf_activities_stat_top25$GeneID) %>%
  column_to_rownames(var = "GeneID") %>%
  as.matrix()
tf_activities_vector <- as.vector(tf_activities_counts_filter)

paletteLength <- 100
myColor <- 
  colorRampPalette(c("darkblue", "whitesmoke","indianred"))(paletteLength)

dorotheaBreaks <- c(seq(min(tf_activities_vector), 0, 
                        length.out=ceiling(paletteLength/2) + 1),
                    seq(max(tf_activities_vector)/paletteLength, 
                        max(tf_activities_vector), 
                        length.out=floor(paletteLength/2)))
dorothea_hmap <- pheatmap(tf_activities_counts_filter,
                          fontsize=14, fontsize_row = 8, fontsize_col = 8, 
                          color=myColor, breaks = dorotheaBreaks,
                          main = "Dorothea ABC", angle_col = 45,
                          treeheight_col = 0,  border_color = NA)
ggsave(
  filename = "hmap_dorothea.png",
  plot = dorothea_hmap,
  width = 8,
  height = 10,
  dpi = 300
)
##### most important outputs here #####
net
sample_acts
contrast_acts
tf_activities_CARNIVALinput

#visualisation
top25_heatmap
p #barplot DEG ASD vs CTRL TFs
p1 #volcano plot SP1
tftop25
SP1
dorothea_hmap


##### support code #####
volcano_nice <- function (df, hAss = 0.05, FCIndex, pValIndex, IDIndex, vAss = NULL,
                          label = FALSE, straight = FALSE, nlabels, manual_labels = NA)
{
  df <- df[complete.cases(df), ]
  names(df)[1] <- "X1"
  hAssOri <- hAss
  hAss <- -log(hAss)
  names(df) <- gsub("adj.P.Val", "FDR", names(df))
  names(df)[FCIndex] <- "logFC"
  names(df)[pValIndex] <- "adj.P.Val"
  if (max(abs(df[, FCIndex])) >= 1) {
    xlimAbs <- ceiling(max(abs(df[, FCIndex])))
    ylimAbs <- ceiling(max(abs(-log(df[, pValIndex]))))
  }
  else {
    xlimAbs <- max(abs(df[, FCIndex]))
    ylimAbs <- max(abs(-log(df[, pValIndex])))
  }
  if (is.null(vAss)) {
    vAss <- xlimAbs/10
  }
  xneg <- function(x) abs(hAss - 1 + x/(x + vAss))
  xpos <- function(x) abs(hAss - 1 + x/(x - vAss))
  test <- function(x, y, vAss) {
    if (x < -vAss) {
      if (xneg(x) < -log(y)) {
        return("1")
      }
      else {
        return("0")
      }
    }
    else {
      if (x > vAss) {
        if (xpos(x) < -log(y)) {
          return("1")
        }
        else {
          return("0")
        }
      }
      else {
        return("0")
      }
    }
  }
  if (straight) {
    df$couleur <- ifelse(abs(df$logFC) >= vAss & df$adj.P.Val <=
                           hAssOri, "1", "0")
  }
  else {
    df$couleur <- "0"
    df$couleur <- apply(df, 1, FUN = function(x) test(as.numeric(x[FCIndex]),
                                                      as.numeric(x[pValIndex]), vAss))
  }
  df <- df[order(df$adj.P.Val, decreasing = F), ]
  df$condLabel <- df[, IDIndex]
  df[df$couleur == "0", "condLabel"] <- NA
  labels_to_keep <- c(df[c(1:nlabels), "condLabel"],manual_labels)
  df[!(df$condLabel %in% labels_to_keep), "condLabel"] <- NA
  df$couleur <- ifelse(df$couleur == "1" & df$logFC < 0, "2",
                       df$couleur)
  if (label) {
    a <- ggplot(df, aes(x = logFC, y = -log(adj.P.Val), color = couleur)) +
      geom_point(alpha = 0.5) + geom_label_repel(aes(label = condLabel)) +
      stat_function(fun = xneg, xlim = c(-xlimAbs, -vAss),
                    color = "black", alpha = 0.7) + ylim(c(0, ylimAbs)) +
      xlim(c(-xlimAbs, xlimAbs)) + stat_function(fun = xpos,
                                                 xlim = c(vAss, xlimAbs), color = "black", alpha = 0.7) +
      scale_colour_manual(values = c("grey30", "red",
                                     "royalblue3")) + theme_minimal() + theme(legend.position = "none")
  }
  else {
    if (straight) {
      a <- ggplot(df, aes(x = logFC, y = -log(adj.P.Val),
                          color = couleur)) + geom_point(alpha = 0.5) +
        geom_vline(xintercept = -vAss, color = "blue") +
        geom_vline(xintercept = vAss, color = "blue") +
        ylim(c(0, ylimAbs)) + xlim(c(-xlimAbs, xlimAbs)) +
        geom_hline(yintercept = hAss, color = "red") +
        scale_colour_manual(values = c("grey30", "red",
                                       "royalblue3")) + theme_minimal() + theme(legend.position = "none")
    }
    else {
      a <- ggplot(df, aes(x = logFC, y = -log(adj.P.Val),
                          color = couleur)) + geom_point(alpha = 0.5) +
        stat_function(fun = xneg, xlim = c(-xlimAbs,
                                           -vAss), color = "black", alpha = 0.7) + ylim(c(0,
                                                                                          ylimAbs)) + xlim(c(-xlimAbs, xlimAbs)) + stat_function(fun = xpos,
                                                                                                                                                 xlim = c(vAss, xlimAbs), color = "black", alpha = 0.7) +
        scale_colour_manual(values = c("grey30", "red",
                                       "royalblue3")) + theme_minimal() + theme(legend.position = "none")
    }
  }
  return(a)
}
