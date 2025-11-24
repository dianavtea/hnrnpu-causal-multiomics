##### 21.11.2025 #####

##### packages #####
library(decoupleR)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(ggrepel)

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

##### most important outputs here #####
net
sample_acts
contrast_acts

#visualisation
top25_heatmap
p #barplot DEG ASD vs CTRL TFs
p1 #volcano plot SP1