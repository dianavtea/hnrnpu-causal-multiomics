##### 21.11.2025 #####
##### library installation #####
BiocManager::install("decoupleR")
BiocManager::install("progeny")
BiocManager::install("dorothea")
library(Biobase)
library(stats)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(decoupleR)
library(reshape2)
library(progeny)
library(dorothea)
library(tibble)
library(pheatmap)
library(readr)

##### trials and tests of PROGENy #####
#This function is designed for getting a model matrix with top significant genes for each pathway
progenymodel <- getModel(organism = "Human", top = 100, decoupleR = T)
head(progenymodel)
#top = top n of genes according to significance
#decoupleR true = compatible w decoupleR

model_human_full
#22479 genes, associated pathways, weight and the p-value
get("model_human_full", envir = .GlobalEnv)

##### Trial 1: Progeny with 'progeny' document as tutorial #####
#hgnc symbols needed in rows, vst counts in columns

pathwayscores <- progeny(vst_progeny, 
        scale = TRUE, 
        organism = "Human", 
        top = 100, 
        perm = 10000, 
        verbose = TRUE, 
        z_scores = TRUE, 
        get_nulldist = TRUE, 
        assay_name = "RNA", 
        return_assay = FALSE)

pathwayscores[[1]] #individual pathway scores for each

#perm = 10000 to ensure stat significance accurately determined, recommended by saez
head(pathwayscores)

##### visualisation #####
# Compute average pathway score per group
ctrl_avg <- rowMeans(pathwayscores[[1]][, c("CD28_146", "CD28_147", "CD28_148")])
asd_avg  <- rowMeans(pathwayscores[[1]][, c("AD28_144", "AD28_145", "AD28_149")])

# Create a data frame for plotting
df <- data.frame(
  pathway = names(ctrl_avg),
  Ctrl = ctrl_avg,
  ASD  = asd_avg
)
head(df)

#scatterplot
scatterplot <- ggplot(df, aes(x = Ctrl, y = ASD, label = pathway)) +
  geom_point(color = "darkred", size = 3) +
  geom_text(vjust = -0.5, size = 3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(
    title = "PROGENY Pathway Activity: Ctrl vs ASD",
    x = "Ctrl (average z-score)",
    y = "ASD (average z-score)"
  )
ggsave(
  filename = "PROGENY_scatterplot.png",
  plot = scatterplot,
  width = 8,
  height = 8,
  dpi = 300
)

#barplot
df_long <- data.frame(
  pathway = rep(names(ctrl_avg), 2),
  value = c(ctrl_avg, asd_avg),
  group = rep(c("Ctrl", "ASD"), each = length(ctrl_avg))
)

activitybarplot <- ggplot(df_long, aes(x = pathway, y = value, fill = group)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(
    title = "PROGENY Pathway Activity",
    y = "Average z-score",
    x = "Pathway"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(
  filename = "PROGENY_pathway_activity.png",
  plot = activitybarplot,
  width = 8,
  height = 8,
  dpi = 300
)


#per sample analysis
# Extract the pathway scores matrix
scores <- as.matrix(pathwayscores[[1]])  # pathways x samples

# Define sample groups
ctrl_samples <- c("CD28_146", "CD28_147", "CD28_148")
asd_samples  <- c("AD28_144", "AD28_145", "AD28_149")

#Convert to long format for ggplot
df_long <- data.frame(
  pathway = rep(rownames(scores), times = ncol(scores)),
  sample = rep(colnames(scores), each = nrow(scores)),
  value  = as.vector(scores),
  group  = rep(c(rep("Ctrl", length(ctrl_samples)), rep("ASD", length(asd_samples))),
               each = nrow(scores))
)

#Plot boxplots + jitter for all pathways
boxplot <- ggplot(df_long, aes(x = pathway, y = value, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.4, position = position_dodge(width = 0.8)) +
  geom_jitter(aes(color = group), size = 2, position = position_dodge(width = 0.8)) +
  theme_minimal() +
  labs(title = "PROGENY Pathway Scores per Sample",
       x = "Pathway",
       y = "PROGENY z-score") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(
  filename = "PROGENY_boxplotpersample.png",
  plot = boxplot,
  width = 8,
  height = 8,
  dpi = 300
)

##### Trial 2: progeny with transcriptutorial #####
PathwayActivity_counts <- progeny(vst_progeny, scale=TRUE, 
                                  organism="Human", top = 100)
Activity_counts <- as.vector(PathwayActivity_counts)

#heatmap
paletteLength <- 100
myColor <- 
  colorRampPalette(c("darkblue", "whitesmoke","indianred"))(paletteLength)

progenyBreaks <- c(seq(min(Activity_counts), 0, 
                       length.out=ceiling(paletteLength/2) + 1),
                   seq(max(Activity_counts)/paletteLength, 
                       max(Activity_counts), 
                       length.out=floor(paletteLength/2)))

progeny_hmap <- pheatmap(t(PathwayActivity_counts),fontsize=14, 
                         fontsize_row = 10, fontsize_col = 10, 
                         color=myColor, breaks = progenyBreaks, 
                         main = "PROGENy (Top 100)", angle_col = 45,
                         treeheight_col = 0,  border_color = NA)
ggsave(
  filename = "PROGENY_heatmap_t2.png",
  plot = progeny_hmap,
  width = 8,
  height = 8,
  dpi = 300
)

#Now, we run an enrichment analysis using a competitive permutation approach to 
#assess the significance of the pathway activity. 
#We end up with Normalised Enrichment Scores (NES) for each pathway.

PathwayActivity_zscore <- progeny(degstat_mat, 
                                  scale=TRUE, organism="Human", top = 100, perm = 10000, z_scores = TRUE) %>%
  t()
colnames(PathwayActivity_zscore) <- "NES"

PathwayActivity_zscore_df <- as.data.frame(PathwayActivity_zscore) %>% 
  rownames_to_column(var = "Pathway") %>%
  dplyr::arrange(NES) %>%
  dplyr::mutate(Pathway = factor(Pathway))

NES_pathway <- ggplot(PathwayActivity_zscore_df,aes(x = reorder(Pathway, NES), y = NES)) + 
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
  xlab("Pathways")
ggsave(
  filename = "PROGENY_NES_pathways.png",
  plot = NES_pathway,
  width = 8,
  height = 8,
  dpi = 300
)

#TGFb most active pathway upon perturbation (HNRNPUdel vs CTRL)
#visualise MAPK most responsive genes (progeny_weights) with t vals to interpret results
#in scatterplot we can see genes contributing most to pathway enrichment
prog_matrix <- getModel("Human", top=100) %>% 
  as.data.frame()  %>%
  tibble::rownames_to_column("GeneID")

degstat_mat_df <- degstat_mat %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("ID")

scat_plots <- progeny::progenyScatter(df = degstat_mat_df, 
                                      weight_matrix = prog_matrix, 
                                      statName = "statistic", verbose = FALSE)
plot(scat_plots[[1]]$`TGFb`) 

##### progeny input for CARNIVAL #####
#CARNIVAL sets weights based on **PROGENy** scores in each pathway-related node in order
#to find more relevant solutions. We therefore run **PROGENy** again with 
#slightly different parameters, setting `z_scores = FALSE` so that **PROGENy** returns pathway activity 
#values between 1 and -1, rather than converting to Z-Scores.

PathwayActivity_CARNIVALinput <- progeny(degstat_mat, 
                                         scale=TRUE, organism="Human", top = 100, perm = 10000, z_scores = FALSE) %>%
  t () %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "Pathway") 
colnames(PathwayActivity_CARNIVALinput)[2] <- "score"

write.csv(PathwayActivity_CARNIVALinput, "PathwayActivity_CARNIVALinput.csv")

##### most important outputs here #####
#trial 1 'progeny' doc
pathwayscores[[1]]
pathwayscores
df #avg CTRL and ASD pathway activity scores

#visualisation
scatterplot
activitybarplot
boxplot #per sample boxplot
progeny_hmap
NES_pathway
plot(scat_plots[[1]]$`TGFb`) 

#trial 2 'transcriptutorial'
PathwayActivity_counts
PathwayActivity_zscore_df
PathwayActivity_CARNIVALinput 