##### D28 only #####
##### PROGENy #####
library(progeny)
library(dorothea)
library(tibble)
library(tidyr)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(readr)

## We read the normalised counts and the experimental design 
nc_D28only <- read_csv("~/hnrnpu-causal-multiomics/processeddata/VST_progD28only.csv") 
#vst_progD28only
head(nc_D28only)

metafiltD28 <- meta_filtered[21:26,]

## We read the results from the differential analysis. 
DED28only <- read_csv("~/hnrnpu-causal-multiomics/DESeq2_res_D28_only.csv") #21262 genes
#res_D28_only.annot

nc_D28only_mat <- as.matrix(
  data.frame(nc_D28only, row.names = nc_D28only$...1)[, -1]
)
head(nc_D28only_mat)
# literally same as just vst_progD28only

#clean so hgnc symbols all recognised properly
DED28only <- DED28only %>%
  mutate(hgnc_symbol = as.character(hgnc_symbol)) %>%
  mutate(hgnc_symbol = trimws(hgnc_symbol))

#format for progeny
DED28only_mat <- as.matrix(DED28only$stat)
rownames(DED28only_mat) <- DED28only$hgnc_symbol
colnames(DED28only_mat) <- "statistic"
head(DED28only_mat)

#top 100 genes
PathwayActivity_countsD28only <- progeny(nc_D28only_mat, scale = TRUE,
                                     organism = "Human", top = 100)

#define heatmap colour breaks
Activity_countsD28only <- as.vector(PathwayActivity_countsD28only) 

paletteLength <- 100
myColor <- 
  colorRampPalette(c("darkblue", "whitesmoke","indianred"))(paletteLength)

progenyBreaksD28only <- c(seq(min(Activity_countsD28only), 0, 
                          length.out=ceiling(paletteLength/2) + 1),
                      seq(max(Activity_countsD28only)/paletteLength, 
                          max(Activity_countsD28only), 
                          length.out=floor(paletteLength/2)))

annotation_col <- metafiltD28 %>%
  dplyr::select(Sample_ID, condition) %>%
  as.data.frame()

rownames(annotation_col) <- annotation_col$Sample_ID
annotation_col$Sample_ID <- NULL

# force condition order
annotation_col$condition <- factor(
  annotation_col$condition,
  levels = c("CTRL", "HNRNPUdel")
)

# order samples
sample_order <- rownames(annotation_col)[order(annotation_col$condition)]

# reorder matrix
matD28only <- t(PathwayActivity_countsD28only)[, sample_order, drop = FALSE]

# reorder annotation (IMPORTANT: keep as data.frame)
annotation_col2 <- annotation_col[sample_order, , drop = FALSE]

# plot
progeny_hmapD28only <- pheatmap(
  matD28only,
  fontsize = 14,
  fontsize_row = 10,
  fontsize_col = 10,
  color = myColor,
  breaks = progenyBreaksD28only,
  main = "PROGENy (100)",
  angle_col = 45,
  treeheight_col = 0,
  border_color = NA,
  annotation_col = annotation_col2,
  cluster_cols = FALSE,
  gaps_col = sum(annotation_col2$condition == "CTRL")
)


ggsave(
  filename = "PROGENY_heatmapD28only.png",
  plot = progeny_hmapD28only,
  width = 8,
  height = 12,
  dpi = 300
)

#enrichment analysis using competitive permutation approach to assess
#significance of pathway activity to end with NES for each pathway
PathwayActivity_zscoreD28only <- progeny(DED28only_mat, 
                                     scale=TRUE, organism="Human", top = 100, perm = 10000, z_scores = TRUE) %>%
  t()
colnames(PathwayActivity_zscoreD28only) <- "NES"

PathwayActivity_zscoreD28only_df <- as.data.frame(PathwayActivity_zscoreD28only) %>% 
  rownames_to_column(var = "Pathway") %>%
  dplyr::arrange(NES) %>%
  dplyr::mutate(Pathway = factor(Pathway))

NES_pathwayD28only <- ggplot(PathwayActivity_zscoreD28only_df,aes(x = reorder(Pathway, NES), y = NES)) + 
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
  filename = "PROGENY_NES_pathwaysD28only.png",
  plot = NES_pathwayD28only,
  width = 10,
  height = 8,
  dpi = 300
)

# pathway most active - visualise most responsive genes (progeny weights)
#along with their statistic values to interpret results
#scatterplot shows which genes contributing most to the pathway enrichment

prog_matrixD28only <- getModel("Human", top=100) %>% 
  as.data.frame()  %>%
  tibble::rownames_to_column("GeneID")

DED28only_df <- DED28only_mat %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("GeneID")

sum(is.na(DED28only_df$statistic))
DED28only_df <- DED28only_df %>%
  dplyr::filter(!is.na(statistic))

scat_plotsD28 <- progeny::progenyScatter(df = DED28only_df, 
                                         weight_matrix = prog_matrixD28only, 
                                         statName = "statistic", verbose = FALSE)
plot(scat_plotsD28[[1]]$`TGFb`) 

#progeny results as input for CARNIVAL
#CARNIVAL sets weights based on PROGENy scores in each pathway-related node
#to find more relevant solutions
#set progeny zscores false to return pathway activity values between 1 and -1

PathwayActivity_CARNIVALinputD28only <- progeny(DED28only_mat, 
                                            scale=TRUE, organism="Human", top = 100, perm = 10000, z_scores = FALSE) %>%
  t () %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "Pathway") 
colnames(PathwayActivity_CARNIVALinputD28only)[2] <- "score"
head(PathwayActivity_CARNIVALinputD28only)
write_csv(PathwayActivity_CARNIVALinputD28only, 
          "~/hnrnpu-causal-multiomics/processeddata/PathwayActivity_CARNIVALinputD28.csv")


##### D28 vs D0  #####
##### D5 only #####
## We read the normalised counts and the experimental design 
nc_D5only <- read_csv("~/hnrnpu-causal-multiomics/processeddata/VST_progD5only_clean.csv") 
#vst_D5only_clean
head(nc_D5only)

metafiltD5 <- meta_filtered[c(6:8, 14:20), ]

## We read the results from the differential analysis. 
DED5only <- read_csv("~/hnrnpu-causal-multiomics/DESeq2_res_D5only_unique.csv") #21262 genes
#res_D5_only.annot

nc_D5only_mat <- as.matrix(
  data.frame(nc_D5only, row.names = nc_D5only$...1)[, -1]
)
head(nc_D5only_mat)
# literally same as just vst_progD5only

#clean so hgnc symbols all recognised properly
DED5only <- DED5only %>%
  mutate(hgnc_symbol = as.character(hgnc_symbol)) %>%
  mutate(hgnc_symbol = trimws(hgnc_symbol))

#format for progeny
DED5only_mat <- as.matrix(DED5only$stat)
rownames(DED5only_mat) <- DED5only$hgnc_symbol
colnames(DED5only_mat) <- "statistic"
head(DED5only_mat)

#top 100 genes
PathwayActivity_countsD5only <- progeny(nc_D5only_mat, scale = TRUE,
                                         organism = "Human", top = 100)

#define heatmap colour breaks
Activity_countsD5only <- as.vector(PathwayActivity_countsD5only) 

paletteLength <- 100
myColor <- 
  colorRampPalette(c("darkblue", "whitesmoke","indianred"))(paletteLength)

progenyBreaksD5only <- c(seq(min(Activity_countsD5only), 0, 
                              length.out=ceiling(paletteLength/2) + 1),
                          seq(max(Activity_countsD5only)/paletteLength, 
                              max(Activity_countsD5only), 
                              length.out=floor(paletteLength/2)))

annotation_col <- metafiltD5 %>%
  dplyr::select(Sample_ID, condition) %>%
  as.data.frame()

rownames(annotation_col) <- annotation_col$Sample_ID
annotation_col$Sample_ID <- NULL

# force condition order
annotation_col$condition <- factor(
  annotation_col$condition,
  levels = c("CTRL", "HNRNPUdel")
)

# order samples
sample_order <- rownames(annotation_col)[order(annotation_col$condition)]

# reorder matrix
matD5only <- t(PathwayActivity_countsD5only)[, sample_order, drop = FALSE]

# reorder annotation (IMPORTANT: keep as data.frame)
annotation_col2 <- annotation_col[sample_order, , drop = FALSE]

# plot
progeny_hmapD5only <- pheatmap(
  matD5only,
  fontsize = 14,
  fontsize_row = 10,
  fontsize_col = 10,
  color = myColor,
  breaks = progenyBreaksD5only,
  main = "PROGENy (100)",
  angle_col = 45,
  treeheight_col = 0,
  border_color = NA,
  annotation_col = annotation_col2,
  cluster_cols = FALSE,
  gaps_col = sum(annotation_col2$condition == "CTRL")
)


ggsave(
  filename = "PROGENY_heatmapD5only.png",
  plot = progeny_hmapD5only,
  width = 8,
  height = 12,
  dpi = 300
)

#enrichment analysis using competitive permutation approach to assess
#significance of pathway activity to end with NES for each pathway
PathwayActivity_zscoreD5only <- progeny(DED5only_mat, 
                                         scale=TRUE, organism="Human", top = 100, perm = 10000, z_scores = TRUE) %>%
  t()
colnames(PathwayActivity_zscoreD5only) <- "NES"

PathwayActivity_zscoreD5only_df <- as.data.frame(PathwayActivity_zscoreD5only) %>% 
  rownames_to_column(var = "Pathway") %>%
  dplyr::arrange(NES) %>%
  dplyr::mutate(Pathway = factor(Pathway))

NES_pathwayD5only <- ggplot(PathwayActivity_zscoreD5only_df,aes(x = reorder(Pathway, NES), y = NES)) + 
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
  filename = "PROGENY_NES_pathwaysD5only.png",
  plot = NES_pathwayD5only,
  width = 10,
  height = 8,
  dpi = 300
)

# pathway most active - visualise most responsive genes (progeny weights)
#along with their statistic values to interpret results
#scatterplot shows which genes contributing most to the pathway enrichment

prog_matrixD5only <- getModel("Human", top=100) %>% 
  as.data.frame()  %>%
  tibble::rownames_to_column("GeneID")

DED5only_df <- DED5only_mat %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("GeneID")

sum(is.na(DED5only_df$statistic))
DED5only_df <- DED5only_df %>%
  dplyr::filter(!is.na(statistic))

scat_plotsD5 <- progeny::progenyScatter(df = DED5only_df, 
                                         weight_matrix = prog_matrixD5only, 
                                         statName = "statistic", verbose = FALSE)
plot(scat_plotsD5[[1]]$`TGFb`) 

#progeny results as input for CARNIVAL
#CARNIVAL sets weights based on PROGENy scores in each pathway-related node
#to find more relevant solutions
#set progeny zscores false to return pathway activity values between 1 and -1

PathwayActivity_CARNIVALinputD5only <- progeny(DED5only_mat, 
                                                scale=TRUE, organism="Human", top = 100, perm = 10000, z_scores = FALSE) %>%
  t () %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "Pathway") 
colnames(PathwayActivity_CARNIVALinputD5only)[2] <- "score"
head(PathwayActivity_CARNIVALinputD5only)
write_csv(PathwayActivity_CARNIVALinputD5only, 
          "~/hnrnpu-causal-multiomics/processeddata/PathwayActivity_CARNIVALinputD5.csv")

##### D5 vs D0 #####
##### D0 #####
## We read the normalised counts and the experimental design 
nc_D0only <- read_csv("~/hnrnpu-causal-multiomics/processeddata/VST_progD0_clean.csv") 
#vst_progD0only
head(nc_D0only)

metafiltD0 <- meta_filtered[c(1:5, 9:13),]

## We read the results from the differential analysis. 
DED0only <- read_csv("~/hnrnpu-causal-multiomics/DESeq2_res_D0only_unique.csv") #17017 genes
#res_D0_only.annot

nc_D0only_mat <- as.matrix(
  data.frame(nc_D0only, row.names = nc_D0only$...1)[, -1]
)
head(nc_D0only_mat)
# literally same as just vst_progD0only

#clean so hgnc symbols all recognised properly
DED0only <- DED0only %>%
  mutate(hgnc_symbol = as.character(hgnc_symbol)) %>%
  mutate(hgnc_symbol = trimws(hgnc_symbol))

#format for progeny
DED0only_mat <- as.matrix(DED0only$stat)
rownames(DED0only_mat) <- DED0only$hgnc_symbol
colnames(DED0only_mat) <- "statistic"
head(DED0only_mat)

#top 100 genes
PathwayActivity_countsD0only <- progeny(nc_D0only_mat, scale = TRUE,
                                         organism = "Human", top = 100)

#define heatmap colour breaks
Activity_countsD0only <- as.vector(PathwayActivity_countsD0only) 

paletteLength <- 100
myColor <- 
  colorRampPalette(c("darkblue", "whitesmoke","indianred"))(paletteLength)

progenyBreaksD0only <- c(seq(min(Activity_countsD0only), 0, 
                              length.out=ceiling(paletteLength/2) + 1),
                          seq(max(Activity_countsD0only)/paletteLength, 
                              max(Activity_countsD0only), 
                              length.out=floor(paletteLength/2)))

annotation_col <- metafiltD0 %>%
  dplyr::select(Sample_ID, condition) %>%
  as.data.frame()

rownames(annotation_col) <- annotation_col$Sample_ID
annotation_col$Sample_ID <- NULL

# force condition order
annotation_col$condition <- factor(
  annotation_col$condition,
  levels = c("CTRL", "HNRNPUdel")
)

# order samples
sample_order <- rownames(annotation_col)[order(annotation_col$condition)]

# reorder matrix
matD0only <- t(PathwayActivity_countsD0only)[, sample_order, drop = FALSE]

# reorder annotation (IMPORTANT: keep as data.frame)
annotation_col2 <- annotation_col[sample_order, , drop = FALSE]

# plot
progeny_hmapD0only <- pheatmap(
  matD0only,
  fontsize = 14,
  fontsize_row = 10,
  fontsize_col = 10,
  color = myColor,
  breaks = progenyBreaksD0only,
  main = "PROGENy (100)",
  angle_col = 45,
  treeheight_col = 0,
  border_color = NA,
  annotation_col = annotation_col2,
  cluster_cols = FALSE,
  gaps_col = sum(annotation_col2$condition == "CTRL")
)


ggsave(
  filename = "PROGENY_heatmapD0only.png",
  plot = progeny_hmapD0only,
  width = 8,
  height = 12,
  dpi = 300
)

#enrichment analysis using competitive permutation approach to assess
#significance of pathway activity to end with NES for each pathway
PathwayActivity_zscoreD0only <- progeny(DED0only_mat, 
                                         scale=TRUE, organism="Human", top = 100, perm = 10000, z_scores = TRUE) %>%
  t()
colnames(PathwayActivity_zscoreD0only) <- "NES"

PathwayActivity_zscoreD0only_df <- as.data.frame(PathwayActivity_zscoreD0only) %>% 
  rownames_to_column(var = "Pathway") %>%
  dplyr::arrange(NES) %>%
  dplyr::mutate(Pathway = factor(Pathway))

NES_pathwayD0only <- ggplot(PathwayActivity_zscoreD0only_df,aes(x = reorder(Pathway, NES), y = NES)) + 
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
  filename = "PROGENY_NES_pathwaysD0only.png",
  plot = NES_pathwayD0only,
  width = 10,
  height = 8,
  dpi = 300
)

# pathway most active - visualise most responsive genes (progeny weights)
#along with their statistic values to interpret results
#scatterplot shows which genes contributing most to the pathway enrichment

prog_matrixD0only <- getModel("Human", top=100) %>% 
  as.data.frame()  %>%
  tibble::rownames_to_column("GeneID")

DED0only_df <- DED0only_mat %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("GeneID")

sum(is.na(DED0only_df$statistic))
DED0only_df <- DED0only_df %>%
  dplyr::filter(!is.na(statistic))

scat_plotsD0 <- progeny::progenyScatter(df = DED0only_df, 
                                         weight_matrix = prog_matrixD0only, 
                                         statName = "statistic", verbose = FALSE)
plot(scat_plotsD0[[1]]$`TGFb`) 

#progeny results as input for CARNIVAL
#CARNIVAL sets weights based on PROGENy scores in each pathway-related node
#to find more relevant solutions
#set progeny zscores false to return pathway activity values between 1 and -1

PathwayActivity_CARNIVALinputD0only <- progeny(DED0only_mat, 
                                                scale=TRUE, organism="Human", top = 100, perm = 10000, z_scores = FALSE) %>%
  t () %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "Pathway") 
colnames(PathwayActivity_CARNIVALinputD0only)[2] <- "score"
head(PathwayActivity_CARNIVALinputD0only)
write_csv(PathwayActivity_CARNIVALinputD0only, 
          "~/hnrnpu-causal-multiomics/processeddata/PathwayActivity_CARNIVALinputD0.csv")


