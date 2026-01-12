#####libraries#####
library(progeny)
library(dorothea)
library(tibble)
library(tidyr)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(readr)
library(viper)

##### D0 #####
###### collectri ######
net <- decoupleR::get_collectri(organism='human', split_complexes=FALSE)
net #43159 - source, target, mor

resD0.mat <- resD0_clean[, c("hgnc_symbol", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
head(resD0.mat)
row.names(resD0.mat) <- resD0.mat$hgnc_symbol
#17017 transcripts

names(resD0.mat)[1] <- "ID"
head(resD0.mat)

#make into numeric factor
resD0.num <- resD0.mat$stat
names(resD0.num) <- resD0.mat$ID
resD0.num <- resD0.num[!is.na(resD0.num)]
any(is.na(resD0.num)) #FALSE
head(resD0.num)

TFscores_D0 <- decoupleR::run_viper(
  mat = resD0.num, #statistic matrix 
  net = net,                 
  .source = 'source', 
  .target = 'target', 
  .mor = 'mor', 
  minsize = 25
)
TFscores_D0$condition <- "HNRNPUdel_vs_CTRL"
#280 TF scores
head(TFscores_D0)

TFscoresD0 <- TFscores_D0[, c('source', 'score')]
names(TFscoresD0) <- c("TF", "score")
head(TFscoresD0)

TFscoresD0 <- TFscoresD0 %>%
  tibble::column_to_rownames("TF") %>%
  as.matrix()
head(TFscoresD0)
write.csv(TFscoresD0, "~/hnrnpu-causal-multiomics/processeddata/TFscoresD0.csv", row.names = TRUE) #280

##### D5 #####
###### collectri ######
net <- decoupleR::get_collectri(organism='human', split_complexes=FALSE)
net #43159 - source, target, mor

resD5.mat <- res_D5only_clean[, c("hgnc_symbol", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
head(resD5.mat)
row.names(resD5.mat) <- resD5.mat$hgnc_symbol
#17017 transcripts

names(resD5.mat)[1] <- "ID"
head(resD5.mat)

#make into numeric factor
resD5.num <- resD5.mat$stat
names(resD5.num) <- resD5.mat$ID
resD5.num <- resD5.num[!is.na(resD5.num)]
any(is.na(resD5.num))
head(resD5.num)

TFscores_D5 <- decoupleR::run_viper(
  mat = resD5.num, #statistic matrix 
  net = net,                 
  .source = 'source', 
  .target = 'target', 
  .mor = 'mor', 
  minsize = 25
)
TFscores_D5$condition <- "HNRNPUdel_vs_CTRL"
#280 TF scores
head(TFscores_D5)

TFscoresD5 <- TFscores_D5[, c('source', 'score')]
names(TFscoresD5) <- c("TF", "score")
head(TFscoresD5)

TFscoresD5 <- TFscoresD5 %>%
  tibble::column_to_rownames("TF") %>%
  as.matrix()
head(TFscoresD5)
write.csv(TFscoresD5, "~/hnrnpu-causal-multiomics/processeddata/TFscoresD5.csv", row.names = TRUE) #280

##### D28 #####
###### collectri ######
net <- decoupleR::get_collectri(organism='human', split_complexes=FALSE)
net #43159 - source, target, mor

resD28.mat <- res_D28only_clean[, c("hgnc_symbol", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
head(resD28.mat)
row.names(resD28.mat) <- resD28.mat$hgnc_symbol
#17017 transcripts

names(resD28.mat)[1] <- "ID"
head(resD28.mat)

#make into numeric factor
resD28.num <- resD28.mat$stat
names(resD28.num) <- resD28.mat$ID
resD28.num <- resD28.num[!is.na(resD28.num)]
any(is.na(resD28.num))
head(resD28.num)

TFscores_D28 <- decoupleR::run_viper(
  mat = resD28.num, #statistic matrix 
  net = net,                 
  .source = 'source', 
  .target = 'target', 
  .mor = 'mor', 
  minsize = 25
)
TFscores_D28$condition <- "HNRNPUdel_vs_CTRL"
#280 TF scores
head(TFscores_D28)

TFscoresD28 <- TFscores_D28[, c('source', 'score')]
names(TFscoresD28) <- c("TF", "score")
head(TFscoresD28)

TFscoresD28 <- TFscoresD28 %>%
  tibble::column_to_rownames("TF") %>%
  as.matrix()
head(TFscoresD28)
write.csv(TFscoresD28, "~/hnrnpu-causal-multiomics/processeddata/TFscoresD28.csv", row.names = TRUE) #280

##### TF dot plot #####
library(tidyverse)

top10_D0 <- TFscores_D0 %>%
  arrange(desc(abs(score))) %>%
  slice_head(n = 10)
head(top10_D0)

top10_D0 <- top10_D0 %>%
  dplyr::select(TF = source, score, p_value) %>%
  dplyr::mutate(timepoint = "D0")

top10_D5 <- TFscores_D5 %>%
  arrange(desc(abs(score))) %>%
  slice_head(n = 10)
head(top10_D5)

top10_D5 <- top10_D5 %>%
  dplyr::select(TF = source, score, p_value) %>%
  dplyr::mutate(timepoint = "D5")

top10_D28 <- TFscores_D28 %>%
  arrange(desc(abs(score))) %>%
  slice_head(n = 10)
head(top10_D28)

top10_D28 <- top10_D28 %>%
  dplyr::select(TF = source, score, p_value) %>%
  dplyr::mutate(timepoint = "D28")

# Combine them into one large dataframe
df_all <- bind_rows(top10_D0, top10_D5, top10_D28)

df_plot <- df_all %>%
  # Create the inverse P-value score for dot size
  # We use -log10 so that smaller P-values become larger numbers
  dplyr::mutate(log_p = -log10(p_value)) %>%
  
  # Ensure the timepoints appear in the correct chronological order (not alphabetical)
  dplyr::mutate(timepoint = factor(timepoint, levels = c("D0", "D5", "D28")))

limit <- max(abs(df_plot$score)) * c(-1, 1)

TFplot <- ggplot(df_plot, aes(x = timepoint, y = reorder(TF, score))) +
  geom_point(aes(color = score, size = log_p)) +
  scale_color_gradientn(
    colors = my_colors,
    limits = limit,    # Forces 0 to be the exact center
    name = "Activity\nScore"
  ) +
  scale_size_continuous(
    name = "-log10(P-val)",
    range = c(5, 10) # Adjust these numbers to make dots bigger/smaller
  ) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Differential TF Activity (decoupleR)",
    x = "Timepoint",
    y = "Transcription Factor"
  ) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 10), 
    panel.grid.major.x = element_blank()  # Clean up vertical grid lines
  )

ggsave("TFdotplot.png", width = 8, height = 10, dpi = 300)
