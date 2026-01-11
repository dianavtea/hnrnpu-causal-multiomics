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
any(is.na(resD0.num))
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

n_tfs <- 25

library(dplyr)
## top 25 TF in general (abs)
# Filter for top 25 most active TFs (by absolute score)
top_tfsD0 <- TFscores_D0 %>%
  arrange(desc(abs(score))) %>%
  slice_head(n = 25)

# Plot them directly 
top25allD0 <- ggplot(top_tfsD0, aes(x = reorder(source, score), y = score, fill = score > 0)) +
  geom_col() +
  coord_flip() + # Makes it horizontal
  theme_minimal() +
  labs(title = "Top 25 TF Activities (ASD vs CTRL)", 
       y = "Normalized Enrichment Score (NES)", 
       x = "Transcription Factor") +
  scale_fill_manual(values = c("TRUE" = "#E41A1C", "FALSE" = "#377EB8"), 
                    labels = c("Activated", "Inhibited"),
                    name = "Status")

ggsave(
  filename = "top25abs_D0.png",      # use PDF or PNG
  plot = top25allD0,           
  width = 8, height = 10,          # size in inches
  dpi = 300                       # resolution for PNG (ignored for PDF)
)

##top 15 for activated and top 15 for inhibited TFs
# Select Top 15 Activated (Highest Positive Scores)
top_up <- TFscores_D0 %>%
  arrange(desc(score)) %>%
  slice_head(n = 15)

# Select Top 15 Inhibited (Lowest Negative Scores)
top_down <- TFscores_D0 %>%
  arrange(score) %>%
  slice_head(n = 15)

# Combine them into one dataframe
top_30_balanced <- bind_rows(top_up, top_down)

# We use reorder() to ensure they appear sorted by score, not alphabetically
top30D0 <- ggplot(top_30_balanced, aes(x = reorder(source, score), y = score, fill = score > 0)) +
  geom_col() +
  coord_flip() + # Flips to horizontal bars for readability
  theme_minimal() +
  labs(title = "Top 15 Activated & Inhibited TFs",
       subtitle = "Differential Analysis (ASD vs CTRL)",
       x = "Transcription Factor",
       y = "NES (Activity Score)",
       fill = "Regulation") +
  scale_fill_manual(values = c("TRUE" = "#E41A1C", "FALSE" = "#377EB8"), 
                    labels = c("Inhibited", "Activated")) +
  theme(
    axis.text.y = element_text(size = 10), # Adjust text size if needed
    legend.position = "top"
  )

ggsave(
  filename = "top30TFs_D0.png",      # use PDF or PNG
  plot = top30D0,           
  width = 10, height = 8,          # size in inches
  dpi = 300                       # resolution for PNG (ignored for PDF)
)

## For the volcano plot (related to support functions)
library(ggrepel)
## We also load the support functions
source("support_functions.R")

#To interpret the results, we can look at the expression of targets of one of the most deregulated TFs
targets_ETS1 <- net$target[net$source == "ETS1"]
DEG_ETS1 <- DED0only[DED0only$hgnc_symbol %in% targets_ETS1, ]
volcano_nice(as.data.frame(DEG_ETS1), 
             FCIndex = 3,     # log2FoldChange column
             pValIndex = 7,   # padj column
             IDIndex = 9,     # gene ID column
             nlabels = 20, 
             label = TRUE, 
             straight = FALSE)

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

n_tfs <- 25

library(dplyr)
## top 25 TF in general (abs)
# Filter for top 25 most active TFs (by absolute score)
top_tfsD5 <- TFscores_D5 %>%
  arrange(desc(abs(score))) %>%
  slice_head(n = 25)

# Plot them directly 
top25allD5 <- ggplot(top_tfsD5, aes(x = reorder(source, score), y = score, fill = score > 0)) +
  geom_col() +
  coord_flip() + # Makes it horizontal
  theme_minimal() +
  labs(title = "Top 25 TF Activities (ASD vs CTRL)", 
       y = "Normalized Enrichment Score (NES)", 
       x = "Transcription Factor") +
  scale_fill_manual(values = c("TRUE" = "#E41A1C", "FALSE" = "#377EB8"), 
                    labels = c("Activated", "Inhibited"),
                    name = "Status")

ggsave(
  filename = "top25abs_D5.png",      # use PDF or PNG
  plot = top25allD5,           
  width = 8, height = 10,          # size in inches
  dpi = 300                       # resolution for PNG (ignored for PDF)
)

##top 15 for activated and top 15 for inhibited TFs
# Select Top 15 Activated (Highest Positive Scores)
top_up <- TFscores_D5 %>%
  arrange(desc(score)) %>%
  slice_head(n = 15)

# Select Top 15 Inhibited (Lowest Negative Scores)
top_down <- TFscores_D5 %>%
  arrange(score) %>%
  slice_head(n = 15)

# Combine them into one dataframe
top_30_balanced <- bind_rows(top_up, top_down)

# We use reorder() to ensure they appear sorted by score, not alphabetically
top30D5 <- ggplot(top_30_balanced, aes(x = reorder(source, score), y = score, fill = score > 0)) +
  geom_col() +
  coord_flip() + # Flips to horizontal bars for readability
  theme_minimal() +
  labs(title = "Top 15 Activated & Inhibited TFs",
       subtitle = "Differential Analysis (ASD vs CTRL)",
       x = "Transcription Factor",
       y = "NES (Activity Score)",
       fill = "Regulation") +
  scale_fill_manual(values = c("TRUE" = "#E41A1C", "FALSE" = "#377EB8"), 
                    labels = c("Inhibited", "Activated")) +
  theme(
    axis.text.y = element_text(size = 10), # Adjust text size if needed
    legend.position = "top"
  )

ggsave(
  filename = "top30TFs_D5.png",      # use PDF or PNG
  plot = top30D5,           
  width = 10, height = 8,          # size in inches
  dpi = 300                       # resolution for PNG (ignored for PDF)
)

## For the volcano plot (related to support functions)
library(ggrepel)
## We also load the support functions
source("support_functions.R")

#To interpret the results, we can look at the expression of targets of one of the most deregulated TFs
targets_SP1 <- net$target[net$source == "SP1"]
DEG_SP1 <- DED5only[DED5only$hgnc_symbol %in% targets_SP1, ]
volcano_nice(as.data.frame(DEG_SP1), 
             FCIndex = 3,     # log2FoldChange column
             pValIndex = 7,   # padj column
             IDIndex = 9,     # gene ID column
             nlabels = 20, 
             label = TRUE, 
             straight = FALSE)

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

n_tfs <- 25

library(dplyr)
## top 25 TF in general (abs)
# Filter for top 25 most active TFs (by absolute score)
top_tfsD28 <- TFscores_D28 %>%
  arrange(desc(abs(score))) %>%
  slice_head(n = 25)

# Plot them directly 
top25allD28 <- ggplot(top_tfsD28, aes(x = reorder(source, score), y = score, fill = score > 0)) +
  geom_col() +
  coord_flip() + # Makes it horizontal
  theme_minimal() +
  labs(title = "Top 25 TF Activities (ASD vs CTRL)", 
       y = "Normalized Enrichment Score (NES)", 
       x = "Transcription Factor") +
  scale_fill_manual(values = c("TRUE" = "#E41A1C", "FALSE" = "#377EB8"), 
                    labels = c("Activated", "Inhibited"),
                    name = "Status")

ggsave(
  filename = "top25abs_D28.png",      # use PDF or PNG
  plot = top25allD28,           
  width = 8, height = 10,          # size in inches
  dpi = 300                       # resolution for PNG (ignored for PDF)
)

##top 15 for activated and top 15 for inhibited TFs
# Select Top 15 Activated (Highest Positive Scores)
top_up <- TFscores_D28 %>%
  arrange(desc(score)) %>%
  slice_head(n = 15)

# Select Top 15 Inhibited (Lowest Negative Scores)
top_down <- TFscores_D28 %>%
  arrange(score) %>%
  slice_head(n = 15)

# Combine them into one dataframe
top_30_balanced <- bind_rows(top_up, top_down)

# We use reorder() to ensure they appear sorted by score, not alphabetically
top30D28 <- ggplot(top_30_balanced, aes(x = reorder(source, score), y = score, fill = score > 0)) +
  geom_col() +
  coord_flip() + # Flips to horizontal bars for readability
  theme_minimal() +
  labs(title = "Top 15 Activated & Inhibited TFs",
       subtitle = "Differential Analysis (ASD vs CTRL)",
       x = "Transcription Factor",
       y = "NES (Activity Score)",
       fill = "Regulation") +
  scale_fill_manual(values = c("TRUE" = "#E41A1C", "FALSE" = "#377EB8"), 
                    labels = c("Inhibited", "Activated")) +
  theme(
    axis.text.y = element_text(size = 10), # Adjust text size if needed
    legend.position = "top"
  )

ggsave(
  filename = "top30TFs_D28.png",      # use PDF or PNG
  plot = top30D28,           
  width = 10, height = 8,          # size in inches
  dpi = 300                       # resolution for PNG (ignored for PDF)
)

## For the volcano plot (related to support functions)
library(ggrepel)
## We also load the support functions
source("support_functions.R")

#To interpret the results, we can look at the expression of targets of one of the most deregulated TFs
targets_SP1 <- net$target[net$source == "SP1"]
DEG_SP1 <- DED28only[DED28only$hgnc_symbol %in% targets_SP1, ]
volcano_nice(as.data.frame(DEG_SP1), 
             FCIndex = 3,     # log2FoldChange column
             pValIndex = 7,   # padj column
             IDIndex = 9,     # gene ID column
             nlabels = 20, 
             label = TRUE, 
             straight = FALSE)


###### OPTIONAL #######
#ALT way of representing TFs
tf_activities_statf_top25 <- TFscores_D0 %>%
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