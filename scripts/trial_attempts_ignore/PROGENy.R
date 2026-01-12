##### Trial with filtered set and transcriptutorial#####
library(progeny)
library(dorothea)
library(tibble)
library(tidyr)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(readr)
## We read the normalised counts and the experimental design 
Normalised_counts <- read_csv("~/hnrnpu-causal-multiomics/VST_progeny_filt.csv") #15116 obs
Experimental_design <- read_csv("~/hnrnpu-causal-multiomics/processeddata/targets.csv")

## We read the results from the differential analysis. 
DEGASDvsCTRL <- read_csv("~/hnrnpu-causal-multiomics/processeddata/DEG_allstatsf.csv") #17627 genes

head(Normalised_counts)
colnames(Normalised_counts)[1] <- "gene"
Normalised_counts_matrix <- as.matrix(Normalised_counts[,-1])  # remove gene column
rownames(Normalised_counts_matrix) <- Normalised_counts$gene     # set rownames

head(Normalised_counts_matrix)

DEGASDvsCTRL <- DEGASDvsCTRL %>%
  mutate(ID = as.character(ID)) %>%
  mutate(ID = trimws(ID))

DEGASDvsCTRL_matrix <- as.matrix(DEGASDvsCTRL$statistic)
rownames(DEGASDvsCTRL_matrix) <- DEGASDvsCTRL$ID
colnames(DEGASDvsCTRL_matrix) <- "statistic"
head(DEGASDvsCTRL_matrix)

PathwayActivity_countsf <- progeny(Normalised_counts_matrix, scale = TRUE,
                                   organism = "Human", top = 100)
Activity_countsf <- as.vector(PathwayActivity_countsf)

paletteLength <- 100
myColor <- 
  colorRampPalette(c("darkblue", "whitesmoke","indianred"))(paletteLength)

progenyBreaksf <- c(seq(min(Activity_countsf), 0, 
                       length.out=ceiling(paletteLength/2) + 1),
                   seq(max(Activity_countsf)/paletteLength, 
                       max(Activity_countsf), 
                       length.out=floor(paletteLength/2)))

progeny_hmapf <- pheatmap(t(PathwayActivity_countsf),fontsize=14, 
                         fontsize_row = 10, fontsize_col = 10, 
                         color=myColor, breaks = progenyBreaksf, 
                         main = "PROGENy (100)", angle_col = 45,
                         treeheight_col = 0,  border_color = NA)
ggsave(
  filename = "PROGENY_heatmap_t2f.png",
  plot = progeny_hmapf,
  width = 8,
  height = 8,
  dpi = 300
)

#enrichment analysis using competitive permutation approach to assess
#significance of pathway activity to end with NES for each pathway
PathwayActivity_zscoref <- progeny(DEGASDvsCTRL_matrix, 
                                  scale=TRUE, organism="Human", top = 100, perm = 10000, z_scores = TRUE) %>%
  t()
colnames(PathwayActivity_zscoref) <- "NES"

PathwayActivity_zscoref_df <- as.data.frame(PathwayActivity_zscoref) %>% 
  rownames_to_column(var = "Pathway") %>%
  dplyr::arrange(NES) %>%
  dplyr::mutate(Pathway = factor(Pathway))

NES_pathwayf <- ggplot(PathwayActivity_zscoref_df,aes(x = reorder(Pathway, NES), y = NES)) + 
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
  filename = "PROGENY_NES_pathwaysf.png",
  plot = NES_pathwayf,
  width = 8,
  height = 8,
  dpi = 300
)

#JAK.STAT pathway most active - visualise most responsive genes (progeny weights)
#along with their statistic values to interpret results
#scatterplot shows which genes contributing most to the pathway enrichment

prog_matrixf <- getModel("Human", top=100) %>% 
  as.data.frame()  %>%
  tibble::rownames_to_column("GeneID")

DEGASDvsCTRL_df <- DEGASDvsCTRL_matrix %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("GeneID")

scat_plotsf <- progeny::progenyScatter(df = DEGASDvsCTRL_df, 
                                      weight_matrix = prog_matrixf, 
                                      statName = "statistic", verbose = FALSE)
plot(scat_plots[[1]]$`JAK-STAT`) 

#progeny results as input for CARNIVAL
#CARNIVAL sets weights based on PROGENy scores in each pathway-related node
#to find more relevant solutions
#set progeny zscores false to return pathway activity values between 1 and -1

PathwayActivity_CARNIVALinputf <- progeny(DEGASDvsCTRL_matrix, 
                                         scale=TRUE, organism="Human", top = 100, perm = 10000, z_scores = FALSE) %>%
  t () %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "Pathway") 
colnames(PathwayActivity_CARNIVALinputf)[2] <- "score"
write_csv(PathwayActivity_CARNIVALinputf, 
          "~/hnrnpu-causal-multiomics/processeddata/PathwayActivity_CARNIVALinputf.csv")
