##### D5 analysis #####
metafiltD5 <- meta %>%
  filter(
    !grepl("^si", Sample_Name),               # removes siNTC & siHNRNPU
    grepl("^CTRL|^HNRNPU", Sample_Name),       # keeps CTRL & HNRNPUdel
    Day %in% "D5"
  )
metafiltD5$condition <- factor(
  metafiltD5$StatusHNRNPU,
  levels = c("HNRNPU_WT", "HNRNPU_KD"),
  labels = c("CTRL", "HNRNPUdel")
)
countsfiltD5 <- merged_gene_counts_all[, colnames(merged_gene_counts_all) %in% metafiltD5$Sample_ID]