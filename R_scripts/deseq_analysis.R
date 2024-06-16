### Differential gene expression analysis of stage 1 vs. stage 4 breast cancer patients from the TCGA-BRCA project ###

## Loading packages ##

library("TCGAbiolinks")
library("tidyverse")
library("pheatmap")
library("DESeq2")
library("ggplot2")
library("clusterProfiler")
library("org.Hs.eg.db")

## Querying clinical data to use for filtering by stage of cancer ##

clinical_data <- GDCquery_clinic("TCGA-BRCA")

stage_1_cases <- clinical_data[clinical_data$ajcc_pathologic_stage %in% c("Stage I", "Stage IA", "Stage IA1", "Stage IA2", "Stage IA3", "Stage IB", "Stage IB1", "Stage IB2", "Stage IC"), ]
stage_4_cases <- clinical_data[clinical_data$ajcc_pathologic_stage %in% c("Stage IV", "Stage IVA", "Stage IVB", "Stage IVC"), ]

stage_1_ids <- stage_1_cases$submitter_id
stage_4_ids <- stage_4_cases$submitter_id

## Querying RNA-seq data for each group ##

query_stage1 <- GDCquery(project = "TCGA-BRCA", data.category = "Transcriptome Profiling", 
                         data.type = "Gene Expression Quantification", 
                         sample.type = c("Primary Tumor"), barcode = stage_1_ids)

GDCdownload(query_stage1)

query_stage4 <- GDCquery(project = "TCGA-BRCA", data.category = "Transcriptome Profiling", 
                           data.type = "Gene Expression Quantification", 
                           sample.type = c("Primary Tumor"), barcode = stage_4_ids)

GDCdownload(query_stage4)

## Data Processing ##

stage1_RNASeq <- GDCprepare(query_stage1, save = TRUE, summarizedExperiment = TRUE, save.filename = "stage1_RNASeq_raw.rda")
stage4_RNASeq <- GDCprepare(query_stage4, save = TRUE, summarizedExperiment = TRUE, save.filename = "stage4_RNASeq_raw.rda")

  ## Preparing gene expression matrix
stage1_removedOutliers <- TCGAanalyze_Preprocessing(object = stage1_RNASeq, cor.cut = 0.85, datatype = "unstranded", filename = "stage1_correlation_RNASeq.png")
stage4_removedOutliers <- TCGAanalyze_Preprocessing(object = stage4_RNASeq, cor.cut = 0.85, datatype = "unstranded", filename = "stage4_correlation_RNASeq.png")


normalized_counts <- TCGAanalyze_Normalization(tabDF = cbind(stage1_removedOutliers, stage4_removedOutliers), geneInfo = geneInfoHT)

filtered_counts <- TCGAanalyze_Filtering(tabDF = normalized_counts, method = "quantile", qnt.cut = 0.05)
save(filtered_counts, file = paste0("combined_normalized_RNASeq_counts.rda"))

## DESeq2 Analysis ##

  ## Formatting data properly for DESeq2
samples <- colnames(filtered_counts)
sample_prefixes <- substr(samples, 1, 12)
stage_info <- data.frame(Sample = samples, Stage = NA, stringsAsFactors = FALSE)
rownames(stage_info) <- samples

stage_info$Stage[sample_prefixes %in% stage_1_ids] <- "Stage 1"
stage_info$Stage[sample_prefixes %in% stage_4_ids] <- "Stage 4"

stage_info <- stage_info[, "Stage", drop = FALSE]

all(colnames(filtered_counts %in% rownames(stage_info)))
all(colnames(filtered_counts) == rownames(stage_info))

des_dataset <- DESeqDataSetFromMatrix(countData = filtered_counts,
                       colData = stage_info,
                       design = ~ Stage)

des_dataset$Stage <- relevel(des_dataset$Stage, ref = "Stage 1")

  ## Running gene expression analysis
des_dataset <- DESeq(des_dataset)
results_0.1 <- results(des_dataset)

results_0.1
summary(results_0.1)

results_0.05 <- results(des_dataset, alpha = 0.05)
summary(results_0.05)

res_df <- as.data.frame(results_0.1)

# Visualization

  # MA plot and Dispersion Estimate plot
png(filename = "ma_plot.png", width = 800, height = 600, res = 100)
plotMA(results_0.05, ylim = c(-7, 7))
dev.off()

png(filename = "disp_est_plot.png", width = 800, height = 600, res = 100)
plotDispEsts(des_dataset)
dev.off()

  # Volcano plot
res_df$logP <- -log10(res_df$pvalue)
res_df$significance <- ifelse(res_df$padj < 0.05, "Significant", "Not Significant")
p <- ggplot(res_df, aes(x = log2FoldChange, y = logP, color = significance)) +
  geom_point(alpha = 0.4, size = 1.75) +
  theme_minimal() +
  xlim(c(-8, 8)) +
  ylim(c(0, 30)) +
  xlab("Log2 Fold Change") +
  ylab("-Log10 P-value") +
  ggtitle("Volcano Plot") +
  scale_color_manual(values = c("Not Significant" = "gray", "Significant" = "red")) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.title = element_blank()
  )
ggsave(filename = "volcano_plot.png", plot = p, width = 8, height = 6, dpi = 100)

# Pathway Analysis

gene_list <- results_0.1$log2FoldChange
names(gene_list) <- rownames(results_0.1)
gene_list <- sort(gene_list, decreasing = TRUE)

  # GO enrichment analysis
ego <- enrichGO(gene = names(gene_list),
                OrgDb = org.Hs.eg.db,
                keyType = "ENSEMBL",
                ont = "BP",
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05)

  # GO plot
png(filename = "go_barplot.png", width = 800, height = 600, res = 100)
barplot(ego, showCategory = 10)
dev.off()
