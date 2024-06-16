### Survival analysis of TCGA-BRCA patients based on expression levels of CT83, the most significant ###
### gene from differential gene expression analysis                                                  ###

## Loading packages ##

library("TCGAbiolinks")
library("tidyverse")
library("survminer")
library("survival")
library("SummarizedExperiment")
library("DESeq2")

## Querying clinical data ##

clinical_data <- GDCquery_clinic("TCGA-BRCA")

clinical_data$deceased <- ifelse(clinical_data$vital_status == "Alive", FALSE, TRUE)

clinical_data$overall_survival <- ifelse(clinical_data$vital_status == "Alive",
                                         clinical_data$days_to_last_follow_up,
                                         clinical_data$days_to_death)

## Querying gene expression data ##

query <- GDCquery(project = "TCGA-BRCA",
                         data.category = "Transcriptome Profiling",
                         experimental.strategy = "RNA-Seq",
                         workflow.type = "STAR - Counts",
                         data.type = "Gene Expression Quantification", 
                         sample.type = c("Primary Tumor"),
                         access = "open")

output <- getResults(query)
primaryTumor <- output[output$sample_type == 'Primary Tumor', "cases"]

query <- GDCquery(project = "TCGA-BRCA",
                  data.category = "Transcriptome Profiling",
                  experimental.strategy = "RNA-Seq",
                  workflow.type = "STAR - Counts",
                  data.type = "Gene Expression Quantification", 
                  sample.type = c("Primary Tumor"),
                  access = "open",
                  barcode = primaryTumor)

GDCdownload(query)

## Pre-processing data and isolating CT83 expression data ##

ge_data = GDCprepare(query, summarizedExperiment = TRUE)

ge_matrix <- assay(ge_data, "unstranded")

gene_metadata <- as.data.frame(rowData(ge_data))
colData <- as.data.frame(colData(ge_data))

dds <- DESeqDataSetFromMatrix(countData = ge_matrix,
                              colData = colData,
                              design = ~ 1)

dds <- dds[rowSums(counts(dds)) >= 10]

matrix_vst <- assay(vst(dds, blind = FALSE))

ct83_data <- matrix_vst %>% 
  as.data.frame() %>%
  rownames_to_column(var = 'gene_id') %>%
  gather(key = 'case_id', value = 'counts', -gene_id) %>%
  left_join(., gene_metadata, by = 'gene_id') %>%
  filter(gene_name == 'CT83')

ct83_data$strata <- ifelse(ct83_data$counts > median(ct83_data$counts), "HIGH", "LOW")

ct83_data$case_id <- gsub('-01.*', '', ct83_data$case_id)

ct83_data <- merge(ct83_data, clinical_data, by.x = 'case_id', by.y = 'submitter_id')

## Survival analysis ##
fit <- survfit(Surv(overall_survival, deceased) ~ strata, data = ct83_data)

png(filename = "survplot.png", width = 800, height = 600, res = 100)
ggsurvplot(fit,
           data = ct83_data,
           pval = T,
           risk.table = T)
dev.off()
