if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# BiocManager::install(c("TCGAbiolinks", "SummarizedExperiment"))

library(TCGAbiolinks)
library(SummarizedExperiment)
library(data.table)

query <- GDCquery(
  project = "TCGA-HNSC",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

GDCdownload(query, method = "api")

hnsc_se <- GDCprepare(query)


gene_info <- as.data.table(rowData(hnsc_se))

names(gene_info)
head(gene_info)
dim(gene_info)


fwrite(
  gene_info,
  "TCGA_HNSC_STAR_Counts_all_genes_rowData.tsv",
  sep = "\t"
)


a <- fread( "C:/Peter/gene_enrichment/code_project/data/HNSC/TCGA-HNSC.TUMOR.pos")
a <- fread( ""//wsl.localhost/Ubuntu-22.04/home/jcc623/fusion_project/TCGA-HNSC.TUMOR/fusion_weights_summary.csv"")
dim(a)
