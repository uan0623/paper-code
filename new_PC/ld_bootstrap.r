library(data.table)
library(dplyr)
library(tidyverse)
library(stringr)
library(ggplot2)
library(ggrepel)
library(matrixStats) # rowMins
library(qvalue)
library(bestNormalize)
library(ggplot2)
library(openxlsx)
library(magrittr)



gt <- fread("D:/Peter/rawData_eQTL/trash/raw_maf_gt_N_pvalue.txt")
gene_pos <- fread("D:/Peter/oral_cancer/expression/expression_data/probe24526infoB.txt", header = T)


# 按基因、pvalue 排序
setkey(gt, gene, `p-value`)
gt <- gt[, .SD[`p-value` == min(`p-value`, na.rm = TRUE)], by = "gene"]
gt_lots_asso <- gt[, if (.N > 1) .SD, by = gene]
uniqueN(gt$gene)
lead_snp <- gt[, if (.N == 1) .SD, by = gene] %>%
  unique(., by = "SNP") %>%
  select(SNP)



lead_vec <- unique(lead_snp$SNP)
gt[, in_lead := fifelse(SNP %in% lead_vec, 1L, 0L)]

# 把 in_lead=1 的，也就是不只出現一次的 probe，snp 與只出現一次的 probe 的 snp 重疊，優先選
gt_one <- gt[
  order(gene, -in_lead, `p-value`)
][
  , .SD[1L],
  by = gene
][
  , in_lead := NULL
]


uniqueN(gt_one$gene)
uniqueN(gt_one$SNP)
gt_one[, .N, by = SNP]$N %>% table()


# liftOver 轉換 pos

gene_pos <- fread("D:/Peter/oral_cancer/expression/expression_data/probe_pos_mulInterval.txt", header = T)
gene_probe <- fread("D:/Peter/oral_cancer/expression/expression_data/probe24526infoB.txt", header = T)
gene_pos <- merge(gene_pos[, .(PROBE_ID = Gene, chr = CHROMOSOME, start, end)], gene_probe[, .(PROBE_ID, TargetID)])

fwrite(
  data.table(
    a = paste0("chr", gene_pos$chr),
    b = gene_pos$start,
    c = gene_pos$end,
    d =  gene_pos$PROBE_ID
  ),
  "D:/Peter/gene_enrichment/code_project/data/ld_bootstrap/gene_hg18.txt",
  row.names = F, col.names = F, sep = "\t"
)


cmd <- paste(
  "/mnt/d/Peter/oral_cancer/liftover/liftOver",
  "/mnt/d/Peter/gene_enrichment/code_project/data/ld_bootstrap/gene_hg18.txt",
  "/mnt/d/Peter/oral_cancer/liftover/hg18ToHg19.over.chain",
  "/mnt/d/Peter/gene_enrichment/code_project/data/ld_bootstrap/gene_hg18TOhg19.bed",
  "/mnt/d/Peter/gene_enrichment/code_project/data/ld_bootstrap/hg18TOhg19_unmapped.bed"
)

system2("wsl", cmd)

a <- fread("D:/Peter/gene_enrichment/code_project/data/ld_bootstrap/gene_hg18TOhg19.bed", header = F)
# 忽略轉失敗 .bed 檔案的 #Deleted in new 該行

unmapped_path <- "D:/Peter/gene_enrichment/code_project/data/ld_bootstrap/hg18TOhg19_unmapped.bed"

if (!file.exists(unmapped_path) || file.info(unmapped_path)$size == 0) {
  a_unmapped <- data.table(V1 = character(), V2 = integer(), V3 = integer())
} else {
  a_unmapped <- fread(
    unmapped_path,
    header = FALSE,
    comment.char = "#"
  )
}


col_order <- names(gene_pos)
names(a_unmapped) <- c("chr", "start", "end")
names(a) <- c("chr_hg19", "start_hg19", "end_hg19")
# aa
a_unmapped[, chr := sub("^chr", "", chr) %>% as.numeric()]
a[, chr_hg19 := sub("^chr", "", chr_hg19) %>% as.numeric()]



# 新增轉換後 pos
a_unmapped[, c("chr_hg19", "start_hg19", "end_hg19") := 0]
final <- merge(gene_pos, a_unmapped, all.x = T)

success_idx <- is.na(final$chr_hg19)
final[
  success_idx,
  c("chr_hg19", "start_hg19", "end_hg19") := a[, .(chr_hg19, start_hg19, end_hg19)]
]

# 轉換紀錄
# fwrite(final,
#   "D:/Peter/gene_enrichment/code_project/data/hg38_hg18.txt",
#   row.names = F, col.names = T, sep = "\t"
# )


# 下載 grch37 annotation ----
library(data.table)

gtf <- fread(
  "D:/Peter/gene_enrichment/code_project/data/ld_bootstrap/gencode.v37lift37.annotation.gtf",
  sep = "\t",
  header = FALSE,
  quote = "",
  data.table = FALSE
)

colnames(gtf) <- c(
  "chr", "source", "feature", "start", "end",
  "score", "strand", "frame", "attribute"
)

gene_gtf <- gtf[gtf$feature == "gene", ]
extract_attr <- function(x, key) {
  sub(
    paste0(".*", key, " \"([^\"]+)\".*"),
    "\\1",
    x
  )
}

gene_gtf$gene_id <- extract_attr(gene_gtf$attribute, "gene_id")
gene_gtf$gene_name <- extract_attr(gene_gtf$attribute, "gene_name")
gene_gtf$gene_type <- extract_attr(gene_gtf$attribute, "gene_type")
gene_gtf$ensembl_id <- sub("\\..*", "", gene_gtf$gene_id)

gene_gtf$tss <- ifelse(
  gene_gtf$strand == "-",
  gene_gtf$end,
  gene_gtf$start
)

gene_gtf$chr <- gsub("^chr", "", gene_gtf$chr)
gene_gtf <- gene_gtf[
  gene_gtf$chr %in% as.character(1:22),
]

gene_gtf <- as.data.table(gene_gtf)




# under 0.8 threshold, 有用 cis snp 的 probe ----
all_eQTL <- fread("D:/Peter/rawData_eQTL/trash/raw_maf_gt_N_pvalue.txt", header = T)[
  ,
  .(pval = `p-value`, SNP, probe = gene)
]
R2_filter <- fread("D:/Peter/oral_cancer/expression/outcome/multiple_nucleotide_variant/cis_snp_R2_0.8.txt", header = T)

setkey(all_eQTL, pval)
all_eQTL <- all_eQTL[SNP %in% R2_filter$hg18_snpID, ]
all_eQTL <- all_eQTL[, .SD[1], by = probe]
all_eQTL[, probe := gsub("^([^_]+_[^_]+).*", "\\1", probe)]

probe_info <- fread("D:/Peter/rawData_eQTL/r2_filter_0.8/outcome/raw_exp_different_r2_0.8.txt")[
  ,
  .(probe = PROBE_ID, Gene)
]
all_eQTL <- probe_info[all_eQTL, on = .(probe), nomatch = 0]
all_eQTL_gene <- unique(all_eQTL$Gene)
all_eQTL_probe <- unique(all_eQTL$probe)


# 挑出我們資料有 cis snp 的 gene
gene_gtf <- gene_gtf[gene_name %in% all_eQTL_gene]



#  hg19 (same as 1000 Genomes phase 1)
ld_blocks <- fread(
  "D:/Peter/gene_enrichment/code_project/data/ld_bootstrap/ldetect-data/ASN/fourier_ls-all.bed"
)
colnames(ld_blocks)[1:3] <- c("chr", "start", "end")
setDT(ld_blocks)

ld_blocks[, chr := gsub("^chr", "", chr)]

ld_blocks[, LD_block := paste0(
  "ASN_chr", chr, "_",
  start, "_",
  end
)]

gene_gtf[, LD_block := NA_character_]

gene_gtf[
  ld_blocks,
  on = .(
    chr,
    tss >= start,
    tss <= end
  ),
  LD_block := i.LD_block
]

png("D:/Peter/gene_enrichment/code_project/outcome/ld_bootstrap/LD_block_box.png",
  width = 1000, height = 1000, res = 200
)
gene_gtf[, .N, by = LD_block]$N %>% boxplot(main = "Gene number in each LD_block")
dev.off()

fwrite(ld_blocks,
  "D:/Peter/gene_enrichment/code_project/outcome/ld_bootstrap/ld_block.txt",
  row.names = F, col.names = T, sep = "\t"
)
fwrite(gene_gtf,
  "D:/Peter/gene_enrichment/code_project/outcome/ld_bootstrap/OC_gene_blocks.txt",
  row.names = F, col.names = T, sep = "\t"
)
