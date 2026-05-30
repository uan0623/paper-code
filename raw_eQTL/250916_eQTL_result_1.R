# 目的
# 生成 chr1-22 的 eQTL result raw_maf_gt_N_pvalue.txt, differentially express 的 raw_exp_different.txt


## package ----
rm(list = ls())
gc()

library(data.table)
library(dplyr)
library(tidyverse)
library(stringr) 
library(ggplot2)
library(ggrepel)
library(matrixStats)  # rowMins
library(qvalue)


# exp split N, T ----
exp <- fread("D:/oral_cancer/expression/expression_data/ExpRes11sv20120601Gene.txt", header = T)

a <- exp[, c(2, 43:82)]
names(a) <-  c("PROBE_ID", sprintf("0%dB", 1:9), sprintf("%dB", 10:40))
fwrite(a, "C:/Peter/rawData_eQTL/raw_Exp_mulInterval_N.txt",
  row.names = F, col.names = T, sep = "\t"
)

a <- exp[, c(2,3:42)] 
names(a) <-  c("PROBE_ID", sprintf("0%dB", 1:9), sprintf("%dB", 10:40))
fwrite(a, "C:/Peter/rawData_eQTL/raw_Exp_mulInterval_T.txt",
  row.names = F, col.names = T, sep = "\t"
)


# differentially express ----
exp <- fread("D:/oral_cancer/expression/expression_data/ExpRes11sv20120601Gene.txt", header = T)
exp_diff <- exp[, 43:82] - exp[, 3:42] 
names(exp_diff) <- c(sprintf("0%dB", 1:9), sprintf("%dB", 10:40))

n <- ncol(exp_diff)
dbar <- rowMeans(exp_diff)
sd_d <- apply(exp_diff, 1, sd)
tval <- dbar / (sd_d / sqrt(n))
pval <- 2 * pt(-abs(tval), df = n - 1)

final <- cbind(exp[, 1:2], pval, tval)
names(final) <- c("Gene", "PROBE_ID", "pval_raw", "t_raw")
probe <- fread("D:/oral_cancer/expression/expression_data/probe24526infoB.txt",
  header = T
)[, .(chr = CHROMOSOME, PROBE_ID, PROBE_COORDINATES)]

final <- probe[final, on=.(PROBE_ID)]
final <- final[chr %in% (c(1:22) %>% as.character())]


final[, pval_raw_BH := p.adjust(pval_raw, method = "BH") %>%
  format(digits = 4, scientific = T) %>%
  as.numeric()]
final[, sig_raw := ifelse(pval_raw_BH < 0.05, 1, 0)]
final[, sig_raw_Bonfi := ifelse(pval_raw < 0.05 / nrow(final), 1, 0)]

# 調整數值
final[, pval_raw := format(pval_raw, digits = 4, scientific = T)]
# differentially express pval 排序
final[, rank_pval := frank(pval_raw, ties.method = "dense")]
fwrite(final, "D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/exp_different.txt",
  row.names = F, col.names = T, sep = "\t"
)







# 建立資料夾 ----
output_dir <- "C:/Peter/rawData_eQTL"
r2_filters <- c("0.6", "0.7", "0.8", "0.9", "no")

dirs_to_create <- c(
  file.path(output_dir, paste0("r2_filter_", r2_filters), "trash", "no_intersect_LD", "EUR"),
  file.path(output_dir, paste0("r2_filter_", r2_filters), "trash", "no_intersect_LD", "FIN"),
  file.path(output_dir, paste0("r2_filter_", r2_filters), "outcome"),
  file.path(output_dir, "outcome"),
  file.path(output_dir, "trash", "Normal"),
  file.path(output_dir, "trash", "Tumor")
)

for (dir_path in dirs_to_create) {
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
}


# 複製 "D:/oral_cancer/expression/trash/maf_gt_N_pvalue.txt" 檔案到 "C:/Peter/rawData_eQTL/trash/raw_maf_gt_N_pvalue.txt"
# 複製 "D:/oral_cancer/expression/trash/maf_gt_T_pvalue.txt" 檔案到 "C:/Peter/rawData_eQTL/trash/raw_maf_gt_T_pvalue.txt"
# 複製 "D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/exp_different.txt" 到 "C:/Peter/rawData_eQTL/outcome/exp_different.txt"

copy_list <- data.frame(
  from = c(
    "D:/oral_cancer/expression/trash/maf_gt_N_pvalue.txt",
    "D:/oral_cancer/expression/trash/maf_gt_T_pvalue.txt",
    "D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/exp_different.txt"
  ),
  to = c(
    "C:/Peter/rawData_eQTL/trash/raw_maf_gt_N_pvalue.txt",
    "C:/Peter/rawData_eQTL/trash/raw_maf_gt_T_pvalue.txt",
    "C:/Peter/rawData_eQTL/outcome/raw_exp_different.txt"
  )
)

file.copy(copy_list$from, copy_list$to, overwrite = TRUE)


# "ILMN_1343291_1" -> "ILMN_1343291"
a <- fread("C:/Peter/rawData_eQTL/trash/raw_maf_gt_N_pvalue.txt")
a[, gene := gsub("^([^_]+_[^_]+).*", "\\1", gene)]
a <- unique(a)
fwrite(a, "C:/Peter/rawData_eQTL/trash/raw_maf_gt_N_pvalue.txt",
  row.names = F, col.names = T, sep = "\t"
)


a <- fread("D:/oral_cancer/expression/trash/maf_gt_T_pvalue.txt")
a[, gene := gsub("^([^_]+_[^_]+).*", "\\1", gene)]
a <- unique(a)
fwrite(a, "C:/Peter/rawData_eQTL/trash/raw_maf_gt_T_pvalue.txt",
  row.names = F, col.names = T, sep = "\t"
)





# each probe select most sig snp, 生成 raw_maf_gt_N_MOSTpvalue.txt ----
gt <- fread("C:/Peter/rawData_eQTL/trash/raw_maf_gt_N_pvalue.txt")
gene_pos <- fread("D:/oral_cancer/expression/expression_data/probe_pos_mulInterval.txt", header = T)


# 按基因、pvalue 排序
setkey(gt, gene, `p-value`)
gt <- merge(gt, gene_pos, by.x = "gene", by.y = "Gene")
gt[, snp_pos := str_extract(SNP, "(?<=\\:)\\d+") %>% as.numeric()]

# 算距離
gt[, dist := pmin(abs(snp_pos - start), abs(snp_pos - end))]
setkey(gt, gene, `p-value`, dist)
gt <- gt[, .SD[1L], by = "gene"]

# 對於 gene ILMN_2353642_1, ILMN_2353642_2，只留下一個 snp
gt[, gene2 := gsub("^([^_]+_[^_]+).*", "\\1", gene)]
setkey(gt, gene2, `p-value`, dist)
gt <- gt[, .SD[1L], by = "gene2"]
gt[, c("CHROMOSOME", "start", "end", "snp_pos", "dist", "gene2") := NULL]

# 這挑出的最近 snp，可能是錯的，因為沒考慮落在區間內的，也許有 snp 不在 probe 範圍內，卻是最近的
fwrite(gt, "C:/Peter/rawData_eQTL/trash/raw_maf_gt_N_MOSTpvalue.txt",
  row.names = F, col.names = T, sep = "\t"
)


gt <- fread("C:/Peter/rawData_eQTL/trash/raw_maf_gt_T_pvalue.txt")
# 同樣pvalue 大小，篩選距離最近的
gene_pos <- fread("D:/oral_cancer/expression/expression_data/probe_pos_mulInterval.txt", header = T)

# 按基因、pvalue 排序
setkey(gt, gene, `p-value`)
gt <- merge(gt, gene_pos, by.x = "gene", by.y = "Gene")
gt[, snp_pos := str_extract(SNP, "(?<=\\:)\\d+") %>% as.numeric()]

# 算距離
gt[, dist := pmin(abs(snp_pos - start), abs(snp_pos - end))]
setkey(gt, gene, `p-value`, dist)
gt <- gt[, .SD[1L], by = "gene"]

# 對於 gene ILMN_2353642_1, ILMN_2353642_2，只留下一個 snp
gt[, gene2 := gsub("^([^_]+_[^_]+).*", "\\1", gene)]
setkey(gt, gene2, `p-value`, dist)
gt <- gt[, .SD[1L], by = "gene2"]
gt[, c("CHROMOSOME", "start", "end", "snp_pos", "dist", "gene2") := NULL]

# 這挑出的最近 snp，可能是錯的，因為沒考慮落在區間內的，也許有 snp 不在 probe 範圍內，卻是最近的
fwrite(gt, "C:/Peter/rawData_eQTL/trash/raw_maf_gt_T_MOSTpvalue.txt",
  row.names = F, col.names = T, sep = "\t"
)



