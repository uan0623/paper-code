# 釋放記憶體
rm(list = ls())
gc()

library(data.table)
library(dplyr)
library(tidyverse)
library(stringr)
library(ggplot2)
library(ggrepel)
library(matrixStats) # rowMins
library(qvalue)

plot_density <- function(dt, plot_title) {
  if (!is.character(plot_title) || length(plot_title) != 1 || !is.data.table(dt)) {
    stop("plot_title 必須是一個字串, my_dt 必須是一個 data.table")
  }

  ds <- dt

  bin_width <- 1 / 100
  h <- hist(ds$`p-value`, breaks = seq(0, 1, by = bin_width), plot = FALSE)

  # 轉成 bin_data畫圖，避免記憶體爆
  bin_data <- data.table(
    mid = h$mids, # 每個 bin 的中點
    N = h$counts
  )

  # step3: 計算密度 (N / (總樣本數 * bin 寬度))
  total_n <- sum(bin_data$N)
  bin_data[, density := N / (total_n * bin_width)]

  # step4: 繪製密度直方圖
  a <- ggplot(bin_data, aes(x = mid, y = density)) +
    geom_col(width = bin_width, color = "black", fill = "skyblue") +
    geom_hline(yintercept = 1, linetype = "solid", color = "red") + # 水平線
    labs(
      x = "P-value", y = "Density",
      title = plot_title
    ) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))

  return(a)
}


# DS ####
pvalue_0.05_N <- 0
pvalue_0.05_T <- 0


# Normal part####

file_paths <- sprintf("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/Normal/ds_N_cis_chr%d_snpNotRepeat.txt", 1:22)

ds <- lapply(file_paths, fread, header = T) %>%
  rbindlist(use.names = TRUE)

pvalue_0.05_N <- pvalue_0.05_N + ds %>%
  filter(`p-value` < 0.05)          %>%
  nrow()


ggsave("D:/oral_cancer/expression/trash/ds_N_pvalue_Autosome.png",
  plot = plot_density(ds, "Normal Part DS: P-value of Cis-SNPs in Autosome"),
  width = 10, height = 5, dpi = 300
)

# ds_N_pvalue.txt 會53GB
# fwrite(ds, "D:/oral_cancer/expression/trash/ds_N_pvalue.txt",
#        row.names = F, col.names = T, sep = "\t")


# 同樣pvalue 大小，篩選距離最近的 ####
gene_pos <- fread("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/probe_pos.txt", header = T)

# 按基因、pvalue 排序
setkey(ds, gene, `p-value`)
ds <- merge(ds, gene_pos, by.x = "gene", by.y = "Gene")
ds[, snp_pos := str_extract(SNP, "(?<=\\:)\\d+") %>% as.numeric()]

# 算距離
ds[, dist := pmin(abs(snp_pos - start), abs(snp_pos - end))]
setkey(ds, gene, `p-value`, dist)
ds <- ds[, .SD[1L], by = "gene"]

# 對於 gene ILMN_2353642_1, ILMN_2353642_2，只留下一個 snp
ds[, gene2 := gsub("^([^_]+_[^_]+).*", "\\1", gene)]
setkey(ds, gene2, `p-value`, dist)
ds <- ds[, .SD[1L], by = "gene2"]



ggsave("D:/oral_cancer/expression/trash/ds_N_MOSTpvalue_Autosome.png",
  plot = plot_density(ds, "Normal Part DS: Each Probe Choose the Most Significant Cis-eQTL in Autosome"),
  width = 10, height = 5, dpi = 300
)

fwrite(ds, "D:/oral_cancer/expression/trash/ds_N_MOSTpvalue.txt",
  row.names = F, col.names = T, sep = "\t"
)

rm(ds)
gc()



# Tumor part####

file_paths <- sprintf("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/Tumor/ds_T_cis_chr%d_snpNotRepeat.txt", 1:22)

ds <- lapply(file_paths, fread, header = T) %>%
  rbindlist(use.names = TRUE)


pvalue_0.05_T <- pvalue_0.05_T + ds %>%
  filter(`p-value` < 0.05)          %>%
  nrow()


ggsave("D:/oral_cancer/expression/trash/ds_T_pvalue_Autosome.png",
  plot = plot_density(ds, "Tumor Part DS: P-value of Cis-SNPs in Autosome"),
  width = 10, height = 5, dpi = 300
)

# ds_T_pvalue.txt 會53GB
# fwrite(ds, "D:/oral_cancer/expression/trash/ds_T_pvalue.txt",
#        row.names = F, col.names = T, sep = "\t")


# 同樣pvalue 大小，篩選距離最近的 ####
gene_pos <- fread("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/probe_pos.txt", header = T)

# 按基因、pvalue 排序
setkey(ds, gene, `p-value`)
ds <- merge(ds, gene_pos, by.x = "gene", by.y = "Gene")
ds[, snp_pos := str_extract(SNP, "(?<=\\:)\\d+") %>% as.numeric()]

# 算距離
ds[, dist := pmin(abs(snp_pos - start), abs(snp_pos - end))]
setkey(ds, gene, `p-value`, dist)
ds <- ds[, .SD[1L], by = "gene"]

# 對於 gene ILMN_2353642_1, ILMN_2353642_2，只留下一個 snp
ds[, gene2 := gsub("^([^_]+_[^_]+).*", "\\1", gene)]
setkey(ds, gene2, `p-value`, dist)
ds <- ds[, .SD[1L], by = "gene2"]




ggsave("D:/oral_cancer/expression/trash/ds_T_MOSTpvalue_Autosome.png",
  plot = plot_density(ds, "Tumor part DS: Each Probe Choose the Most Significant Cis-eQTL in Autosome"),
  width = 10, height = 5, dpi = 300
)

fwrite(ds, "D:/oral_cancer/expression/trash/ds_T_MOSTpvalue.txt",
  row.names = F, col.names = T, sep = "\t"
)


# 輸出 pvalue < 0.05 association 個數
cat("Normal part DS autosome, pvalue < 0.05 association number: ", pvalue_0.05_N)
cat("Tumor part DS autosome, pvalue < 0.05 association number: ", pvalue_0.05_T)



rm(ds)
gc()


# GT ####
pvalue_0.05_N <- 0
pvalue_0.05_T <- 0



# Normal part####
file_paths <- sprintf("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/Normal/gt_N_cis_chr%d_snpNotRepeat.txt", 1:22)

gt <- lapply(file_paths, fread, header = T) %>%
  rbindlist(use.names = TRUE)

pvalue_0.05_N <- pvalue_0.05_N + gt %>%
  filter(`p-value` < 0.05)          %>%
  nrow()

ggsave("D:/oral_cancer/expression/trash/gt_N_pvalue_Autosome.png",
  plot = plot_density(gt, "Normal Part GT: P-value of Cis-SNPs in Autosome"),
  width = 10, height = 5, dpi = 300
)

fwrite(gt, "D:/oral_cancer/expression/trash/gt_N_pvalue.txt",
  row.names = F, col.names = T, sep = "\t"
)


# 同樣pvalue 大小，篩選距離最近的 ####
gene_pos <- fread("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/probe_pos.txt", header = T)

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




ggsave("D:/oral_cancer/expression/trash/gt_N_MOSTpvalue_Autosome.png",
  plot = plot_density(gt, "Normal Part GT: Each Probe Choose the Most Significant Cis-eQTL in Autosome"),
  width = 10, height = 5, dpi = 300
)

fwrite(gt, "D:/oral_cancer/expression/trash/gt_N_MOSTpvalue.txt",
  row.names = F, col.names = T, sep = "\t"
)


rm(gt)
gc()

# Tumor part####

file_paths <- sprintf("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/Tumor/gt_T_cis_chr%d_snpNotRepeat.txt", 1:22)

gt <- lapply(file_paths, fread, header = T) %>%
  rbindlist(use.names = TRUE)

pvalue_0.05_T <- pvalue_0.05_T + gt %>%
  filter(`p-value` < 0.05)          %>%
  nrow()

ggsave("D:/oral_cancer/expression/trash/gt_T_pvalue_Autosome.png",
  plot = plot_density(gt, "Tumor Part GT: P-value of Cis-SNPs in Autosome"),
  width = 10, height = 5, dpi = 300
)

fwrite(gt, "D:/oral_cancer/expression/trash/gt_T_pvalue.txt",
  row.names = F, col.names = T, sep = "\t"
)


# 同樣pvalue 大小，篩選距離最近的 ####
gene_pos <- fread("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/probe_pos.txt", header = T)

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



ggsave("D:/oral_cancer/expression/trash/gt_T_MOSTpvalue_Autosome.png",
  plot = plot_density(gt, "Tumor Part GT: Each Probe Choose the Most Significant Cis-eQTL in Autosome"),
  width = 10, height = 5, dpi = 300
)

fwrite(gt, "D:/oral_cancer/expression/trash/gt_T_MOSTpvalue.txt",
  row.names = F, col.names = T, sep = "\t"
)


# 輸出 pvalue < 0.05 association 個數
cat("Normal part GT autosome, pvalue < 0.05 association number: ", pvalue_0.05_N)
cat("Tumor part GT autosome, pvalue < 0.05 association number: ", pvalue_0.05_T)


rm(gt)
gc()

# GT ####
pvalue_0.05_N <- 0
pvalue_0.05_T <- 0


# Normal part####
file_paths <- sprintf("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/Normal/maf_gt_N_cis_chr%d_snpNotRepeat.txt", 1:22)

gt <- lapply(file_paths, fread, header = T) %>%
  rbindlist(use.names = TRUE)

pvalue_0.05_N <- pvalue_0.05_N + gt %>%
  filter(`p-value` < 0.05)          %>%
  nrow()


ggsave("D:/oral_cancer/expression/trash/maf_gt_N_pvalue_Autosome.png",
  plot = plot_density(gt, "Normal Part GT & 0.05<MAF: P-value of Cis-SNPs in Autosome"),
  width = 10, height = 5, dpi = 300
)

fwrite(gt, "D:/oral_cancer/expression/trash/maf_gt_N_pvalue.txt",
  row.names = F, col.names = T, sep = "\t"
)


# 同樣pvalue 大小，篩選距離最近的 ####
gene_pos <- fread("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/probe_pos.txt", header = T)

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


ggsave("D:/oral_cancer/expression/trash/maf_gt_N_MOSTpvalue_Autosome.png",
  plot = plot_density(gt, "Normal Part GT & 0.05<MAF: Each Probe Choose the Most Significant Cis-eQTL in Autosome"),
  width = 10, height = 5, dpi = 300
)

fwrite(gt, "D:/oral_cancer/expression/trash/maf_gt_N_MOSTpvalue.txt",
  row.names = F, col.names = T, sep = "\t"
)


# Tumor part####

file_paths <- sprintf("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/Tumor/maf_gt_T_cis_chr%d_snpNotRepeat.txt", 1:22)

gt <- lapply(file_paths, fread, header = T) %>%
  rbindlist(use.names = TRUE)

pvalue_0.05_T <- pvalue_0.05_T + gt %>%
  filter(`p-value` < 0.05)          %>%
  nrow()


ggsave("D:/oral_cancer/expression/trash/maf_gt_T_pvalue_Autosome.png",
  plot = plot_density(gt, "Tumor Part GT & 0.05<MAF: P-value of Cis-SNPs in Autosome"),
  width = 10, height = 5, dpi = 300
)

fwrite(gt, "D:/oral_cancer/expression/trash/maf_gt_T_pvalue.txt",
  row.names = F, col.names = T, sep = "\t"
)


# 同樣pvalue 大小，篩選距離最近的 ####
gene_pos <- fread("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/probe_pos.txt", header = T)

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

ggsave("D:/oral_cancer/expression/trash/maf_gt_T_MOSTpvalue_Autosome.png",
  plot = plot_density(gt, "Tumor Part GT & 0.05<MAF: Each Probe Choose the Most Significant Cis-eQTL in Autosome"),
  width = 10, height = 5, dpi = 300
)

fwrite(gt, "D:/oral_cancer/expression/trash/maf_gt_T_MOSTpvalue.txt",
  row.names = F, col.names = T, sep = "\t"
)


# 輸出 pvalue < 0.05 association 個數
cat("Normal part GT autosome & 0.05<MAF, pvalue < 0.05 association number: ", pvalue_0.05_N)
cat("Tumor part GT autosome & 0.05<MAF, pvalue < 0.05 association number: ", pvalue_0.05_T)


rm(gt)
gc()


# DS ####
pvalue_0.05_N <- 0
pvalue_0.05_T <- 0


# Normal part####
file_paths <- sprintf("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/Normal/maf_ds_N_cis_chr%d_snpNotRepeat.txt", 1:22)

ds <- lapply(file_paths, fread, header = T) %>%
  rbindlist(use.names = TRUE)

pvalue_0.05_N <- pvalue_0.05_N + ds %>%
  filter(`p-value` < 0.05)          %>%
  nrow()


ggsave("D:/oral_cancer/expression/trash/maf_ds_N_pvalue_Autosome.png",
  plot = plot_density(ds, "Normal Part DS & 0.05<MAF: P-value of Cis-SNPs in Autosome"),
  width = 10, height = 5, dpi = 300
)

fwrite(ds, "D:/oral_cancer/expression/trash/maf_ds_N_pvalue.txt",
  row.names = F, col.names = T, sep = "\t"
)


# 同樣pvalue 大小，篩選距離最近的 ####
gene_pos <- fread("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/probe_pos.txt", header = T)

# 按基因、pvalue 排序
setkey(ds, gene, `p-value`)
ds <- merge(ds, gene_pos, by.x = "gene", by.y = "Gene")
ds[, snp_pos := str_extract(SNP, "(?<=\\:)\\d+") %>% as.numeric()]

# 算距離
ds[, dist := pmin(abs(snp_pos - start), abs(snp_pos - end))]
setkey(ds, gene, `p-value`, dist)
ds <- ds[, .SD[1L], by = "gene"]

# 對於 gene ILMN_2353642_1, ILMN_2353642_2，只留下一個 snp
ds[, gene2 := gsub("^([^_]+_[^_]+).*", "\\1", gene)]
setkey(ds, gene2, `p-value`, dist)
ds <- ds[, .SD[1L], by = "gene2"]
ds[, c("CHROMOSOME", "start", "end", "snp_pos", "dist", "gene2") := NULL]


ggsave("D:/oral_cancer/expression/trash/maf_ds_N_MOSTpvalue_Autosome.png",
  plot = plot_density(ds, "Normal Part DS & 0.05<MAF: Each Probe Choose the Most Significant Cis-eQTL in Autosome"),
  width = 10, height = 5, dpi = 300
)
fwrite(ds, "D:/oral_cancer/expression/trash/maf_ds_N_MOSTpvalue.txt",
  row.names = F, col.names = T, sep = "\t"
)



rm(ds)





# Tumor part####

file_paths <- sprintf("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/Tumor/maf_ds_T_cis_chr%d_snpNotRepeat.txt", 1:22)

ds <- lapply(file_paths, fread, header = T) %>%
  rbindlist(use.names = TRUE)

pvalue_0.05_T <- pvalue_0.05_T + ds %>%
  filter(`p-value` < 0.05)          %>%
  nrow()


ggsave("D:/oral_cancer/expression/trash/maf_ds_T_pvalue_Autosome.png",
  plot = plot_density(ds, "Tumor Part DS & 0.05<MAF: P-value of Cis-SNPs in Autosome"),
  width = 10, height = 5, dpi = 300
)

fwrite(ds, "D:/oral_cancer/expression/trash/maf_ds_T_pvalue.txt",
  row.names = F, col.names = T, sep = "\t"
)


# 同樣pvalue 大小，篩選距離最近的 ####
gene_pos <- fread("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/probe_pos.txt", header = T)

# 按基因、pvalue 排序
setkey(ds, gene, `p-value`)
ds <- merge(ds, gene_pos, by.x = "gene", by.y = "Gene")
ds[, snp_pos := str_extract(SNP, "(?<=\\:)\\d+") %>% as.numeric()]

# 算距離
ds[, dist := pmin(abs(snp_pos - start), abs(snp_pos - end))]
setkey(ds, gene, `p-value`, dist)
ds <- ds[, .SD[1L], by = "gene"]

# 對於 gene ILMN_2353642_1, ILMN_2353642_2，只留下一個 snp
ds[, gene2 := gsub("^([^_]+_[^_]+).*", "\\1", gene)]
setkey(ds, gene2, `p-value`, dist)
ds <- ds[, .SD[1L], by = "gene2"]
ds[, c("CHROMOSOME", "start", "end", "snp_pos", "dist", "gene2") := NULL]


ggsave("D:/oral_cancer/expression/trash/maf_ds_T_MOSTpvalue_Autosome.png",
  plot = plot_density(ds, "Tumor part DS & 0.05<MAF: Each Probe Choose the Most Significant Cis-eQTL in Autosome"),
  width = 10, height = 5, dpi = 300
)

fwrite(ds, "D:/oral_cancer/expression/trash/maf_ds_T_MOSTpvalue.txt",
  row.names = F, col.names = T, sep = "\t"
)


# 輸出 pvalue < 0.05 association 個數
cat("Normal part DS autosome & 0.05<MAF, pvalue < 0.05 association number: ", pvalue_0.05_N)
cat("Tumor part DS autosome & 0.05<MAF, pvalue < 0.05 association number: ", pvalue_0.05_T)


# 計算並挑出FDR<0.05 的row
compute_FDR <- function(inputname, outputname, qvalue_type) {
  if (!is.character(inputname) || length(inputname) != 1 || !is.character(outputname) || length(outputname) != 1) {
    stop("inputname, outputname 必須是一個字串")
  }

  x <- fread(inputname, header = T)

  x[, sig_pval_Bonfi := ifelse(`p-value` < 0.05 / nrow(x), 1, 0)]
  x[, FDR := p.adjust(`p-value`, method = "BH") %>%
    format(digits = 10, scientific = T)         %>%
    as.numeric()]

  if (qvalue_type == 1) {
    pi0_hat <- min(max(
      sum(x$`p-value` > 0.7) / ((1 - 0.7) * length(x$`p-value`)),
      0
    ), 1)

    x[, qvalue := qvalue(`p-value`, pi0 = pi0_hat)$qvalues]
  }

  if (qvalue_type == 2) {
    x[, qvalue := qvalue(`p-value`)$qvalues]
  }


  x <- x %>% filter(FDR < 0.05)
  x <- x[order(FDR), ]
  x[, gene2 := gsub("^([^_]+_[^_]+).*", "\\1", gene)]


  fwrite(x, outputname, row.names = F, col.names = T, sep = "\t")
}


file_name <- c(
  "maf_ds_N_MOSTpvalue",
  "maf_ds_T_MOSTpvalue",
  "maf_gt_N_MOSTpvalue",
  "maf_gt_T_MOSTpvalue"
)

file_name_all <- c(
  "maf_gt_N_pvalue",
  "maf_gt_T_pvalue"
)


# 對這些檔案生成有FDR 檔
# 如果合併ds_T_pvalue.txt 會53GB
for (i in file_name) {
  compute_FDR(
    paste0("D:/oral_cancer/expression/trash/", i, ".txt"),
    paste0("D:/oral_cancer/expression/trash/", i, "_FDR.txt"), 1
  )
}

for (i in file_name_all) {
  compute_FDR(
    paste0("D:/oral_cancer/expression/trash/", i, ".txt"),
    paste0("D:/oral_cancer/expression/trash/", i, "_FDR.txt"), 2
  )
}


Ordered_normal_quantile <- function(df, by.row = T) {
  if (!is.matrix(df)) {
    stop("df 必須是一個 metrix")
  }

  if (!by.row) {
    # 對每col 做排序####
    i <- 2
    n <- nrow(df)

    # 按照row 個數n，把(0,1) 切n+1 刀，再減掉0.5/n。這樣就能取n等分的值
    avg_sorted <- qnorm((1:n - 0.5) / n, mean = 0, sd = 1)

    # i=2，表示根據col運算，找出df 每個col不同row 的值大小順位，並賦予相應順位的值
    mat_qn <- apply(df, i, function(col) {
      # ties.method= "min" 代表 同值取最小的 rank，兩個expression 同值，normal_quantile 的結果會相同
      ranks <- rank(col, ties.method = "min")
      avg_sorted[ranks]
    })
  } else {
    # 對每row 做排序####
    i <- 1
    n <- ncol(df)
    avg_sorted <- qnorm((1:n - 0.5) / n, mean = 0, sd = 1)

    # 找出df 每個row不同row 的值大小順位，並賦予相應順位的值
    mat_qn <- apply(df, i, function(row) {
      # i=2，表示根據row 運算，ties.method= "min" 代表 同值取最小的 rank，兩個expression 同值，normal_quantile 的結果會相同
      ranks <- rank(row, ties.method = "min")
      avg_sorted[ranks]
    }) %>% t()
  }


  return(mat_qn)
}


normal_quantile <- function(df, by.row = T) {
  if (!is.matrix(df)) {
    stop("df 必須是一個 metrix")
  }

  if (!by.row) {
    # 對每col 做排序####
    i <- 2
    n <- nrow(df)
    avg_sorted <- apply(df, i, sort)
    avg_sorted <- rowMeans(avg_sorted)

    mat_qn <- apply(df, i, function(col) {
      ranks <- rank(col, ties.method = "min")
      avg_sorted[ranks]
    })
  } else {
    # 對每row 做排序####
    i <- 1
    n <- ncol(df)
    avg_sorted <- apply(df, i, sort)
    avg_sorted <- colMeans(avg_sorted)

    mat_qn <- apply(df, i, function(row) {
      ranks <- rank(row, ties.method = "min")
      avg_sorted[ranks]
    }) %>% t()
  }

  return(mat_qn)
}

# 不能用probe 重複出現的資料做 QN，probe值相同會變不同
exp <- fread("D:/oral_cancer/expression/expression_data/ExpRes11sv20120601Gene.txt", header = T)



# Ordered_normal_quantile ####
exp_diff_afQN <- Ordered_normal_quantile(as.matrix(exp[, 43:82]), by.row = F)   %>%
  as.data.table() - Ordered_normal_quantile(as.matrix(exp[, 3:42]), by.row = F) %>%
  as.data.table()


n <- ncol(exp_diff_afQN)
dbar <- rowMeans(exp_diff_afQN)

sd_d <- apply(exp_diff_afQN, 1, sd)
tval <- dbar / (sd_d / sqrt(n))
pval <- 2 * pt(-abs(tval), df = n - 1)

pvalue_oqn <- cbind(exp[, 1:2], pval, tval)
names(pvalue_oqn) <- c("Gene", "PROBE_ID", "pval_OQN", "t_OQN")


# normal_quantile ####
exp_diff_afQN <- normal_quantile(as.matrix(exp[, 43:82]), by.row = F)   %>%
  as.data.table() - normal_quantile(as.matrix(exp[, 3:42]), by.row = F) %>%
  as.data.table()


n <- ncol(exp_diff_afQN)
dbar <- rowMeans(exp_diff_afQN)

sd_d <- apply(exp_diff_afQN, 1, sd)
tval <- dbar / (sd_d / sqrt(n))
pval <- 2 * pt(-abs(tval), df = n - 1)

pvalue_oqn <- cbind(pvalue_oqn, pval, tval)
names(pvalue_oqn) <- c("Gene", "PROBE_ID", "pval_OQN", "t_OQN", "pval_QN", "t_QN")


# 計算
pvalue_oqn[, pval_OQN_BH := p.adjust(pval_OQN, method = "BH") %>%
  format(digits = 4, scientific = T)                          %>%
  as.numeric()]

pvalue_oqn[, pval_QN_BH := p.adjust(pval_QN, method = "BH") %>%
  format(digits = 4, scientific = T)                        %>%
  as.numeric()]


# 判斷pvalue 是否小於0.05
pvalue_oqn[, sig_OQN := ifelse(pval_OQN_BH < 0.05, 1, 0)]
pvalue_oqn[, sig_QN := ifelse(pval_QN_BH < 0.05, 1, 0)]
pvalue_oqn[, sig_OQN_Bonfi := ifelse(pval_OQN < 0.05 / nrow(pvalue_oqn), 1, 0)]
pvalue_oqn[, sig_QN_Bonfi := ifelse(pval_QN < 0.05 / nrow(pvalue_oqn), 1, 0)]


fwrite(pvalue_oqn, "D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/exp_different.txt",
  row.names = F, col.names = T, sep = "\t"
)


df <- fread("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/exp_different.txt", header = T)

file_name <- c(
  "maf_gt_N_pvalue_FDR",
  "maf_gt_T_pvalue_FDR",
  "maf_gt_N_MOSTpvalue_FDR",
  "maf_gt_T_MOSTpvalue_FDR",
  "maf_ds_N_pvalue_FDR",
  "maf_ds_T_pvalue_FDR",
  "maf_ds_N_MOSTpvalue_FDR",
  "maf_ds_T_MOSTpvalue_FDR"
)


# 幫 exp_different.txt 增加各檔案 sigCis-SNP_number
for (i in file_name) {
  a <- fread(paste0("D:/oral_cancer/expression/trash/", i, ".txt"), header = T)
  a[, gene2 := gsub("^([^_]+_[^_]+).*", "\\1", gene)]

  sig <- a[, .N, by = gene2]
  df <- merge(df, sig, by.x = "PROBE_ID", by.y = "gene2", all.x = T)
  df[is.na(N), N := 0]

  setnames(df, old = "N", new = paste0(i, "_sigCis-SNP_number"))
}

# 調整數值
df[, pval_OQN := format(pval_OQN, digits = 4, scientific = T)]
df[, pval_QN := format(pval_QN, digits = 4, scientific = T)]
df[, pval_OQN_BH := format(pval_OQN_BH, digits = 4, scientific = T)]
df[, pval_QN_BH := format(pval_QN_BH, digits = 4, scientific = T)]


df[, (9:20) := lapply(.SD, function(x) format(x, scientific = FALSE)), .SDcols = 9:20]



fwrite(df,
  "D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/exp_different.txt",
  row.names = F, col.names = T, sep = "\t"
)

rm(list = ls())
gc()

mix_info <- function(inputname, outputname) {
  if (!is.character(inputname) || length(inputname) != 1 || !is.character(outputname) || length(outputname) != 1) {
    stop("inputname, outputname 必須是一個字串")
  }

  fdr <- fread(inputname, header = T)

  # 對於多區間的probe 名稱"ILMN_2215025_2" 取出 "ILMN_2215025"
  fdr[, gene2 := gsub("^([^_]+_[^_]+).*", "\\1", gene)]

  # merge, due to some snp MAF<0.05,
  fdr <- merge(fdr, pvalue_oqn, by.x = "gene2", by.y = "PROBE_ID")
  fdr <- merge(fdr, probe_pos, by.x = "gene", by.y = "Gene")
  fdr <- merge(fdr, maf, by.x = "SNP", by.y = "hg18_snpID")
  fdr <- merge(fdr, probe_info, by.x = "gene2", by.y = "PROBE_ID")


  # 調整
  fdr[, PROBE_COORDINATES := paste0(start, "-", end)]
  fdr[, c("start", "end", "gene", "INFO", "ER2", "chr") := NULL]
  setnames(fdr, old = c("gene2", "CHROMOSOME"), new = c("Probe", "CHR"))


  setcolorder(fdr, c(
    "CHR", "ProbeID", "Probe", "Gene", "PROBE_COORDINATES",
    "SNP", "rsID", "R2", "impute_type", "MAF", "REF", "ALT", "beta", "t-stat",
    "p-value", "FDR", "qvalue", "sig_pval_Bonfi", "pval_OQN", "t_OQN",
    "pval_OQN_BH", "sig_OQN", "sig_OQN_Bonfi",
    "pval_QN", "t_QN", "pval_QN_BH", "sig_QN", "sig_QN_Bonfi"
  ))

  fdr <- fdr[order(FDR), ]

  # 調整數值
  fdr[, MAF := round(MAF, digits = 4)]
  fdr[, beta := round(beta, digits = 4)]
  fdr[, `t-stat` := round(`t-stat`, digits = 4)]
  fdr[, `p-value` := format(`p-value`, digits = 4, scientific = T)]
  fdr[, FDR := format(FDR, digits = 4, scientific = T)]
  fdr[, pval_OQN := format(pval_OQN, digits = 4, scientific = T)]
  fdr[, pval_QN := format(pval_QN, digits = 4, scientific = T)]
  fdr[, t_OQN := format(t_OQN, digits = 4, scientific = T)]
  fdr[, t_QN := format(t_QN, digits = 4, scientific = T)]
  fdr[, R2 := round(R2, digits = 4)]
  fdr[, qvalue := format(qvalue, digits = 4, scientific = T)]


  fwrite(fdr, outputname,
    row.names = F, col.names = T, sep = "\t"
  )
}


pvalue_oqn <- fread("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/exp_different.txt", header = T)

# 把各檔案cis snp number 資料刪掉，info 檔案用不到
pvalue_oqn <- pvalue_oqn[, c(
  "maf_gt_N_pvalue_FDR_sigCis-SNP_number", "maf_gt_T_pvalue_FDR_sigCis-SNP_number",
  "maf_gt_N_MOSTpvalue_FDR_sigCis-SNP_number", "maf_gt_T_MOSTpvalue_FDR_sigCis-SNP_number",
  "maf_ds_N_pvalue_FDR_sigCis-SNP_number", "maf_ds_T_pvalue_FDR_sigCis-SNP_number",
  "maf_ds_N_MOSTpvalue_FDR_sigCis-SNP_number",
  "maf_ds_T_MOSTpvalue_FDR_sigCis-SNP_number"
) := NULL]


probe_pos <- fread("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/probe_pos.txt", header = T)
maf <- fread("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/cis_hg18_19_snp_MAF.txt", header = T)
maf[, hg19_snpID := NULL]

probe_info <- fread("D:/oral_cancer/expression/expression_data/probe24526infoB.txt", header = T)
probe_info[, c("TargetID", "CHROMOSOME", "PROBE_COORDINATES") := NULL]


file_name <- c(
  "D:/oral_cancer/expression/trash/maf_gt_N_pvalue_FDR",
  "D:/oral_cancer/expression/trash/maf_gt_T_pvalue_FDR",
  "D:/oral_cancer/expression/trash/maf_gt_N_MOSTpvalue_FDR",
  "D:/oral_cancer/expression/trash/maf_gt_T_MOSTpvalue_FDR",
  "D:/oral_cancer/expression/trash/maf_ds_N_MOSTpvalue_FDR",
  "D:/oral_cancer/expression/trash/maf_ds_T_MOSTpvalue_FDR"
)

file_name_new <- c(
  "D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/maf_gt_N_pvalue_FDR",
  "D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/maf_gt_T_pvalue_FDR",
  "D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/maf_gt_N_MOSTpvalue_FDR",
  "D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/maf_gt_T_MOSTpvalue_FDR",
  "D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/maf_ds_N_MOSTpvalue_FDR",
  "D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/maf_ds_T_MOSTpvalue_FDR"
)

sig_number <- fread("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/exp_different.txt", header = T)
sig_number <- sig_number %>%
  select(c(
    "PROBE_ID",
    "maf_gt_N_pvalue_FDR_sigCis-SNP_number",
    "maf_gt_T_pvalue_FDR_sigCis-SNP_number",
    "maf_gt_N_MOSTpvalue_FDR_sigCis-SNP_number",
    "maf_gt_T_MOSTpvalue_FDR_sigCis-SNP_number",
    "maf_ds_N_pvalue_FDR_sigCis-SNP_number",
    "maf_ds_T_pvalue_FDR_sigCis-SNP_number",
    "maf_ds_N_MOSTpvalue_FDR_sigCis-SNP_number",
    "maf_ds_T_MOSTpvalue_FDR_sigCis-SNP_number"
  ))


# 對這些檔案生成info 檔
for (idx in 1:length(file_name)) {
  i <- file_name[idx]
  j <- file_name_new[idx]

  mix_info(
    paste0(i, ".txt"),
    paste0(j, "_info.txt")
  )

  # 新增 cis sig snp
  a <- paste0(j, "_info.txt") %>%
    fread()

  # 挑出1, idx+1 col
  b <- sig_number[, .SD, .SDcols = c(1, idx + 1)]
  names(b) <- c("PROBE_ID", "sig_Cis-SNP_number")
  a <- merge(a, b, by.x = "Probe", by.y = "PROBE_ID")


  setcolorder(a, c("ProbeID", "Probe", "Gene", "CHR", "PROBE_COORDINATES", "sig_Cis-SNP_number", "SNP", "rsID", "impute_type", "R2", "MAF", "REF", "ALT", "beta", "t-stat", "p-value", "FDR", "qvalue", "sig_pval_Bonfi", "t_OQN", "pval_OQN", "pval_OQN_BH", "sig_OQN", "sig_OQN_Bonfi", "t_QN", "pval_QN", "pval_QN_BH", "sig_QN", "sig_QN_Bonfi"))

  fwrite(a, paste0(j, "_info.txt"),
    row.names = F, col.names = T, sep = "\t"
  )
}



gt_N <- fread("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/maf_gt_N_pvalue_FDR_info.txt")
gt_N[
  !gt_N$REF %in% c("A", "T", "C", "G") | !gt_N$ALT %in% c("A", "T", "C", "G"),
  rsID := NA_character_
]

fwrite(gt_N, "D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/maf_gt_N_pvalue_FDR_info.txt",
  row.names = F, col.names = T, sep = "\t"
)

gt_T <- fread("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/maf_gt_T_pvalue_FDR_info.txt")
gt_T[
  !gt_T$REF %in% c("A", "T", "C", "G") | !gt_T$ALT %in% c("A", "T", "C", "G"),
  rsID := NA_character_
]

fwrite(gt_T, "D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/maf_gt_T_pvalue_FDR_info.txt",
  row.names = F, col.names = T, sep = "\t"
)

file_name <- c(
  "D:/oral_cancer/expression/trash/maf_gt_N_pvalue.txt",
  "D:/oral_cancer/expression/trash/maf_gt_T_pvalue.txt"
)


for (i in file_name) {
  a <- fread(i)
  a[, gene2 := gsub("^([^_]+_[^_]+).*", "\\1", gene)]
  a[, gene := NULL]

  sig_fdr_snp <- unique(a, by = "SNP") %>% nrow()
  sig_fdr_asso <- a %>% nrow()
  sig_fdr_probe <- unique(a, by = "gene2") %>% nrow()


  cat(i, "\n")

  cat("eQTL :", "\n")
  cat("unique Probe number: ")
  print(sig_fdr_probe)
  cat("unique SNPs number: ")
  print(sig_fdr_snp)
  cat("association number: ")
  print(sig_fdr_asso)

  cat(rep("\n", 3))
}

a <- fread("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/exp_different.txt", header = T)
zero_index <- (rowSums(a[, 13:20]) == 0) %>%
  which()


# 刪掉association 都是 MAF<0.05 snp 的 probe
a <- a[-zero_index, ]
names(a)

# OQN 顯著不同，eQTL FDR<0.05 的probe 數目
b <- a %>% filter(sig_OQN == 1)
for (i in 13:20) {
  which(b[[i]] != 0) %>%
    length()         %>%
    print()
}
cat("\n")


# QN 顯著不同，eQTL FDR<0.05 的probe 數目
b <- a %>% filter(sig_QN == 1)
for (i in 13:20) {
  which(b[[i]] != 0) %>%
    length()         %>%
    print()
}
cat("\n")





# bonfi ####
# OQN 顯著不同，eQTL FDR<0.05 的probe 數目
b <- a %>% filter(sig_OQN_Bonfi == 1)
for (i in 13:20) {
  which(b[[i]] != 0) %>%
    length()         %>%
    print()
}
cat("\n")

# QN 顯著不同，eQTL FDR<0.05 的probe 數目
b <- a %>% filter(sig_QN_Bonfi == 1)
for (i in 13:20) {
  which(b[[i]] != 0) %>%
    length()         %>%
    print()
}
cat("\n")





file_name <- c(
  "D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/maf_gt_N_pvalue_FDR_info.txt",
  "D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/maf_gt_T_pvalue_FDR_info.txt"
)

total <- fread("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/exp_different.txt", header = T)

# 抓出22217 probe 中有cis-SNP 的20913個，藉此找出FDR>0.05 的 probe
a <- fread("D:/oral_cancer/expression/trash/maf_gt_N_MOSTpvalue.txt")
a[, gene2 := gsub("^([^_]+_[^_]+).*", "\\1", gene)]
a[, gene := NULL]

total <- total[PROBE_ID %in% a$gene2]


# 對這些檔案生成N,T 顯著snp, probe, association number
for (i in file_name) {
  a <- fread(i)

  # 挑出FDR>0.05 probe
  probe <- setdiff(total$PROBE_ID, a$Probe)

  # N,T different FDR sig (OQN)
  sig_OQN_probe <- total %>% filter(sig_OQN == 1)
  sig_OQN_probe <- sig_OQN_probe[PROBE_ID %in% probe] %>% nrow()


  # N,T different FDR sig (QN)
  sig_QN_probe <- total %>% filter(sig_QN == 1)
  sig_QN_probe <- sig_QN_probe[PROBE_ID %in% probe] %>% nrow()

  # N,T different Bonfi sig (OQN)
  Bonfi_OQN_probe <- total %>% filter(sig_OQN_Bonfi == 1)
  Bonfi_OQN_probe <- Bonfi_OQN_probe[PROBE_ID %in% probe] %>% nrow()

  # N,T different Bonfi sig (QN)
  Bonfi_QN_probe <- total %>% filter(sig_QN_Bonfi == 1)
  Bonfi_QN_probe <- Bonfi_QN_probe[PROBE_ID %in% probe] %>% nrow()


  cat(i, "\n")

  cat("eQTL q-value>0.05:", "\n")
  cat("unique Probe number: ", probe %>% length(), "\n")

  cat("N,T different q-value<0.05 (OQN):", "\n")
  cat("unique Probe number: ", sig_OQN_probe, "\n")

  cat("N,T different Bonfi (OQN):", "\n")
  cat("unique Probe number: ", Bonfi_OQN_probe, "\n")

  cat("N,T different q-value<0.05 (QN):", "\n")
  cat("unique Probe number: ", sig_QN_probe, "\n")

  cat("N,T different Bonfi (QN):", "\n")
  cat("unique Probe number: ", Bonfi_QN_probe, "\n")

  cat(rep("\n", 3))
}


total <- fread("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/exp_different.txt", header = T)
zero_index <- (rowSums(total[, 13:20]) == 0) %>%
  which()

# 刪掉association 都是 MAF<0.05 snp 的 probe
total <- total[-zero_index, ]

# gt pvalue ####
a <- fread("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/maf_gt_N_pvalue_FDR_info.txt")
b <- fread("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/maf_gt_T_pvalue_FDR_info.txt")
both_probe <- intersect(a$Probe, b$Probe)
cat("gt pvalue", "\n")

both_probe %>% length()
# OQN
sig_OQN_probe <- total %>% filter(sig_OQN == 1)
sig_OQN_probe <- sig_OQN_probe[PROBE_ID %in% both_probe] %>% nrow()
sig_OQN_probe

# OQN_bon
sig_OQN_probe <- total %>% filter(sig_OQN_Bonfi == 1)
sig_OQN_probe <- sig_OQN_probe[PROBE_ID %in% both_probe] %>% nrow()
sig_OQN_probe

# QN
sig_QN_probe <- total %>% filter(sig_QN == 1)
sig_QN_probe <- sig_QN_probe[PROBE_ID %in% both_probe] %>% nrow()
sig_QN_probe

# QN_bon
sig_QN_probe <- total %>% filter(sig_QN_Bonfi == 1)
sig_QN_probe <- sig_QN_probe[PROBE_ID %in% both_probe] %>% nrow()
sig_QN_probe

cat(rep("\n", 2))


# gt MOSTpvalue ####
a <- fread("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/maf_gt_N_MOSTpvalue_FDR_info.txt")
b <- fread("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/maf_gt_T_MOSTpvalue_FDR_info.txt")
both_probe <- intersect(a$Probe, b$Probe)
cat("gt MOSTpvalue", "\n")

both_probe %>% length()
# OQN
sig_OQN_probe <- total %>% filter(sig_OQN == 1)
sig_OQN_probe <- sig_OQN_probe[PROBE_ID %in% both_probe] %>% nrow()
sig_OQN_probe

# OQN_bon
sig_OQN_probe <- total %>% filter(sig_OQN_Bonfi == 1)
sig_OQN_probe <- sig_OQN_probe[PROBE_ID %in% both_probe] %>% nrow()
sig_OQN_probe

# QN
sig_QN_probe <- total %>% filter(sig_QN == 1)
sig_QN_probe <- sig_QN_probe[PROBE_ID %in% both_probe] %>% nrow()
sig_QN_probe

# QN_bon
sig_QN_probe <- total %>% filter(sig_QN_Bonfi == 1)
sig_QN_probe <- sig_QN_probe[PROBE_ID %in% both_probe] %>% nrow()
sig_QN_probe




# ds pvalue ####
a <- fread("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/maf_ds_N_pvalue_FDR_info.txt")
b <- fread("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/maf_ds_T_pvalue_FDR_info.txt")
both_probe <- intersect(a$Probe, b$Probe)
cat("ds pvalue", "\n")

both_probe %>% length()
# OQN
sig_OQN_probe <- total %>% filter(sig_OQN == 1)
sig_OQN_probe <- sig_OQN_probe[PROBE_ID %in% both_probe] %>% nrow()
sig_OQN_probe

# OQN_bon
sig_OQN_probe <- total %>% filter(sig_OQN_Bonfi == 1)
sig_OQN_probe <- sig_OQN_probe[PROBE_ID %in% both_probe] %>% nrow()
sig_OQN_probe

# QN
sig_QN_probe <- total %>% filter(sig_QN == 1)
sig_QN_probe <- sig_QN_probe[PROBE_ID %in% both_probe] %>% nrow()
sig_QN_probe

# QN_bon
sig_QN_probe <- total %>% filter(sig_QN_Bonfi == 1)
sig_QN_probe <- sig_QN_probe[PROBE_ID %in% both_probe] %>% nrow()
sig_QN_probe

cat(rep("\n", 2))


# ds MOSTpvalue ####
a <- fread("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/maf_ds_N_MOSTpvalue_FDR_info.txt")
b <- fread("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/maf_ds_T_MOSTpvalue_FDR_info.txt")
both_probe <- intersect(a$Probe, b$Probe)
cat("ds MOSTpvalue", "\n")

both_probe %>% length()
# OQN
sig_OQN_probe <- total %>% filter(sig_OQN == 1)
sig_OQN_probe <- sig_OQN_probe[PROBE_ID %in% both_probe] %>% nrow()
sig_OQN_probe

# OQN_bon
sig_OQN_probe <- total %>% filter(sig_OQN_Bonfi == 1)
sig_OQN_probe <- sig_OQN_probe[PROBE_ID %in% both_probe] %>% nrow()
sig_OQN_probe

# QN
sig_QN_probe <- total %>% filter(sig_QN == 1)
sig_QN_probe <- sig_QN_probe[PROBE_ID %in% both_probe] %>% nrow()
sig_QN_probe

# QN_bon
sig_QN_probe <- total %>% filter(sig_QN_Bonfi == 1)
sig_QN_probe <- sig_QN_probe[PROBE_ID %in% both_probe] %>% nrow()
sig_QN_probe


file_name <- c(
  "maf_gt_N_pvalue",
  "maf_gt_T_pvalue",
  "maf_ds_N_pvalue",
  "maf_ds_T_pvalue"
)

for (i in file_name) {
  a <- paste0("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/", i, "_FDR_info.txt") %>%
    fread()
  k <- a[, .N, by = "Probe"]
  plot(k$N, ylab = "SNP Number", main = i)
  hist(k$N, ylab = "Frequency", xlab = "SNP Number", main = i, breaks = 100, col = "skyblue")
  boxplot(k$N, main = i)
  cat(i, "\n")
  cat("mean of SNP Number in each probe: ", mean(k$N), "\n")
  cat("median of SNP Number in each probe: ", median(k$N), "\n")
  cat("\n")
}


rm(list = ls())
gc()

file_name <- c(
  "D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/maf_gt_N_pvalue_FDR_info.txt",
  "D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/maf_gt_T_pvalue_FDR_info.txt",
  "D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/maf_gt_N_MOSTpvalue_FDR_info.txt",
  "D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/maf_gt_T_MOSTpvalue_FDR_info.txt"
)

# 對這些檔案生成N,T 顯著snp, probe, association number
for (i in file_name) {
  a <- fread(i)

  # FDR sig
  sig_fdr_snp <- unique(a, by = "SNP") %>% nrow()
  sig_fdr_asso <- a %>% nrow()
  sig_fdr_probe <- unique(a, by = "ProbeID") %>% nrow()

  # N,T different FDR sig (OQN)
  sig_OQN_snp <- a       %>%
    filter(sig_OQN == 1) %>%
    unique(by = "SNP")   %>%
    nrow()
  sig_OQN_asso <- a      %>%
    filter(sig_OQN == 1) %>%
    nrow()
  sig_OQN_probe <- a       %>%
    filter(sig_OQN == 1)   %>%
    unique(by = "ProbeID") %>%
    nrow()

  # N,T different FDR sig (QN)
  sig_QN_snp <- a       %>%
    filter(sig_QN == 1) %>%
    unique(by = "SNP")  %>%
    nrow()
  sig_QN_asso <- a      %>%
    filter(sig_QN == 1) %>%
    nrow()
  sig_QN_probe <- a        %>%
    filter(sig_QN == 1)    %>%
    unique(by = "ProbeID") %>%
    nrow()

  # N,T different Bonfi sig (OQN)
  Bonfi_OQN_snp <- a           %>%
    filter(sig_OQN_Bonfi == 1) %>%
    unique(by = "SNP")         %>%
    nrow()
  Bonfi_OQN_asso <- a          %>%
    filter(sig_OQN_Bonfi == 1) %>%
    nrow()
  Bonfi_OQN_probe <- a         %>%
    filter(sig_OQN_Bonfi == 1) %>%
    unique(by = "ProbeID")     %>%
    nrow()

  # N,T different Bonfi sig (QN)
  Bonfi_QN_snp <- a           %>%
    filter(sig_QN_Bonfi == 1) %>%
    unique(by = "SNP")        %>%
    nrow()
  Bonfi_QN_asso <- a          %>%
    filter(sig_QN_Bonfi == 1) %>%
    nrow()
  Bonfi_QN_probe <- a         %>%
    filter(sig_QN_Bonfi == 1) %>%
    unique(by = "ProbeID")    %>%
    nrow()


  cat(i, "\n")

  cat("eQTL q-value<0.05:", "\n")
  cat("unique Probe number: ")
  print(sig_fdr_probe)
  cat("unique SNPs number: ")
  print(sig_fdr_snp)
  cat("association number: ")
  print(sig_fdr_asso)

  cat("N,T different q-value<0.05 (OQN):", "\n")
  cat("unique Probe number: ")
  print(sig_OQN_probe)
  cat("unique SNPs number: ")
  print(sig_OQN_snp)
  cat("association number: ")
  print(sig_OQN_asso)

  cat("N,T different Bonfi (OQN):", "\n")
  cat("unique Probe number: ")
  print(Bonfi_OQN_probe)
  cat("unique SNPs number: ")
  print(Bonfi_OQN_snp)
  cat("association number: ")
  print(Bonfi_OQN_asso)

  cat("N,T different q-value<0.05 (QN):", "\n")
  cat("unique Probe number: ")
  print(sig_QN_probe)
  cat("unique SNPs number: ")
  print(sig_QN_snp)
  cat("association number: ")
  print(sig_QN_asso)


  cat("N,T different Bonfi (QN):", "\n")
  cat("unique Probe number: ")
  print(Bonfi_QN_probe)
  cat("unique SNPs number: ")
  print(Bonfi_QN_snp)
  cat("association number: ")
  print(Bonfi_QN_asso)

  cat(rep("\n", 3))
}


rm(list = ls())
gc()

file_name <- c(
  "D:/oral_cancer/expression/trash/maf_gt_N_MOSTpvalue.txt",
  "D:/oral_cancer/expression/trash/maf_gt_T_MOSTpvalue.txt",
  "D:/oral_cancer/expression/trash/maf_gt_N_MOSTpvalue_FDR.txt",
  "D:/oral_cancer/expression/trash/maf_gt_T_MOSTpvalue_FDR.txt",
  "D:/oral_cancer/expression/trash/maf_ds_N_MOSTpvalue.txt",
  "D:/oral_cancer/expression/trash/maf_ds_T_MOSTpvalue.txt",
  "D:/oral_cancer/expression/trash/maf_ds_N_MOSTpvalue_FDR.txt",
  "D:/oral_cancer/expression/trash/maf_ds_T_MOSTpvalue_FDR.txt"
)

# 對這些檔案生成N,T 顯著snp, probe, association number
for (i in file_name) {
  a <- fread(i)

  # FDR sig
  sig_fdr_snp <- str_extract(a$SNP, ".*(?=:)") %>%
    as.numeric()                               %>%
    table()

  cat(i, "\n")

  print(sig_fdr_snp)
  cat(rep("\n", 3))
}

maf_hg18_19 <- fread("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/cis_hg18_19_snp_MAF.txt", header = T)

maf_hg18_19 <- maf_hg18_19 %>%
  filter(MAF > 0.05)

fwrite(
  data.table(
    chr = str_extract(maf_hg18_19$hg19_snpID, ".*(?=:)"),
    pos_start = str_extract(maf_hg18_19$hg19_snpID, "(?<=\\:)(\\d+)"),
    pos_end = str_extract(maf_hg18_19$hg19_snpID, "(?<=\\:)(\\d+)")
  ),
  "D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/maf_cis_snpID_hg19.txt",
  row.names = F, col.names = F, sep = "\t"
)


# maf_cis_snpID_hg19.txt 檔案放 22   16050115   16050115，挑出 maf_cis_snp
cmd <- paste0(
  "wsl bash -c 'set -e; for i in {1..22}; do ",
  "echo Processing chr${i}...; ",
  "tabix -R /mnt/d/oral_cancer/expression/outcome/multiple_nucleotide_variant/maf_cis_snpID_hg19.txt /mnt/d/oral_cancer/expression/trash/chr${i}_cis_snp.vcf.gz > /mnt/d/oral_cancer/expression/trash/maf_chr${i}.vcf; ",
  "done'"
)
system(cmd)


# 對maf_cis_snp 轉成 rsID
cmd <- paste0(
  "wsl bash -c 'set -e; for i in {1..22}; do ",
  "echo Running VEP on chr${i}...; ",
  "/root/ensembl-vep/vep --cache --offline ",
  "--dir_cache ~/.vep ",
  "--assembly GRCh37 ",
  "--fasta ~/.vep/homo_sapiens/115_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa ",
  "--input_file /mnt/d/oral_cancer/expression/trash/maf_chr${i}.vcf ",
  "--output_file /mnt/d/oral_cancer/expression/outcome/multiple_nucleotide_variant/chr${i}_vep.vcf ",
  "--vcf --force_overwrite --check_existing ",
  "--custom ~/.vep/homo_sapiens-chr${i}.vcf.gz,dbSNP,vcf,exact,0,ID; ",
  "echo Done chr${i}; ",
  "done'"
)
system(cmd)




# 合併
file_paths <- sprintf("D:/oral_cancer/expression/trash/chr%d_vep.txt", 1:22)

ds <- lapply(file_paths, fread, header = F) %>%
  rbindlist(use.names = F)

names(ds) <- c("rsID", "pos", "REF", "ALT")

fwrite(ds, "D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/chr1-22_vep.txt",
  row.names = F, col.names = T, sep = "\t"
)

snp_pos <- fread("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/cis_hg18_19_snp_MAF.txt", header = T)
rsID <- fread("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/chr1-22_vep.txt", header = T)
mix_rs_pos <- merge(snp_pos, rsID, by.x = "hg19_snpID", by.y = "pos", all.x = T)

mix_rs_pos <- mix_rs_pos %>%
  unique()

setcolorder(mix_rs_pos, c("hg18_snpID", "hg19_snpID", "rsID", "MAF", "REF", "ALT"))
fwrite(mix_rs_pos, "D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/cis_hg18_19_snp_MAF.txt",
  row.names = F, col.names = T, sep = "\t"
)



df <- fread("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/cis_hg18_19_snp_MAF.txt")


# 去除重複
rs_list <- unique(df$rsID)

get_grch37_coord <- function(rsid) {
  numeric_id <- gsub("rs", "", rsid)
  url <- paste0("https://api.ncbi.nlm.nih.gov/variation/v0/beta/refsnp/", numeric_id)

  resp <- GET(url)
  if (resp$status_code != 200) {
    return(data.frame(rsid = rsid, chr = NA, pos = NA))
  }

  json <- tryCatch(
    fromJSON(content(resp, "text", encoding = "UTF-8")),
    error = function(e) {
      return(NULL)
    }
  )
  if (is.null(json)) {
    return(data.frame(rsid = rsid, chr = NA, pos = NA))
  }

  placements <- json$primary_snapshot_data$placements_with_allele
  if (length(placements) == 0) {
    return(data.frame(rsid = rsid, chr = NA, pos = NA))
  }

  # ----------- 找 GRCh37 使用 seq_id_traits_by_assembly -----------

  grch37_index <- NULL

  for (i in 1:nrow(placements)) {
    pa <- placements$placement_annot$seq_id_traits_by_assembly[[i]]

    if (length(pa) != 0) {
      assembly_names <- pa$assembly_name
      if (assembly_names == "GRCh37.p13") {
        grch37_index <- i
        break
      }
    }
  }

  if (is.null(grch37_index)) {
    return(data.frame(rsid = rsid, chr = NA, pos = NA))
  }

  # ----------- 抽取 GRCh37 座標 -----------

  target <- placements[grch37_index, ]

  # seq_id 例如 "NC_000020.10" → chr 20
  chr_raw <- target$seq_id
  chr <- gsub("NC_0+|\\..*$", "", chr_raw) # 最穩定抽法

  # 取第一個 allele 的位置（GRCh37）
  allele <- target$alleles[[1]]
  pos <- allele$allele$spdi$position + 1 # NCBI 是 0-based，要轉為 1-based

  data.frame(
    rsid = rsid,
    chr = chr,
    pos = pos
  )
}

results <- list()
for (i in rs_list) {
  results[[i]] <- get_grch37_coord(i)
  Sys.sleep(0.5) # 每次呼叫間隔 0.5 秒 → 每秒最多 2 次
}

final_df <- bind_rows(results) %>%
  unique()
rownames(final_df) <- NULL


# 把hg19/GRCH37 id 併進資料
final_df <- as.data.table(final_df)
final_df[, SNP_GRCH37 := paste0(chr, ":", pos)]
final_df[, c("chr", "pos") := NULL]


{ # aa
  data_T <- fread("D:/oral_cancer/expression/trash/maf_gt_T_pvalue.txt")
  data_N <- fread("D:/oral_cancer/expression/trash/maf_gt_N_pvalue.txt")
  m_total_T <- nrow(data_T)
  m_total_N <- nrow(data_N)


  # T ####
  # 篩選 p < 0.05
  p_T <- sort(data_T[, `p-value`])
  rank_T <- 1:length(p_T)

  # 找Bonferroni cutoff
  bonfi_threshold_T <- 0.05 / m_total_T

  # 找每個pvalue 的BH cutoff
  bh_threshold_T <- rank_T / m_total_T * 0.05
  # 找小於BH cutoff 的pvalue 當中最大的(顯著pvalue 最大值)
  sig_cut_T <- max(p_T[p_T <= bh_threshold_T])

  # 顯著pvalue 最大值的 1.5 倍
  ylim_max_T <- sig_cut_T * 1.5

  # 篩選出在「縮放範圍內」的pvalue
  p_T_zoomed <- p_T[p_T <= ylim_max_T]
  # 找縮放範圍內的 BH cutoff
  bh_threshold_T_zoomed <- bh_threshold_T[which(p_T <= ylim_max_T)]

  # X 軸是 rank/m，找縮放範圍內的
  plot_x_T_zoomed <- which(p_T <= ylim_max_T)


  # N####
  p_N <- sort(data_N[, `p-value`])
  m_subset_N <- length(p_N)
  rank_N <- 1:m_subset_N

  bonfi_threshold_N <- 0.05 / m_total_N
  bh_threshold_N <- rank_N / m_total_N * 0.05
  sig_cut_N <- max(p_N[p_N <= bh_threshold_N])

  ylim_max_N <- sig_cut_N * 1.5

  p_N_zoomed <- p_N[which(p_N <= ylim_max_N)]
  bh_threshold_N_zoomed <- bh_threshold_N[which(p_N <= ylim_max_N)]
  plot_x_N_zoomed <- which(p_N <= ylim_max_N)


  # --- 繪圖 (合併到同一頁) ---####
  png("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/BHplot_maf_gt_pvalue_original.png", width = 2000, height = 1000, res = 200)
  par(mfrow = c(1, 2))

  # --- 開始畫 T 組的圖 ---
  plot(rank_T[c(seq(1, length(rank_T), by = 4), length(rank_T))],
    p_T[c(seq(1, length(p_T), by = 4), length(p_T))],
    pch = ".",
    col = "black",
    xlab = "Rank/m (m=88867137)",
    ylab = "P-value",
    main = "Tumor Part: All Associaiton"
  )

  lines(rank_T[c(seq(1, length(rank_T), by = 4), length(rank_T))],
    bh_threshold_T[c(seq(1, length(bh_threshold_T), by = 4), length(bh_threshold_T))],
    col = "red",
    lwd = 2
  )

  # 標示橘色點 (在縮放範圍內，且 p <= sig_cut 的點)
  sig_indices_T_zoomed <- (p_T_zoomed <= sig_cut_T)
  points(plot_x_T_zoomed[sig_indices_T_zoomed], p_T_zoomed[sig_indices_T_zoomed],
    col = "orange",
    pch = "."
  )
  abline(h = bonfi_threshold_T, col = "blue", lty = 2)
  abline(a = 0, b = 1, col = "green", lty = 2)

  legend("topleft",
    legend = c(
      "BH (q-value< 0.05)",
      paste("Bonferroni (p-value <=", format(bonfi_threshold_T, digits = 2), ")"),
      "slope 1"
    ),
    col = c("red", "blue", "green"),
    lty = c(1, 2, 2),
    lwd = c(2, 2, 2)
  )

  # --- 開始畫 N 組的圖 (只畫縮放後的資料) ---
  plot(rank_N[c(seq(1, length(rank_N), by = 4), length(rank_N))],
    p_N[c(seq(1, length(p_N), by = 4), length(p_N))],
    pch = ".",
    col = "black",
    xlab = "Rank/m (m=88867137)",
    ylab = "P-value",
    main = "Normal Part: All Associaiton"
  )

  lines(rank_N[c(seq(1, length(rank_N), by = 4), length(rank_N))],
    bh_threshold_N[c(seq(1, length(bh_threshold_N), by = 4), length(bh_threshold_N))],
    col = "red",
    lwd = 2
  )

  sig_indices_N_zoomed <- (p_N_zoomed <= sig_cut_N)
  points(plot_x_N_zoomed[sig_indices_N_zoomed], p_N_zoomed[sig_indices_N_zoomed],
    col = "orange",
    pch = "."
  )
  abline(h = bonfi_threshold_N, col = "blue", lty = 2)
  abline(a = 0, b = 1, col = "green", lty = 2)

  legend("topleft",
    legend = c(
      "BH (q-value< 0.05)",
      paste("Bonferroni (p-value <=", format(bonfi_threshold_N, digits = 2), ")"),
      "slope 1"
    ),
    col = c("red", "blue", "green"),
    lty = c(1, 2, 2),
    lwd = c(2, 2, 2)
  )

  # 關閉檔案
  dev.off()

  data_N[, FDR := p.adjust(`p-value`, method = "BH") %>%
    format(digits = 10, scientific = T)              %>%
    as.numeric()]

  data_T[, FDR := p.adjust(`p-value`, method = "BH") %>%
    format(digits = 10, scientific = T)              %>%
    as.numeric()]

  data_N[, qvalue := qvalue(`p-value`)$qvalues]
  data_scale_N <- data_N %>%
    filter(FDR < 0.05)

  data_T[, qvalue := qvalue(`p-value`)$qvalues]
  data_scale_T <- data_T %>%
    filter(FDR < 0.05)

  setkey(data_T, FDR, qvalue)
  setkey(data_N, FDR, qvalue)

  # 4個取一個
  data_sub_T <- data_T[c(seq(1, nrow(data_T), by = 4), nrow(data_T)), ]
  data_sub_N <- data_N[c(seq(1, nrow(data_N), by = 4), nrow(data_N)), ]


  # FDR-Qvalue 長條圖 ####
  png("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/FDR_qvalue_maf_gt_pvalue_differnet.png", width = 2000, height = 1000, res = 200)
  par(mfrow = c(1, 2))

  plot(1:nrow(data_sub_N),
    data_sub_N$FDR - data_sub_N$qvalue,
    pch = ".",
    col = "red",
    main = "Normal Part Difference of BH and Qvalue",
    xlab = "Index",
    ylab = "Freq"
  )

  plot(1:nrow(data_sub_T),
    data_sub_T$FDR - data_sub_T$qvalue,
    pch = ".",
    col = "red",
    main = "Tumor Part Difference of BH and Qvalue",
    xlab = "Index",
    ylab = "Freq"
  )

  dev.off()


  png("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/FDR_qvalue_maf_gt_pvalue_differnet_2.png", width = 2000, height = 1000, res = 200)
  par(mfrow = c(1, 2))

  hist(data_sub_N$FDR - data_sub_N$qvalue,
    main = "Normal Part Difference of BH and Qvalue",
    xlab = "BH-Qvalue",
    ylab = "Freq",
    breaks = 100
  )
  hist(data_sub_T$FDR - data_sub_T$qvalue,
    main = "Tumor Part Difference of BH and Qvalue",
    xlab = "BH-Qvalue",
    ylab = "Freq",
    breaks = 100
  )
  dev.off()


  # 畫Qvalue vs FDR圖 ####
  png("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/FDR_qvalue_maf_gt_pvalue_original.png", width = 2000, height = 1000, res = 200)
  par(mfrow = c(1, 2))


  # --- 開始畫 T 組的圖---
  plot(data_sub_T$FDR,
    data_sub_T$qvalue,
    pch = ".",
    col = "black",
    main = "Tumor Part BH V.S. Qvalue",
    xlab = "BH",
    ylab = "Qvalue"
  )
  abline(a = 0, b = 1, col = "green", lty = 2)

  legend("topleft",
    legend = c("slope 1"),
    col = c("green"),
    lty = c(2),
    lwd = c(2)
  )

  # --- 開始畫 N 組的圖---
  plot(data_sub_N$FDR,
    data_sub_N$qvalue,
    pch = ".",
    col = "black",
    main = "Normal Part BH V.S. Qvalue",
    xlab = "BH",
    ylab = "Qvalue"
  )
  abline(a = 0, b = 1, col = "green", lty = 2)

  legend("topleft",
    legend = c("slope 1"),
    col = c("green"),
    lty = c(2),
    lwd = c(2)
  )


  # 關閉檔案
  dev.off()

  # --- 繪圖 (合併到同一頁) ---####
  png("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/BHplot_maf_gt_pvalue_scale.png", width = 2000, height = 1000, res = 200)
  par(mfrow = c(1, 2))

  # --- 開始畫 T 組的圖 (只畫縮放後的資料) ----
  plot(plot_x_T_zoomed, p_T_zoomed,
    pch = ".",
    col = "black",
    xlab = "Rank/m (m=88867137)",
    ylab = "P-value",
    main = "Tumor Part: All Associaiton",
    ylim = c(0, ylim_max_T)
  )

  lines(plot_x_T_zoomed, bh_threshold_T_zoomed, col = "red", lwd = 2)

  # 標示橘色點 (在縮放範圍內，且 p <= sig_cut 的點)
  sig_indices_T_zoomed <- (p_T_zoomed <= sig_cut_T)
  points(plot_x_T_zoomed[sig_indices_T_zoomed], p_T_zoomed[sig_indices_T_zoomed],
    col = "orange",
    pch = "."
  )
  abline(h = bonfi_threshold_T, col = "blue", lty = 2)
  abline(a = 0, b = 1, col = "green", lty = 2)

  legend("topleft",
    legend = c(
      "BH (q-value< 0.05)",
      paste("Bonferroni (p-value <=", format(bonfi_threshold_T, digits = 2), ")"),
      "slope 1"
    ),
    col = c("red", "blue", "green"),
    lty = c(1, 2, 2),
    lwd = c(2, 2, 2)
  )

  # --- 開始畫 N 組的圖 (只畫縮放後的資料) ---
  plot(plot_x_N_zoomed, p_N_zoomed,
    pch = ".",
    col = "black",
    xlab = "Rank/m (m=88867137)",
    ylab = "P-value",
    main = "Normal Part: All Associaiton",
    ylim = c(0, ylim_max_N)
  ) # 【關鍵】強制設定 Y 軸範圍

  lines(plot_x_N_zoomed, bh_threshold_N_zoomed, col = "red", lwd = 2)

  sig_indices_N_zoomed <- (p_N_zoomed <= sig_cut_N)
  points(plot_x_N_zoomed[sig_indices_N_zoomed], p_N_zoomed[sig_indices_N_zoomed],
    col = "orange",
    pch = "."
  )
  abline(h = bonfi_threshold_N, col = "blue", lty = 2)
  abline(a = 0, b = 1, col = "green", lty = 2)

  legend("topleft",
    legend = c(
      "BH (q-value< 0.05)",
      paste("Bonferroni (p-value <=", format(bonfi_threshold_N, digits = 2), ")"),
      "slope 1"
    ),
    col = c("red", "blue", "green"),
    lty = c(1, 2, 2),
    lwd = c(2, 2, 2)
  )

  # 關閉檔案
  dev.off()


  cat("How many qvalue<0.05 in tumor part: ", plot_x_T_zoomed[sig_indices_T_zoomed] %>% length())
  cat("How many qvalue<0.05 in normal part: ", plot_x_N_zoomed[sig_indices_N_zoomed] %>% length())


  # 設定 x,y 軸範圍
  cutoff <- 0.05
  file_name <- paste0("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/FDR_qvalue_maf_gt_pvalue_scale_", cutoff, ".png")

  data_N[, FDR := p.adjust(`p-value`, method = "BH") %>%
    format(digits = 10, scientific = T)              %>%
    as.numeric()]

  data_T[, FDR := p.adjust(`p-value`, method = "BH") %>%
    format(digits = 10, scientific = T)              %>%
    as.numeric()]

  data_N[, qvalue := qvalue(`p-value`)$qvalues]
  data_scale_N <- data_N %>%
    filter(FDR < 0.05)

  data_T[, qvalue := qvalue(`p-value`)$qvalues]
  data_scale_T <- data_T %>%
    filter(FDR < 0.05)

  setkey(data_T, FDR, qvalue)
  setkey(data_N, FDR, qvalue)

  # 4個取一個
  data_sub_T <- data_T[c(seq(1, nrow(data_T), by = 4), nrow(data_T)), ]
  data_sub_N <- data_N[c(seq(1, nrow(data_N), by = 4), nrow(data_N)), ]


  png(file_name, width = 2000, height = 1000, res = 200)
  par(mfrow = c(1, 2))

  # --- 開始畫 T 組的圖---
  plot(data_scale_T$FDR,
    data_scale_T$qvalue,
    pch = ".",
    col = "black",
    main = "Tumor Part BH V.S. Qvalue",
    xlab = "BH",
    ylab = "Qvalue",
    ylim = c(0, cutoff),
    xlim = c(0, cutoff)
  )
  abline(a = 0, b = 1, col = "green", lty = 2)

  legend("topleft",
    legend = c("slope 1"),
    col = c("green"),
    lty = c(2),
    lwd = c(2)
  )

  # --- 開始畫 N 組的圖---
  plot(data_scale_N$FDR,
    data_scale_N$qvalue,
    pch = ".",
    col = "black",
    main = "Normal Part BH V.S. Qvalue",
    xlab = "BH",
    ylab = "Qvalue",
    ylim = c(0, cutoff),
    xlim = c(0, cutoff)
  )
  abline(a = 0, b = 1, col = "green", lty = 2)

  legend("topleft",
    legend = c("slope 1"),
    col = c("green"),
    lty = c(2),
    lwd = c(2)
  )


  # 關閉檔案
  dev.off()
}

data_T <- fread("D:/oral_cancer/expression/trash/maf_gt_T_MOSTpvalue.txt")
data_N <- fread("D:/oral_cancer/expression/trash/maf_gt_N_MOSTpvalue.txt")
m_total_T <- nrow(data_T)
m_total_N <- nrow(data_N)


# T ####
p_T <- sort(data_T[, `p-value`])
rank_T <- 1:length(p_T)

bh_threshold_T <- rank_T / m_total_T * 0.05
sig_cut_T <- max(p_T[p_T <= bh_threshold_T])
bonfi_threshold_T <- 0.05 / m_total_T
ylim_max_T <- sig_cut_T * 2

p_T_zoomed <- p_T[which(p_T <= ylim_max_T)]
bh_threshold_T_zoomed <- bh_threshold_T[which(p_T <= ylim_max_T)]
plot_x_T_zoomed <- which(p_T <= ylim_max_T)





# N####
p_N <- sort(data_N[, `p-value`])
rank_N <- 1:length(p_N)

bh_threshold_N <- rank_N / m_total_N * 0.05
sig_cut_N <- max(p_N[p_N <= bh_threshold_N])
bonfi_threshold_N <- 0.05 / m_total_N
ylim_max_N <- sig_cut_N * 2

p_N_zoomed <- p_N[which(p_N <= ylim_max_N)]
bh_threshold_N_zoomed <- bh_threshold_N[which(p_N <= ylim_max_N)]
plot_x_N_zoomed <- which(p_N <= ylim_max_N)


# --- 繪圖 (合併到同一頁) ---####
png("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/BHplot_maf_gt_mostpvalue_original.png", width = 2000, height = 1000, res = 200)
par(mfrow = c(1, 2))

# --- 開始畫 T 組的圖 (只畫縮放後的資料) ---
plot(rank_T, p_T,
  pch = ".",
  col = "black",
  xlab = "Rank",
  ylab = "P-value",
  main = "Tumor Part: Each Probe with One Sig-SNP"
)

lines(rank_T, bh_threshold_T, col = "red", lwd = 2)

# 標示橘色點 (在縮放範圍內，且 p <= sig_cut 的點)
sig_indices_T_zoomed <- (p_T_zoomed <= sig_cut_T)
points(plot_x_T_zoomed[sig_indices_T_zoomed], p_T_zoomed[sig_indices_T_zoomed],
  col = "orange",
  pch = "."
)
abline(h = bonfi_threshold_T, col = "blue", lty = 2)
abline(a = 0, b = 1, col = "green", lty = 2)

legend("topleft",
  legend = c(
    "BH (q-value< 0.05)",
    paste("Bonferroni (p-value <=", format(bonfi_threshold_T, digits = 2), ")"),
    "slope 1"
  ),
  col = c("red", "blue", "green"),
  lty = c(1, 2, 2),
  lwd = c(2, 2, 2)
)


# --- 開始畫 N 組的圖 (只畫縮放後的資料) ---
plot(rank_N, p_N,
  pch = ".",
  col = "black",
  xlab = "Rank",
  ylab = "P-value",
  main = "Normal Part: Each Probe with One Sig-SNP"
)

lines(rank_N, bh_threshold_N, col = "red", lwd = 2)

sig_indices_N_zoomed <- (p_N_zoomed <= sig_cut_N)
points(plot_x_N_zoomed[sig_indices_N_zoomed], p_N_zoomed[sig_indices_N_zoomed],
  col = "orange",
  pch = "."
)
abline(h = bonfi_threshold_N, col = "blue", lty = 2)
abline(a = 0, b = 1, col = "green", lty = 2)

legend("topleft",
  legend = c(
    "BH (q-value< 0.05)",
    paste("Bonferroni (p-value <=", format(bonfi_threshold_N, digits = 2), ")"),
    "slope 1"
  ),
  col = c("red", "blue", "green"),
  lty = c(1, 2, 2),
  lwd = c(2, 2, 2)
)

# 關閉檔案
dev.off()




library(qvalue)

data_N[, FDR := p.adjust(`p-value`, method = "BH") %>%
  format(digits = 10, scientific = T)              %>%
  as.numeric()]

data_T[, FDR := p.adjust(`p-value`, method = "BH") %>%
  format(digits = 10, scientific = T)              %>%
  as.numeric()]


pi0_hat <- min(max(
  sum(data_N$`p-value` > 0.7) / ((1 - 0.7) * length(data_N$`p-value`)),
  0
), 1)

data_N[, qvalue := qvalue(`p-value`, pi0 = pi0_hat)$qvalues]
data_scale_N <- data_N %>%
  filter(FDR < 0.05)

pi0_hat <- min(max(
  sum(data_T$`p-value` > 0.7) / ((1 - 0.7) * length(data_T$`p-value`)),
  0
), 1)

data_T[, qvalue := qvalue(`p-value`, pi0 = pi0_hat)$qvalues]
data_scale_T <- data_T %>%
  filter(FDR < 0.05)

setkey(data_T, FDR, qvalue)
setkey(data_N, FDR, qvalue)


# 畫圖 ####
png("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/FDR_qvalue_maf_gt_mostpvalue_original.png", width = 2000, height = 1000, res = 200)
par(mfrow = c(1, 2))

plot(data_N$FDR,
  data_N$qvalue,
  col = "black",
  main = "Normal Part: Each Probe with One Sig-SNP",
  xlab = "BH",
  ylab = "Qvalue"
)
abline(a = 0, b = 1, col = "green", lty = 2)

legend("topleft",
  legend = c("slope 1"),
  col = c("green"),
  lty = c(2),
  lwd = c(2)
)


plot(data_T$FDR,
  data_T$qvalue,
  col = "black",
  main = "Tumor Part: Each Probe with One Sig-SNP",
  xlab = "BH",
  ylab = "Qvalue"
)
abline(a = 0, b = 1, col = "green", lty = 2)

legend("topleft",
  legend = c("slope 1"),
  col = c("green"),
  lty = c(2),
  lwd = c(2)
)

# 關閉檔案
dev.off()


# T ####
ylim_max_T <- sig_cut_T * 2

p_T_zoomed <- p_T[which(p_T <= ylim_max_T)]
bh_threshold_T_zoomed <- bh_threshold_T[which(p_T <= ylim_max_T)]
plot_x_T_zoomed <- which(p_T <= ylim_max_T)


# N####
ylim_max_N <- sig_cut_N * 2

p_N_zoomed <- p_N[which(p_N <= ylim_max_N)]
bh_threshold_N_zoomed <- bh_threshold_N[which(p_N <= ylim_max_N)]
plot_x_N_zoomed <- which(p_N <= ylim_max_N)


# --- 繪圖 (合併到同一頁) ---####
png("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/BHplot_maf_gt_mostpvalue_scale.png", width = 2000, height = 1000, res = 200)
par(mfrow = c(1, 2))

# --- 開始畫 T 組的圖 (只畫縮放後的資料) ---
plot(plot_x_T_zoomed, p_T_zoomed,
  pch = ".",
  col = "black",
  xlab = "Rank",
  ylab = "P-value",
  main = "Tumor Part: Each Probe with One Sig-SNP",
  ylim = c(0, ylim_max_T)
) # 【關鍵】強制設定 Y 軸範圍

lines(plot_x_T_zoomed, bh_threshold_T_zoomed, col = "red", lwd = 2)

# 標示橘色點 (在縮放範圍內，且 p <= sig_cut 的點)
sig_indices_T_zoomed <- (p_T_zoomed <= sig_cut_T)
points(plot_x_T_zoomed[sig_indices_T_zoomed], p_T_zoomed[sig_indices_T_zoomed],
  col = "orange",
  pch = "."
)
abline(h = bonfi_threshold_T, col = "blue", lty = 2)
abline(a = 0, b = 1, col = "green", lty = 2)

legend("topleft",
  legend = c(
    "BH (q-value< 0.05)",
    paste("Bonferroni (p-value <=", format(bonfi_threshold_T, digits = 2), ")"),
    "slope 1"
  ),
  col = c("red", "blue", "green"),
  lty = c(1, 2, 2),
  lwd = c(2, 2, 2)
)


# --- 開始畫 N 組的圖 (只畫縮放後的資料) ---
plot(plot_x_N_zoomed, p_N_zoomed,
  pch = ".",
  col = "black",
  xlab = "Rank",
  ylab = "P-value",
  main = "Normal Part: Each Probe with One Sig-SNP",
  ylim = c(0, ylim_max_N)
) # 【關鍵】強制設定 Y 軸範圍

lines(plot_x_N_zoomed, bh_threshold_N_zoomed, col = "red", lwd = 2)

sig_indices_N_zoomed <- (p_N_zoomed <= sig_cut_N)
points(plot_x_N_zoomed[sig_indices_N_zoomed], p_N_zoomed[sig_indices_N_zoomed],
  col = "orange",
  pch = "."
)
abline(h = bonfi_threshold_N, col = "blue", lty = 2)
abline(a = 0, b = 1, col = "green", lty = 2)

legend("topleft",
  legend = c(
    "BH (q-value< 0.05)",
    paste("Bonferroni (p-value <=", format(bonfi_threshold_N, digits = 2), ")"),
    "slope 1"
  ),
  col = c("red", "blue", "green"),
  lty = c(1, 2, 2),
  lwd = c(2, 2, 2)
)

# 關閉檔案
dev.off()



png("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/FDR_qvalue_maf_gt_mostpvalue_scale.png", width = 2000, height = 1000, res = 200)
par(mfrow = c(1, 2))

plot(data_scale_N$FDR,
  data_scale_N$qvalue,
  pch = ".",
  col = "black",
  main = "Normal Part: Each Probe with One Sig-SNP",
  xlab = "BH",
  ylab = "Qvalue",
  xlim = c(0, 0.02)
)
abline(a = 0, b = 1, col = "green", lty = 2)

legend("topleft",
  legend = c("slope 1"),
  col = c("green"),
  lty = c(2),
  lwd = c(2)
)


plot(data_scale_T$FDR,
  data_scale_T$qvalue,
  pch = ".",
  col = "black",
  main = "Tumor Part: Each Probe with One Sig-SNP",
  xlab = "BH",
  ylab = "Qvalue",
  xlim = c(0, 0.02)
)
abline(a = 0, b = 1, col = "green", lty = 2)

legend("topleft",
  legend = c("slope 1"),
  col = c("green"),
  lty = c(2),
  lwd = c(2)
)

# 關閉檔案
dev.off()

t.test(gt_N$FDR, qvalue1$qvalues, paired = TRUE, alternative = "two.sided")
t.test(x$FDR, x_q$qvalues, paired = TRUE, alternative = "two.sided")





file_name <- c(
  "maf_ds_N_MOSTpvalue_FDR.txt",
  "maf_ds_T_MOSTpvalue_FDR.txt",
  "maf_gt_N_MOSTpvalue_FDR.txt",
  "maf_gt_T_MOSTpvalue_FDR.txt",
  "maf_ds_T_pvalue_FDR.txt",
  "maf_ds_N_pvalue_FDR.txt",
  "maf_gt_N_pvalue_FDR.txt",
  "maf_gt_T_pvalue_FDR.txt"
)


for (i in file_name) {
  a <- paste0("D:/oral_cancer/expression/trash/", i) %>%
    fread(header = T)

  qvalue <- qvalue(a$`p-value`)
  t.test(a$FDR, qvalue$qvalues, paired = TRUE, alternative = "two.sided") %>%
    print()
}

inputname <- "D:/oral_cancer/expression/trash/maf_gt_N_pvalue.txt"
x <- fread(inputname, header = T)
qvalue_1 <- qvalue(x$`p-value`)
# pi0= 0.996

inputname <- "D:/oral_cancer/expression/trash/maf_gt_T_pvalue.txt"
x <- fread(inputname, header = T)
qvalue_1 <- qvalue(x$`p-value`)
# pi0= 0.993