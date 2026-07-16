## package ----
library(data.table)
library(dplyr)
library(tidyverse)
library(stringr)
library(ggplot2)
library(ggrepel)
library(matrixStats) # rowMins
library(qvalue)
library(bestNormalize)
library(rtracklayer)


# ft ----
plot_density <- function(dt, plot_title, part) {
  if (!is.character(plot_title) || length(plot_title) != 1 || !is.data.table(dt)) {
    stop("plot_title 必須是一個字串, my_dt 必須是一個 data.table")
  }

  ds <- dt

  bin_width <- 1 / part
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


plot_gt_exp_box <- function(snp_id, probe_id, gt_dt, allele_dt, exp_dt_list, title_plot, exp_dt_name = names(exp_dt_list)[1]) {
  if (!is.list(exp_dt_list) || is.null(names(exp_dt_list)) || any(names(exp_dt_list) == "")) {
    stop("exp_dt_list must be a named list, e.g. list(`Raw eQTL` = exp_N, `OQN FIXpeople eQTL` = OQN_exp_N)")
  }
  if (!exp_dt_name %in% names(exp_dt_list)) {
    stop("exp_dt_name must be one of: ", paste(names(exp_dt_list), collapse = ", "))
  }

  exp_dt <- exp_dt_list[[exp_dt_name]]

  # 對list 的 exp_N, OQN_exp_N 資料砍掉 "Gene", "PROBE_ID" 變數後，取交集
  sample_cols <- setdiff(names(exp_dt), c("Gene", "PROBE_ID"))
  sample_cols <- intersect(sample_cols, setdiff(names(gt_dt), "ID"))

  # measure.vars 拉下來當 variable.name 欄位的值
  gt_long <- gt_dt[ID == snp_id] %>%
    melt(
      id.vars = "ID",
      measure.vars = sample_cols,
      variable.name = "sample",
      value.name = "gt_code"
    ) %>%
    as.data.table()

  # get REF / ALT
  allele_info <- allele_dt[SNP == snp_id]
  ref <- allele_info$REF[1]
  alt <- allele_info$ALT[1]

  # 轉 gt
  gt_long[, genotype := fifelse(
    gt_code == 0, paste0(ref, ref),
    fifelse(
      gt_code == 1, paste0(ref, alt),
      fifelse(gt_code == 2, paste0(alt, alt), NA_character_)
    )
  )]


  gt_long[, genotype := factor(
    genotype,
    levels = c(paste0(ref, ref), paste0(ref, alt), paste0(alt, alt))
  )]

  genotype_levels <- levels(gt_long$genotype)
  gt_count <- gt_long[!is.na(genotype), .N, by = genotype]

  # match 找第一個向量的元素在第二個向量中首次出現的位置
  # 把 gt_count 強制順序從 gt_code=0,1,2
  gt_count <- gt_count[match(genotype_levels, genotype)]
  x_labels <- paste0(genotype_levels, "\n(n = ", fifelse(is.na(gt_count$N), 0L, gt_count$N), ")")
  names(x_labels) <- genotype_levels

  # expression: wide -> long for selected data source
  exp_long <- exp_dt[PROBE_ID == probe_id] %>%
    melt(
      id.vars = c("Gene", "PROBE_ID"),
      measure.vars = sample_cols,
      variable.name = "sample",
      value.name = "expression"
    ) %>%
    as.data.table()

  plot_dt <- merge(exp_long, gt_long, by = "sample")


  ggplot(plot_dt, aes(x = genotype, y = expression)) +
    geom_boxplot(
      width = 0.55,
      outlier.shape = NA,
      fill = "#27b30e",
      alpha = 0.7
    ) +
    geom_point(
      position = position_jitter(width = 0.12, height = 0),
      size = 2,
      alpha = 0.8
    ) +
    scale_x_discrete(labels = x_labels, drop = FALSE) +
    labs(
      title = title_plot,
      x = "Genotype",
      y = "Expression"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 11),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
    )
}


# all pval dist 圖 ----
for (i in c("N", "T")) {
  for (bin_number in c(20)) {
    eqtl <- fread(sprintf("C:/Peter/rawData_eQTL/trash/raw_maf_gt_%s_pvalue.txt", i))

    R2_filter <- fread("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/cis_snp_R2_0.8.txt", header = T)
    eqtl <- eqtl[SNP %in% R2_filter$hg18_snpID, ]

    ggsave(sprintf("C:/Peter/why_T_MOREthan_N/maf_gt_%s_pvalue_Autosome_%s.png", i, bin_number),
      plot = plot_density(
        eqtl, if (i == "T") {
          "Tumor Part eQTL P-value of Cis-SNPs in Autosome"
        } else {
          "Normal Part eQTL P-value of Cis-SNPs in Autosome"
        },
        bin_number
      ),
      width = 10, height = 5, dpi = 300
    )
  }
}


# 過 FDR dist 圖 ----
a <- fread("C:/Peter/rawData_eQTL/r2_filter_0.8/outcome/raw_maf_gt_N_pvalue_FDR_R2_0.8.txt")
b <- fread("C:/Peter/rawData_eQTL/r2_filter_0.8/outcome/raw_maf_gt_T_pvalue_FDR_R2_0.8.txt")

png("C:/Peter/why_T_MOREthan_N/gt_pvalue_passFDR.png", width = 1000, height = 2000, res = 200)
par(mfrow = c(1, 1))
hist(a$`p-value`,
  col = rgb(0.2, 0.5, 1, 0.4),
  seq(0, max(a$`p-value`, b$`p-value`), length.out = 30),
  main = "eQTL pval: pass FDR",
  xlab = "pval"
)
hist(b$`p-value`,
  col = rgb(1, 0.2, 0.2, 0.4),
  seq(0, max(a$`p-value`, b$`p-value`), length.out = 30),
  add = TRUE
)

legend(
  "topright",
  legend = c(
    sprintf("Normal part (%s pairs)", nrow(a)),
    sprintf("Tumor part (%s pairs)", nrow(b))
  ),
  fill = c(
    rgb(0.2, 0.5, 1, 0.4),
    rgb(1, 0.2, 0.2, 0.4)
  )
)

# 關閉檔案
dev.off()



# 過 bon pval dist 圖 ----
bon_cutoff <- 0.05 / 40666601
png("C:/Peter/why_T_MOREthan_N/gt_pvalue_passBON.png", width = 1000, height = 2000, res = 200)
par(mfrow = c(1, 1))
hist(a$`p-value`[a$`p-value` < bon_cutoff],
  col = rgb(0.2, 0.5, 1, 0.4),
  seq(0, max(a$`p-value`[a$`p-value` < bon_cutoff], b$`p-value`[b$`p-value` < bon_cutoff]), length.out = 20),
  main = "eQTL pval: pass Bon",
  xlab = "pval"
)
hist(b$`p-value`[b$`p-value` < bon_cutoff],
  col = rgb(1, 0.2, 0.2, 0.4),
  seq(0, max(a$`p-value`[a$`p-value` < bon_cutoff], b$`p-value`[b$`p-value` < bon_cutoff]), length.out = 20),
  add = TRUE
)

legend(
  "topright",
  legend = c(
    sprintf("Normal part (%s pairs)", length(a$`p-value`[a$`p-value` < bon_cutoff])),
    sprintf("Tumor part (%s pairs)", length(b$`p-value`[b$`p-value` < bon_cutoff]))
  ),
  fill = c(
    rgb(0.2, 0.5, 1, 0.4),
    rgb(1, 0.2, 0.2, 0.4)
  )
)

# 關閉檔案
dev.off()






a <- fread("C:/Peter/rawData_eQTL/r2_filter_0.8/outcome/raw_maf_gt_N_pvalue_FDR_R2_0.8.txt")
b <- fread("C:/Peter/rawData_eQTL/r2_filter_0.8/outcome/raw_maf_gt_T_pvalue_FDR_R2_0.8.txt")
bon_cutoff <- 0.05 / 40666601
b <- b[`p-value` < bon_cutoff]
setkey(b, gene, `p-value`)

# 挑出只有在 Tumor 過 FDR 的 BON association，選最顯著的
t_only <- setdiff(b[, .(SNP, gene)], a[, .(SNP, gene)])
t_only <- b[t_only, on = .(SNP, gene)]
t_only <- t_only[`p-value` < bon_cutoff]
setkey(t_only, gene, `p-value`)
t_only <- t_only[, .SD[1], by = gene]




# box plot ----
# gt data for rawData, FIXpeople eQTL pass bon asso. snp
snp_bon <- fread("C:/Peter/rawData_eQTL/r2_filter_0.8/outcome/raw_T_FDR_R2_0.8.txt") %>%
  filter(sig_pval_Bonfi == 1)                                                        %>%
  setkey(`p-value`)                                                                  %>%
  .[, .SD[1], by = Gene]                                                             %>%
  dplyr::select(CHR, SNP)

snp_bon[, CHR := sub(":.*", "", SNP) %>% as.numeric()]
setkey(snp_bon, CHR)

snp_list <- list()
for (chr in sort(unique(snp_bon$CHR))) {
  chr_snp <- snp_bon[CHR == chr, unique(SNP)]
  chr_file <- sprintf(
    "C:/Peter/OQN_FIXpeople_before_eQTL/cis_snp_gt_maf_chr%d.txt",
    chr
  )

  snp_list[[as.character(chr)]] <- fread(chr_file)[ID %in% chr_snp]
}

final <- rbindlist(snp_list, use.names = TRUE, fill = TRUE)



# get alt, ref
gt_T <- fread("C:/Peter/rawData_eQTL/r2_filter_0.8/outcome/raw_T_FDR_R2_0.8.txt")
gt_T <- gt_T[SNP %in% final$ID, c("SNP", "ALT", "REF")] %>% setkey(SNP)


# exp data
exp <- fread("C:/Peter/oral_cancer/expression/expression_data/ExpRes11sv20120601Gene.txt", header = T)
exp_N <- exp[, c(1, 2, 43:82)]
names(exp_N) <- c("Gene", "PROBE_ID", sprintf("0%dB", 1:9), sprintf("%dB", 10:40))

exp <- fread("C:/Peter/oral_cancer/expression/expression_data/ExpRes11sv20120601Gene.txt", header = T)
exp_T <- exp[, c(1, 2, 3:42)]
names(exp_T) <- c("Gene", "PROBE_ID", sprintf("0%dB", 1:9), sprintf("%dB", 10:40))


for (i in nrow(t_only)) {
  p <- plot_gt_exp_box(
    snp_id = t_only$SNP[i],
    probe_id = t_only$gene[i],
    gt_dt = final,
    allele_dt = gt_T,
    exp_dt_list = list(
      `Normal Part` = exp_N,
      `Tumor Part` = exp_T
    ),
    title_plot = sprintf("Expression of Probe %s by SNP %s", t_only$gene[i], t_only$SNP[i])
  )

  ggsave(
    filename = file.path("C:/Peter/why_T_MOREthan_N", sprintf("BON_NT_compare_%s-1.png", i)),
    plot = p,
    width = 8,
    height = 5,
    dpi = 300
  )
}