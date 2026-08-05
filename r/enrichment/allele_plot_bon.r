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




# plot_gt_exp_box ----


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


# n ----
# gt data for rawData, FIXpeople eQTL pass bon asso. snp
a <- fread("C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_0.8/outcome/OQN_N_FDR_R2_0.8.txt") %>%
  filter(sig_pval_Bonfi == 1) %>%
  setkey(`p-value`)

snp_bon <- a[, .SD[1], by = Gene] %>%
  dplyr::select(CHR, SNP)

a <- fread("C:/Peter/rawData_eQTL/r2_filter_0.8/outcome/raw_N_FDR_R2_0.8.txt") %>%
  filter(sig_pval_Bonfi == 1) %>%
  setkey(`p-value`) %>%
  .[, .SD[1], by = Gene] %>%
  dplyr::select(CHR, SNP)

snp_bon <- rbind(snp_bon, a) %>% unique()
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
gt_N <- fread("C:/Peter/rawData_eQTL/r2_filter_0.8/outcome/raw_N_FDR_R2_0.8.txt")
gt_N <- gt_N[SNP %in% final$ID, c("SNP", "ALT", "REF")] %>% setkey(SNP)


# exp data
exp <- fread("C:/Peter/oral_cancer/expression/expression_data/ExpRes11sv20120601Gene.txt", header = T)
exp_N <- exp[, c(1, 2, 43:82)]
names(exp_N) <- c("Gene", "PROBE_ID", sprintf("0%dB", 1:9), sprintf("%dB", 10:40))

OQN_exp_N <- fread("C:/Peter/OQN_FIXpeople_before_eQTL/trash/OQN_ExpRes11sv20120601Gene_N.txt")


# we interest
probe_diff <- c("ILMN_1798177", "ILMN_1809147", "ILMN_2130441")
snp_diff <- c("14:64465196", "22:44101417", "6:30034171")

out_dir <- "C:/Peter/rawData_eQTL/exp_different"

for (i in seq_along(snp_diff)) {
  p <- plot_gt_exp_box(
    snp_id = snp_diff[i],
    probe_id = probe_diff[i],
    gt_dt = final,
    allele_dt = gt_N,
    exp_dt_list = list(
      `Raw eQTL` = exp_N,
      `OQN FIXpeople eQTL` = OQN_exp_N
    ),
    exp_dt_name = "Raw eQTL",
    title_plot = sprintf("Normal part eQTL %s %s",snp_diff[i],probe_diff[i]  )
  )

  ggsave(
    filename = file.path(out_dir, sprintf("BON_N_compare_%s.png", i)),
    plot = p,
    width = 8,
    height = 5,
    dpi = 300
  )
}


# t ----
# get alt, ref
# gt data for rawData, FIXpeople eQTL pass bon asso. snp
a <- fread("C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_0.8/outcome/OQN_T_FDR_R2_0.8.txt") %>%
  filter(sig_pval_Bonfi == 1) %>%
  setkey(`p-value`)

snp_bon <- a[, .SD[1], by = Gene] %>%
  dplyr::select(CHR, SNP)

a <- fread("C:/Peter/rawData_eQTL/r2_filter_0.8/outcome/raw_T_FDR_R2_0.8.txt") %>%
  filter(sig_pval_Bonfi == 1) %>%
  setkey(`p-value`) %>%
  .[, .SD[1], by = Gene] %>%
  dplyr::select(CHR, SNP)

snp_bon <- rbind(snp_bon, a) %>% unique()
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
exp_T <- exp[, c(1, 2, 3:42)]
names(exp_T) <- c("Gene", "PROBE_ID", sprintf("0%dB", 1:9), sprintf("%dB", 10:40))

OQN_exp_T <- fread("C:/Peter/OQN_FIXpeople_before_eQTL/trash/OQN_ExpRes11sv20120601Gene_T.txt")


# we interest
probe_diff <- c("ILMN_1798177")
snp_diff <- c("14:64444933")

out_dir <- "C:/Peter/rawData_eQTL/exp_different"

for (i in seq_along(snp_diff)) {
  p <- plot_gt_exp_box(
    snp_id = snp_diff[i],
    probe_id = probe_diff[i],
    gt_dt = final,
    allele_dt = gt_T,
    exp_dt_list = list(
      `Raw eQTL` = exp_T,
      `OQN FIXpeople eQTL` = OQN_exp_T
    ),
    exp_dt_name = "Raw eQTL",
    title_plot = sprintf("Tumor part eQTL %s %s",snp_diff[i],probe_diff[i]  )
  )

  ggsave(
    filename = file.path(out_dir, sprintf("BON_T_compare_%s.png", i)),
    plot = p,
    width = 8,
    height = 5,
    dpi = 300
  )
}























# N, T 同時跑

{

plot_gt_exp_box <- function(snp_id, probe_id, gt_dt, allele_dt, exp_dt_list, title_plot) {
  if (!is.list(exp_dt_list) || is.null(names(exp_dt_list)) || any(names(exp_dt_list) == "")) {
    stop("exp_dt_list must be a named list, e.g. list(`Raw eQTL` = exp_N, `OQN FIXpeople eQTL` = OQN_exp_N)")
  }

  # 對list 的 exp_N, OQN_exp_N 資料砍掉 "Gene", "PROBE_ID" 變數後，取交集
  sample_cols <- Reduce(intersect, lapply(exp_dt_list, function(x) {
    setdiff(names(x), c("Gene", "PROBE_ID"))
  }))
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

  # expression: wide -> long for each data source
  exp_long <- rbindlist(lapply(names(exp_dt_list), function(data_type) {
    exp_dt_list[[data_type]][PROBE_ID == probe_id] %>%
      melt(
        id.vars = c("Gene", "PROBE_ID"),
        measure.vars = sample_cols,
        variable.name = "sample",
        value.name = "expression"
      ) %>%
      as.data.table() %>%
      .[, data_type := data_type]
  }), use.names = TRUE, fill = TRUE)

  exp_long[, data_type := factor(data_type, levels = names(exp_dt_list))]

  plot_dt <- merge(exp_long, gt_long, by = "sample")

  datatype_levels <- levels(plot_dt$data_type)

  data_type_colors <- c("#f7500e", "#27b30e")
  names(data_type_colors) <- names(exp_dt_list)


  ggplot(plot_dt, aes(x = genotype, y = expression, fill = data_type)) +
    geom_boxplot(
      aes(group = interaction(genotype, data_type), alpha = data_type),
      width = 0.55,
      outlier.shape = NA,
      position = position_dodge(width = 0.75)
    ) +
    geom_point(
      position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.75),
      size = 2,
      alpha = 0.8
    ) +
    scale_x_discrete(labels = x_labels, drop = FALSE) +
    scale_fill_manual(values = data_type_colors, drop = FALSE) +
    scale_alpha_manual(values = seq(0.55, 0.85, length.out = length(exp_dt_list))) +
    labs(
      title = title_plot,
      x = "Genotype",
      y = "Expression",
      fill = "Data",
      alpha = "Data"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 11),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
    )
}




# n ----
# gt data for rawData, FIXpeople eQTL pass bon asso. snp
a <- fread("C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_0.8/outcome/OQN_N_FDR_R2_0.8.txt") %>%
  filter(sig_pval_Bonfi == 1) %>%
  setkey(`p-value`)

snp_bon <- a[, .SD[1], by = Gene] %>%
  dplyr::select(CHR, SNP)

a <- fread("C:/Peter/rawData_eQTL/r2_filter_0.8/outcome/raw_N_FDR_R2_0.8.txt") %>%
  filter(sig_pval_Bonfi == 1) %>%
  setkey(`p-value`) %>%
  .[, .SD[1], by = Gene] %>%
  dplyr::select(CHR, SNP)

snp_bon <- rbind(snp_bon, a) %>% unique()
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
gt_N <- fread("C:/Peter/rawData_eQTL/r2_filter_0.8/outcome/raw_N_FDR_R2_0.8.txt")
gt_N <- gt_N[SNP %in% final$ID, c("SNP", "ALT", "REF")] %>% setkey(SNP)


# exp data
exp <- fread("C:/Peter/oral_cancer/expression/expression_data/ExpRes11sv20120601Gene.txt", header = T)
exp_N <- exp[, c(1, 2, 43:82)]
names(exp_N) <- c("Gene", "PROBE_ID", sprintf("0%dB", 1:9), sprintf("%dB", 10:40))

exp <- fread("C:/Peter/oral_cancer/expression/expression_data/ExpRes11sv20120601Gene.txt", header = T)
exp_T <- exp[, c(1, 2, 3:42)]
names(exp_T) <- c("Gene", "PROBE_ID", sprintf("0%dB", 1:9), sprintf("%dB", 10:40))


# we interest
probe_diff <- c("ILMN_1798177", "ILMN_1809147", "ILMN_2130441")
snp_diff <- c("14:64465196", "22:44101417", "6:30034171")


out_dir <- "C:/Peter/rawData_eQTL/exp_different"

for (i in seq_along(snp_diff)) {
  p <- plot_gt_exp_box(
    snp_id = snp_diff[i],
    probe_id = probe_diff[i],
    gt_dt = final,
    allele_dt = gt_N,
    exp_dt_list = list(
      `Normal Part` = exp_N,
      `Tumor Part` = exp_T
    ),
    title_plot = sprintf("Expression of Probe %s by SNP %s",probe_diff[i],snp_diff[i])
  )

  ggsave(
    filename = file.path(out_dir, sprintf("BON_NT_compare_%s.png", i)),
    plot = p,
    width = 8,
    height = 5,
    dpi = 300
  )
}







# t ----
# gt data for rawData, FIXpeople eQTL pass bon asso. snp
a <- fread("C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_0.8/outcome/OQN_T_FDR_R2_0.8.txt") %>%
  filter(sig_pval_Bonfi == 1) %>%
  setkey(`p-value`)

snp_bon <- a[, .SD[1], by = Gene] %>%
  dplyr::select(CHR, SNP)

a <- fread("C:/Peter/rawData_eQTL/r2_filter_0.8/outcome/raw_T_FDR_R2_0.8.txt") %>%
  filter(sig_pval_Bonfi == 1) %>%
  setkey(`p-value`) %>%
  .[, .SD[1], by = Gene] %>%
  dplyr::select(CHR, SNP)

snp_bon <- rbind(snp_bon, a) %>% unique()
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


# we interest
probe_diff <- c("ILMN_1798177")
snp_diff <- c("14:64444933")

out_dir <- "C:/Peter/rawData_eQTL/exp_different"

for (i in seq_along(snp_diff)) {
  p <- plot_gt_exp_box(
    snp_id = snp_diff[i],
    probe_id = probe_diff[i],
    gt_dt = final,
    allele_dt = gt_T,
    exp_dt_list = list(
      `Normal Part` = exp_N,
      `Tumor Part` = exp_T
    ),
    title_plot = sprintf("Expression of Probe %s by SNP %s",probe_diff[i],snp_diff[i])
  )

  ggsave(
    filename = file.path(out_dir, sprintf("BON_NT_compare_%s-1.png", i)),
    plot = p,
    width = 8,
    height = 5,
    dpi = 300
  )
}



}
