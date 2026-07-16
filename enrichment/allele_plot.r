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



# d.e. in OQN raw data plot ----


d.e.raw <- fread("C:/Peter/rawData_eQTL/outcome/raw_exp_different.txt")
d.e.OQN <- fread("C:/Peter/OQN_FIXpeople_before_eQTL/outcome/OQN_exp_different.txt")
# setdiff(d.e.OQN[sig_OQN_Bonfi==1,PROBE_ID], d.e.raw[sig_raw_Bonfi==1,PROBE_ID])
# setdiff(d.e.raw[sig_raw_Bonfi==1,PROBE_ID], d.e.OQN[sig_OQN_Bonfi==1,PROBE_ID])



exp <- fread("C:/Peter/oral_cancer/expression/expression_data/ExpRes11sv20120601Gene.txt", header = T)
exp_N <- exp[, c(1, 2, 43:82)]
names(exp_N) <- c("Gene", "PROBE_ID", sprintf("0%dB", 1:9), sprintf("%dB", 10:40))
exp_T <- exp[, c(1, 2, 3:42)]
names(exp_T) <- c("Gene", "PROBE_ID", sprintf("0%dB", 1:9), sprintf("%dB", 10:40))

sample_cols <- setdiff(names(exp_N), c("Gene", "PROBE_ID"))

n_probe <- c("ILMN_2186216", "ILMN_2404850", "ILMN_2130441", "ILMN_1798177", "ILMN_1733103")
n_probe <- setdiff(d.e.raw[sig_raw_Bonfi == 1, PROBE_ID], d.e.OQN[sig_OQN_Bonfi == 1, PROBE_ID])[1:10]
diff_dt_raw <- exp_N[, -c("Gene", "PROBE_ID")] - exp_T[, -c("Gene", "PROBE_ID")]
diff_dt_raw <- cbind(exp_N[, c("Gene", "PROBE_ID")], diff_dt_raw) %>%
  filter(PROBE_ID %in% n_probe)

OQN_exp_N <- fread("C:/Peter/OQN_FIXpeople_before_eQTL/trash/OQN_ExpRes11sv20120601Gene_N.txt")
OQN_exp_T <- fread("C:/Peter/OQN_FIXpeople_before_eQTL/trash/OQN_ExpRes11sv20120601Gene_T.txt")



diff_dt <- OQN_exp_N[, -c("Gene", "PROBE_ID")] - OQN_exp_T[, -c("Gene", "PROBE_ID")]
diff_dt <- cbind(OQN_exp_N[, c("Gene", "PROBE_ID")], diff_dt) %>%
  filter(PROBE_ID %in% n_probe)


diff_dt_raw[, source := "raw"]
diff_dt[, source := "OQN_FIXpeople"]

combine_dt <- rbind(diff_dt, diff_dt_raw)

# measure.vars 拉下來當 variable.name 欄位的值
plot_dt <- melt(
  combine_dt,
  id.vars = c("Gene", "PROBE_ID", "source"),
  measure.vars = sample_cols,
  variable.name = "sample",
  value.name = "diff_N_minus_T"
)


# match 找第一個向量的元素在第二個向量中首次出現的位置
plot_dt[, sample_index := match(sample, sample_cols)]

for (i in seq_along(n_probe)) {
  png(sprintf("C:/Peter/rawData_eQTL/exp_different/probe_diff_%s.png", i),
    width = 1600, height = 1200, res = 200
  )

  p <- ggplot(plot_dt %>% filter(PROBE_ID %in% n_probe[i]), aes(
    x = sample_index,
    y = diff_N_minus_T,
    color = source,
    group = source
  )) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.8) +
    scale_x_continuous(
      breaks = seq_along(sample_cols),
      labels = sample_cols
    ) +
    labs(
      x = "sample index",
      y = "N-T",
      color = "source",
      title = n_probe[i]
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1)
    )

  print(p)
  dev.off()
}


# plot_gt_exp_box ----


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
      )               %>%
      as.data.table() %>%
      .[, data_type := data_type]
  }), use.names = TRUE, fill = TRUE)

  exp_long[, data_type := factor(data_type, levels = names(exp_dt_list))]

  plot_dt <- merge(exp_long, gt_long, by = "sample")

  datatype_levels <- levels(plot_dt$data_type)

  data_type_colors <- c("#e95c0a", "#59A14F")
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
      title = paste0(title_plot, "   ", snp_id, " / ", probe_id),
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
  filter(sig_pval_Bonfi == 1)                                                               %>%
  setkey(`p-value`)

snp_bon <- a[, .SD[1], by = Gene] %>%
  select(CHR, SNP)

a <- fread("C:/Peter/rawData_eQTL/r2_filter_0.8/outcome/raw_N_FDR_R2_0.8.txt") %>%
  filter(sig_pval_Bonfi == 1)                                                  %>%
  setkey(`p-value`)                                                            %>%
  .[, .SD[1], by = Gene]                                                       %>%
  select(CHR, SNP)

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
probe_diff <- c("ILMN_2323979", "ILMN_2099745", "ILMN_2130441", "ILMN_2186216", "ILMN_2215025", "ILMN_1693323", "ILMN_2404850")
snp_diff <- c("1:119254023", "11:128310767", "6:30034171", "3:169211728", "1:149016537", "6:88083637", "3:40443900")

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
    title_plot = "N eQTL"
  )

  ggsave(
    filename = file.path(out_dir, sprintf("box_N_compare_%s.png", i)),
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
  filter(sig_pval_Bonfi == 1)                                                               %>%
  setkey(`p-value`)

snp_bon <- a[, .SD[1], by = Gene] %>%
  select(CHR, SNP)

a <- fread("C:/Peter/rawData_eQTL/r2_filter_0.8/outcome/raw_T_FDR_R2_0.8.txt") %>%
  filter(sig_pval_Bonfi == 1)                                                  %>%
  setkey(`p-value`)                                                            %>%
  .[, .SD[1], by = Gene]                                                       %>%
  select(CHR, SNP)

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
probe_diff <- c("ILMN_2323979", "ILMN_1798177", "ILMN_1794767", "ILMN_2130441", "ILMN_2323979", "ILMN_1733103")
snp_diff <- c("1:119346137", "14:64444933", "11:5812440", "6:29889028", "1:119254023", "7:65452244")

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
    title_plot = "T eQTL"
  )

  ggsave(
    filename = file.path(out_dir, sprintf("box_T_compare_%s.png", i)),
    plot = p,
    width = 8,
    height = 5,
    dpi = 300
  )
}









# OQN, raw exp 不同? ----

# exp data
exp <- fread("C:/Peter/oral_cancer/expression/expression_data/ExpRes11sv20120601Gene.txt", header = T)
exp_N <- exp[, c(1, 2, 43:82)]
names(exp_N) <- c("Gene", "PROBE_ID", sprintf("0%dB", 1:9), sprintf("%dB", 10:40))

OQN_exp_N <- fread("C:/Peter/OQN_FIXpeople_before_eQTL/trash/OQN_ExpRes11sv20120601Gene_N.txt")

sample_cols <- c(sprintf("0%dB", 1:9), sprintf("%dB", 10:40))


probe_diff <- c(
  "ILMN_2323979", "ILMN_2099745", "ILMN_2130441", "ILMN_2186216",
  "ILMN_2215025", "ILMN_1693323", "ILMN_2404850"
)

for (i in seq_along(probe_diff)) {
  exp_sub <- exp_N[PROBE_ID %in% probe_diff[i]]
  exp_sub_raw <- melt(
    exp_sub,
    id.vars = c("PROBE_ID"),
    measure.vars = sample_cols,
    variable.name = "sample",
    value.name = "raw_N"
  )

  exp_sub <- OQN_exp_N[PROBE_ID %in% probe_diff[i]]
  exp_sub_oqn <- melt(
    exp_sub,
    id.vars = c("PROBE_ID"),
    measure.vars = sample_cols,
    variable.name = "sample",
    value.name = "OQN_N"
  )
  exp_sub <- exp_sub_oqn[exp_sub_raw, on = .(PROBE_ID, sample)]
  x_cutoff <- c(exp_sub$OQN_N, exp_sub$raw_N) %>%
    abs()                                     %>%
    max()

  png(sprintf("C:/Peter/rawData_eQTL/exp_different/exp_compare_%s.png", i), width = 1000, height = 1000, res = 200)
  plot(exp_sub$OQN_N %>% abs(),
    exp_sub$raw_N %>% abs(),
    pch = 20,
    col = "black",
    main = paste0(probe_diff[i], " N Part Expression"),
    xlab = "abs OQN exp",
    ylab = "abs Raw exp",
    xlim = c(0, x_cutoff),
    ylim = c(0, x_cutoff)
  )
  abline(a = 0, b = 1, col = "green", lty = 2)

  legend("topleft",
    legend = c("slope 1"),
    col = c("green"),
    lty = c(2),
    lwd = c(2)
  )

  dev.off()
}


# pval data ----
d.e.oqn <- fread("C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_0.8/outcome/OQN_exp_different_r2_0.8.txt")
d.e.raw <- fread("C:/Peter/rawData_eQTL/r2_filter_0.8/outcome/raw_exp_different_r2_0.8.txt")

probe_diff <- c(
  "ILMN_2323979", "ILMN_2099745", "ILMN_2130441", "ILMN_2186216",
  "ILMN_2215025", "ILMN_1693323", "ILMN_2404850"
)


a <- d.e.oqn[PROBE_ID %in% probe_diff, c("PROBE_ID", "pval_OQN")]
b <- d.e.raw[PROBE_ID %in% probe_diff, c("PROBE_ID", "pval_raw")]
pval_compare <- a[b, on = .(PROBE_ID)]
pval_compare[, pval_OQN := -log10(pval_OQN)]
pval_compare[, pval_raw := -log10(pval_raw)]

x_max <- c(pval_compare$pval_OQN, pval_compare$pval_raw) %>%
  abs()                                                  %>%
  max()
x_min <- c(pval_compare$pval_OQN, pval_compare$pval_raw) %>%
  abs()                                                  %>%
  min()


png("C:/Peter/rawData_eQTL/exp_different/pval_compare_N.png", width = 1000, height = 1000, res = 200)

plot(pval_compare$pval_OQN,
  pval_compare$pval_raw,
  pch = 20,
  xaxt = "n",
  col = "black",
  main = "pval Compare",
  xlab = "- log10 OQN pval",
  ylab = "- log10 Raw pval",
  xlim = c(x_min, x_max),
  ylim = c(x_min, x_max)
)
axis(1, at = pval_compare$pval_OQN, labels = FALSE)
abline(a = 0, b = 1, col = "green", lty = 2)

text(
  x = pval_compare$pval_OQN,
  y = par("usr")[3],
  labels = pval_compare$PROBE_ID,
  srt = 45,
  adj = 1,
  xpd = TRUE,
  cex = 0.7
)

legend("topleft",
  legend = c("slope 1"),
  col = c("green"),
  lty = c(2),
  lwd = c(2)
)

dev.off()


# all pval ----
d.e.oqn <- fread("C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_0.8/outcome/OQN_exp_different_r2_0.8.txt")
d.e.raw <- fread("C:/Peter/rawData_eQTL/r2_filter_0.8/outcome/raw_exp_different_r2_0.8.txt")
names(d.e.oqn)

a <- d.e.oqn[sig_OQN_Bonfi == 1, c("PROBE_ID", "pval_OQN")]
b <- d.e.raw[sig_raw_Bonfi == 1, c("PROBE_ID", "pval_raw")]
pval_compare <- a[b, on = .(PROBE_ID)]
pval_compare[, pval_OQN := -log10(pval_OQN)]
pval_compare[, pval_raw := -log10(pval_raw)]

all_pval <- c(pval_compare$pval_OQN, pval_compare$pval_raw)

x_max <- all_pval[!is.na(all_pval)] %>%
  abs()                             %>%
  max()
x_min <- all_pval[!is.na(all_pval)] %>%
  abs()                             %>%
  min()

uniqueN(a$PROBE_ID)
uniqueN(b$PROBE_ID)
uniqueN(intersect(a$PROBE_ID, b$PROBE_ID))


png("C:/Peter/rawData_eQTL/exp_different/pval_compare_all_log.png", width = 1000, height = 1000, res = 200)

plot(pval_compare$pval_OQN,
  pval_compare$pval_raw,
  pch = ".",
  col = "black",
  main = "pass bon probe pval Compare",
  xlab = "- log10 OQN pval",
  ylab = "- log10 Raw pval",
  xlim = c(x_min, x_max),
  ylim = c(x_min, x_max)
)
abline(a = 0, b = 1, col = "green", lty = 2)


legend("topleft",
  legend = c("slope 1"),
  col = c("green"),
  lty = c(2),
  lwd = c(2)
)

dev.off()



# all  hist ----
a <- d.e.oqn[, c("PROBE_ID", "pval_OQN")]
b <- d.e.raw[, c("PROBE_ID", "pval_raw")]
pval_compare <- a[b, on = .(PROBE_ID)]
pval_compare[, diff := pval_OQN - pval_raw]
png("C:/Peter/rawData_eQTL/exp_different/pval_compare.png", width = 1000, height = 1000, res = 200)
ggplot(pval_compare, aes(x = diff)) +
  geom_histogram(bins = 1000, fill = "#4993b1", color = "white") +
  coord_cartesian(xlim = c(-0.02, 0.02)) +
  labs(
    title = "Distribution of Difference in p-value Significance",
    x = expression(-log[10](p[OQN]) - -log[10](p[Raw])),
    y = "Count"
  ) +
  theme_minimal()
dev.off()

range(pval_compare$diff)



x_range <- c(-0.2, 0.2)
png("C:/Peter/rawData_eQTL/exp_different/pval_compare_scale.png", width = 1000, height = 1000, res = 200)
ggplot(pval_compare, aes(x = diff)) +
  geom_histogram(bins = 300, fill = "#4993b1", color = "white") +
  coord_cartesian(xlim = x_range) +
  labs(
    title = "Distribution of Difference in p-value Significance",
    x = expression(-log[10](p[OQN]) - -log[10](p[Raw])),
    y = "Count"
  ) +
  theme_minimal()
dev.off()

#


# eQTL pval ----
# 1. 每個 probe 選最顯著的 snp eQTL pval 畫
gt_N <- fread("C:/Peter/rawData_eQTL/r2_filter_0.8/outcome/raw_N_FDR_R2_0.8.txt")[
  ,
  .(Probe, Gene, SNP_raw = SNP, pval_raw = `p-value`)
]
gt_N_oqn <- fread("C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_0.8/outcome/OQN_N_FDR_R2_0.8.txt")[
  ,
  .(Probe, Gene, SNP_OQN = SNP, pval_OQN = `p-value`)
]

setkey(gt_N, pval_raw)
gt_N <- gt_N[, .SD[1], by = Gene]
setkey(gt_N_oqn, pval_OQN)
gt_N_oqn <- gt_N_oqn[, .SD[1], by = Gene]
gt_N_combine <- gt_N_oqn[gt_N, on = .(Probe, Gene), nomatch = 0]

all_pval <- c(
  gt_N_combine$pval_OQN,
  gt_N_combine$pval_raw
)
x_max <- all_pval[!is.na(all_pval)] %>%
  abs()                             %>%
  max()
x_min <- all_pval[!is.na(all_pval)] %>%
  abs()                             %>%
  min()

png("C:/Peter/rawData_eQTL/exp_different/eQTLpval_compare_fdr.png", width = 1000, height = 1000, res = 200)


plot(gt_N_combine$pval_OQN,
  gt_N_combine$pval_raw,
  pch = 20,
  col = "black",
  main = "pass FDR probe pval Compare",
  xlab = "eQTL OQN pval",
  ylab = "eQTL Raw pval",
  xlim = c(x_min, x_max),
  ylim = c(x_min, x_max)
)
abline(a = 0, b = 1, col = "green", lty = 2)


legend("topleft",
  legend = c("slope 1"),
  col = c("green"),
  lty = c(2),
  lwd = c(2)
)

dev.off()

# log 10 transform ----

gt_N_combine[, pval_OQN := -log10(pval_OQN)]
gt_N_combine[, pval_raw := -log10(pval_raw)]

all_pval <- c(
  gt_N_combine$pval_OQN,
  gt_N_combine$pval_raw
)
x_max <- all_pval[!is.na(all_pval)] %>%
  abs()                             %>%
  max()
x_min <- all_pval[!is.na(all_pval)] %>%
  abs()                             %>%
  min()

png("C:/Peter/rawData_eQTL/exp_different/eQTLpval_compare_fdr_log.png", width = 1000, height = 1000, res = 200)


plot(gt_N_combine$pval_OQN,
  gt_N_combine$pval_raw,
  pch = 20,
  col = "black",
  main = "pass FDR probe pval Compare",
  xlab = "eQTL OQN -log10(pval)",
  ylab = "eQTL Raw -log10(pval)",
  xlim = c(x_min, x_max),
  ylim = c(x_min, x_max)
)
abline(a = 0, b = 1, col = "green", lty = 2)


legend("topleft",
  legend = c("slope 1"),
  col = c("green"),
  lty = c(2),
  lwd = c(2)
)

dev.off()

# 2. 每個過 bon probe 選最顯著的 snp eQTL pval 畫 ----
gt_N <- fread("C:/Peter/rawData_eQTL/r2_filter_0.8/outcome/raw_N_FDR_R2_0.8.txt") %>%
  filter(sig_pval_Bonfi == 1)
gt_N_oqn <- fread("C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_0.8/outcome/OQN_N_FDR_R2_0.8.txt") %>%
  filter(sig_pval_Bonfi == 1)


gt_N <- gt_N[
  ,
  .(Probe, Gene, SNP_raw = SNP, pval_raw = `p-value`)
]

gt_N_oqn <- gt_N_oqn[
  ,
  .(Probe, Gene, SNP_OQN = SNP, pval_OQN = `p-value`)
]

setkey(gt_N, pval_raw)
gt_N <- gt_N[, .SD[1], by = Gene]
setkey(gt_N_oqn, pval_OQN)
gt_N_oqn <- gt_N_oqn[, .SD[1], by = Gene]
gt_N_combine <- gt_N_oqn[gt_N, on = .(Probe, Gene), nomatch = 0]

all_pval <- c(
  gt_N_combine$pval_OQN,
  gt_N_combine$pval_raw
)
x_max <- all_pval[!is.na(all_pval)] %>%
  abs()                             %>%
  max()
x_min <- all_pval[!is.na(all_pval)] %>%
  abs()                             %>%
  min()

png("C:/Peter/rawData_eQTL/exp_different/eQTLpval_compare_bon.png", width = 1000, height = 1000, res = 200)


plot(gt_N_combine$pval_OQN,
  gt_N_combine$pval_raw,
  pch = 20,
  col = "black",
  main = "pass bon probe pval Compare",
  xlab = "eQTL OQN pval",
  ylab = "eQTL Raw pval",
  xlim = c(x_min, x_max),
  ylim = c(x_min, x_max)
)
abline(a = 0, b = 1, col = "green", lty = 2)


legend("topleft",
  legend = c("slope 1"),
  col = c("green"),
  lty = c(2),
  lwd = c(2)
)

dev.off()







# 跑出回歸線 ----

snp_id <- snp_diff[i]
probe_id <- probe_diff[i]
gt_dt <- final
allele_dt <- gt_N
exp_dt_list <- list(
  `Raw eQTL` = exp_N,
  `OQN FIXpeople eQTL` = OQN_exp_N
)
title_plot <- "N eQTL"

exp_long <- rbindlist(lapply(names(exp_dt_list), function(data_type) {
  exp_dt_list[[data_type]][PROBE_ID == probe_id] %>%
    melt(
      id.vars = c("Gene", "PROBE_ID"),
      measure.vars = sample_cols,
      variable.name = "sample",
      value.name = "expression"
    )               %>%
    as.data.table() %>%
    .[, data_type := data_type]
}), use.names = TRUE, fill = TRUE)

exp_long[, data_type := factor(data_type, levels = names(exp_dt_list))]



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

plot_dt <- merge(exp_long, gt_long, by = "sample")
plot_dt_sub <- plot_dt[data_type %in% c("Raw eQTL", "OQN FIXpeople eQTL")]
data_type_cols <- c(
  `Raw eQTL` = "#e95c0a",
  `OQN FIXpeople eQTL` = "#59A14F"
)

windows()
plot(
  plot_dt_sub$gt_code,
  plot_dt_sub$expression,
  col = data_type_cols[as.character(plot_dt_sub$data_type)],
  pch = 19,
  main = paste0(title_plot, "   ", snp_id, " / ", probe_id),
  xlab = "genotype code",
  ylab = "expression"
)
for (data_type_i in names(data_type_cols)) {
  plot_dt_i <- plot_dt_sub[data_type == data_type_i]
  abline(lm(expression ~ gt_code, data = plot_dt_i), col = data_type_cols[data_type_i], lwd = 2)
}
legend(
  "topright",
  legend = names(data_type_cols),
  col = data_type_cols,
  pch = 19,
  lwd = 2,
  bty = "n"
)

plot_gt_exp_lm <- function(
  snp_id,
  probe_id,
  gt_dt,
  allele_dt,
  exp_dt_list,
  title_plot = "N eQTL",
  data_type_cols = c(`Raw eQTL` = "#e95c0a", `OQN FIXpeople eQTL` = "#59A14F"),
  open_window = TRUE
) {
  if (!is.list(exp_dt_list) || is.null(names(exp_dt_list)) || any(names(exp_dt_list) == "")) {
    stop("exp_dt_list must be a named list, e.g. list(`Raw eQTL` = exp_N, `OQN FIXpeople eQTL` = OQN_exp_N)")
  }

  sample_cols <- Reduce(intersect, lapply(exp_dt_list, function(x) {
    setdiff(names(x), c("Gene", "PROBE_ID"))
  }))
  sample_cols <- intersect(sample_cols, setdiff(names(gt_dt), "ID"))

  exp_long <- rbindlist(lapply(names(exp_dt_list), function(data_type) {
    exp_dt_list[[data_type]][PROBE_ID == probe_id] %>%
      melt(
        id.vars = c("Gene", "PROBE_ID"),
        measure.vars = sample_cols,
        variable.name = "sample",
        value.name = "expression"
      )               %>%
      as.data.table() %>%
      .[, data_type := data_type]
  }), use.names = TRUE, fill = TRUE)

  exp_long[, data_type := factor(data_type, levels = names(exp_dt_list))]

  gt_long <- gt_dt[ID == snp_id] %>%
    melt(
      id.vars = "ID",
      measure.vars = sample_cols,
      variable.name = "sample",
      value.name = "gt_code"
    ) %>%
    as.data.table()

  allele_info <- allele_dt[SNP == snp_id]
  ref <- allele_info$REF[1]
  alt <- allele_info$ALT[1]

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

  plot_dt <- merge(exp_long, gt_long, by = "sample")
  plot_dt_sub <- plot_dt[data_type %in% names(data_type_cols)]

  if (open_window) {
    windows()
  }

  plot(
    plot_dt_sub$gt_code,
    plot_dt_sub$expression,
    col = data_type_cols[as.character(plot_dt_sub$data_type)],
    pch = 19,
    main = paste0(title_plot, "   ", snp_id, " / ", probe_id),
    xlab = "genotype code",
    ylab = "expression"
  )

  for (data_type_i in names(data_type_cols)) {
    plot_dt_i <- plot_dt_sub[data_type == data_type_i]
    abline(lm(expression ~ gt_code, data = plot_dt_i), col = data_type_cols[data_type_i], lwd = 2)
  }

  legend(
    "topright",
    legend = names(data_type_cols),
    col = data_type_cols,
    pch = 19,
    lwd = 2,
    bty = "n"
  )

  invisible(plot_dt_sub)
}

plot_gt_exp_lm_vec <- function(
  snp_id,
  probe_id,
  gt_dt,
  allele_dt,
  exp_dt_list,
  title_plot = "N eQTL",
  data_type_cols = c(`Raw eQTL` = "#e95c0a", `OQN FIXpeople eQTL` = "#59A14F")
) {
  if (length(snp_id) != length(probe_id)) {
    stop("snp_id and probe_id must have the same length.")
  }

  invisible(lapply(seq_along(snp_id), function(i) {
    plot_gt_exp_lm(
      snp_id = snp_id[i],
      probe_id = probe_id[i],
      gt_dt = gt_dt,
      allele_dt = allele_dt,
      exp_dt_list = exp_dt_list,
      title_plot = title_plot,
      data_type_cols = data_type_cols,
      open_window = TRUE
    )
  }))
}