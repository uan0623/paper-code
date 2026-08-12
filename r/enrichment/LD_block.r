library(data.table)
library(dplyr)
library(tidyverse)
library(stringr)
library(ggrepel)
library(matrixStats) # rowMins
library(qvalue)
library(bestNormalize)
library(ggplot2)
library(openxlsx)
library(magrittr)


# 不同相關定義，gene 過fdr 門檻的數量 v.s. gene 高度相關數量 兩個值的相關性 ----


  make_group_id <- function(gap_prev, threshold = 1e6) {
  acc <- 0
  gid <- integer(length(gap_prev))
  current_group <- 1L

  for (i in seq_along(gap_prev)) {
    gap <- gap_prev[i]
    if (is.na(gap)) {
      acc <- 0
    } else {
      acc <- acc + gap
      if (acc > threshold) {
        current_group <- current_group + 1L
        acc <- 0
      }
    }
    gid[i] <- current_group
  }
  gid
}



cor_summary <- vector("list", length(seq(0,1,by=0.05)[-1]))
  for (j in seq(0,1,by=0.05)[-c(1,21)] %>% seq_along()) {
   
    cor_cutoff <- seq(0,1,by=0.05)[-1][j]
exp_N <- fread("C:/Peter/rawData_eQTL/raw_Exp_mulInterval_N.txt", header = T)
exp_N <- t(exp_N) %>% as.data.table()
names(exp_N) <- as.character(exp_N[1, ])
exp_N <- exp_N[-c(1), ]
exp_N[, (names(exp_N)) := lapply(.SD, as.numeric)]

gene_pos <- fread("D:/oral_cancer/expression/expression_data/probe_pos_mulInterval.txt", header = T)
gene_pos <- gene_pos[Gene %in% names(exp_N)]
setDT(gene_pos)

gene_pos <- gene_pos[order(CHROMOSOME, start, end)]
gene_pos[, gap_prev := start - data.table::shift(end),
  by = CHROMOSOME
]

# 距離累計 > 1e6，更新 group_id
gene_pos[, group_id := make_group_id(gap_prev), by = CHROMOSOME]
gene_pos[, group_id := paste0("chr", CHROMOSOME, "_group", group_id)]
gene_pos[, gap_prev := NULL]


ld_type <- gene_pos$group_id %>% unique()
cor_group <- vector("list", length(ld_type))
for (i in seq_along(ld_type)) {
  df <- gene_pos[group_id %in% ld_type[i]]
  group_probe <- df$Gene
  cor_dt <- cor(exp_N[, ..group_probe]) %>% 
  as.table() %>% 
  as.data.table()

  setnames(cor_dt, c("probe_1", "probe_2", "cor"))

  # 留下對角線下的 cor 且 cor 夠大的
  cor_dt <- cor_dt[probe_1< probe_2 & abs(cor) > cor_cutoff]
  cor_group[[i]] <- cor_dt
}

cor_dt <- rbindlist(cor_group)
cor_dt <- gene_pos[, .(probe_1=Gene, chr=CHROMOSOME,group_id)][cor_dt, on=.(probe_1)] 

gt_N <- fread("C:/Peter/rawData_eQTL/r2_filter_0.8/outcome/raw_N_FDR_R2_0.8.txt")
eQTL_gene <- unique(gt_N$Gene)
eQTL_probe <- unique(gt_N$Probe)


# 計算每個 group，有顯著 eQTL 的數量
gene_pos[,sig_FDR := 0]
gene_pos[Gene %in% eQTL_probe, sig_FDR := 1]
group_count <- gene_pos[, .(sig_FDR_sum = sum(sig_FDR, na.rm = TRUE)),
  by = group_id
]

# 挑出 cor>cor_cutoff 的 probe
cor_pass_probe <- gene_pos[Gene %in% union(cor_dt$probe_1,cor_dt$probe_2)]
record_sig <- cor_pass_probe[,.(High_Cor_probe_number=.N),by=group_id][group_count, on=.(group_id)]
record_sig[is.na(High_Cor_probe_number), High_Cor_probe_number := 0]
fit <- lm(High_Cor_probe_number ~ sig_FDR_sum, data = record_sig)
cor_result <- cor(record_sig$sig_FDR_sum, record_sig$High_Cor_probe_number)

# 算 Spearman cor
spearman <- cor.test(
    record_sig$High_Cor_probe_number,
    record_sig$sig_FDR_sum,
    method = "spearman"
)


# poisson reg
fit_poisson <- glm(
  sig_FDR_sum ~ High_Cor_probe_number,
  family = poisson(link = "log"),
  data = record_sig
)

summary(fit_poisson)$coefficients["High_Cor_probe_number", "Pr(>|z|)"]
sum(residuals(fit_poisson, type = "pearson")^2) / fit_poisson$df.residual


cor_summary[[j]] <- data.table(
  probe_cor_cutoff = cor_cutoff,
  cor = cor_result,
  reg_pval = summary(fit)$coefficients[2, 4],
  spearman_pval = spearman$p.value,
  spearman_rho = spearman$estimate %>% as.numeric(),
  poisson_pval = summary(fit_poisson)$coefficients["High_Cor_probe_number", "Pr(>|z|)"],
  poisson_beta = summary(fit_poisson)$coefficients["High_Cor_probe_number", "Estimate"],
  dispersion = sum(residuals(fit_poisson, type = "pearson")^2) / fit_poisson$df.residual
)
cor_summary[[j]][, overdispersion := ifelse(dispersion>1,1,0)]






}

cor_final <- rbindlist(cor_summary)
fwrite(cor_final,
  "C:/Peter/gene_enrichment/code_project/outcome/ld_bootstrap/cor_final.txt",
  row.names = F, col.names = T, sep = "\t"
)









# 紀錄 0.7 的 ----

 
cor_cutoff <- 0.7
exp_N <- fread("C:/Peter/rawData_eQTL/raw_Exp_mulInterval_N.txt", header = T)
exp_N <- t(exp_N) %>% as.data.table()
names(exp_N) <- as.character(exp_N[1, ])
exp_N <- exp_N[-c(1), ]
exp_N[, (names(exp_N)) := lapply(.SD, as.numeric)]

gene_pos <- fread("D:/oral_cancer/expression/expression_data/probe_pos_mulInterval.txt", header = T)
gene_pos <- gene_pos[Gene %in% names(exp_N)]
setDT(gene_pos)

gene_pos <- gene_pos[order(CHROMOSOME, start, end)]
gene_pos[, gap_prev := start - data.table::shift(end),
  by = CHROMOSOME
]

# 距離累計 > 1e6，更新 group_id
gene_pos[, group_id := make_group_id(gap_prev), by = CHROMOSOME]
gene_pos[, group_id := paste0("chr", CHROMOSOME, "_group", group_id)]
gene_pos[, gap_prev := NULL]


ld_type <- gene_pos$group_id %>% unique()
cor_group <- vector("list", length(ld_type))
for (i in seq_along(ld_type)) {
  df <- gene_pos[group_id %in% ld_type[i]]
  group_probe <- df$Gene
  cor_dt <- cor(exp_N[, ..group_probe]) %>% 
  as.table() %>% 
  as.data.table()

  setnames(cor_dt, c("probe_1", "probe_2", "cor"))

  # 留下對角線下的 cor 且 cor 夠大的
  cor_dt <- cor_dt[probe_1< probe_2 & abs(cor) > cor_cutoff]
  cor_group[[i]] <- cor_dt
}

cor_dt <- rbindlist(cor_group)
cor_dt <- gene_pos[, .(probe_1=Gene, chr=CHROMOSOME,group_id)][cor_dt, on=.(probe_1)] 

gt_N <- fread("C:/Peter/rawData_eQTL/r2_filter_0.8/outcome/raw_N_FDR_R2_0.8.txt")
eQTL_gene <- unique(gt_N$Gene)
eQTL_probe <- unique(gt_N$Probe)


# 計算每個 group，有顯著 eQTL 的數量 ----
gene_pos[,sig_FDR := 0]
gene_pos[Gene %in% eQTL_probe, sig_FDR := 1]
group_count <- gene_pos[, .(sig_FDR_sum = sum(sig_FDR, na.rm = TRUE)),
  by = group_id
]

# 挑出 cor>cor_cutoff 的 probe
cor_pass_probe <- gene_pos[Gene %in% union(cor_dt$probe_1,cor_dt$probe_2)]
record_sig <- cor_pass_probe[,.(High_Cor_probe_number=.N),by=group_id][group_count, on=.(group_id)]
record_sig[is.na(High_Cor_probe_number), High_Cor_probe_number := 0]



  fwrite(gene_pos,
    "C:/Peter/gene_enrichment/code_project/outcome/ld_bootstrap/probe_group.txt",
    row.names = F, col.names = T, sep = "\t"
  )

  fwrite(record_sig,
    "C:/Peter/gene_enrichment/code_project/outcome/ld_bootstrap/0.7number_record.txt",
    row.names = F, col.names = T, sep = "\t"
  )


