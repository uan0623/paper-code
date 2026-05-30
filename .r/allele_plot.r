library(data.table)
library(dplyr)
library(tidyverse)
library(stringr) 
library(ggplot2)
library(ggrepel)
library(matrixStats)  # rowMins
library(qvalue)
library(bestNormalize)
library(ggplot2)
library(openxlsx)
library(magrittr)



# 畫 d.e. in OQN raw data plot ----


d.e.raw <- fread("C:/Peter/rawData_eQTL/outcome/raw_exp_different.txt")
d.e.OQN <- fread("C:/Peter/OQN_FIXpeople_before_eQTL/outcome/OQN_exp_different.txt")
# setdiff(d.e.OQN[sig_OQN_Bonfi==1,PROBE_ID], d.e.raw[sig_raw_Bonfi==1,PROBE_ID])
# setdiff(d.e.raw[sig_raw_Bonfi==1,PROBE_ID], d.e.OQN[sig_OQN_Bonfi==1,PROBE_ID])


exp <- fread("D:/oral_cancer/expression/expression_data/ExpRes11sv20120601Gene.txt", header = T)

exp_N <- exp[, c(1,2, 43:82)]
names(exp_N) <-  c("Gene","PROBE_ID", sprintf("0%dB", 1:9), sprintf("%dB", 10:40))
exp_T <- exp[, c(1,2, 3:42)]
names(exp_N) <-  c("Gene","PROBE_ID", sprintf("0%dB", 1:9), sprintf("%dB", 10:40))

n_probe <- c("ILMN_2186216", "ILMN_2404850", "ILMN_2130441", "ILMN_1798177", "ILMN_1733103")
n_probe <- setdiff(d.e.raw[sig_raw_Bonfi==1,PROBE_ID], d.e.OQN[sig_OQN_Bonfi==1,PROBE_ID])[1:10]
diff_dt_raw <- exp_N[,-c("Gene","PROBE_ID")]-exp_T[,-c("Gene","PROBE_ID")]
diff_dt_raw <- cbind(exp_N[, c("Gene", "PROBE_ID")], diff_dt_raw) %>%
  filter(PROBE_ID %in% n_probe)

OQN_exp_N <- fread("C:/Peter/OQN_FIXpeople_before_eQTL/trash/OQN_ExpRes11sv20120601Gene_N.txt")
OQN_exp_T <- fread("C:/Peter/OQN_FIXpeople_before_eQTL/trash/OQN_ExpRes11sv20120601Gene_T.txt")

diff_dt <- OQN_exp_N[,-c("Gene","PROBE_ID")]-OQN_exp_T[,-c("Gene","PROBE_ID")]
diff_dt <- cbind(OQN_exp_N[, c("Gene", "PROBE_ID")], diff_dt) %>%
  filter(PROBE_ID %in% n_probe)


diff_dt_raw[, source:= "raw"]
diff_dt[, source:= "OQN_FIXpeople"]

combine_dt <- rbind(diff_dt,diff_dt_raw)

# id.vars = 解釋這個數值是誰的欄位; measure.vars = 真正要被攤開的數值欄位
plot_dt <- melt(
  combine_dt,
  id.vars = c("Gene", "PROBE_ID","source"),
  measure.vars = sample_cols,
  variable.name = "sample",
  value.name = "diff_N_minus_T"
)



plot_dt[, sample_index := match(sample, sample_cols)]

for (i in seq_along(n_probe)) {
   
  png(sprintf("C:/Peter/rawData_eQTL/outcome/probe_diff_%s.png", i),
      width = 1600, height = 1200, res = 200)
    
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






# 畫 genotype plot ----

# gt data ----
a <- fread("C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_0.8/outcome/OQN_N_FDR_R2_0.8.txt") %>%
  filter(sig_pval_Bonfi == 1) %>%
  setkey(`p-value`)

snp_bon <- a[, .SD[1], by = Gene] %>%
  select(CHR, SNP)

a <- fread("C:/Peter/rawData_eQTL/r2_filter_0.8/outcome/raw_N_FDR_R2_0.8.txt") %>%
  filter(sig_pval_Bonfi == 1) %>%
  setkey(`p-value`) %>%
  .[, .SD[1], by=Gene] %>% 
  select(CHR,SNP)

snp_bon <- rbind(snp_bon, a) %>% unique()
snp_bon[,CHR := sub(":.*", "", SNP) %>% as.numeric()]
setkey(snp_bon,CHR)

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

# alt, ref ----
gt_N <- fread("C:/Peter/rawData_eQTL/r2_filter_0.8/outcome/raw_N_FDR_R2_0.8.txt") 
test_2 <- gt_N[SNP %in% final$ID, c("SNP","ALT","REF")] %>% setkey(SNP)

# 轉 gt
geno_sub <- geno_dt[, .(
  sample,
  genotype = factor(
    get(target_snp),
    levels = c(0, 1, 2),
    labels = c("CC", "CT", "TT")
  )
)]












# plot_gt_exp_box ----


plot_gt_exp_box <- function(snp_id, probe_id, gt_dt, allele_dt, exp_dt) {
  
  sample_cols <- setdiff(names(exp_dt), c("Gene", "PROBE_ID"))
  
  # genotype: wide -> long
  gt_long <- gt_dt[ID == snp_id] %>%
    melt(
      id.vars = "ID",
      measure.vars = sample_cols,
      variable.name = "sample",
      value.name = "gt_code"
    ) %>%
    as.data.table()
  
  # expression: wide -> long
  exp_long <- exp_dt[PROBE_ID == probe_id] %>%
    melt(
      id.vars = c("Gene", "PROBE_ID"),
      measure.vars = sample_cols,
      variable.name = "sample",
      value.name = "expression"
    ) %>%
    as.data.table()
  
  # get REF / ALT
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
  
  ggplot(plot_dt, aes(x = genotype, y = expression, fill = genotype)) +
    geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.75) +
    geom_jitter(width = 0.12, size = 2, alpha = 0.8) +
    labs(
      title = paste0(snp_id, " / ", probe_id),
      x = "Genotype",
      y = "Expression"
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
    )
}



p <- plot_gt_exp_box(
  snp_id = "10:105086011",
  probe_id = "ILMN_1762337",
  gt_dt = final,
  allele_dt = test_2,
  exp_dt = exp_N
)