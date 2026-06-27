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

library(data.table)
library(openxlsx)
library(magrittr)



# 要填 excel，所需檔案路徑
# BHplot
# C:/Peter/repeatSNP_clumping_FIXpeople/r2_filter_0.8/outcome/BHplot_r2_0.8_original.png
# C:/Peter/repeatSNP_clumping_FIXpeople/r2_filter_0.8/outcome/BHplot_r2_0.8_scale.png

# FDR_qvalue
# C:/Peter/repeatSNP_clumping_FIXpeople/r2_filter_0.8/outcome/FDR_qvalue_r2_0.8_original.png
# C:/Peter/repeatSNP_clumping_FIXpeople/r2_filter_0.8/outcome/FDR_qvalue_r2_0.8_scale_0.05.png

# N,T Comparison
# C:/Peter/rawData_eQTL/outcome/paper_table_r2_0.8_FDR.xlsx
# C:/Peter/OQN_FIXpeople_before_eQTL/outcome/paper_table_r2_0.8_bon.xlsx

# probe 對到幾個 snp 圖
# C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_0.8/outcome/OQN_N_SNP_number_plot.png

# probe_info
# C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_0.8/outcome/OQN_exp_different_r2_0.8.txt

# maf_gt_N_FDR
# C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_0.8/outcome/OQN_N_FDR_R2_0.8.txt

# Finngen_1000gEUR
# C:/Peter/OQN_FIXpeople_before_eQTL/outcome/OQN_FDR_EUR_r2_filter_summary.xlsx
# C:/Peter/OQN_FIXpeople_before_eQTL/outcome/OQN_FDR_FIN_r2_filter_summary.xlsx
# C:/Peter/OQN_FIXpeople_before_eQTL/outcome/OQN_bon_EUR_r2_filter_summary.xlsx
# C:/Peter/OQN_FIXpeople_before_eQTL/outcome/OQN_bon_FIN_r2_filter_summary.xlsx

# clumping_pass_FDR
# C:/Peter/repeatSNP_clumping_FIXpeople/af_clump_info.xlsx
# C:/Peter/repeatSNP_clumping_FIXpeople/af_clumping_SNPID.txt
# C:/Peter/repeatSNP_clumping_FIXpeople/af_clumping_SNPID_1.txt

# 5783_snp_3_tech_FIN
# C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_0.8/outcome/OQN_MixFinngenPval_7_FIN.txt

# rs4975538_eQTL
# C:/Peter/rs4975538_permutation/outcome/rawData_eQTL/rs4975538_N_eQTL.txt
# C:/Peter/rs4975538_permutation/outcome/OQN_FIXpeople_before_eQTL/rs4975538_asso_probe_correlation.png


# 要上 paper table
# C:/Peter/rawData_eQTL/outcome/paper_table_r2_0.8_FDR.xlsx

1 + 1
names(a)

setdiff(names(a),c(
  "CHR", "Probe", "PROBE_COORDINATES", "Gene", "SNP_hg18", "rsID",
  "impute_type", "R2", "beta", "p-value", "FDR", "pval_raw", "FDR_raw (BH)",
  "rank_pval", "rank_FIXgene_pval", "rank_FIXpeople_pval"
))
