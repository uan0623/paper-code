## Setup ----
# 釋放記憶體
rm(list=ls())
gc()


library(data.table)
library(dplyr)
library(tidyverse)
library(stringr) 
library(ggplot2)
library(ggrepel)
library(matrixStats)  # rowMins
library(qvalue)


# 存成 .r 檔執行 (跑)


## Run Source ----
finngen_raw <- fread("D:/oral_cancer/snp_repeat_Finngen/outcome/C3_ORALCAVITY_EXALLC_2_hg19SNPunique.txt")[ ,.(chr = `#chrom`, hg19_snpID, P = pval, ref, alt, beta,mlogp,sebeta)]

# 跑完迴圈花費 40 mins
for (r2_threshold in c("0.6", "0.7", "0.9", "0.8", "no")) {
  for (eqtl_threshold in c("bon","FDR")) {
    cat("\n開始執行 r2_threshold =", r2_threshold, "eqtl_threshold =", eqtl_threshold, "\n")
    env <- new.env(parent = globalenv())
    env$r2_threshold <- r2_threshold
    env$eqtl_threshold <- eqtl_threshold
    env$finngen <- finngen_raw

    source("C:/Peter/vscode/.r/clumping.R", local = env)
    rm(env)
    gc(full = TRUE)
    cat("\n執行完:", eqtl_threshold, " ", r2_threshold, "\n")
  }
}

## Record ----
# C:\Peter\vscode\.r\3_tech_r2Filter_tableMaker.r 可生成 excel




