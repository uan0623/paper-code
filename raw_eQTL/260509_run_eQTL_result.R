
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


## Run Source ----
maf_raw <- fread("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/cis_hg18_19_snp_MAF.txt")
probe_pos_raw <- fread("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/probe_pos.txt")
probe_info_raw <- fread("D:/oral_cancer/expression/expression_data/probe24526infoB.txt")


for (r2_threshold in c(seq(0.6, 0.9, length.out = 4) %>% as.character(), "no")) {
  
  cat("\n開始執行 r2_threshold =", r2_threshold, "\n")
  env <- new.env(parent = globalenv())
  env$r2_threshold <- r2_threshold
  env$maf_raw <- maf_raw
  env$probe_pos_raw <- probe_pos_raw
  env$probe_info_raw <- probe_info_raw

  source("C:/Peter/vscode/raw_eQTL/250916_eQTL_result_2.R", local = env)
  rm(env)
  gc(full = TRUE)
  cat("\n執行完:", r2_threshold, "\n")
}


