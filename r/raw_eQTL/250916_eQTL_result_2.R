# 目的 ----
# 1. 對 r2 = c(seq(0.6, 0.9, length.out = 4) %>% as.character(), "no") 的 eQTL 算 FDR, qval
# 2. 每個 probe 加上 cisSNP number，放進 raw_%s_FDR_R2_%s.txt, raw_%s_bon_R2_%s.txt 等不同 imputed quality r2 紀錄 differentialy express的檔案 
# 用 .r/260509_run_eQTL_result.R 執行



## package ----
# 釋放記憶體
options(error = stop)


library(data.table)
library(dplyr)
library(tidyverse)
library(stringr) 
library(ggplot2)
library(ggrepel)
library(matrixStats)  # rowMins
library(qvalue)



# 計算所有FDR
compute_FDR <- function(inputname,R2filter ,outputname, qvalue_type){
  
  if (!is.character(inputname) || length(inputname) != 1 || !is.character(outputname) || length(outputname) != 1) {
    stop("inputname, outputname 必須是一個字串")
  }
  
  x <- fread(inputname,header = T)
  R2_filter <- fread(R2filter,header = T)
  
  # R2 filter
  x <- x[SNP %in% R2_filter$hg18_snpID,]
  
  x[,sig_pval_Bonfi := ifelse(`p-value`< 0.05/nrow(x), 1,0)]
  x[,FDR :=p.adjust(`p-value`, method = "BH") %>% 
     format(digits = 10,scientific = T) %>% 
     as.numeric()]
  
  if(qvalue_type==1){
    pi0_hat <- min(max(
    sum(x$`p-value` > 0.7) / ((1-0.7) * length(x$`p-value`)),
    0), 1)
    
    x[,qvalue := qvalue(`p-value`,pi0=pi0_hat)$qvalues]
  }
  
  if(qvalue_type==2){
    x[,qvalue := qvalue(`p-value`)$qvalues]
  }
  
  # FDR filter
  x <- x %>% filter(FDR< 0.05)
  x <- x[order(FDR),]
  x[, gene2 := gsub("^([^_]+_[^_]+).*", "\\1",gene)]
  
  
  fwrite(x, outputname,row.names = F, col.names = T, sep = "\t")
}





  if(r2_threshold =="no"){
    for (tissue in c("N","T")) {
          compute_FDR(sprintf("C:/Peter/rawData_eQTL/trash/raw_maf_gt_%s_pvalue.txt", tissue),
                    "D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/cis_hg18_19_snp_MAF.txt",
                    sprintf("C:/Peter/rawData_eQTL/r2_filter_%s/outcome/raw_maf_gt_%s_pvalue_FDR_R2_%s.txt",
                            r2_threshold, tissue, r2_threshold),
                    2)


      }
  } else{
    for (tissue in c("N","T")) {
        compute_FDR(sprintf("C:/Peter/rawData_eQTL/trash/raw_maf_gt_%s_pvalue.txt", tissue),
                  sprintf("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/cis_snp_R2_%s.txt",
                          r2_threshold),
                  sprintf("C:/Peter/rawData_eQTL/r2_filter_%s/outcome/raw_maf_gt_%s_pvalue_FDR_R2_%s.txt",
                          r2_threshold, tissue, r2_threshold),
                  2)


    }
  }
  

### choose most pval ----
# 對各種門檻的結果，找到有算出 associaiton 的 probe

gt <- fread("C:/Peter/rawData_eQTL/trash/raw_maf_gt_N_pvalue.txt")

  if(r2_threshold =="no"){
    a <- fread("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/cis_hg18_19_snp_MAF.txt")

  } else{
    a <- sprintf("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/cis_snp_R2_%s.txt",
    r2_threshold) %>% 
    fread()
  }
  
  gt_cisSNP <- gt[SNP %in% a$hg18_snpID,]
  probe_cisSNP <- unique(gt_cisSNP$gene)
  
  # 取兩個底線以前的字元，像是 "ILMN_1343291_1" -> "ILMN_1343291"
  probe_cisSNP <- gsub("^([^_]+_[^_]+).*", "\\1",probe_cisSNP)
  probe_cisSNP <- unique(probe_cisSNP)
  
  
  fwrite(data.table(Has_association_probe = probe_cisSNP),
         sprintf("C:/Peter/rawData_eQTL/r2_filter_%s/outcome/Has_association_probe_R2_%s.txt",
          r2_threshold,r2_threshold),
         row.names = F, col.names = T, sep = "\t")
  







### function  
rm(list = setdiff(ls(), c("probe_pos", "maf", "probe_info","r2_threshold")))
gc()



### Add sig Cis-SNP number ----
df <- fread("C:/Peter/rawData_eQTL/outcome/raw_exp_different.txt",header=T)
file_name <- c("N","T")

  for (file in file_name) {
  
  a <- sprintf("C:/Peter/rawData_eQTL/r2_filter_%s/outcome/raw_maf_gt_%s_pvalue_FDR_R2_%s.txt", 
  r2_threshold, file, r2_threshold ) %>% 
    fread(.,header=T)
  
  a[, gene2 := gsub("^([^_]+_[^_]+).*", "\\1",gene)]
  sig <- a[, .N, by = gene2]
  df <- merge(df,sig,by.x="PROBE_ID",by.y="gene2", all.x = T)
  df[is.na(N),N := 0]
  
  setnames(df,
    old = "N",
    new = paste0(r2_threshold, "_", file, "_sigCis-SNP_number")
  )

  # 加上 eQTL 過 bon 門檻的 cisSNP 數量 
  a <- a[sig_pval_Bonfi == 1, ]
  sig <- a[, .N, by = gene2]
  df <- merge(df,sig,by.x="PROBE_ID",by.y="gene2", all.x = T)
  df[is.na(N),N := 0]
  
  setnames(df,
    old = "N",
    new = paste0(r2_threshold, "_", file, "_sigCis-SNP_number (bon)")
  )
  }


  
    # 調整數值
  df[,pval_raw := format(pval_raw, digits = 4, scientific = T)]
  fwrite(df,
         sprintf("C:/Peter/rawData_eQTL/r2_filter_%s/outcome/raw_exp_different_r2_%s.txt",
         r2_threshold, r2_threshold),
         row.names = F, col.names = T, sep = "\t")




### use ----

final_cols <- c(
  "ProbeID", "Probe", "Gene", "CHR", "PROBE_COORDINATES",
  "sig_Cis-SNP_number", "sig_Cis-SNP_number_bon",
  "SNP", "rsID", "impute_type", "R2", "MAF", "REF", "ALT", "beta", "t-stat",
  "p-value", "FDR", "qvalue", "sig_pval_Bonfi", "t_raw", "pval_raw",
  "pval_raw_BH", "sig_raw", "sig_raw_Bonfi"
)

base_cols <- setdiff(final_cols, c("sig_Cis-SNP_number", "sig_Cis-SNP_number_bon"))

mix_info <- function(inputname) {
  fdr <- fread(inputname, header = TRUE)
  fdr[, gene2 := gsub("^([^_]+_[^_]+).*", "\\1", gene)]

  fdr <- pvalue_raw[fdr, on = .(gene2)]
  fdr <- probe_info[fdr, on = .(gene2)]
  fdr <- probe_pos[fdr, on = .(gene)]
  fdr <- maf[fdr, on = .(SNP)]

  fdr[, PROBE_COORDINATES := paste0(start, "-", end)]
  fdr[, c("start", "end", "gene", "INFO", "ER2", "chr") := NULL]
  setnames(fdr, old = c("gene2", "CHROMOSOME"), new = c("Probe", "CHR"))

  setcolorder(fdr, base_cols)
  fdr <- fdr[order(FDR)]

  fdr[, `:=`(
    MAF = round(MAF, 4),
    beta = round(beta, 4),
    `t-stat` = round(`t-stat`, 4),
    `p-value` = format(`p-value`, digits = 4, scientific = TRUE),
    FDR = format(FDR, digits = 4, scientific = TRUE),
    pval_raw = format(pval_raw, digits = 4, scientific = TRUE),
    t_raw = format(t_raw, digits = 4, scientific = TRUE),
    R2 = round(R2, 4),
    qvalue = format(qvalue, digits = 4, scientific = TRUE)
  )]

  fdr
}

pvalue_raw <- fread(sprintf(
  "C:/Peter/rawData_eQTL/r2_filter_%s/outcome/raw_exp_different_r2_%s.txt",
  r2_threshold, r2_threshold
), header = TRUE)

# 接外層讀的檔案
maf <- copy(maf_raw)
probe_pos <- copy(probe_pos_raw)
probe_info_raw <- copy(probe_info_raw)


maf[, hg19_snpID := NULL]



probe_info <- copy(probe_info_raw)
probe_info[, c("TargetID", "CHROMOSOME") := NULL]

probe_coord <- probe_info_raw[, .(
  Probe = PROBE_ID,
  PROBE_COORDINATES_new = PROBE_COORDINATES
)]

setnames(pvalue_raw, "PROBE_ID", "gene2", skip_absent = TRUE)
setnames(probe_info, "PROBE_ID", "gene2", skip_absent = TRUE)
setnames(probe_pos, "Gene", "gene", skip_absent = TRUE)
setnames(maf, "hg18_snpID", "SNP", skip_absent = TRUE)

cols_remove <- grep("_sigCis-SNP_number", names(pvalue_raw), value = TRUE)
pvalue_raw[, (cols_remove) := NULL]

setkey(pvalue_raw, gene2)
setkey(probe_info, gene2)
setkey(probe_pos, gene)
setkey(maf, SNP)
setkey(probe_coord, Probe)

suffix <- c("_sigCis-SNP_number", "_sigCis-SNP_number (bon)")

grid <- data.table(type = c("N", "T"), r2 = r2_threshold)
grid[, input := sprintf(
  "C:/Peter/rawData_eQTL/r2_filter_%s/outcome/raw_maf_gt_%s_pvalue_FDR_R2_%s.txt",
  r2, type, r2
)]
grid[, output := sprintf(
  "C:/Peter/rawData_eQTL/r2_filter_%s/outcome/raw_%s_FDR_R2_%s.txt",
  r2, type, r2
)]
grid[, bon_output := sprintf(
  "C:/Peter/rawData_eQTL/r2_filter_%s/outcome/raw_%s_bon_R2_%s.txt",
  r2, type, r2
)]

sig_number_path <- sprintf(
  "C:/Peter/rawData_eQTL/r2_filter_%s/outcome/raw_exp_different_r2_%s.txt",
  r2_threshold, r2_threshold
)

sig_number <- fread(sig_number_path, header = TRUE)

sig_number_cols <- as.vector(outer(
  paste0(grid$r2, "_", grid$type),
  suffix,
  paste0
))

missing_sig_number_cols <- setdiff(c("PROBE_ID", sig_number_cols), names(sig_number))
if (length(missing_sig_number_cols) > 0) {
  stop("Missing columns in ", sig_number_path, ": ",
       paste(missing_sig_number_cols, collapse = ", "))
}

sig_number <- sig_number[, c("PROBE_ID", sig_number_cols), with = FALSE]

for (idx in seq_len(nrow(grid))) {
  cols_now <- paste0(grid$r2[idx], "_", grid$type[idx], suffix)

  sig_now <- sig_number[, .SD, .SDcols = c("PROBE_ID", cols_now)]
  setnames(sig_now, cols_now, c("sig_Cis-SNP_number", "sig_Cis-SNP_number_bon"))

  gt <- mix_info(grid$input[idx])
  gt <- sig_now[gt, on = .(PROBE_ID = Probe)]
  setnames(gt, "PROBE_ID", "Probe")

  gt[
    !REF %in% c("A", "T", "C", "G") | !ALT %in% c("A", "T", "C", "G"),
    rsID := NA_character_
  ]

  gt <- probe_coord[gt, on = .(Probe)]
  gt[, PROBE_COORDINATES := NULL]
  setnames(gt, "PROBE_COORDINATES_new", "PROBE_COORDINATES")

  setcolorder(gt, final_cols)

  fwrite(gt, grid$output[idx], row.names = FALSE, col.names = TRUE, sep = "\t")

  gt_bon <- gt[sig_pval_Bonfi==1,]
  fwrite(gt_bon, grid$bon_output[idx], row.names = FALSE, col.names = TRUE, sep = "\t")
}


# 多出i.PROBE_COORDINATES, i.chr 欄位
for (idx in seq_len(nrow(grid))) {
  a <- fread(grid$output[idx])
  b <- fread(grid$bon_output[idx])
  if ("i.PROBE_COORDINATES" %in% names(a) &
    "i.PROBE_COORDINATES" %in% names(b)) {
    a[, "i.PROBE_COORDINATES" := NULL]
    b[, "i.PROBE_COORDINATES" := NULL]
  }
  
  if ("i.chr" %in% names(a) &
    "i.chr" %in% names(b)) {
    a[, "i.chr" := NULL]
    b[, "i.chr" := NULL]
  }
  
  fwrite(a, grid$output[idx], row.names = FALSE, col.names = TRUE, sep = "\t")
  fwrite(b, grid$bon_output[idx], row.names = FALSE, col.names = TRUE, sep = "\t")
}






