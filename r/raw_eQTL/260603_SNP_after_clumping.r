library(data.table)
library(dplyr)
library(tidyverse)
library(stringr) 
library(ggplot2)
library(ggrepel)
library(matrixStats)  # rowMins
library(qvalue)
library(bestNormalize)

library(data.table)
library(openxlsx)
library(magrittr)






# 挑出不同 r2_threshold 做 LD clumping 後剩下的的 snp ID ----
maf_hg18_19 <- fread("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/cis_hg18_19_snp_MAF.txt",
  header = T
)
iteration <- 1
data_list <- vector("list")

for (pop in c("EUR", "FIN")) {
  for (eq in c("bon", "FDR")) {
    cat("\n開始整理 population =", pop, "\n")
    for (r2_threshold in c("0.6" ,"0.7","0.9", "0.8","no")) {
    base_dir <- sprintf("C:/Peter/repeatSNP_clumping_raw/r2_filter_%s/eQTL_%s/%s", r2_threshold, eq, pop)
    outcome_dir <- sprintf("%s/outcome", base_dir)
    trash_dir <- sprintf("%s/trash", base_dir)

    if (!dir.exists(outcome_dir) || !dir.exists(trash_dir)) {
      message(sprintf("Skip %s: missing outcome/trash directory under %s", pop, base_dir))
      next
    }

    outcome_prefix <- sprintf("%s/", outcome_dir)
    files <- c("2_fdr.txt", "2_qval.txt", "1_qval.txt", "1_fdr.txt")
    file_name <- paste0(outcome_prefix, files)
    
    for (i in seq_along(file_name)) {
      base <- basename(file_name[i])
      prefix <- sub("\\.txt$", "", base) # 取得 1_qval, 1_fdr 等前綴

      for (para in c(0.2)) {
        # 1000kb
        final_out <- sprintf("%s/1000kb_LD%s/%s.clumped", outcome_dir, para, prefix)
        if (!file.exists(final_out)) {
          next
        }
        a <- fread(final_out)
        k <- data.table(
          hg19_snpID = a$SNP,
          r2_threshold = r2_threshold,
          file = prefix,
          eqtl_threshold = eq,
          population = pop
        )
        data_list[[iteration]] <- k
        iteration <- iteration+1
      }
    }
    if (!any(lengths(data_list) > 0)) {
      message(sprintf("Skip %s %s: no clumped files found", pop, eq))
      next
    }
  }
  }
  }

result <- rbindlist(data_list, use.names = TRUE, fill = TRUE)
result_record <- maf_hg18_19[
  ,
  c("rsID", "hg18_snpID", "hg19_snpID")
][result, on = .(hg19_snpID)]
setkey(result_record, population, r2_threshold, eqtl_threshold, file)


# 加上是否 FDR<0.05，被紀錄再 eQTL 檔案
result_record_sub <- unique(result_record, by=c("hg18_snpID","r2_threshold"))
result_record[, FDR_sig := 0L]

for (r2 in c("0.6", "0.7", "0.9", "0.8", "no")) {

  eqtl_file <- sprintf(
    "C:/Peter/rawData_eQTL/r2_filter_%s/outcome/raw_N_FDR_R2_%s.txt",
    r2, r2
  ) %>% fread()
  
  af_clump_snphg18 <- result_record_sub[r2_threshold == r2, unique(na.omit(hg18_snpID))]
  matched_snp <- af_clump_snphg18[af_clump_snphg18 %in% eqtl_file$SNP]
  if (length(matched_snp) > 0) {
    result_record[r2_threshold == r2 & hg18_snpID %in% matched_snp, FDR_sig := 1L]
  }
}

fwrite(result_record,
  "C:/Peter/repeatSNP_clumping_raw/af_clumping_SNPID.txt",
  row.names = F, col.names = T, sep = "\t"
)

result_record[,c("hg18_snpID", "hg19_snpID"):=NULL]
result_record <- result_record[, .(rsID = paste(rsID,collapse = ", ")),
 by= .(r2_threshold	,population,eqtl_threshold, file)]

fwrite(result_record,
  "C:/Peter/repeatSNP_clumping_raw/af_clumping_SNPID_1.txt",
  row.names = F, col.names = T, sep = "\t"
)







# 挑出這些 clumping 後的 snp，結合 eQTL info ----

result_record <- fread("C:/Peter/repeatSNP_clumping_raw/af_clumping_SNPID.txt")
result_record_sub <- unique(result_record, by=c("hg18_snpID","r2_threshold"))

wb <- createWorkbook()
has_sheet <- FALSE

for (r2 in c("0.6", "0.7", "0.9", "0.8", "no")) {
  eqtl_file <- sprintf(
    "C:/Peter/rawData_eQTL/r2_filter_%s/outcome/raw_N_FDR_R2_%s.txt",
    r2, r2
  ) %>% fread()
  
  af_clump_snphg18 <- result_record_sub[r2_threshold == r2, unique(na.omit(hg18_snpID))]
  af_clump_info <- eqtl_file[SNP %in% af_clump_snphg18]

  sheet_name <- sprintf("eQTL_N_%s", r2)
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, af_clump_info)
  has_sheet <- TRUE
}

if (has_sheet) {
  saveWorkbook(
    wb,
    "C:/Peter/repeatSNP_clumping_raw/af_clump_info.xlsx",
    overwrite = TRUE
  )
}
