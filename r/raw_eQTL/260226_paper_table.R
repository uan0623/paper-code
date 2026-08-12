# 目的: 生成放論文的 excel


# package ----

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
library(bestNormalize)



# ft ----
write_table_sheet <- function(dt, sheet_name, xlsx_path) {
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("Package 'openxlsx' is required. Please install it with install.packages('openxlsx').")
  }

  sheet_name <- substr(gsub("[][*?/\\\\:]", "_", sheet_name), 1, 31)
  
  wb <- if (file.exists(xlsx_path)) {
    openxlsx::loadWorkbook(xlsx_path)
  } else {
    openxlsx::createWorkbook()
  }
  
  if (sheet_name %in% openxlsx::sheets(wb)) {
    openxlsx::removeWorksheet(wb, sheet_name)
  }
  
  openxlsx::addWorksheet(wb, sheet_name)
  openxlsx::writeData(wb, sheet_name, as.data.frame(dt))
  if (ncol(dt) > 0) {
    openxlsx::setColWidths(wb, sheet_name, cols = 1:ncol(dt), widths = "auto")
  }
  openxlsx::freezePane(wb, sheet_name, firstRow = TRUE)
  openxlsx::saveWorkbook(wb, xlsx_path, overwrite = TRUE)
}

count_unique_col <- function(dt, cols) {
  col <- intersect(cols, names(dt))[1]
  if (is.na(col)) {
    return(NA_integer_)
  }
  uniqueN(dt[[col]])
}

build_paper_extra_records <- function(eqtl_threshold,
                                      r2_threshold = "0.8",
                                      xlsx_path) {
  eqtl_files <- sprintf(
    "C:/Peter/rawData_eQTL/r2_filter_%s/outcome/raw_%s_%s_R2_%s.txt",
    r2_threshold, c("N", "T"), eqtl_threshold, r2_threshold
  )
  
  exp_r2_path <- sprintf(
    "C:/Peter/rawData_eQTL/r2_filter_%s/outcome/raw_exp_different_r2_%s.txt",
    r2_threshold, r2_threshold
  )
  
  exp_dt <- fread(exp_r2_path, header = TRUE)
  
  nt_comparison <- rbindlist(lapply(seq_along(eqtl_files), function(idx) {
    tissue <- c("N", "T")[idx]
    dt <- fread(eqtl_files[idx], header = TRUE)
    
    rbind(
      data.table(
        tissue = tissue,
        r2_threshold = r2_threshold,
        condition = "sig_pval_Bonfi == 1 & sig_raw_Bonfi == 1",
        Total_Rows = dt[sig_pval_Bonfi == 1 & sig_raw_Bonfi == 1, .N],
        Probe = dt[sig_pval_Bonfi == 1 & sig_raw_Bonfi == 1, uniqueN(Probe)]
      ),
      data.table(
        tissue = tissue,
        r2_threshold = r2_threshold,
        condition = "sig_raw == 1 & FDR < 0.05",
        Total_Rows = dt[sig_raw == 1 & FDR < 0.05, .N],
        Probe = dt[sig_raw == 1 & FDR < 0.05, uniqueN(Probe)]
      )
    )
  }))
  
  write_table_sheet(nt_comparison, sheet_name = "probe_ratio", xlsx_path = xlsx_path)
  

  n_dt <- fread(eqtl_files[1], header = TRUE)
  t_dt <- fread(eqtl_files[2], header = TRUE)
  
  both_probe <- intersect(n_dt$Probe, t_dt$Probe)
  n_bon <- n_dt[sig_pval_Bonfi == 1]
  t_bon <- t_dt[sig_pval_Bonfi == 1]
  both_probe_bon <- intersect(n_bon$Probe, t_bon$Probe)
  
  sig_raw_probe <- exp_dt %>% filter(sig_raw == 1)
  sig_raw_probe <- sig_raw_probe[PROBE_ID %in% both_probe] %>% nrow()
  
  probe_bon <- exp_dt %>% filter(sig_raw_Bonfi == 1)
  probe_bon <- probe_bon[PROBE_ID %in% both_probe_bon] %>% nrow()
  
  common_probe <- data.table(
    eQTL_threshold = eqtl_threshold,
    r2_threshold = r2_threshold,
    record = c(
      "common probe in N,T",
      "common probe in N,T and N,T different FDR < 0.05 (raw)",
      "common Bonferroni probe in N,T",
      "common Bonferroni probe in N,T and N,T different Bonfi (raw)"
    ),
    unique_probe_number = c(
      length(both_probe),
      sig_raw_probe,
      length(both_probe_bon),
      probe_bon
    )
  )
  
  write_table_sheet(common_probe, sheet_name = "sig_ratio", xlsx_path = xlsx_path)


  
  sig_record <- rbindlist(lapply(eqtl_files, function(path) {
    dt <- fread(path, header = TRUE)
    tissue <- str_extract(basename(path), "(?<=raw_)\\w(?=_)")
    
    sig_summary <- data.table(
      tissue = tissue,
      eQTL_threshold = eqtl_threshold,
      r2_threshold = r2_threshold,
      record = "eQTL sig",
      unique_probe_number = count_unique_col(dt, c("ProbeID", "Probe")),
      unique_gene_number = count_unique_col(dt, "Gene"),
      unique_snp_number = count_unique_col(dt, "SNP"),
      association_number = nrow(dt)
    )
    
    oqn_fdr <- dt %>% filter(sig_raw == 1)
    oqn_fdr_summary <- data.table(
      tissue = tissue,
      eQTL_threshold = eqtl_threshold,
      r2_threshold = r2_threshold,
      record = "N,T different FDR < 0.05 (raw)",
      unique_probe_number = count_unique_col(oqn_fdr, c("ProbeID", "Probe")),
      unique_gene_number = count_unique_col(oqn_fdr, "Gene"),
      unique_snp_number = count_unique_col(oqn_fdr, "SNP"),
      association_number = nrow(oqn_fdr)
    )
    
    oqn_bon <- dt %>% filter(sig_raw_Bonfi == 1)
    oqn_bon_summary <- data.table(
      tissue = tissue,
      eQTL_threshold = eqtl_threshold,
      r2_threshold = r2_threshold,
      record = "N,T different Bonfi (raw)",
      unique_probe_number = count_unique_col(oqn_bon, c("ProbeID", "Probe")),
      unique_gene_number = count_unique_col(oqn_bon, "Gene"),
      unique_snp_number = count_unique_col(oqn_bon, "SNP"),
      association_number = nrow(oqn_bon)
    )
    
    rbind(sig_summary, oqn_fdr_summary, oqn_bon_summary)
  }))
  
  write_table_sheet(sig_record, sheet_name = "eQTL_DE_ratio", xlsx_path = xlsx_path)
  
  probe_source <- fread("C:/Peter/rawData_eQTL/trash/raw_maf_gt_N_MOSTpvalue.txt")
  probe_source[, gene2 := gsub("^([^_]+_[^_]+).*", "\\1", gene)]
  total_has_cis <- exp_dt[PROBE_ID %in% probe_source$gene2]
  
  not_sig_record <- rbindlist(lapply(eqtl_files, function(path) {
    dt <- fread(path, header = TRUE)
    tissue <- str_extract(basename(path), "(?<=raw_)\\w(?=_)")
    probe <- setdiff(total_has_cis$PROBE_ID, dt$Probe)
    
    sig_raw_probe <- total_has_cis %>% filter(sig_raw == 1)
    sig_raw_probe <- sig_raw_probe[PROBE_ID %in% probe] %>% nrow()
    
    Bonfi_raw_probe <- total_has_cis %>% filter(sig_raw_Bonfi == 1)
    Bonfi_raw_probe <- Bonfi_raw_probe[PROBE_ID %in% probe] %>% nrow()
    
    data.table(
      tissue = tissue,
      eQTL_threshold = eqtl_threshold,
      r2_threshold = r2_threshold,
      record = c(
        "eQTL not sig",
        "Among eQTL not sig, N,T different FDR < 0.05 (raw)",
        "Among eQTL not sig, N,T different Bonfi (raw)"
      ),
      unique_probe_number = c(length(probe), sig_raw_probe, Bonfi_raw_probe)
    )
  }))
  
  write_table_sheet(not_sig_record, sheet_name = "eQTL_DE_ratio_part2", xlsx_path = xlsx_path)
}







## get_info ft. ----

# 顯著 differentially express 的 probe (FDR<0.05), 有eQTL的 probe 比例

# 定義函數
get_info <- function(df, filter_exprs, target_cols) {
  # df: 輸入的 data.table
  # filter_exprs: 篩選條件的字串向量，例如 c("pval_raw < 0.05", "sig_raw_Bonfi == 1")
  # target_cols: 想要計算 uniqueN 的欄位名稱向量
  
# 1. 組合成篩選字串
  combined_condition <- paste(filter_exprs, collapse = " & ")
  
  # 2. 執行過濾與計算
  # 在 j 區塊使用 c(list(Total_Rows = .N), lapply(.SD, uniqueN))
  result <- df[eval(parse(text = combined_condition)), 
               c(list(Total_Rows = .N), lapply(.SD, uniqueN)), 
               .SDcols = target_cols]
  
  return(result)
}















# Build N/T comparison Excel sheet ----
# 顯著 differentially express 的 probe (FDR<0.05), 有eQTL的 probe 比例
df <- fread("C:/Peter/rawData_eQTL/r2_filter_0.8/outcome/raw_N_FDR_R2_0.8.txt",header=T)
# eQTL, differentially express 都用 Bonferroni 篩選顯著
get_info(df,
         c("sig_pval_Bonfi==1","sig_raw_Bonfi==1"), c("Probe"))



## 製作 excel 分頁 maf_gt_N, maf_gt_T ----
fdr_table_xlsx <- "C:/Peter/rawData_eQTL/outcome/paper_table_r2_0.8_FDR.xlsx"
a <- fread("C:/Peter/rawData_eQTL/r2_filter_0.8/outcome/raw_N_FDR_R2_0.8.txt") 

a <- a %>% 
  filter(sig_pval_Bonfi==1)
setkey(a, `p-value`)
uniqueN(a$Probe)
a <- a[, .SD[1], by=Gene] %>% 
  select(CHR,Probe,PROBE_COORDINATES,Gene,SNP,rsID,impute_type,R2,beta,`t-stat`,`p-value`,FDR,t_raw,
  pval_raw,pval_raw_BH,sig_raw_Bonfi)


# 補上不同 OQN 的rank_pval
exp <- fread("C:/Peter/rawData_eQTL/outcome/raw_exp_different.txt")[
  ,
  .(Probe = PROBE_ID, rank_pval)
]
exp_other1 <- fread("C:/Peter/OQN_FIXpeople_before_eQTL/outcome/OQN_exp_different.txt")[
  ,
  .(Probe = PROBE_ID, rank_FIXpeople_pval = rank_pval)
]
exp_other2 <- fread("C:/Peter/OQN_FIXgene_before_eQTL/outcome/OQN_exp_different.txt")[
  ,
  .(Probe = PROBE_ID, rank_FIXgene_pval = rank_pval)
]
a <- exp[a, on=.(Probe)]
a <- exp_other1[a, on=.(Probe)]
a <- exp_other2[a, on = .(Probe)]
setnames(a, old = c("SNP", "pval_raw_BH"), new = c("SNP_hg18", "FDR_raw (BH)"))

setcolorder(a,c(
  "CHR", "Probe", "PROBE_COORDINATES", "Gene", "SNP_hg18", "rsID",
  "impute_type", "R2", "beta","t-stat", "p-value", "FDR", "t_raw","pval_raw", "FDR_raw (BH)","sig_raw_Bonfi",
  "rank_pval", "rank_FIXgene_pval", "rank_FIXpeople_pval"
))
write_table_sheet(a, sheet_name = "maf_gt_N", xlsx_path = fdr_table_xlsx)




a <- fread("C:/Peter/rawData_eQTL/r2_filter_0.8/outcome/raw_T_FDR_R2_0.8.txt")

a <- a %>% 
  filter(sig_pval_Bonfi==1)
setkey(a, `p-value`)

a <- a[, .SD[1], by=Gene] %>% 
  select(CHR,Probe,PROBE_COORDINATES,Gene,SNP,rsID,impute_type,R2,beta,`t-stat`,`p-value`,FDR,t_raw,
  pval_raw,pval_raw_BH,sig_raw_Bonfi)


exp <- fread("C:/Peter/rawData_eQTL/outcome/raw_exp_different.txt")[
  ,
  .(Probe = PROBE_ID, rank_pval)
]
exp_other1 <- fread("C:/Peter/OQN_FIXpeople_before_eQTL/outcome/OQN_exp_different.txt")[
  ,
  .(Probe = PROBE_ID, rank_FIXpeople_pval = rank_pval)
]
exp_other2 <- fread("C:/Peter/OQN_FIXgene_before_eQTL/outcome/OQN_exp_different.txt")[
  ,
  .(Probe = PROBE_ID, rank_FIXgene_pval = rank_pval)
]
a <- exp[a, on=.(Probe)]
a <- exp_other1[a, on=.(Probe)]
a <- exp_other2[a, on = .(Probe)]
setnames(a, old = c("SNP", "pval_raw_BH"), new = c("SNP_hg18", "FDR_raw (BH)"))

setcolorder(a,c(
  "CHR", "Probe", "PROBE_COORDINATES", "Gene", "SNP_hg18", "rsID",
  "impute_type", "R2", "beta","t-stat", "p-value", "FDR", "t_raw","pval_raw", "FDR_raw (BH)","sig_raw_Bonfi",
  "rank_pval", "rank_FIXgene_pval", "rank_FIXpeople_pval"
))
write_table_sheet(a, sheet_name = "maf_gt_T", xlsx_path = fdr_table_xlsx)





# diferentially express 分頁 ----
a <- fread("C:/Peter/rawData_eQTL/r2_filter_0.8/outcome/raw_exp_different_r2_0.8.txt")

b <- a[
  sig_raw_Bonfi == 1,
  c("chr", "PROBE_ID", "PROBE_COORDINATES", "Gene", "pval_raw", "pval_raw_BH")
]
all(b$pval_raw_BH==1)

setkey(b, pval_raw)



a <- b[1:10, ]
setnames(a, old = "PROBE_ID", new = "Probe")
exp <- fread("C:/Peter/rawData_eQTL/outcome/raw_exp_different.txt")[
  ,
  .(Probe = PROBE_ID, rank_pval)
]
exp_other1 <- fread("C:/Peter/OQN_FIXpeople_before_eQTL/outcome/OQN_exp_different.txt")[
  ,
  .(Probe = PROBE_ID, rank_FIXpeople_pval = rank_pval)
]
exp_other2 <- fread("C:/Peter/OQN_FIXgene_before_eQTL/outcome/OQN_exp_different.txt")[
  ,
  .(Probe = PROBE_ID, rank_FIXgene_pval = rank_pval)
]
a <- exp[a, on=.(Probe)]
a <- exp_other1[a, on=.(Probe)]
a <- exp_other2[a, on = .(Probe)]

setcolorder(
  a,
  c(
    "chr", "Probe", "PROBE_COORDINATES", "Gene", "pval_raw", "pval_raw_BH",
    "rank_pval", "rank_FIXpeople_pval", "rank_FIXgene_pval"
  )
)


write_table_sheet(a, sheet_name = "diferentially_express", xlsx_path = fdr_table_xlsx)

build_paper_extra_records(
  eqtl_threshold = "FDR",
  r2_threshold = "0.8",
  xlsx_path = fdr_table_xlsx
)



## bon paper table sheets 
bon_table_xlsx <- "C:/Peter/rawData_eQTL/outcome/paper_table_r2_0.8_bon.xlsx"

a <- fread("C:/Peter/rawData_eQTL/r2_filter_0.8/outcome/raw_N_bon_R2_0.8.txt")

a <- a %>% 
  filter(sig_pval_Bonfi==1) %>% 
  select(CHR,Probe,PROBE_COORDINATES,Gene,SNP,rsID,impute_type,R2,beta,`p-value`,FDR,qvalue)

setorder(a, CHR,Gene,FDR,-R2)
write_table_sheet(a, sheet_name = "N_bon", xlsx_path = bon_table_xlsx)


a <- fread("C:/Peter/rawData_eQTL/r2_filter_0.8/outcome/raw_T_bon_R2_0.8.txt")

a <- a %>% 
  filter(sig_pval_Bonfi==1) %>% 
  select(CHR,Probe,PROBE_COORDINATES,Gene,SNP,rsID,impute_type,R2,beta,`p-value`,FDR,qvalue)

setorder(a, CHR,Gene,FDR,-R2)
write_table_sheet(a, sheet_name = "T_bon", xlsx_path = bon_table_xlsx)


a <- fread("C:/Peter/rawData_eQTL/r2_filter_0.8/outcome/raw_N_bon_R2_0.8.txt")

a <- a %>% 
  filter(sig_pval_Bonfi==1)
setkey(a, `p-value`)
uniqueN(a$Probe)
a <- a[, .SD[1], by=Gene] %>% 
  select(CHR,Probe,PROBE_COORDINATES,Gene,SNP,rsID,impute_type,R2,beta,`p-value`,FDR,qvalue,pval_raw,pval_raw_BH)

write_table_sheet(a, sheet_name = "maf_gt_N", xlsx_path = bon_table_xlsx)


a <- fread("C:/Peter/rawData_eQTL/r2_filter_0.8/outcome/raw_T_bon_R2_0.8.txt")

a <- a %>% 
  filter(sig_pval_Bonfi==1)
setkey(a, `p-value`)

a <- a[, .SD[1], by=Gene] %>% 
  select(CHR,Probe,PROBE_COORDINATES,Gene,SNP,rsID,impute_type,R2,beta,`p-value`,FDR,qvalue,pval_raw,pval_raw_BH)

write_table_sheet(a, sheet_name = "maf_gt_T", xlsx_path = bon_table_xlsx)



a <- fread("C:/Peter/rawData_eQTL/r2_filter_0.8/outcome/raw_exp_different_r2_0.8.txt")
b <- a[
  sig_raw_Bonfi == 1,
  c("chr", "PROBE_ID", "PROBE_COORDINATES", "Gene", "pval_raw", "pval_raw_BH")
]

setkey(b,pval_raw)

write_table_sheet(b[1:10,], sheet_name = "diferentially_express", xlsx_path = bon_table_xlsx)

build_paper_extra_records(
  eqtl_threshold = "bon",
  r2_threshold = "0.8",
  xlsx_path = bon_table_xlsx
)



































# 顯著的probe, 有eQTL的比例 ----
# FDR (BH)<0.05且有eQTL的gene，是否比較容易有N,T expression 不同？ All snp 的eQTL，反而是 FDR (BH)>0.05且有eQTL的gene，容易有exp 不同?


# N, T EXPRESSION 顯著不同，eQTL FDR<0.05 的probe 數目
for (eqtl_threshold in c("bon","FDR")) {

message("==========================================")
message(sprintf("eQTL threshold %s", eqtl_threshold))
message("==========================================")

a <- fread("C:/Peter/rawData_eQTL/r2_filter_0.8/outcome/raw_exp_different_r2_0.8.txt",header=T)

# expression FDR 顯著不同，且至少一個 associatioin eQTL FDR<0.05 的probe 數目
b <- a %>% filter(sig_raw == 1)
which(b$`0.8_N_sigCis-SNP_number` != 0) %>% length() %>% print()
which(b$`0.8_T_sigCis-SNP_number` != 0) %>% length() %>% print()
cat("\n")



# expression bonfi 顯著不同，且至少一個 associatioin eQTL bon 過門檻 的probe 數目
b <- a %>% filter(sig_raw_Bonfi==1)
which(b$`0.8_N_sigCis-SNP_number (bon)` != 0) %>%
  length() %>%
  print()
which(b$`0.8_T_sigCis-SNP_number (bon)` != 0) %>%
  length() %>%
  print()
cat("\n")




# FDR>0.05 result ----
file_name <- sprintf(
  "C:/Peter/rawData_eQTL/r2_filter_0.8/outcome/raw_%s_%s_R2_0.8.txt",
   c("N", "T"),
  eqtl_threshold
)

total <- fread("C:/Peter/rawData_eQTL/r2_filter_0.8/outcome/raw_exp_different_r2_0.8.txt",header=T)

file_name <- sprintf(
  "C:/Peter/rawData_eQTL/r2_filter_0.8/outcome/raw_%s_%s_R2_0.8.txt",
   c("N", "T"),
  eqtl_threshold
)

  
# 對這些檔案生成N,T 顯著snp, probe, association number
for (i in file_name) {
  a <- fread(i)
  
  # 挑出FDR>0.05 probe
  probe <- setdiff(total$PROBE_ID, a$Probe)
  
  # N,T different FDR sig (raw)
  sig_raw_probe <- total %>% filter(sig_raw ==1)
  sig_raw_probe <- sig_raw_probe[PROBE_ID %in% probe] %>% nrow()

  # N,T different Bonfi sig (raw)
  Bonfi_raw_probe <- total %>% filter(sig_raw_Bonfi ==1) 
  Bonfi_raw_probe <- Bonfi_raw_probe[PROBE_ID %in% probe] %>% nrow()
  
  
  
  cat(i,"\n")
  cat("eQTL not sig: ","\n")
  cat("unique Probe number: ",probe %>% length(),"\n")
  
  cat("Among eQTL not sig,  N,T different FDR<0.05 (raw): ","\n")
  cat("unique Probe number: ",sig_raw_probe,"\n")
  
  cat("Among eQTL not sig,  N,T different Bonfi (raw):","\n")
  cat("unique Probe number: ",Bonfi_raw_probe,"\n")
  
  cat(rep("\n", 3))
}



# FDR<0.05 result ----


file_name <- sprintf(
  "C:/Peter/rawData_eQTL/r2_filter_0.8/outcome/raw_%s_%s_R2_0.8.txt",
   c("N", "T"),
  eqtl_threshold
)

# 對這些檔案生成N,T 顯著snp, probe, association number
for (i in file_name) {
  a <- fread(i)
  
  # FDR sig
  sig_fdr_snp <- unique(a, by="SNP") %>% nrow() 
  sig_fdr_asso <- a %>% nrow()
  sig_fdr_probe <- unique(a, by="ProbeID") %>% nrow()
  sig_fdr_gene <- unique(a, by="Gene") %>% nrow()
  
  # N,T different FDR sig (raw)
  sig_raw_snp <- a %>% filter(sig_raw ==1) %>% unique( by="SNP") %>% nrow()
  sig_raw_asso <- a %>% filter(sig_raw ==1) %>% nrow()
  sig_raw_probe <- a %>% filter(sig_raw ==1) %>% unique( by="ProbeID") %>% nrow()
  sig_raw_gene<- a %>% filter(sig_raw ==1) %>% unique( by="Gene") %>% nrow()

  # N,T different Bonfi sig (raw)
  Bonfi_raw_snp <- a %>% filter(sig_raw_Bonfi ==1) %>% unique( by="SNP") %>% nrow()
  Bonfi_raw_asso <- a %>% filter(sig_raw_Bonfi ==1) %>% nrow()
  Bonfi_raw_probe <- a %>% filter(sig_raw_Bonfi ==1) %>% unique( by="ProbeID") %>% nrow()
  Bonfi_raw_gene <- a %>% filter(sig_raw_Bonfi ==1) %>% unique( by="Gene") %>% nrow()
  
  
  
  cat(i,"\n")
  
  cat("eQTL sig: ","\n")
  cat("unique Probe number: ")
  print(sig_fdr_probe)
  cat("unique gene number: ")
  print(sig_fdr_gene)
  cat("unique SNPs number: ")
  print(sig_fdr_snp)
  cat("association number: ")
  print(sig_fdr_asso)
  
  

  
  cat("N,T different FDR<0.05 (raw):","\n")
  cat("unique Probe number: ")
  print(sig_raw_probe)
  cat("unique gene number: ")
  print(sig_raw_gene)
  cat("unique SNPs number: ")
  print(sig_raw_snp)
  cat("association number: ")
  print(sig_raw_asso)
  
  
  cat("N,T different Bonfi (raw):","\n")
  cat("unique Probe number: ")
  print(Bonfi_raw_probe)
  cat("unique gene number: ")
  print(Bonfi_raw_gene)
  cat("unique SNPs number: ")
  print(Bonfi_raw_snp)
  cat("association number: ")
  print(Bonfi_raw_asso)
  
  cat(rep("\n", 3))
}



# common probe in N,T ----
# N,T 都有cis eQTL 的probe，在N,T expression 有顯著不同的比例很高?


total <- fread("C:/Peter/rawData_eQTL/r2_filter_0.8/outcome/raw_exp_different_r2_0.8.txt",header=T)


# gt pvalue

a <- sprintf(
  "C:/Peter/rawData_eQTL/r2_filter_0.8/outcome/raw_N_%s_R2_0.8.txt",
  eqtl_threshold
) %>% fread()
b <- sprintf(
  "C:/Peter/rawData_eQTL/r2_filter_0.8/outcome/raw_T_%s_R2_0.8.txt",
  eqtl_threshold
) %>% fread()

both_probe <- intersect(a$Probe, b$Probe)

a_bon <- a[sig_pval_Bonfi==1,]
b_bon <- b[sig_pval_Bonfi==1,]
both_probe_bon <- intersect(a_bon$Probe,b_bon$Probe)

names(a)
# 都用 FDR(BH) 
both_probe %>% length()
sig_raw_probe <- total %>% filter(sig_raw ==1)
sig_raw_probe <- sig_raw_probe[PROBE_ID %in% both_probe] %>% nrow()
sig_raw_probe


# 都用 Bonferroni
both_probe_bon %>% length()
probe_bon <- total %>% filter(sig_raw_Bonfi ==1)
probe_bon <- probe_bon[PROBE_ID %in% both_probe_bon] %>% nrow()
probe_bon

cat(rep("\n",2))



# probe 分別有幾個snp ----
# 每個probe 選一個snp，很明顯每個probe 分別有1個snp

file_name <- c("raw_N",
               "raw_T")

plot_output_dir <- "C:/Peter/rawData_eQTL/r2_filter_0.8/outcome"

for (i in file_name) {
  a <- sprintf(
    "C:/Peter/rawData_eQTL/r2_filter_0.8/outcome/%s_%s_R2_0.8.txt",
    i, eqtl_threshold
  ) %>% fread()
  
  k <- a[, .N, by = "Probe"]
  if (eqtl_threshold == "FDR"){

plot_title <- ifelse( i == "raw_N",
    sprintf("Number of significant cis-SNPs for %d probes in normal tissue (FDR < 0.05)", 
    nrow(k)),
    sprintf("Number of significant cis-SNPs for %d probes in tumor tissue (FDR < 0.05)",
     nrow(k))
  )
  } else if (eqtl_threshold == "bon"){
     plot_title <- ifelse( i == "raw_N",
    sprintf("Number of significant cis-SNPs for %d probes in normal tissue (bon pass)", nrow(k)),
    sprintf("Number of significant cis-SNPs for %d probes in tumor tissue (bon pass)", nrow(k))
  )
  }
  
  
  png(sprintf("%s/%s_SNP_number_plot.png", plot_output_dir, i),
      width = 1600, height = 1200, res = 200)
  plot(k$N,ylab = "SNP Number", main = plot_title)
  dev.off()
  
  png(sprintf("%s/%s_SNP_number_hist.png", plot_output_dir, i),
      width = 1600, height = 1200, res = 200)
  hist(k$N,ylab = "Frequency",xlab = "SNP Number", main = plot_title, breaks=50, col = "skyblue") 
  dev.off()
  
  
  cat(i,"\n")
  cat("mean of SNP Number in each probe: ",mean(k$N), "\n") 
  cat("median of SNP Number in each probe: ",median(k$N), "\n")  
  cat("\n")
}


}