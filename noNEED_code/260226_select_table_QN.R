# ============================================================
# TITLE: Select Table QN
# SUBTITLE: Build selected quantile-normalized summary tables for downstream review.
# SEARCH TAGS: QN, select table, summary table, eQTL, SNP
# NOTE: Code below is unchanged; only this navigation header was added.
# ============================================================

# OUTLINE: Load required packages ----
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

# *Ordered_normal_quantile 函數*

# - Ordered Quantile normalization 公式來自 <https://petersonr.github.io/bestNormalize/reference/orderNorm.html#references-1>

# - Ordered_normal_quantile()跟內建的OQN函數 orderNorm() 落差絕對值最大 1.44329e-13，算出來是差不多的，所以可以信賴

Ordered_normal_quantile <- function(df, by.row = T){
  
  if (!is.matrix(df) ) {
    stop("df 必須是一個 metrix")
  }
  
  if(!by.row){
# 對每col 做排序####
    i=2
    n <- nrow(df)
    
    # 按照row 個數n，把(0,1) 切n+1 刀，再減掉0.5/n。這樣就能取n等分的值
    avg_sorted <- qnorm((1:n - 0.5)/n, mean = 0, sd = 1)
    
    # i=2，表示根據col運算，找出df 每個col不同row 的值大小順位，並賦予相應順位的值
    mat_qn <- apply(df, i, function(col) {
      
    # ties.method= "min" 代表 同值取最小的 rank，兩個expression 同值，normal_quantile 的結果會相同
    ranks <- rank(col, ties.method="min") 
    avg_sorted[ranks]})
    
    } else{
# 對每row 做排序####
    i=1
    n <- ncol(df)
    avg_sorted <- qnorm((1:n - 0.5)/n, mean = 0, sd = 1)
    
    # 找出df 每個row不同row 的值大小順位，並賦予相應順位的值
    mat_qn <- apply(df, i, function(row) {
      
    # i=2，表示根據row 運算，ties.method= "min" 代表 同值取最小的 rank，兩個expression 同值，normal_quantile 的結果會相同
    ranks <- rank(row, ties.method="min") 
    avg_sorted[ranks]
    
    }) %>% t()
    
    
    }
  
  
  return(mat_qn)
}






# *normal_quantile 函數*


normal_quantile <- function(df, by.row = T){
  
  if (!is.matrix(df) ) {stop("df 必須是一個 metrix")}
  
  if(!by.row){
# 對每col 做排序####
    i=2
    n <- nrow(df)
    avg_sorted <- apply(df,i,sort)
    avg_sorted <- rowMeans(avg_sorted)
    
    mat_qn <- apply(df, i, function(col) {
      ranks <- rank(col, ties.method="min") 
      avg_sorted[ranks]})
    } else{
# 對每row 做排序####
      i=1
      n <- ncol(df)
      avg_sorted <- apply(df,i,sort)
      avg_sorted <- colMeans(avg_sorted)
      
      mat_qn <- apply(df, i, function(row) {
        ranks <- rank(row, ties.method="min") 
        avg_sorted[ranks]
      }) %>% t()
    
    }
  
    return(mat_qn)
}



# OUTLINE: Define table information helper ----
## get_info ft. ----

# 顯著 differentially express 的 probe (FDR<0.05), 有eQTL的 probe 比例


# 定義函數
get_info <- function(df, filter_exprs, target_cols) {
  # df: 輸入的 data.table
  # filter_exprs: 篩選條件的字串向量，例如 c("pval_QN < 0.05", "sig_QN_Bonfi == 1")
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









# 挑出顯著的

# OUTLINE: Build N/T comparison Excel sheet ----
## 製作 excel 分頁 N,T Comparison ----
# 顯著 differentially express 的 probe (FDR<0.05), 有eQTL的 probe 比例

df <- fread("C:/Peter/QN_before_eQTL/trash/QN_maf_gt_N_pvalue_FDR_info.txt",header=T)

# eQTL, differentially express 都用 Bonferroni 篩選顯著
get_info(df,
         c("sig_pval_Bonfi==1","sig_QN_Bonfi==1"), c("Probe"))


# eQTL, differentially express 都用 FDR(BH) 篩選顯著
get_info(df,
         c("sig_QN==1","FDR<0.05"), c("Probe"))
get_info(df,
         c("sig_QN==1","FDR<0.05"), c("Probe"))



# OUTLINE: Build Bonferroni significant sheets ----
## 製作 excel 分頁 T_bon, N_bon ----
a <- fread("C:/Peter/QN_before_eQTL/trash/QN_maf_gt_N_pvalue_FDR_info.txt") 

a <- a %>% 
  filter(sig_pval_Bonfi==1) %>% 
  select(CHR,Probe,PROBE_COORDINATES,Gene,SNP,rsID,impute_type,R2,beta,`p-value`,FDR,qvalue)

setorder(a, CHR,Gene,FDR,-R2)
fwrite(a,"C:/Users/user/Desktop/N_bon.txt",
         row.names = F, col.names = T, sep = "\t")


a <- fread("C:/Peter/QN_before_eQTL/trash/QN_maf_gt_T_pvalue_FDR_info.txt") 

a <- a %>% 
  filter(sig_pval_Bonfi==1) %>% 
  select(CHR,Probe,PROBE_COORDINATES,Gene,SNP,rsID,impute_type,R2,beta,`p-value`,FDR,qvalue)

setorder(a, CHR,Gene,FDR,-R2)
fwrite(a,"C:/Users/user/Desktop/T_bon.txt",
         row.names = F, col.names = T, sep = "\t")


# OUTLINE: Build MAF genotype sheets ----
## 製作 excel 分頁 maf_gt_N, maf_gt_T ----
a <- fread("C:/Peter/QN_before_eQTL/trash/QN_maf_gt_N_pvalue_FDR_info.txt") 

a <- a %>% 
  filter(sig_pval_Bonfi==1)
setorder(a, FDR,-R2)

a <- a[, .SD[1], by=Gene] %>% 
  select(CHR,Probe,PROBE_COORDINATES,Gene,SNP,rsID,impute_type,R2,beta,`p-value`,FDR,qvalue,pval_QN,pval_QN_BH)

fwrite(a,"C:/Users/user/Desktop/maf_gt_N.txt",
         row.names = F, col.names = T, sep = "\t")


a <- fread("C:/Peter/QN_before_eQTL/trash/QN_maf_gt_T_pvalue_FDR_info.txt")
a <- a %>% 
  filter(sig_pval_Bonfi==1)
setorder(a, FDR,-R2)

a <- a[, .SD[1], by=Gene] %>% 
  select(CHR,Probe,PROBE_COORDINATES,Gene,SNP,rsID,impute_type,R2,beta,`p-value`,FDR,qvalue,pval_QN,pval_QN_BH)

fwrite(a,"C:/Users/user/Desktop/maf_gt_T.txt",
         row.names = F, col.names = T, sep = "\t")


# *顯著的probe, 有eQTL的比例*
# N, T EXPRESSION 顯著不同，eQTL FDR<0.05 的probe 數目

a <- fread("C:/Peter/QN_before_eQTL/trash/QN_exp_different.txt",header=T)
ncolumn <- ncol(a)
zero_index <- (rowSums(a[,(ncolumn-4):ncolumn])==0) %>% 
    which()


# 刪掉association 都是 MAF<0.05 snp 的 probe
a <- a[-zero_index,]
names(a)

cat("顯著的probe, 有eQTL的比例","\n")
# QN 顯著不同，eQTL FDR<0.05 的probe 數目
b <- a %>% filter(sig_QN==1)
for (i in (ncolumn-4):ncolumn) {
  which(b[[i]]!=0) %>% length() %>% print()
}
cat("\n")





# bonfi ----

# QN 顯著不同，eQTL FDR<0.05 的probe 數目
b <- a %>% filter(sig_QN_Bonfi==1)
for (i in (ncolumn-4):ncolumn) {
  which(b[[i]]!=0) %>% length() %>% print()
}
cat("\n")






# *FDR>0.05 result*
# 為了知道
# N,T 都有cis eQTL 的probe，在N,T expression 有顯著不同的比例很高?

file_name <- c("C:/Peter/QN_before_eQTL/trash/QN_maf_gt_N_pvalue_FDR_info.txt",
               "C:/Peter/QN_before_eQTL/trash/QN_maf_gt_T_pvalue_FDR_info.txt")

total <- fread("C:/Peter/QN_before_eQTL/trash/QN_exp_different.txt",header=T)

# 抓出22217 probe 中有cis-SNP 的20913個，藉此找出FDR>0.05 的 probe
a <- fread("C:/Peter/QN_before_eQTL/trash/QN_maf_gt_N_MOSTpvalue.txt")
a[, gene2 := gsub("^([^_]+_[^_]+).*", "\\1",gene)]
a[,gene := NULL]

total <- total[PROBE_ID %in% a$gene2] 
  
cat("FDR>0.05 result, N,T 都有cis eQTL 的probe，在N,T expression 有顯著不同的比例很高?","\n")
# 對這些檔案生成N,T 顯著snp, probe, association number
for (i in file_name) {
  a <- fread(i)
  
  # 挑出FDR>0.05 probe
  probe <- setdiff(total$PROBE_ID, a$Probe)
  
  # N,T different FDR sig (QN)
  sig_QN_probe <- total %>% filter(sig_QN ==1)
  sig_QN_probe <- sig_QN_probe[PROBE_ID %in% probe] %>% nrow()

  # N,T different Bonfi sig (QN)
  Bonfi_QN_probe <- total %>% filter(sig_QN_Bonfi ==1) 
  Bonfi_QN_probe <- Bonfi_QN_probe[PROBE_ID %in% probe] %>% nrow()
  
  
  
  cat(i,"\n")
  
  cat("eQTL q-value>0.05:","\n")
  cat("unique Probe number: ",probe %>% length(),"\n")
  
  cat("N,T different q-value<0.05 (QN):","\n")
  cat("unique Probe number: ",sig_QN_probe,"\n")
  
  cat("N,T different Bonfi (QN):","\n")
  cat("unique Probe number: ",Bonfi_QN_probe,"\n")
  
  cat(rep("\n", 3))
}



# OUTLINE: Build significant ratio sheet ----
## 製作 excel 分頁 sig_ratio ----
# *common probe in N,T*
# N,T 都有cis eQTL 的probe，在N,T expression 有顯著不同的比例很高?


total <- fread("C:/Peter/QN_before_eQTL/trash/QN_exp_different.txt",header=T)
ncolumn <- ncol(total)

zero_index <- (rowSums(total[,(ncolumn-4):ncolumn])==0) %>% 
    which()

# 刪掉association 都是 MAF<0.05 snp 的 probe
total <- total[-zero_index,]

# gt pvalue ----
a <- fread("C:/Peter/QN_before_eQTL/trash/QN_maf_gt_N_pvalue_FDR_info.txt")
b <- fread("C:/Peter/QN_before_eQTL/trash/QN_maf_gt_T_pvalue_FDR_info.txt")
both_probe <- intersect(a$Probe,b$Probe)
cat("N,T 都有cis eQTL 的probe，在N,T expression 有顯著不同的比例很高?","\n")
cat("gt pvalue","\n")

both_probe %>% length()

# QN
sig_QN_probe <- total %>% filter(sig_QN ==1)
sig_QN_probe <- sig_QN_probe[PROBE_ID %in% both_probe] %>% nrow()
sig_QN_probe

# QN_bon
sig_QN_probe <- total %>% filter(sig_QN_Bonfi ==1)
sig_QN_probe <- sig_QN_probe[PROBE_ID %in% both_probe] %>% nrow()
sig_QN_probe

cat(rep("\n",2))


# gt MOSTpvalue ----
a <- fread("C:/Peter/QN_before_eQTL/trash/QN_maf_gt_N_MOSTpvalue_FDR_info.txt")
b <- fread("C:/Peter/QN_before_eQTL/trash/QN_maf_gt_T_MOSTpvalue_FDR_info.txt")
both_probe <- intersect(a$Probe,b$Probe)
cat("gt MOSTpvalue","\n")

both_probe %>% length()

# QN
sig_QN_probe <- total %>% filter(sig_QN ==1)
sig_QN_probe <- sig_QN_probe[PROBE_ID %in% both_probe] %>% nrow()
sig_QN_probe

# QN_bon
sig_QN_probe <- total %>% filter(sig_QN_Bonfi ==1)
sig_QN_probe <- sig_QN_probe[PROBE_ID %in% both_probe] %>% nrow()
sig_QN_probe




##  probe 分別有幾個snp ----
# 每個probe 選一個snp，很明顯每個probe 分別有1個snp


# "C:/Peter/QN_before_eQTL/trash/QN_maf_gt_N_pvalue_FDR_info.txt"

file_name <- c("QN_maf_gt_N_pvalue",
               "QN_maf_gt_T_pvalue")

for (i in file_name) {
  a <- paste0("C:/Peter/QN_before_eQTL/trash/",i,"_FDR_info.txt") %>% 
    fread()
  k <- a[,.N ,by="Probe"]
  plot(k$N,ylab = "SNP Number", main = i)
  hist(k$N,ylab = "Frequency",xlab = "SNP Number", main = i, breaks=50, col = "skyblue") 
  boxplot(k$N, main = i)
  cat(i,"\n")
  cat("mean of SNP Number in each probe: ",mean(k$N), "\n") 
  cat("median of SNP Number in each probe: ",median(k$N), "\n")  
  cat("\n")
}


## 製作 excel 分頁 N,T Comparison 表格 ----

# FDR (BH)<0.05且有eQTL的gene，是否比較容易有N,T expression 不同？ All snp 的eQTL，反而是 FDR (BH)>0.05且有eQTL的gene，容易有exp 不同

rm(list=ls())
gc()

file_name <- c("C:/Peter/QN_before_eQTL/trash/QN_maf_gt_N_pvalue_FDR_info.txt",
               "C:/Peter/QN_before_eQTL/trash/QN_maf_gt_T_pvalue_FDR_info.txt")

# 對這些檔案生成N,T 顯著snp, probe, association number
for (i in file_name) {
  a <- fread(i)
  
  # FDR sig
  sig_fdr_snp <- unique(a, by="SNP") %>% nrow() 
  sig_fdr_asso <- a %>% nrow()
  sig_fdr_probe <- unique(a, by="ProbeID") %>% nrow()
  
  # N,T different FDR sig (QN)
  sig_QN_snp <- a %>% filter(sig_QN ==1) %>% unique( by="SNP") %>% nrow()
  sig_QN_asso <- a %>% filter(sig_QN ==1) %>% nrow()
  sig_QN_probe <- a %>% filter(sig_QN ==1) %>% unique( by="ProbeID") %>% nrow()

  # N,T different Bonfi sig (QN)
  Bonfi_QN_snp <- a %>% filter(sig_QN_Bonfi ==1) %>% unique( by="SNP") %>% nrow()
  Bonfi_QN_asso <- a %>% filter(sig_QN_Bonfi ==1) %>% nrow()
  Bonfi_QN_probe <- a %>% filter(sig_QN_Bonfi ==1) %>% unique( by="ProbeID") %>% nrow()
  
  
  
  cat(i,"\n")
  
  cat("eQTL q-value<0.05:","\n")
  cat("unique Probe number: ")
  print(sig_fdr_probe)
  cat("unique SNPs number: ")
  print(sig_fdr_snp)
  cat("association number: ")
  print(sig_fdr_asso)
  

  
  cat("N,T different q-value<0.05 (QN):","\n")
  cat("unique Probe number: ")
  print(sig_QN_probe)
  cat("unique SNPs number: ")
  print(sig_QN_snp)
  cat("association number: ")
  print(sig_QN_asso)
  
  
  cat("N,T different Bonfi (QN):","\n")
  cat("unique Probe number: ")
  print(Bonfi_QN_probe)
  cat("unique SNPs number: ")
  print(Bonfi_QN_snp)
  cat("association number: ")
  print(Bonfi_QN_asso)
  
  cat(rep("\n", 3))
}




##  FDR>0.05 result ----
# 為了知道
# N,T 都有cis eQTL 的probe，在N,T expression 有顯著不同的比例很高?

file_name <- c("C:/Peter/QN_before_eQTL/trash/QN_maf_gt_N_pvalue_FDR_info.txt",
               "C:/Peter/QN_before_eQTL/trash/QN_maf_gt_T_pvalue_FDR_info.txt")

total <- fread("C:/Peter/QN_before_eQTL/trash/QN_exp_different.txt",header=T)

# 抓出22217 probe 中有cis-SNP 的20913個，藉此找出FDR>0.05 的 probe
a <- fread("C:/Peter/QN_before_eQTL/trash/QN_maf_gt_N_MOSTpvalue.txt")
a[, gene2 := gsub("^([^_]+_[^_]+).*", "\\1",gene)]
a[,gene := NULL]

total <- total[PROBE_ID %in% a$gene2] 
  
  
# 對這些檔案生成N,T 顯著snp, probe, association number
for (i in file_name) {
  a <- fread(i)
  
  # 挑出FDR>0.05 probe
  probe <- setdiff(total$PROBE_ID, a$Probe)
  
  # N,T different FDR sig (QN)
  sig_QN_probe <- total %>% filter(sig_QN ==1)
  sig_QN_probe <- sig_QN_probe[PROBE_ID %in% probe] %>% nrow()

  # N,T different Bonfi sig (QN)
  Bonfi_QN_probe <- total %>% filter(sig_QN_Bonfi ==1) 
  Bonfi_QN_probe <- Bonfi_QN_probe[PROBE_ID %in% probe] %>% nrow()
  
  
  
  cat(i,"\n")
  
  cat("eQTL q-value>0.05:","\n")
  cat("unique Probe number: ",probe %>% length(),"\n")
  
  cat("N,T different q-value<0.05 (QN):","\n")
  cat("unique Probe number: ",sig_QN_probe,"\n")
  
  cat("N,T different Bonfi (QN):","\n")
  cat("unique Probe number: ",Bonfi_QN_probe,"\n")
  
  cat(rep("\n", 3))
}




# for diff R2 filter
# 不同的 r2 filter 做表格

## 製作 excel 分頁 N,T Comparison ----
# 顯著 differentially express 的 probe (FDR<0.05), 有eQTL的 probe 比例

for (i in seq(0.6,0.9, length.out=4) %>% as.character()) {

  df <- sprintf("C:/Peter/QN_before_eQTL/r2_filter_%s/outcome/QN_N_FDR_R2_%s.txt",i,i) %>% 
    fread(.,header=T)
    
get_info(df,
         c("sig_pval_Bonfi==1","sig_QN_Bonfi==1"), c("Probe")) %>% print()
get_info(df,
         c("sig_QN==1","FDR<0.05"), c("Probe")) %>% print()

df <- sprintf("C:/Peter/QN_before_eQTL/r2_filter_%s/outcome/QN_T_FDR_R2_%s.txt",i,i) %>% 
    fread(.,header=T)


get_info(df,
         c("sig_pval_Bonfi==1","sig_QN_Bonfi==1"), c("Probe")) %>% print()
get_info(df,
         c("sig_QN==1","FDR<0.05"), c("Probe")) %>% print()

  
}



# FDR (BH)<0.05且有eQTL的gene，是否比較容易有N,T expression 不同？ All snp 的eQTL，反而是 FDR (BH)>0.05且有eQTL的gene，容易有exp 不同


file_name <- c("QN_N_FDR_R2",
               "QN_T_FDR_R2")


cat("FDR (BH)<0.05且有eQTL的gene，是否比較容易有N,T expression 不同")
for (i in seq(0.6, 0.9, length.out=4) %>% as.character() ) {
  
  for (j in file_name) {
    
  a <- sprintf("C:/Peter/QN_before_eQTL/r2_filter_%s/outcome/%s_%s.txt",i,j,i) %>% 
    fread(., header = T)
    
  # FDR sig
  sig_fdr_snp <- unique(a, by="SNP") %>% nrow() 
  sig_fdr_asso <- a %>% nrow()
  sig_fdr_probe <- unique(a, by="ProbeID") %>% nrow()
  
  # N,T different FDR sig (QN)
  sig_QN_snp <- a %>% filter(sig_QN ==1) %>% unique( by="SNP") %>% nrow()
  sig_QN_asso <- a %>% filter(sig_QN ==1) %>% nrow()
  sig_QN_probe <- a %>% filter(sig_QN ==1) %>% unique( by="ProbeID") %>% nrow()

  # N,T different Bonfi sig (QN)
  Bonfi_QN_snp <- a %>% filter(sig_QN_Bonfi ==1) %>% unique( by="SNP") %>% nrow()
  Bonfi_QN_asso <- a %>% filter(sig_QN_Bonfi ==1) %>% nrow()
  Bonfi_QN_probe <- a %>% filter(sig_QN_Bonfi ==1) %>% unique( by="ProbeID") %>% nrow()
  
  
  
  cat(i,"\n")
  
  cat("eQTL q-value<0.05:","\n")
  cat("unique Probe number: ")
  print(sig_fdr_probe)
  cat("unique SNPs number: ")
  print(sig_fdr_snp)
  cat("association number: ")
  print(sig_fdr_asso)
  

  
  cat("N,T different q-value<0.05 (QN):","\n")
  cat("unique Probe number: ")
  print(sig_QN_probe)
  cat("unique SNPs number: ")
  print(sig_QN_snp)
  cat("association number: ")
  print(sig_QN_asso)
  
  
  cat("N,T different Bonfi (QN):","\n")
  cat("unique Probe number: ")
  print(Bonfi_QN_probe)
  cat("unique SNPs number: ")
  print(Bonfi_QN_snp)
  cat("association number: ")
  print(Bonfi_QN_asso)
  
  cat(rep("\n", 3))
    
  }
  
  
}





# OUTLINE: Find common probes in normal and tumor ----
## common probe in N,T ----
# N,T 都有cis eQTL 的probe，在N,T expression 有顯著不同的比例很高?


total <- fread("C:/Peter/QN_before_eQTL/trash/QN_exp_different_r2_filter.txt",header=T)
ncolumn <- ncol(total)

for (i in seq(0.6,0.9, length.out=4) %>% as.character()) {

  
a <- sprintf("C:/Peter/QN_before_eQTL/r2_filter_%s/outcome/QN_N_FDR_R2_%s.txt",i,i) %>% 
    fread(.,header=T)
b <- sprintf("C:/Peter/QN_before_eQTL/r2_filter_%s/outcome/QN_T_FDR_R2_%s.txt",i,i) %>% 
    fread(.,header=T)

both_probe <- intersect(a$Probe,b$Probe)

cat("N,T 都有cis eQTL 的probe，在N,T expression 有顯著不同的比例很高?","\n")
cat(i,"\n")

both_probe %>% length() %>% print()

# QN
sig_QN_probe <- total %>% filter(sig_QN ==1)
sig_QN_probe <- sig_QN_probe[PROBE_ID %in% both_probe] %>% nrow()
sig_QN_probe %>% print()

cat(rep("\n",2))
}




## FDR>0.05 result ----
# 為了知道
# N,T 都有cis eQTL 的probe，在N,T expression 有顯著不同的比例很高?



total_all <- fread("C:/Peter/QN_before_eQTL/trash/QN_exp_different.txt",header=T)

# 該 filter 下，有算出 associaiton 的 probe


for (j in seq(0.6,0.9, length.out=4) %>% as.character()) {

    
  file_name <-  sprintf("C:/Peter/QN_before_eQTL/r2_filter_%s/outcome/QN_%s_FDR_R2_%s.txt",j, c("N","T"),j)
  k <- sprintf("C:/Peter/QN_before_eQTL/r2_filter_%s/outcome/Has_association_probe_R2_%s.txt",
                        j,j) %>% fread()
  
  total <- total_all[PROBE_ID %in% k$Has_association_probe]
  cat("FDR>0.05 result")

  # 對這些檔案生成N,T 顯著snp, probe, association number
for (i in file_name) {
  a <- fread(i)
  
  # 挑出FDR>0.05 probe
  probe <- setdiff(total$PROBE_ID, a$Probe)
  
  # N,T different FDR sig (QN)
  sig_QN_probe <- total %>% filter(sig_QN ==1)
  sig_QN_probe <- sig_QN_probe[PROBE_ID %in% probe] %>% nrow()

  # N,T different Bonfi sig (QN)
  Bonfi_QN_probe <- total %>% filter(sig_QN_Bonfi ==1) 
  Bonfi_QN_probe <- Bonfi_QN_probe[PROBE_ID %in% probe] %>% nrow()
  
  
  
  cat(i,"\n")
  cat("r2 ", j, "has association probe number: ", nrow(k),"\n")
  cat("eQTL q-value>0.05:","\n")
  cat("unique Probe number: ",probe %>% length(),"\n")
  
  cat("N,T different q-value<0.05 (QN):","\n")
  cat("unique Probe number: ",sig_QN_probe,"\n")
  
  cat("N,T different Bonfi (QN):","\n")
  cat("unique Probe number: ",Bonfi_QN_probe,"\n")
  
  cat(rep("\n", 1))
}
  cat(rep("\n", 3))
}




# 不放表格

# OUTLINE: Compare expression distributions ----
## t test ----


### 比較有QN 的差別 ----


exp <- fread("D:/oral_cancer/expression/expression_data/ExpRes11sv20120601Gene.txt", header = T)


# raw data ----
# N-T
exp_diff <- as.matrix(exp[,43:82]) %>% as.data.table() - as.matrix(exp[,3:42]) %>% as.data.table()

n <- ncol(exp_diff)
dbar <- rowMeans(exp_diff)
sd_d <- apply(exp_diff, 1, sd)
tval <- dbar / (sd_d / sqrt(n))
pval <- 2 * pt(-abs(tval), df = n-1)

pvalue_raw <- cbind(exp[,1:2],pval,tval)
names(pvalue_raw) <- c("Gene", "PROBE_ID", "pval_raw","t_raw")


# 計算
pvalue_raw[,pval_raw_BH :=p.adjust(pval_raw, method = "BH") %>% 
             format(digits = 4,scientific = T) %>% 
             as.numeric()] 
  
# 判斷pvalue 是否小於0.05
pvalue_raw[,sig_raw :=ifelse(pval_raw_BH<0.05, 1,0)]
pvalue_raw[,sig_raw_Bonfi :=ifelse(pval_raw < 0.05/nrow(pvalue_raw), 1,0)]




# normal_quantile ----
exp_diff_afQN <- normal_quantile(as.matrix(exp[,43:82]), by.row = F) %>% as.data.table() - normal_quantile(as.matrix(exp[,3:42] ), by.row = F) %>% as.data.table()

n <- ncol(exp_diff_afQN)
dbar <- rowMeans(exp_diff_afQN)
sd_d <- apply(exp_diff_afQN, 1, sd)
tval <- dbar / (sd_d / sqrt(n))
pval <- 2 * pt(-abs(tval), df = n-1)

pvalue_qn <- cbind(exp[,1:2],pval,tval)
names(pvalue_qn) <- c("Gene", "PROBE_ID", "pval_QN","t_QN")


# 計算
pvalue_qn[,pval_QN_BH :=p.adjust(pval_QN, method = "BH") %>% 
             format(digits = 4,scientific = T) %>% 
             as.numeric()] 
  
# 判斷pvalue 是否小於0.05
pvalue_qn[,sig_QN :=ifelse(pval_QN_BH<0.05, 1,0)]
pvalue_qn[,sig_QN_Bonfi :=ifelse(pval_QN < 0.05/nrow(pvalue_qn), 1,0)]




# ordered_normal_quantile ----
exp_diff_afOQN <- Ordered_normal_quantile(as.matrix(exp[,43:82]), by.row = F) %>% as.data.table() - Ordered_normal_quantile(as.matrix(exp[,3:42] ), by.row = F) %>% as.data.table()

n <- ncol(exp_diff_afOQN)
dbar <- rowMeans(exp_diff_afOQN)
sd_d <- apply(exp_diff_afOQN, 1, sd)
tval <- dbar / (sd_d / sqrt(n))
pval <- 2 * pt(-abs(tval), df = n-1)

pvalue_oqn <- cbind(exp[,1:2],pval,tval)
names(pvalue_oqn) <- c("Gene", "PROBE_ID", "pval_OQN","t_OQN")


# 計算
pvalue_oqn[,pval_BH :=p.adjust(pval_OQN, method = "BH") %>% 
             format(digits = 4,scientific = T) %>% 
             as.numeric()] 
  
# 判斷pvalue 是否小於0.05
pvalue_oqn[,sig_OQN_BH :=ifelse(pval_BH<0.05, 1,0)]
pvalue_oqn[,sig_OQN_Bonfi :=ifelse(pval_OQN < 0.05/nrow(pvalue_oqn), 1,0)]


df <- merge(pvalue_raw,pvalue_qn,
            by=c("Gene","PROBE_ID"), all.x = T)
df <- merge(df,pvalue_oqn,
            by=c("Gene","PROBE_ID"), all.x = T)
df[,c("t_raw","t_QN","t_OQN"):=NULL]
fwrite(df,"C:/Users/user/Desktop/meeting/QN_OQN_compare.txt",
         row.names = F, col.names = T, sep = "\t")





# *OQN 比 QN 顯著的更多，但 raw data 顯著的較多出現在 QN 顯著 *

# QN 後
# 1. pval<0.05 probe number 從 15581 -> 15300，15192/15300 個 probe 在 15581 顯著probe 中
# 2. bon pval<0.05 probe number 從 6854 -> 6487， 6422/6487 個 probe 在 6854 顯著probe 中
# 3. FDR(BH) pval<0.05 probe number 從 15044 -> 14681， 14602/14681 個 probe 在 15044 顯著probe 中

# OQN 後
# 1. pval<0.05 probe number 從 15581 -> 15513，14094/15513 個 probe 在 15581 顯著probe 中
# 2. bon pval<0.05 probe number 從 6854 -> 7046， 5992/7046 個 probe 在 6854 顯著probe 中
# 3. FDR(BH) pval<0.05 probe number 從 15044 -> 14973， 13577/14973 個 probe 在 15044 顯著probe 中

sum(pvalue_raw$pval_raw<0.05); sum(pvalue_qn$pval_QN<0.05); sum(pvalue_oqn$pval_OQN<0.05)
sum(pvalue_raw$sig_raw_Bonfi); sum(pvalue_qn$sig_QN_Bonfi); sum(pvalue_oqn$sig_OQN_Bonfi)
sum(pvalue_raw$pval_raw_BH<0.05); sum(pvalue_qn$pval_QN_BH<0.05); sum(pvalue_oqn$pval_BH<0.05)

# how many probe QN in raw
pvalue_qn[pval_QN<0.05,Gene] %in% pvalue_raw[pval_raw<0.05,Gene] %>% 
  table()
pvalue_qn[sig_QN_Bonfi==1,Gene] %in% pvalue_raw[sig_raw_Bonfi==1,Gene] %>% 
  table()
pvalue_qn[pval_QN_BH<0.05,Gene] %in% pvalue_raw[pval_raw_BH<0.05,Gene] %>% 
  table()

# how many probe OQN in raw
pvalue_oqn[pval_OQN<0.05,Gene] %in% pvalue_raw[pval_raw<0.05,Gene] %>% 
  table()
pvalue_oqn[sig_OQN_Bonfi==1,Gene] %in% pvalue_raw[sig_raw_Bonfi==1,Gene] %>% 
  table()
pvalue_oqn[pval_BH<0.05,Gene] %in% pvalue_raw[pval_raw_BH<0.05,Gene] %>% 
  table()



# OUTLINE: Compare normal and tumor expression distributions ----
## N,T exp data 分佈不同? ----

# *長條圖*
# 不能用probe 重複出現的資料做 QN，probe值相同會變不同
exp <- fread("D:/oral_cancer/expression/expression_data/ExpRes11sv20120601Gene.txt", header = T)
hist(as.matrix(exp[,43]), breaks=40, main = "01B Normal part exp dist")
hist(as.matrix(exp[,3]), breaks=40, main = "01B Tumor part exp dist")

Ordered_normal_quantile(as.matrix(exp[,43]), by.row = F) %>% 
  hist( breaks=40, main = "After OQN, 01B Normal part exp dist")
Ordered_normal_quantile(as.matrix(exp[,3] ), by.row = F) %>% 
  hist( breaks=40, main = "After OQN, 01B Tumor part exp dist")

normal_quantile(as.matrix(exp[,43]), by.row = F) %>% 
  hist( breaks=40, main = "After QN, 01B Normal part exp dist")
normal_quantile(as.matrix(exp[,3] ), by.row = F) %>% 
  hist( breaks=40, main = "After QN, 01B Tumor part exp dist")

normal_quantile(as.matrix(exp[,43]), by.row = F) %>% 
  summary()

normal_quantile(as.matrix(exp[,3]), by.row = F) %>% 
  summary()


# *原始分佈*

df$y1_median <- ifelse(df$y1 > median(df$y1), 1,0)
# normal part NT=1
sample_1 <- data.table(exp = exp[,43] %>% 
                         as.matrix(),
                       NT = 1)
sample_2 <- data.table(exp = exp[,3] %>% 
                         as.matrix(),
                       NT = 0)
names(sample_1) <- c("exp","part")
names(sample_2) <- c("exp","part")
sample_1 <- rbind(sample_1,sample_2)


print(
  sample_1 %>% 
    ggplot(aes(x = sample_1[[1]], color = factor(part, labels = c("Normal", "Tumor")))) +
    geom_density(alpha = 0.5,  adjust = 0.3) +  # 使用透明度讓不同的密度曲線區分開
    labs(title = "Density Plot", 
         x = names(sample_1)[1], 
         y = "Density", 
         fill = "type") +  # 修改圖例標題
    theme_minimal()+
guides(color = guide_legend(title = NULL))  # 拿掉圖例標題
)


# normal_quantile ----
# normal part NT=1

sample_1 <- data.table(exp = normal_quantile(as.matrix(exp[,43:82] ), by.row = F)[,1] %>% 
                         as.matrix(),
                       NT = 1)
sample_2 <- data.table(exp = normal_quantile(as.matrix(exp[,3:42] ), by.row = F)[,1] %>% 
                         as.matrix(),
                       NT = 0)

names(sample_1) <- c("exp","part")
names(sample_2) <- c("exp","part")
sample_1 <- rbind(sample_1,sample_2)

print(
  sample_1 %>% 
    ggplot(aes(x = sample_1[[1]], color = factor(part, labels = c("Normal", "Tumor")))) +
    geom_density(alpha = 0.5,  adjust = 0.3) +  # 使用透明度讓不同的密度曲線區分開
    labs(title = "Density Plot After QN", 
         x = names(sample_1)[1], 
         y = "Density", 
         fill = "type") +  # 修改圖例標題
    theme_minimal()+
guides(color = guide_legend(title = NULL))  # 拿掉圖例標題
)



# Ordered_normal_quantile ----
# normal part NT=1
sample_1 <- data.table(exp = Ordered_normal_quantile(as.matrix(exp[,43] ), by.row = F) %>% 
                         as.matrix(),
                       NT = 1)
sample_2 <- data.table(exp = Ordered_normal_quantile(as.matrix(exp[,3] ), by.row = F) %>% 
                         as.matrix(),
                       NT = 0)
names(sample_1) <- c("exp","part")
names(sample_2) <- c("exp","part")
sample_1 <- rbind(sample_1,sample_2)

print(
  sample_1 %>% 
    ggplot(aes(x = sample_1[[1]], color = factor(part, labels = c("Normal", "Tumor")))) +
    geom_density(alpha = 0.5,  adjust = 0.3) +  # 使用透明度讓不同的密度曲線區分開
    labs(title = "Density Plot After OQN", 
         x = names(sample_1)[1], 
         y = "Density", 
         fill = "type") +  # 修改圖例標題
    theme_minimal()+
guides(color = guide_legend(title = NULL))  # 拿掉圖例標題
)




# *scale*
df$y1_median <- ifelse(df$y1 > median(df$y1), 1,0)
# normal part NT=1
sample_1 <- data.table(exp = exp[,43] %>% 
                           as.matrix(),
                       NT = 1)
sample_2 <- data.table(exp = exp[,3] %>% 
                           as.matrix(),
                       NT = 0)
names(sample_1) <- c("exp","part")
names(sample_2) <- c("exp","part")
sample_1 <- rbind(sample_1,sample_2)


print(
    sample_1 %>% 
        ggplot(aes(x = sample_1[[1]], color = factor(part, labels = c("Normal", "Tumor")))) +
        geom_density(alpha = 0.5,  adjust = 0.3) +  # 使用透明度讓不同的密度曲線區分開
        labs(title = "Density Plot", 
             x = names(sample_1)[1], 
             y = "Density", 
             fill = "type") +  # 修改圖例標題
        theme_minimal()+
        guides(color = guide_legend(title = NULL))+  # 拿掉圖例標題
        coord_cartesian(xlim = c(-1, 1))  # 限制 x 軸範圍
)


# normal_quantile ----
# normal part NT=1

sample_1 <- data.table(exp = normal_quantile(as.matrix(exp[,43:82] ), by.row = F)[,1] %>% 
                           as.matrix(),
                       NT = 1)
sample_2 <- data.table(exp = normal_quantile(as.matrix(exp[,3:42] ), by.row = F)[,1] %>% 
                           as.matrix(),
                       NT = 0)

names(sample_1) <- c("exp","part")
names(sample_2) <- c("exp","part")
sample_1 <- rbind(sample_1,sample_2)

print(
    sample_1 %>% 
        ggplot(aes(x = sample_1[[1]], color = factor(part, labels = c("Normal", "Tumor")))) +
        geom_density(alpha = 0.5,  adjust = 0.3) +  # 使用透明度讓不同的密度曲線區分開
        labs(title = "Density Plot After QN", 
             x = names(sample_1)[1], 
             y = "Density", 
             fill = "type") +  # 修改圖例標題
        theme_minimal()+
        guides(color = guide_legend(title = NULL))+  # 拿掉圖例標題
        coord_cartesian(xlim = c(-1, 1))  # 限制 x 軸範圍
)


## 40 N,T dist  ----
# 檢查 40樣本 N,T exp dist 是否顯著差異
a <- fread("D:/oral_cancer/expression/expression_data/ExpRes11sv20120601Gene.txt")
a <- a[,-c(1:2)]
new_dt <- melt(a, 
               measure.vars = colnames(a), 
               variable.name = "Sample_Name", 
               value.name = "Value")

setcolorder(new_dt, c("Value", "Sample_Name"))


for (i in 1:16) {

sub_dt <- new_dt[(1+22217*5*(i-1)):(22217*5*i), ]

print(
  sub_dt %>% 
    ggplot(aes(x = sub_dt[[1]],
               color = factor(Sample_Name,
                              labels = unique(sub_dt$Sample_Name)))
           ) +
    geom_density(alpha = 0.5,  adjust = 0.3) +  # 使用透明度讓不同的密度曲線區分開
    labs(title = "Density Plot", 
         x = names(sub_dt)[1], 
         y = "Density", 
         fill = "type") +  # 修改圖例標題
    theme_minimal()+
guides(color = guide_legend(title = NULL))  # 拿掉圖例標題
)
}









