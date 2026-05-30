
## package ----
# 釋放記憶體
rm(list=ls())
gc()


library(data.table)
library(dplyr)
library(tidyverse)
library(stringr)
library(ggplot2)
library(ggrepel)
library(matrixStats) # rowMins
library(qvalue)



## 存成 Rdata，輕量化

# first finngen data
# finngen <- fread(
#   "D:/oral_cancer/snp_repeat_Finngen/outcome/C3_ORALCAVITY_EXALLC_2_hg19SNPunique.txt",
#   select = c("hg38_snpID", "hg19_snpID", "pval", "alt", "ref"),
#   sep = "\t",
#   header = TRUE,
#   nThread = 1
# )
# saveRDS(
#   finngen,
#   "D:/oral_cancer/snp_repeat_Finngen/outcome/C3_ORALCAVITY_EXALLC_2_hg19SNPunique.rds"
# )



# # second finngen data
# finngen_previous <- fread(
#   "D:/oral_cancer/TWAS/summary_stat/FinnGen/C3_ORALCAVITY_EXALLC",
#   sep = "\t",
#   header = TRUE,
#   nThread = 1
# )
# saveRDS(
#   finngen_previous,
#   "D:/oral_cancer/TWAS/summary_stat/FinnGen/C3_ORALCAVITY_EXALLC.rds"
# )






# 存成 .r 檔執行 (跑)
# - 每次 source 用獨立環境

## Run Source ----
# 執行完 c( "0.6","0.7", "0.8","0.9","no")  的 repeatSNP_FinnGen.R 花費 3小時

finngen_previous <- readRDS("D:/oral_cancer/TWAS/summary_stat/FinnGen/C3_ORALCAVITY_EXALLC.rds")
finngen <- readRDS("D:/oral_cancer/snp_repeat_Finngen/outcome/C3_ORALCAVITY_EXALLC_2_hg19SNPunique.rds")
setnames(finngen,old = "hg38_snpID",new = "snp_hg38")


for (r2_threshold in c( "0.6","0.7", "0.8","0.9","no")) {
  
  cat("\n開始執行 r2_threshold =", r2_threshold, "\n")
  
  env <- new.env(parent = globalenv())
  env$r2_threshold <- r2_threshold
  env$finngen_previous <- finngen_previous
  env$finngen <- finngen
  
  source("C:/Peter/vscode/raw_eQTL/repeatSNP_FinnGen.R", local = env)
  
  rm(env)
  gc(full = TRUE)
  
  cat("\n執行完:", r2_threshold, "\n")
}



## Record ----
# 用 C:\Peter\vscode\.r\3_tech_r2Filter_tableMaker.r 生成 excel

r2_threshold <- "0.8"
# 畫 BH plot
## BH Plot ----
for (r2_threshold in c("no","0.6" ,"0.7", "0.8", "0.9")) {
  

data_1 <- fread(
  sprintf("C:/Peter/rawData_eQTL/r2_filter_%s/outcome/raw_MixFinngenPval_7_EUR.txt",
          r2_threshold) 
)

data_1 <- data_1[!is.na(gt_N_hg38_finngen_pval)]
m_total_1 <- nrow(data_1)

plot_dt_1 <- data_1[
  order(gt_N_hg38_finngen_pval),
  .(
    pval = gt_N_hg38_finngen_pval,
    FDR  = gt_N_hg38_finngen_FDR
  )
]

plot_dt_1[, rank := .I]
plot_dt_1[, bh_threshold := rank / m_total_1 * 0.05]
bonfi_threshold_1 <- 0.05 / m_total_1




# 圖降採樣，避免點太多
# idx_full <- unique(c(
#   seq(1, nrow(plot_dt_1), by = 4),
#   nrow(plot_dt_1)
# ))

idx_full <- 1:nrow(plot_dt_1)
plot_dt_1_full_sample <- plot_dt_1[idx_full]


tmp_dir <- sprintf("C:/Peter/repeatSNP_clumping_raw/r2_filter_%s/outcome", r2_threshold)
if(!dir.exists(tmp_dir)) dir.create(tmp_dir, recursive = TRUE)

png(
  sprintf("C:/Peter/repeatSNP_clumping_raw/r2_filter_%s/outcome/BHplot_EUR.png", r2_threshold),
  width = 2000,
  height = 1000,
  res = 200
)

par(mfrow = c(1, 2))

main_title <- if (r2_threshold == "no") {
  "FIN Association (No R2 Filter)"
} else {
  sprintf("FIN R2 > %s Association", r2_threshold)
}

# 完整圖 ####
plot(
  plot_dt_1_full_sample$rank,
  plot_dt_1_full_sample$pval,
  pch = ".",
  col = "black",
  xlab = sprintf("Rank (m=%d)", m_total_1),
  ylab = "P-value",
  main = main_title
)

lines(
  plot_dt_1_full_sample$rank,
  plot_dt_1_full_sample$bh_threshold,
  col = "red",
  lwd = 2
)

points(
  plot_dt_1[FDR <= 0.05, rank],
  plot_dt_1[FDR <= 0.05, pval],
  col = "orange",
  pch = "."
)

abline(h = bonfi_threshold_1, col = "blue", lty = 2)

legend(
  "topleft",
  legend = c(
    "BH cutoff",
    "FDR <= 0.05",
    paste("Bonferroni p <=", format(bonfi_threshold_1, digits = 2))
  ),
  col = c("red", "orange", "blue"),
  lty = c(1, NA, 2),
  pch = c(NA, 20, NA),
  lwd = c(2, NA, 2),
  bty = "n"
)





# 放大圖 ####
if (!any(plot_dt_1$FDR <= 0.05, na.rm = TRUE)) {
  next("沒有任何 gt_N_hg38_finngen_FDR <= 0.05 的 association，不能用顯著 p-value 範圍縮放")
}
sig_cut_1 <- max(plot_dt_1[FDR <= 0.05, pval], na.rm = TRUE)
ylim_max_1 <- sig_cut_1 * 1.5
plot_dt_1_zoomed <- plot_dt_1[pval <= ylim_max_1]


plot(
  plot_dt_1_zoomed$rank,
  plot_dt_1_zoomed$pval,
  pch = 20,
  col = "black",
  xlab = sprintf("Rank (m=%d)", m_total_1),
  ylab = "P-value",
  main = "FDR Sig",
  ylim = c(0, ylim_max_1)
)

lines(
  plot_dt_1_zoomed$rank,
  plot_dt_1_zoomed$bh_threshold,
  col = "red",
  lwd = 2
)

points(
  plot_dt_1_zoomed[FDR <= 0.05, rank],
  plot_dt_1_zoomed[FDR <= 0.05, pval],
  col = "orange",
  pch = 20
)

abline(h = bonfi_threshold_1, col = "blue", lty = 2)



dev.off()

}


## eQTL ----
# 對算出eQTL 的畫圖
R2_filter <- fread("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/cis_snp_R2_0.8.txt",header = T)

data_N_r2  <- fread("C:/Peter/rawData_eQTL/trash/raw_maf_gt_N_pvalue.txt",header = T)
# R2 filter
data_N_r2  <- data_N_r2[SNP %in% R2_filter$hg18_snpID,]


data_T_r2  <- fread("C:/Peter/rawData_eQTL/trash/raw_maf_gt_T_pvalue.txt",header = T)
data_T_r2  <- data_T_r2[SNP %in% R2_filter$hg18_snpID,]




  
m_total_T <- nrow(data_T_r2)
m_total_N <- nrow(data_N_r2)


# T 
# 篩選 p < 0.05
p_T <- sort(data_T_r2[, `p-value`])
rank_T <- 1:length(p_T)

# 找Bonferroni cutoff
bonfi_threshold_T <- 0.05/m_total_T

# 找每個pvalue 的BH cutoff
bh_threshold_T <- rank_T / m_total_T * 0.05
# 找小於BH cutoff 的pvalue 當中最大的(顯著pvalue 最大值)
sig_cut_T <- max(p_T[p_T <= bh_threshold_T])

# 顯著pvalue 最大值的 1.5 倍
ylim_max_T <- sig_cut_T * 1.5

# 篩選出在「縮放範圍內」的pvalue
p_T_zoomed <- p_T[p_T <= ylim_max_T]
# 找縮放範圍內的 BH cutoff
bh_threshold_T_zoomed <- bh_threshold_T[which(p_T <= ylim_max_T)]

# X 軸是 rank/m，找縮放範圍內的
plot_x_T_zoomed <- which(p_T <= ylim_max_T) 



# N
p_N <- sort(data_N_r2[, `p-value`])
m_subset_N <- length(p_N)
rank_N <- 1:m_subset_N

bonfi_threshold_N <- 0.05/m_total_N
bh_threshold_N <- rank_N / m_total_N * 0.05
sig_cut_N <- max(p_N[p_N <= bh_threshold_N])

ylim_max_N <- sig_cut_N * 1.5

p_N_zoomed <- p_N[which(p_N <= ylim_max_N)]
bh_threshold_N_zoomed <- bh_threshold_N[which(p_N <= ylim_max_N)]
plot_x_N_zoomed <- which(p_N <= ylim_max_N) 




# 0.8 BH plot ----
png("C:/Peter/repeatSNP_clumping_raw/r2_filter_0.8/outcome/BHplot_r2_0.8_original.png", width = 2000, height = 1000, res = 200)
par(mfrow = c(1, 2))

# --- 開始畫 T 組的圖 ---
plot(rank_T[c(seq(1,length(rank_T), by = 4), length(rank_T))],
     p_T[c(seq(1,length(p_T), by = 4), length(p_T))],
      pch = ".",
      col = "black",
      xlab = sprintf("Rank (m=%d)",m_total_T),
      ylab = "P-value",
      main = "Tumor Part: All Associaiton") 

lines(rank_T[c(seq(1,length(rank_T), by = 4), length(rank_T))],
      bh_threshold_T[c(seq(1,length(bh_threshold_T), by = 4), length(bh_threshold_T))],
      col = "red",
      lwd = 2)

# 標示橘色點 (在縮放範圍內，且 p <= sig_cut 的點)
sig_indices_T_zoomed <- (p_T_zoomed <= sig_cut_T)
points(plot_x_T_zoomed[sig_indices_T_zoomed], p_T_zoomed[sig_indices_T_zoomed],
       col = "orange",
       pch = ".")
abline(h = bonfi_threshold_T, col = "blue", lty = 2)


legend("topleft",
       legend = c("FDR < 0.05(BH)",
                  paste("Bonferroni (p-value <=", format(bonfi_threshold_T, digits=2), ")")
                  ),
       col = c( "red", "blue"),
       lty = c(1,2),
       lwd = c(2,2))

# --- 開始畫 N 組的圖 (只畫縮放後的資料) ---
plot(rank_N[c(seq(1,length(rank_N), by = 4), length(rank_N))],
     p_N[c(seq(1,length(p_N), by = 4), length(p_N))],
      pch = ".",
      col = "black",
      xlab = sprintf("Rank (m=%d)",m_total_N),
      ylab = "P-value",
      main = "Normal Part: All Associaiton") 

lines(rank_N[c(seq(1,length(rank_N), by = 4), length(rank_N))],
      bh_threshold_N[c(seq(1,length(bh_threshold_N), by = 4), length(bh_threshold_N))],
      col = "red",
      lwd = 2)

sig_indices_N_zoomed <- (p_N_zoomed <= sig_cut_N)
points(plot_x_N_zoomed[sig_indices_N_zoomed], p_N_zoomed[sig_indices_N_zoomed],
       col = "orange",
       pch = ".")
abline(h = bonfi_threshold_N, col = "blue", lty = 2)

legend("topleft",
       legend = c("FDR < 0.05(BH)",
                  paste("Bonferroni (p-value <=", format(bonfi_threshold_N, digits=2), ")")
                  ),
       col = c( "red", "blue"),
       lty = c(1,2),
       lwd = c(2,2))

# 關閉檔案
dev.off()


# *FDR vs Qvalue for all pvalue*

data_N_r2[,FDR :=p.adjust(`p-value`, method = "BH") %>% 
     format(digits = 10,scientific = T) %>% 
     as.numeric()]

data_T_r2[,FDR :=p.adjust(`p-value`, method = "BH") %>% 
     format(digits = 10,scientific = T) %>% 
     as.numeric()]

data_N_r2[,qvalue := qvalue(`p-value`)$qvalues]
data_scale_N <- data_N_r2 %>% 
  filter(FDR<0.05)

data_T_r2[,qvalue := qvalue(`p-value`)$qvalues]
data_scale_T <- data_T_r2 %>% 
  filter(FDR<0.05)

setkey(data_T_r2, FDR, qvalue)
setkey(data_N_r2, FDR, qvalue)

# 4個取一個
data_sub_T <- data_T_r2[c(seq(1,nrow(data_T_r2), by=4), nrow(data_T_r2)),]
data_sub_N <- data_N_r2[c(seq(1,nrow(data_N_r2), by=4), nrow(data_N_r2)),]



# 0.8 Qvalue vs FDR圖 ####
png("C:/Peter/repeatSNP_clumping_raw/r2_filter_0.8/outcome/FDR_qvalue_r2_0.8_original.png", width = 2000, height = 1000, res = 200)
par(mfrow = c(1, 2))


# --- 開始畫 T 組的圖---
plot(data_sub_T$FDR,
     data_sub_T$qvalue,
     pch = ".",
     col = "black",
     main = "Tumor Part FDR V.S. Qvalue",
     xlab = "FDR",
     ylab = "Qvalue")
abline(a=0,b=1, col = "green", lty = 2)

legend("topleft",
       legend = c("slope 1"),
       col = c("green"),
       lty = c(2),
       lwd = c(2))

# --- 開始畫 N 組的圖---
plot(data_sub_N$FDR,
     data_sub_N$qvalue,
     pch = ".",
     col = "black",
     main = "Normal Part FDR V.S. Qvalue",
     xlab = "FDR",
     ylab = "Qvalue")
abline(a=0,b=1, col = "green", lty = 2)

legend("topleft",
       legend = c("slope 1"),
       col = c("green"),
       lty = c(2),
       lwd = c(2))




# 關閉檔案
dev.off()



## 放大圖 ----


# --- 繪圖 (合併到同一頁)
png("C:/Peter/repeatSNP_clumping_raw/r2_filter_0.8/outcome/BHplot_r2_0.8_scale.png", width = 2000, height = 1000, res = 200)
par(mfrow = c(1, 2))

# --- 開始畫 T 組的圖 (只畫縮放後的資料) ---
plot(plot_x_T_zoomed, p_T_zoomed,
      pch = ".",
      col = "black",
      xlab = sprintf("Rank (m=%d)",m_total_T),
      ylab = "P-value",
      main = "Tumor Part: All Associaiton",
      ylim = c(0, ylim_max_T)) 

lines(plot_x_T_zoomed, bh_threshold_T_zoomed, col = "red", lwd = 2)

# 標示橘色點 (在縮放範圍內，且 p <= sig_cut 的點)
sig_indices_T_zoomed <- (p_T_zoomed <= sig_cut_T)
points(plot_x_T_zoomed[sig_indices_T_zoomed], p_T_zoomed[sig_indices_T_zoomed],
       col = "orange",
       pch = ".")
abline(h = bonfi_threshold_T, col = "blue", lty = 2)


legend("topleft",
       legend = c("FDR < 0.05(BH)",
                  paste("Bonferroni (p-value <=", format(bonfi_threshold_T, digits=2), ")")
                  ),
       col = c( "red", "blue"),
       lty = c(1,2),
       lwd = c(2,2))

# --- 開始畫 N 組的圖 (只畫縮放後的資料) ---
plot(plot_x_N_zoomed, p_N_zoomed,
      pch = ".",
      col = "black",
      xlab = sprintf("Rank (m=%d)",m_total_N),
      ylab = "P-value",
      main = "Normal Part: All Associaiton",
      ylim = c(0, ylim_max_N)) # 【關鍵】強制設定 Y 軸範圍

lines(plot_x_N_zoomed, bh_threshold_N_zoomed, col = "red", lwd = 2)

sig_indices_N_zoomed <- (p_N_zoomed <= sig_cut_N)
points(plot_x_N_zoomed[sig_indices_N_zoomed], p_N_zoomed[sig_indices_N_zoomed],
       col = "orange",
       pch = ".")
abline(h = bonfi_threshold_N, col = "blue", lty = 2)

legend("topleft",
       legend = c("FDR < 0.05(BH)",
                  paste("Bonferroni (p-value <=", format(bonfi_threshold_N, digits=2), ")")
                  ),
       col = c( "red", "blue"),
       lty = c(1,2),
       lwd = c(2,2))

# 關閉檔案
dev.off()


cat("How many qvalue<0.05 in tumor part: ",plot_x_T_zoomed[sig_indices_T_zoomed] %>% length())
cat("How many qvalue<0.05 in normal part: ",plot_x_N_zoomed[sig_indices_N_zoomed] %>% length())


# *FDR vs Qvalue for all pvalue*


# 設定 x,y 軸範圍
cutoff <- 0.05
file_name <- paste0("C:/Peter/repeatSNP_clumping_raw/r2_filter_0.8/outcome/FDR_qvalue_r2_0.8_scale_", cutoff, ".png")


data_N_r2[,FDR :=p.adjust(`p-value`, method = "BH") %>% 
     format(digits = 10,scientific = T) %>% 
     as.numeric()]

data_T_r2[,FDR :=p.adjust(`p-value`, method = "BH") %>% 
     format(digits = 10,scientific = T) %>% 
     as.numeric()]

data_N_r2[,qvalue := qvalue(`p-value`)$qvalues]
data_scale_N <- data_N_r2 %>% 
  filter(FDR<0.05)

data_T_r2[,qvalue := qvalue(`p-value`)$qvalues]
data_scale_T <- data_T_r2 %>% 
  filter(FDR<0.05)

setkey(data_T_r2, FDR, qvalue)
setkey(data_N_r2, FDR, qvalue)

# 4個取一個
data_sub_T <- data_T_r2[c(seq(1,nrow(data_T_r2), by=4), nrow(data_T_r2)),]
data_sub_N <- data_N_r2[c(seq(1,nrow(data_N_r2), by=4), nrow(data_N_r2)),]



png(file_name, width = 2000, height = 1000, res = 200)
par(mfrow = c(1, 2))

# --- 開始畫 T 組的圖---
plot(data_scale_T$FDR,
     data_scale_T$qvalue,
     pch = ".",
     col = "black",
     main = "Tumor Part FDR V.S. Qvalue",
     xlab = "FDR",
     ylab = "Qvalue",
     ylim = c(0,cutoff),
     xlim = c(0,cutoff))
abline(a=0,b=1, col = "green", lty = 2)

legend("topleft",
       legend = c("slope 1"),
       col = c("green"),
       lty = c(2),
       lwd = c(2))

# --- 開始畫 N 組的圖---
plot(data_scale_N$FDR,
     data_scale_N$qvalue,
     pch = ".",
     col = "black",
     main = "Normal Part FDR V.S. Qvalue",
     xlab = "FDR",
     ylab = "Qvalue",
     ylim = c(0,cutoff),
     xlim = c(0,cutoff))
abline(a=0,b=1, col = "green", lty = 2)

legend("topleft",
       legend = c("slope 1"),
       col = c("green"),
       lty = c(2),
       lwd = c(2))



# 關閉檔案
dev.off()



