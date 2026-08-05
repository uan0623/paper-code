
options(error = stop)

## Setup ----
library(data.table)
library(dplyr)
library(tidyverse)
library(stringr) 
library(ggplot2)
library(ggrepel)
library(matrixStats)  # rowMins
library(qvalue)




### SNP Lift ----

gt_N <- sprintf(
  "C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_%s/outcome/OQN_maf_gt_N_pvalue_FDR_R2_%s.txt",
  r2_threshold, r2_threshold
) %>% fread()
gt_N[, CHR := str_extract(SNP, ".*(?=:)")]


fwrite(data.table(a=paste0("chr",gt_N$CHR),
                  b=str_extract(gt_N$SNP, "(?<=\\:)\\d+"),
                  c=str_extract(gt_N$SNP, "(?<=\\:)\\d+")),
       sprintf("C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_%s/trash/snp_r2filter_%s.txt",r2_threshold,r2_threshold),
       row.names = F, col.names = F, sep = "\t")


# liftOver 轉換
cmd <- paste(
  "/mnt/d/oral_cancer/leftover/liftOver",
  sprintf("/mnt/c/Peter/OQN_FIXpeople_before_eQTL/r2_filter_%s/trash/snp_r2filter_%s.txt",r2_threshold,r2_threshold),
  "/mnt/d/oral_cancer/leftover/hg18ToHg38.over.chain",
  sprintf("/mnt/c/Peter/OQN_FIXpeople_before_eQTL/r2_filter_%s/trash/snp_r2filter_hg18TOhg38.bed",r2_threshold) ,
  sprintf("/mnt/c/Peter/OQN_FIXpeople_before_eQTL/r2_filter_%s/trash/hg18TOhg38_unmapped.bed",r2_threshold)
)

system2("wsl", cmd)

a <- sprintf("C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_%s/trash/snp_r2filter_hg18TOhg38.bed",r2_threshold) %>% 
  fread(header = F)
# 忽略轉失敗 .bed 檔案的 #Deleted in new 該行

unmapped_path <- sprintf(
  "C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_%s/trash/hg18TOhg38_unmapped.bed",
  r2_threshold
)

if (!file.exists(unmapped_path) || file.info(unmapped_path)$size == 0) {
  a_unmapped <- data.table(V1 = character(), V2 = integer(), V3 = integer())
} else {
  a_unmapped <- fread(
    unmapped_path,
    header = FALSE,
    comment.char = "#"
  )
}

# "chr22:123-124" -> "22:123"
trans_success <-  paste0(a$V1,":",a$V2,"-",a$V2) %>% 
  str_extract( "(\\d+:\\d+)(?=\\-)")

trans_fail <-  paste0(a_unmapped$V1,":",a_unmapped$V2,"-",a_unmapped$V2) %>% 
  str_extract( "(\\d+:\\d+)(?=\\-)")

# 轉失敗的是 c("1:148546405", "1:148567630")，把除了這些以外的 row 依序放入轉換成功的snp in hg38
gt_N[,snp_hg38 := NA_character_]
index_hg38_suc <- ! gt_N$SNP %in% trans_fail
gt_N[index_hg38_suc, snp_hg38 := trans_success]


finngen_previous[, snp:= paste0(`#chrom`,":",pos)]

# finngen_previous 裡有一個snp 不同 alt,pval 的情況，對重複snp 取 pval最小的 
setorder(finngen_previous, pval)
finngen_unique <- finngen_previous[, .SD[1], by = snp]
fwrite(finngen_unique,
       "D:/oral_cancer/TWAS/summary_stat/FinnGen/C3_ORALCAVITY_EXALLC_SNPunique.txt",
       row.names = F, col.names = T, sep = "\t")


# 找出 gt_N$snp_hg38  ==  finngen_unique$snp 的 row，把這些row 的 finngen_unique$pval 寫進 gt_N 新變數 finngen_pval，沒對上的 finngen_pval= NA，使用前確定 finngen_unique$snp 元素唯一出現
gt_N[
  finngen_unique,
  on = .(snp_hg38 = snp),
  finngen_pval := i.pval
]

# 2個snp 沒轉成功，snp number: 13518 -> 13516，其中 9421 個有出現在 finngen OC GWAS
uniqueN(gt_N$SNP, na.rm = T)
uniqueN(gt_N$snp_hg38, na.rm = T)
(unique(gt_N$snp_hg38) %in% finngen_unique$snp) %>% which() %>% length()



# 確認 點 "1:148546405", "1:148567630" snp_hg38=na not null，再存檔
fwrite(gt_N,
       sprintf("C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_%s/outcome/OQN_maf_gt_N_pvalue_FDR_MixFinngenPval.txt",r2_threshold),
       row.names = F, col.names = T, sep = "\t")






### Overlap ----
gt_N <- sprintf("C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_%s/outcome/OQN_maf_gt_N_pvalue_FDR_MixFinngenPval.txt",r2_threshold) %>% fread()


gt_N[, pos:= str_extract(snp_hg38 , "(?<=\\:)\\d+")]
gt_N[, start:= as.numeric(pos)]
gt_N[, end:= as.numeric(pos)]

gt_N <- gt_N %>% 
  select(c("snp_hg38","CHR","start","end"))

gt_N <- gt_N[!is.na(start) & !is.na(end)]
setnames(gt_N,"snp_hg38", "snp_gt_N")

# 調整範圍正負1000KB (10^6 bp)
gt_N[, start := ifelse(start-1e6>0, start-1e6, 0)]
gt_N[, end := end+1e6]

setnames(finngen_unique,"#chrom", "CHR")
finngen_unique[, start:= pos]
finngen_unique[, end:= pos]

finngen_unique <- finngen_unique %>% 
  select(c("snp","CHR","start","end","ref","alt","pval","mlogp","beta","sebeta"))

finngen_unique <- finngen_unique[!is.na(start) & !is.na(end)]

# 要先排序
setkey(finngen_unique, CHR, start, end)
setkey(gt_N, CHR, start, end)

# 對每個 finngen_unique 資料，找他落在哪些gt_N$snp_hg38 的區間裡，finngen_unique$CHR == gt_N$CHR AND finngen_unique$pos >= gt_N$start AND finngen_unique$pos <= gt_N$end
overlap <- finngen_unique[gt_N, 
                          on = .(CHR, start >= start, start <= end), 
                          nomatch = 0L,
                          allow.cartesian = TRUE]


### Chr Save ----
overlap[, CHR := as.character(CHR)]

# 取得所有染色體
chr_list <- unique(overlap$CHR)

# 刪掉舊檔，避免數量錯誤
old_files <- sprintf(
  "C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_%s/outcome/finngen_cisSNP_gt_N_chr%d.txt",
  r2_threshold, 1:22
)

file.remove(old_files[file.exists(old_files)])


# 逐染色體處理
for(chr_i in chr_list){
  
  # 只取當前染色體的資料
  finngen_chr <- overlap[CHR == chr_i]
  
  # 如果有 start/end 欄位，用 foverlaps 前先刪掉 NA
  # finngen_chr <- finngen_chr[!is.na(start) & !is.na(end)]
  
  # 存檔
  fwrite(finngen_chr,
         sprintf("C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_%s/outcome/finngen_cisSNP_gt_N_chr%s.txt",r2_threshold, chr_i),
         sep = "\t",
         row.names = FALSE,
         col.names = TRUE)
  
  # 釋放記憶體
  rm(finngen_chr)
  gc()
  
  cat("chr", chr_i, "done\n")
}




# 範圍調回來 ####

for (i in 1:22) {
  file_path <- sprintf("C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_%s/outcome/finngen_cisSNP_gt_N_chr%d.txt",r2_threshold, i)
  
  # --- 檢查檔案是否存在，不存在則跳過 ---
  if (!file.exists(file_path)) {
    message(sprintf(" 第 %d 號染色體檔案不存在，已跳過。", i))
    next
  }
  
  overlaps <-  fread(file_path,header = T)
  
  overlaps[, start:= str_extract(snp_gt_N , "(?<=\\:)\\d+") %>% as.numeric()]
  overlaps[, end:= start]
  setkey(overlaps, CHR, start)
  
  fwrite(overlaps,
         sprintf("C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_%s/outcome/finngen_cisSNP_gt_N_chr%d.txt",r2_threshold, i),
         sep = "\t",
         row.names = FALSE,
         col.names = TRUE)
}





# 每個 snp +-1000KB 選個pval 最小的 snp ####
data_list <- vector("list", 22)

for (chr_i in 1:22) {
  file_path <- sprintf("C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_%s/outcome/finngen_cisSNP_gt_N_chr%d.txt",r2_threshold, chr_i)
  
  # --- 檢查檔案是否存在，不存在則跳過 ---
  if (!file.exists(file_path)) {
    message(sprintf(" 第 %d 號染色體檔案不存在，已跳過。", chr_i))
    next
  }
  
  a <- fread(file_path,header = T)
  
  setkey(a, snp_gt_N, pval)
  data_list[[chr_i]] <- a[,  .SD[1], by = snp_gt_N] 
}

df <- rbindlist(data_list, use.names = TRUE)

fwrite(df,
       sprintf("C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_%s/outcome/finngen_sigSNP_chr1-22.txt",r2_threshold),
       row.names = F, col.names = T, sep = "\t")






### Nearest SNP ----
data_list <- vector("list", 22)
for (chr_i in 1:22) {
  
  file_path <- sprintf("C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_%s/outcome/finngen_cisSNP_gt_N_chr%d.txt",r2_threshold, chr_i)
  
  # --- 檢查檔案是否存在，不存在則跳過 ---
  if (!file.exists(file_path)) {
    message(sprintf(" 第 %d 號染色體檔案不存在，已跳過。", chr_i))
    next
  }
  
  a <- fread(file_path,header = T)
  
  setkey(a, snp_gt_N, pval)
  data_list[[chr_i]] <- a[,  head(.SD, 10), by = snp_gt_N]
}
df <- rbindlist(data_list, use.names = TRUE)
df[, start:= str_extract(snp_gt_N , "(?<=\\:)\\d+") %>% as.numeric()]
df[, end:= start]
df[, dist := abs(start- start.1)]

setkey(df, CHR, start,pval, dist)

# seq_len(.N) 產生 1: nrow(df)，grp 是一個「每 10 row 一組」的 group index，grp= 0,...,0,1,...,1,...
# 每個 grp 中，比較是否等於第一個 pval
df_1 <- df[
  , grp := (seq_len(.N) - 1) %/% 10
][
  , keep := pval == pval[1], by = grp
][
  keep == TRUE
]




# 最多出現4 個一樣都是最顯著的pval
# MOST_snp_hg38 is chr:pos:ref:alt in hg38
if(df_1[, .N, by = snp_gt_N][, max(N)]>10){
  stop("head(.SD, 10) 不該只選前10 個")
}


df_1[,"MOST_snp_hg38" := paste0(snp,":",ref,":",alt)]
df_1[,"MOST_snp_nearest_hg38" := MOST_snp_hg38]
df_1[,c("start.1" ,"snp","ref","alt","mlogp","beta","sebeta","dist","keep","grp") := NULL]


# 分開
df_most <- df_1[, .SD[1], by=snp_gt_N]
df_others <- setdiff(df_1,df_most)

# df_others 中將相同 snp_gt_N 和 pval 的 MOST_snp_hg38 用逗號串接起來
df_others_collapsed <- df_others[, .(
  others_combined = paste(MOST_snp_hg38, collapse = ", ")
), by = .(snp_gt_N, pval)]

# 將 df_others_collapsed 的內容併入 df_most
df_most[df_others_collapsed, 
        on = .(snp_gt_N, pval), 
        MOST_snp_hg38 := paste0(MOST_snp_hg38, ", ", i.others_combined)]

setnames(df_most, old=c("pval","snp_gt_N"), new=c("MOST_snp_hg38_finngen_pval","snp_hg38"))

fwrite(df_most,
       sprintf("C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_%s/outcome/OQN_maf_gt_N_pvalue_FDR_MixFinngenPval_5.txt",r2_threshold),
       row.names = F, col.names = T, sep = "\t")





### FDR ----




df_most <- sprintf("C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_%s/outcome/OQN_maf_gt_N_pvalue_FDR_MixFinngenPval_5.txt",r2_threshold) %>% fread()

finngen[,"snp_record_finngen" := paste0(snp_hg38,":",ref,":",alt)]

# 將 finngen 的內容併入 df_most
df_most <- finngen[df_most, on = .(snp_hg38)]


setnames(df_most,
         old=c("pval", "snp_record_finngen"),
         new=c("snp_hg38_finngen_pval", "snp_hg38_finngen"))

df_most_1 <- df_most[which(!is.na(snp_hg38_finngen_pval)),]
setkey(df_most_1,CHR, start)



# 算 OQN_gt_N 9421 snp FDR ####
# 要用 1:123 合併，才抓的到不同 ref,alt 的snp
common_snp <- df_most_1[, .SD[1], by=snp_hg38_finngen]
common_snp <- df_most_1 %>% 
  select("snp_hg38_finngen","snp_hg38_finngen_pval")



common_snp[,snp_finngen_FDR :=p.adjust(snp_hg38_finngen_pval, method = "BH") %>% 
             format(digits = 10,scientific = T) %>% 
             as.numeric()]
common_snp[,snp_finngen_qvalue := qvalue(snp_hg38_finngen_pval)$qvalues]

common_snp[,snp_finngen_Bonfi := 
             ifelse(snp_hg38_finngen_pval<(0.05/nrow(common_snp)),
                    1,0)]
df_most <- common_snp[df_most, on = .(snp_hg38_finngen, snp_hg38_finngen_pval)]



rm(df_most_1)


# 算 most pvalue 261 snp FDR，將 finngen 的內容併入 ####
# 這不用抓不同 ref,alt 的snp，因為我們只要 pval 最小的，確認過了ref or alt 不同，pval 就不同
MOSTpval_info <- df_most %>% select(MOST_snp_nearest_hg38)
MOSTpval_info <- unique(MOSTpval_info)
setnames(MOSTpval_info, old = "MOST_snp_nearest_hg38",new = "snp_record_finngen")
MOSTpval_info <- finngen[MOSTpval_info, on = .(snp_record_finngen)]

setnames(MOSTpval_info,
         old=c("pval", "snp_record_finngen"),
         new=c("MOST_snp_nearest_pvalue", "MOST_snp_nearest_hg38"))

MOSTpval_info[,c("snp_hg38","hg19_snpID","alt","ref"):=NULL]
MOSTpval_info[,MOST_snp_hg38_finngen_FDR :=
                p.adjust(MOST_snp_nearest_pvalue, method = "BH") %>% 
                format(digits = 10,scientific = T) %>% 
                as.numeric()]
MOSTpval_info[,MOST_snp_hg38_finngen_Bonfi := 
                ifelse(MOST_snp_nearest_pvalue<  (0.05/nrow(MOSTpval_info)), 1,0)]


# record 紀錄 excel 數量，已經確認過沒有空值, NA, snp重複 等問題
print(nrow(MOSTpval_info[MOST_snp_hg38_finngen_FDR <0.05,]))
print(nrow(MOSTpval_info[MOST_snp_hg38_finngen_Bonfi==1,]))




MOSTpval_info <- MOSTpval_info[df_most, on = .(MOST_snp_nearest_hg38)]
fwrite(data.table(a=paste0("chr",MOSTpval_info$CHR),
                  b=str_extract(MOSTpval_info$snp_hg38, "(?<=\\:)\\d+"), 
                  c=str_extract(MOSTpval_info$snp_hg38, "(?<=\\:)\\d+")) ,
       sprintf("C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_%s/trash/OQN_MixFinngenPval_5_SNPhg38.txt",r2_threshold),
       row.names = F, col.names = F, sep = "\t")

# liftOver 轉換
# 檔案 snp_r2filter.txt 放進 leftover 轉換成 snp_r2filter_success.bed
cmd <- paste(
  "/mnt/d/oral_cancer/leftover/liftOver",
  sprintf("/mnt/c/Peter/OQN_FIXpeople_before_eQTL/r2_filter_%s/trash/OQN_MixFinngenPval_5_SNPhg38.txt",r2_threshold),
  "/mnt/d/oral_cancer/leftover/hg38ToHg19.over.chain",
  sprintf("/mnt/c/Peter/OQN_FIXpeople_before_eQTL/r2_filter_%s/trash/OQN_MixFinngenPval_5_SNPhg38TOhg19_success.bed",r2_threshold) ,
  sprintf("/mnt/c/Peter/OQN_FIXpeople_before_eQTL/r2_filter_%s/trash/SNPhg38TOhg19_unmapped.bed",r2_threshold)
)

system2("wsl", cmd)

a <- sprintf("C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_%s/trash/OQN_MixFinngenPval_5_SNPhg38TOhg19_success.bed",r2_threshold) %>% 
  fread(header = F)


unmapped_path <- sprintf(
  "C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_%s/trash/SNPhg38TOhg19_unmapped.bed",
  r2_threshold
)

if (!file.exists(unmapped_path) || file.info(unmapped_path)$size == 0) {
  a_unmapped <- data.table(V1 = character(), V2 = integer(), V3 = integer())
} else {
  a_unmapped <- fread(
    unmapped_path,
    header = FALSE,
    comment.char = "#"
  )
}

names(a) <- c("chr","start","end")
a[,chr_nochr := str_extract(chr, "\\d+")]
MOSTpval_info_hg19 <- paste0(a$chr_nochr,":",a$start)

trans_fail <- paste0(a_unmapped$V1, ":", a_unmapped$V2, "-", a_unmapped$V2) %>%
  str_extract("(\\d+:\\d+)(?=\\-)")

MOSTpval_info[, hg19_snpID := NA_character_]
index_hg19_suc <- !MOSTpval_info$snp_hg38 %in% trans_fail
MOSTpval_info[index_hg19_suc, hg19_snpID := MOSTpval_info_hg19]



# 刪變數 ####

# MOST_snp_nearest_hg38 比較近
MOSTpval_info[,c( "MOST_snp_hg38_finngen_pval","snp_hg38_finngen") := NULL]
setnames(MOSTpval_info,
         old=c("snp_hg38", "hg19_snpID", "snp_hg38_finngen_pval"),
         new=c("gt_N_hg38", "gt_N_hg19", "gt_N_finngen_pval"))


# 把想要的 "`#chrom`","pos","snp" 移前面
setcolorder(MOSTpval_info,
            c("MOST_snp_nearest_pvalue","MOST_snp_nearest_hg38","MOST_snp_hg38_finngen_FDR","MOST_snp_hg38_finngen_Bonfi","gt_N_finngen_pval","snp_finngen_FDR","snp_finngen_qvalue","snp_finngen_Bonfi","gt_N_hg38","gt_N_hg19","alt","ref","CHR","start","end","MOST_snp_hg38")
)

fwrite(MOSTpval_info,
       sprintf("C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_%s/outcome/OQN_maf_gt_N_pvalue_FDR_MixFinngenPval_6.txt",r2_threshold),
       row.names = F, col.names = T, sep = "\t")


# 重新讀取，才會把 MOSTpval_info$alt 從 NA 讀成 ""
MOSTpval_info <- sprintf("C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_%s/outcome/OQN_maf_gt_N_pvalue_FDR_MixFinngenPval_6.txt",r2_threshold) %>% fread()

# 紀錄不在 finngen snp
not_in_finngen <- MOSTpval_info[MOSTpval_info$alt=="",c("gt_N_hg19","gt_N_hg38")]
fwrite(data.table(k = not_in_finngen$gt_N_hg19),
       sprintf("C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_%s/outcome/not_in_finngen_hg19SNP.txt",r2_threshold),
       row.names = F, col.names = F, sep = "\n") 





# 刪掉舊檔案
OUT_DIR <- sprintf(
  "C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_%s/trash/no_intersect_LD/EUR",
  r2_threshold
)

old_files <- list.files(
  OUT_DIR,
  pattern = "^EUR_LDchr[0-9]+\\.",
  full.names = TRUE
)

if (length(old_files) > 0) {
  file.remove(old_files)
}

## EUR LD ----

plink <- "C:/Program Files/plink/plink.exe"
EQTL_LIST <- sprintf(
  "C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_%s/outcome/not_in_finngen_hg19SNP.txt",
  r2_threshold
)
BFILE_DIR <- "C:/Peter/PCA_1000G_20130502/trash/deal_repeatSNP"

OUT_DIR <- sprintf(
  "C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_%s/trash/no_intersect_LD/EUR",
  r2_threshold
)

for (chr in 1:22) {
  args <- c(
    "--bfile", paste0(BFILE_DIR, "/chr_", chr, "_rename"),
    "--ld-snp-list", EQTL_LIST,
    "--r2",
    "--ld-window-kb", "1000",
    "--ld-window", "999999",
    "--ld-window-r2", "0",
    "--list-all",
    "--memory", "12000",
    "--out", paste0(OUT_DIR, "/EUR_LDchr", chr)
  )
  
  system2(plink, args = args)
  message("第 ", chr, " 號染色體完成")
  gc()
}


# 刪掉舊檔案
OUT_DIR <- sprintf(
  "C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_%s/trash/no_intersect_LD/FIN",
  r2_threshold
)

old_files <- list.files(
  OUT_DIR,
  pattern = "^FIN_LDchr[0-9]+\\.",
  full.names = TRUE
)

if (length(old_files) > 0) {
  file.remove(old_files)
}

## FIN LD ----

plink <- "C:/Program Files/plink/plink.exe"
EQTL_LIST <- sprintf(
  "C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_%s/outcome/not_in_finngen_hg19SNP.txt",
  r2_threshold
)
BFILE_DIR <- "C:/Peter/PCA_1000G_20130502/FIN_sample/trash/deal_repeatSNP"

OUT_DIR <- sprintf(
  "C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_%s/trash/no_intersect_LD/FIN",
  r2_threshold
)

for (chr in 1:22) {
  args <- c(
    "--bfile", paste0(BFILE_DIR, "/chr_", chr, "_rename"),
    "--ld-snp-list", EQTL_LIST,
    "--r2",
    "--ld-window-kb", "1000",
    "--ld-window", "999999",
    "--ld-window-r2", "0",
    "--list-all",
    "--memory", "12000",
    "--out", paste0(OUT_DIR, "/FIN_LDchr", chr)
  )
  
  system2(plink, args = args)
  message("第 ", chr, " 號染色體完成")
  gc()
}







### Combine ----

gt_N_MOSTsnp <- sprintf("C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_%s/outcome/finngen_sigSNP_chr1-22.txt",r2_threshold) %>% fread()

# 併進 gt_N，gt_N_MOSTsnp$snp_gt_N 對應到 gt_N$snp_hg38
gt_N <- sprintf("C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_%s/outcome/OQN_maf_gt_N_pvalue_FDR_MixFinngenPval.txt",r2_threshold) %>% fread()

gt_N_MOSTsnp <- gt_N_MOSTsnp %>%
  select(snp_gt_N,snp)

names(gt_N_MOSTsnp) <- c("snp_hg38","1000kb_MOST_snp")
a <- merge(gt_N, gt_N_MOSTsnp, by = "snp_hg38",  all.x = TRUE)
fwrite(a,
       sprintf("C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_%s/outcome/OQN_maf_gt_N_pvalue_FDR_MixFinngenPval_2.txt",r2_threshold),
       row.names = F, col.names = T, sep = "\t")








### Genome Lift ----


maf_hg1819 <- fread("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/cis_hg18_19_snp_MAF.txt")

gt_N <- sprintf("C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_%s/outcome/OQN_maf_gt_N_pvalue_FDR_MixFinngenPval_2.txt",r2_threshold) %>% fread()

# 確認所有 gt_N snp(hg18) in maf_hg1819$hg18_snpID
all(gt_N$SNP %in% maf_hg1819$hg18_snpID)

a <- merge(gt_N, maf_hg1819 %>%
             select("hg18_snpID", "hg19_snpID"),
           by.x = "SNP",by.y ="hg18_snpID",   all.x = TRUE)


fwrite(data.table(kk = a$hg19_snpID %>% 
                    unique()),
       sprintf("C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_%s/outcome/OQN_eQTL_SNPhg19.txt",r2_threshold),
       row.names = F, col.names = T, sep = "\t")




## EUR ----
### Merge ----

ld_list <- vector("list", 22)

for (i in 1:22) {
  file_path <- sprintf("C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_%s/trash/no_intersect_LD/EUR/EUR_LDchr%d.ld",r2_threshold, i)
  
  # --- 檢查檔案是否存在，不存在則跳過 ---
  if (!file.exists(file_path)) {
    message(sprintf(" 第 %d 號染色體檔案不存在，已跳過。", i))
    next
  }
  
  a <-  fread(file_path,header = T)
  
  # 把 1:123_2 轉成 1:123
  a[,SNP_B_noRepeat:= str_extract(SNP_B, "^[^_]+")]
  # 挑出SNP_A 是 not in finngen，SNP_B 是 in finngen 的組合 (以確認所有SNP_A 都是 not in finngen)
  a <- a[which(SNP_B_noRepeat %in% finngen$hg19_snpID),]
  a[,SNP_B_noRepeat:= NULL]
  
  # 按snp_a, R2從大排到小，用row index 的方式選sub data，省記憶體
  # 挑前180個LD高的，因為有 23 snp 附近有 148 個 LD 一樣大的 snp
  setorder(a, SNP_A, BP_A, -R2)
  ld_list[[i]] <- a[, head(.SD, 180), by = SNP_A]
}

df <- rbindlist(ld_list, use.names = TRUE)

# plink1.9 cmd 不會計算自己跟自己的 LD
df[, dist := abs(BP_A- BP_B)]
setorder(df, CHR_A, BP_A, -R2, dist)

# seq_len(.N) 產生 1: nrow(df)，grp 是一個「每 180 row 一組」的 group index，grp= 0,...,0,1,...,1,...
# 每個 grp 中，比較是否等於第一個 pval
df_1 <- df[
  , grp := (seq_len(.N) - 1) %/% 180
][
  , keep := R2 == R2[1], by = grp
][
  keep == TRUE
]

# check 前面選 180 合理，最多附近有 148個snp 相同 LD (不一定是LD=1)
if(df_1[,.N, by=SNP_A][, max(N)]>180){
  stop("head(.SD, 180) 不該只選前 180 個")
}


# 跟finngen 合併 ####

# 把 16:2045121_2 改成 16:2045121，才能在 finngen 找到 info
setnames(df_1, old = "SNP_B", new = "hg19_snpID")
df_1[, SNP_A := str_extract(SNP_A, "^\\d+:\\d+")]
df_1[, hg19_snpID := str_extract(hg19_snpID, "^\\d+:\\d+")]
df_2 <- finngen[df_1, on= .(hg19_snpID)]

# "mostLD_snp_hg19" := hg19_snpID:ref:alt
df_2[,"mostLD_snp_hg19" := paste0(hg19_snpID,":",ref,":",alt)]


setnames(df_2,
         old = c("snp_hg38","CHR_A","SNP_A","BP_A","R2","pval"),
         new = c("mostLD_snp_hg38", "CHR","snp_gt_N","start","LD","mostLD_snp_hg19_finngen_pval"))

setorder(df_2, CHR, start, -LD, dist, mostLD_snp_hg19_finngen_pval)

df_2[,"mostLD_snp_nearest_hg19" := mostLD_snp_hg19]
df_2[,c("ref","alt","CHR_B","BP_B","dist","keep","grp","hg19_snpID") := NULL]


# 每個gt_N 取最近的高LD snp，再把同樣高LD snp 但距離遠的併入 #### 

# 每個gt_N 取最近的高LD snp
df_most <- df_2[, .SD[1], by=snp_gt_N]
# 其他同樣高LD snp
df_others <- setdiff(df_2,df_most)

# df_others 中將相同 snp_gt_N 和 pval 的 MOST_snp_hg38 用逗號串接起來
df_others_collapsed <- df_others[, .(
  others_combined = paste(mostLD_snp_hg19, collapse = ", ")
), by = .(snp_gt_N, LD)]

# 將 df_most 的內容併入 df_others_collapsed
df_most[df_others_collapsed, 
        on = .(snp_gt_N, LD), 
        mostLD_snp_hg19 := paste0(mostLD_snp_hg19, ", ", i.others_combined)]



rm(a,i,df_others_collapsed,df_others,ld_list,df)



### FDR ----


most_pval_info <- sprintf("C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_%s/outcome/OQN_maf_gt_N_pvalue_FDR_MixFinngenPval_6.txt",r2_threshold) %>% fread()

# 紀錄不在 finngen snp
not_in_finngen <- most_pval_info[most_pval_info$alt=="",c("gt_N_hg19","gt_N_hg38")]

# df_most$snp_gt_N 是 hg19
setnames(df_most,
         old = c("snp_gt_N","mostLD_snp_hg38"),
         new = c("gt_N_hg19","mostLD_snp_nearest_hg38"))


# 這11 snp 不在finngen，且附近1000 kb 沒有 finngen 的 snp (前面算R2就連R2=0都會算出)
setdiff(not_in_finngen$gt_N_hg19, df_most$gt_N_hg19 )

fwrite(data.table(k = setdiff(not_in_finngen$gt_N_hg19, df_most$gt_N_hg19 )),
       sprintf("C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_%s/outcome/not_in_finngen_1000gEUR_hg19SNP.txt",r2_threshold),
       row.names = F, col.names = F, sep = "\n") 



# 跟 maf_gt_N_pvalue_FDR_MixFinngenPval_6.txt 紀錄的gt_N, most_pval 結果合併 ####
most_pval_info <- sprintf("C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_%s/outcome/OQN_maf_gt_N_pvalue_FDR_MixFinngenPval_6.txt",r2_threshold) %>% fread()

most_pval_info[,c("CHR","start"):= NULL]
a <- df_most[most_pval_info, on=.(gt_N_hg19)]


# 確認了，只有23 snp 不在finngen，且附近1000 kb 沒有 finngen 的 snp，才會 mostLD_snp_hg19_finngen_pval, gt_N_finngen_pval 兩個都是 NA，其他snp 都只有一個是 NA
nofin_LD <-  setdiff(not_in_finngen$gt_N_hg19, df_most$gt_N_hg19 )
test <- intersect(which(is.na(a$mostLD_snp_hg19_finngen_pval)), which(is.na(a$gt_N_finngen_pval))) %>%
  a[.,gt_N_hg19]
all(test==nofin_LD)






# 把gt_N_finngen_pval , mostLD_snp_hg19_finngen_pval 合併算 FDR ####
most_LD <- df_most[,.SD[1],by=mostLD_snp_nearest_hg19] %>% 
  select(gt_N_hg19,mostLD_snp_hg19_finngen_pval,mostLD_snp_nearest_hg19)

setnames(most_LD,old = "mostLD_snp_hg19_finngen_pval", new = 'mostLD_or_finngen_pval')

most_pval_info_unique <- most_pval_info[,.SD[1],by=gt_N_hg19] %>%
  select(gt_N_hg19,gt_N_finngen_pval)
setnames(most_pval_info_unique,old = 'gt_N_finngen_pval', new = 'mostLD_or_finngen_pval')
most_pval_info_unique[,mostLD_snp_nearest_hg19 := NA_character_ ]




# 扣掉 4095 snp 不在finngen, pval =NA 的
most_pval_info_unique <- most_pval_info_unique[which(!is.na(mostLD_or_finngen_pval)),]


most_LD[,mostLD_nearest_hg19 := str_extract(mostLD_snp_nearest_hg19, "^[^:]+:[^:]+")]
test <- most_LD[!mostLD_nearest_hg19 %in% most_pval_info_unique$gt_N_hg19]
test[,mostLD_nearest_hg19 := NULL]

record_in_and_LD_pval <- rbind(most_pval_info_unique,test)





record_in_and_LD_pval[,mostLD_or_finngen_FDR :=p.adjust(mostLD_or_finngen_pval,
                                                        method = "BH") %>% 
                        format(digits = 10,scientific = T) %>% 
                        as.numeric()]

record_in_and_LD_pval[,mostLD_or_finngen_qvalue := qvalue(mostLD_or_finngen_pval)$qvalues]
record_in_and_LD_pval[,mostLD_or_finngen_Bonfi := 
                        ifelse(mostLD_or_finngen_pval<(0.05/nrow(record_in_and_LD_pval)),
                               1,0)]






# 不同 gt_N snp 對到同個 LD snp，併入 ####

target_cols <- c("mostLD_or_finngen_pval", 
                 "mostLD_or_finngen_FDR", 
                 "mostLD_or_finngen_qvalue", 
                 "mostLD_or_finngen_Bonfi")

# record_in_and_LD_pval 併入 df_most 
df_result <- record_in_and_LD_pval[, c("mostLD_snp_nearest_hg19", target_cols), with = FALSE][
  df_most, on = "mostLD_snp_nearest_hg19"
]

# 選出紀錄在 finngen的 8599 snp  FDR, qvalue 數值，跟 不在finngen 但有 LD snp 的 3366 snp FDR 數值合併
test <- rbind(record_in_and_LD_pval[which(is.na(mostLD_snp_nearest_hg19))],
              df_result %>% 
                select(names(record_in_and_LD_pval)))
# 更新資料
a <- test[a,on=.(gt_N_hg19,mostLD_snp_nearest_hg19)]





a[,start:= NULL]
setnames(a,
         old=c("end","gt_N_finngen_pval", "snp_finngen_FDR",
               "snp_finngen_qvalue","snp_finngen_Bonfi"),
         new=c("pos","gt_N_hg38_finngen_pval", "gt_N_hg38_finngen_FDR",
               "gt_N_hg38_finngen_qvalue","gt_N_hg38_finngen_Bonfi"))

setcolorder(a,
            c("CHR","pos","gt_N_hg38","gt_N_hg19","alt", "ref",
              "gt_N_hg38_finngen_pval", "gt_N_hg38_finngen_FDR",
              "gt_N_hg38_finngen_qvalue","gt_N_hg38_finngen_Bonfi","MOST_snp_hg38",
              "MOST_snp_nearest_hg38","MOST_snp_nearest_pvalue","MOST_snp_hg38_finngen_FDR",
              "MOST_snp_hg38_finngen_Bonfi","LD","mostLD_snp_hg19","mostLD_snp_nearest_hg19",
              "mostLD_snp_nearest_hg38","mostLD_snp_hg19_finngen_pval",
              "mostLD_or_finngen_pval","mostLD_or_finngen_FDR",
              "mostLD_or_finngen_qvalue","mostLD_or_finngen_Bonfi")
)
setkey(a,CHR,pos)




fwrite(a,
       sprintf("C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_%s/outcome/OQN_MixFinngenPval_7_EUR.txt",r2_threshold),
       row.names = F, col.names = T, sep = "\t")





# 新增 MOST_snp_hg19 ####
a <-sprintf("C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_%s/outcome/OQN_MixFinngenPval_7_EUR.txt",r2_threshold) %>% fread()
all(str_extract(a$gt_N_hg38, ".*(?=:)")==str_extract(a$gt_N_hg19, ".*(?=:)"))
a[,CHR:= str_extract(gt_N_hg38, ".*(?=:)") %>% as.numeric()]

setcolorder(a,
            c("CHR","pos","gt_N_hg38","gt_N_hg19","alt","ref",
              "gt_N_hg38_finngen_pval", "gt_N_hg38_finngen_FDR",
              "gt_N_hg38_finngen_qvalue","gt_N_hg38_finngen_Bonfi","MOST_snp_hg38",
              "MOST_snp_nearest_hg38","MOST_snp_nearest_pvalue","MOST_snp_hg38_finngen_FDR",
              "MOST_snp_hg38_finngen_Bonfi","LD","mostLD_snp_hg19","mostLD_snp_nearest_hg19",
              "mostLD_snp_nearest_hg38","mostLD_snp_hg19_finngen_pval",
              "mostLD_or_finngen_pval","mostLD_or_finngen_FDR",
              "mostLD_or_finngen_qvalue","mostLD_or_finngen_Bonfi")
)
setkey(a,CHR,pos)



# 新增 MOST_snp_nearest_hg19 欄位 ####
most_snp <- a %>% 
  select(MOST_snp_nearest_hg38) %>% 
  unique() %>%
  separate(MOST_snp_nearest_hg38,
           into = c("chr", "pos", "ref", "alt"),
           sep = ":", remove = FALSE) %>% 
  as.data.table()

most_snp[,chr := as.numeric(chr)]
most_snp[,pos := as.numeric(pos)]

setorder(most_snp,chr,pos)
fwrite(data.table(a=paste0("chr",most_snp$chr),
                  b=most_snp$pos, 
                  c=most_snp$pos ) ,
       sprintf("C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_%s/trash/most_snp.txt",r2_threshold),
       row.names = F, col.names = F, sep = "\t")



# liftOver 轉換
cmd <- paste(
  "/mnt/d/oral_cancer/leftover/liftOver",
  sprintf("/mnt/c/Peter/OQN_FIXpeople_before_eQTL/r2_filter_%s/trash/most_snp.txt",r2_threshold),
  "/mnt/d/oral_cancer/leftover/hg38ToHg19.over.chain",
  sprintf("/mnt/c/Peter/OQN_FIXpeople_before_eQTL/r2_filter_%s/trash/most_SNPhg38TOhg19.bed",r2_threshold) ,
  sprintf("/mnt/c/Peter/OQN_FIXpeople_before_eQTL/r2_filter_%s/trash/hg38TOhg19_unmapped.bed",r2_threshold)
)

system2("wsl", cmd)


snp_hg19 <- sprintf("C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_%s/trash/most_SNPhg38TOhg19.bed",r2_threshold) %>% 
  fread(header=F)
snp_hg19[,chr := str_remove(V1, "chr")]


unmapped_path <- sprintf(
  "C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_%s/trash/hg38TOhg19_unmapped.bed",
  r2_threshold
)

if (!file.exists(unmapped_path) || file.info(unmapped_path)$size == 0) {
  a_unmapped <- data.table(V1 = character(), V2 = integer(), V3 = integer())
} else {
  a_unmapped <- fread(
    unmapped_path,
    header = FALSE,
    comment.char = "#"
  )
}

trans_success <- paste0(snp_hg19$chr, ":", snp_hg19$V2)
trans_fail <- paste0(a_unmapped$V1, ":", a_unmapped$V2, "-", a_unmapped$V2) %>%
  str_extract("(\\d+:\\d+)(?=\\-)")

most_snp[, MOST_snp_nearest_hg19 := NA_character_]
index_hg19_suc <- !paste0(most_snp$chr, ":", most_snp$pos) %in% trans_fail
most_snp[index_hg19_suc, MOST_snp_nearest_hg19 := trans_success]

most_snp <- as.data.table(most_snp)  
most_snp[,c("chr","pos","ref","alt" ) := NULL]
a <- most_snp[a, on= .(MOST_snp_nearest_hg38)]

fwrite(a,
       sprintf("C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_%s/outcome/OQN_MixFinngenPval_7_EUR.txt",r2_threshold),
       row.names = F, col.names = T, sep = "\t")



### Final FDR ----

k <- sprintf("C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_%s/outcome/OQN_MixFinngenPval_7_EUR.txt",r2_threshold) %>% fread()

# 挑出 LD snp
ld_snp <- k[mostLD_snp_nearest_hg19!="",
            c("mostLD_snp_nearest_hg38","mostLD_snp_hg19_finngen_pval")] %>% 
  unique()

# 挑出 gt_N record in finngen GWAS snp
gt_N <- k[mostLD_snp_nearest_hg19=="",
          c("gt_N_hg38","gt_N_hg38_finngen_pval")] %>% 
  unique()

# snp 6:31235502 重複出現在 finngen, 挑最顯著的
setorder(gt_N,gt_N_hg38_finngen_pval)
gt_N <- gt_N[, .SD[1], by="gt_N_hg38"]

# 刪掉 11 個沒有 LD snp 且沒紀錄在 GWAS 的
gt_N <- gt_N[!is.na(gt_N_hg38_finngen_pval),]


# LD snp 中有 321 snp 出現在 gt_N，刪掉
combine_ldsnp_gtN <- merge(ld_snp,gt_N,
                           by.x = c("mostLD_snp_nearest_hg38","mostLD_snp_hg19_finngen_pval"),
                           by.y = c("gt_N_hg38","gt_N_hg38_finngen_pval"))

ld_snp <- ld_snp[!mostLD_snp_nearest_hg38 %in% combine_ldsnp_gtN$mostLD_snp_nearest_hg38,]



# 合併不重複的 184 LD snp, 9421 原本就在 finngen 的 snp
setnames(ld_snp, old = c("mostLD_snp_nearest_hg38","mostLD_snp_hg19_finngen_pval"),
         new = c("snp","pval"))
setnames(gt_N, old = c("gt_N_hg38","gt_N_hg38_finngen_pval"),
         new = c("snp","pval"))
final <- rbind(gt_N, ld_snp)

final[,mostLD_or_finngen_FDR :=p.adjust(pval,method = "BH") %>% 
        format(digits = 10,scientific = T) %>% 
        as.numeric()]

final[,mostLD_or_finngen_qvalue := qvalue(pval)$qvalues]
final[,mostLD_or_finngen_Bonfi := 
        ifelse(pval<(0.05/nrow(final)),
               1,0)]


cols_to_fill <- c("mostLD_or_finngen_FDR", "mostLD_or_finngen_qvalue", "mostLD_or_finngen_Bonfi")
k[, (cols_to_fill) := NULL]


# 對 final$snp 對應 k$gt_N_hg38 的 row，final的 cols_to_fill 變數填入 k
k[final, on = .(gt_N_hg38 = snp), 
  (cols_to_fill) := .(i.mostLD_or_finngen_FDR, 
                      i.mostLD_or_finngen_qvalue, 
                      i.mostLD_or_finngen_Bonfi)]


# 第二次填入，對 final$snp 對應 k$mostLD_snp_nearest_hg38 的 row，final的 cols_to_fill 變數填入 k
k[final, on = .(mostLD_snp_nearest_hg38 = snp), 
  (cols_to_fill) := .(i.mostLD_or_finngen_FDR, 
                      i.mostLD_or_finngen_qvalue, 
                      i.mostLD_or_finngen_Bonfi)]



# 加 R2 info
cis_snp <- fread("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/cis_hg18_19_snp_MAF.txt",header=T)
cis_snp[!cis_snp$REF %in% c("A", "T", "C", "G") | !cis_snp$ALT %in% c("A", "T", "C", "G"),
        rsID := NA_character_]

cis_snp <- cis_snp[,c("hg19_snpID","rsID","R2","ER2","impute_type")]
setnames(cis_snp, old = "hg19_snpID", new = "gt_N_hg19")
k <- cis_snp[k, on= .(gt_N_hg19)]

setcolorder(k,c("CHR",	"pos",	"gt_N_hg38",	"gt_N_hg19","rsID","R2","ER2","impute_type",	"alt",	"ref",
                "gt_N_hg38_finngen_pval","gt_N_hg38_finngen_FDR",	"gt_N_hg38_finngen_qvalue",
                "gt_N_hg38_finngen_Bonfi","MOST_snp_hg38",	"MOST_snp_nearest_hg38",	"MOST_snp_nearest_hg19",
                "MOST_snp_nearest_pvalue","MOST_snp_hg38_finngen_FDR",	"MOST_snp_hg38_finngen_Bonfi",
                "LD",	"mostLD_snp_hg19",	"mostLD_snp_nearest_hg19",	"mostLD_snp_nearest_hg38",
                "mostLD_snp_hg19_finngen_pval",	"mostLD_or_finngen_pval",	"mostLD_or_finngen_FDR",
                "mostLD_or_finngen_qvalue",	"mostLD_or_finngen_Bonfi"))

fwrite(k,
       sprintf("C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_%s/outcome/OQN_MixFinngenPval_7_EUR.txt",r2_threshold),
       row.names = F, col.names = T, sep = "\t")


# 造過 eQTL bon 門檻的 3 手法結果 ----
gt_N <- sprintf(
  "C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_%s/outcome/OQN_N_bon_R2_%s.txt",
  r2_threshold, r2_threshold
) %>% fread()
setnames(gt_N, old="SNP", new="hg18_snpID")
snp_number <- nrow(gt_N)

# 用 maf_hg18_19 對照，抓 hg19 id
maf_hg18_19 <- fread("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/cis_hg18_19_snp_MAF.txt", header = T)[
  ,
  .(hg18_snpID, hg19_snpID)
]

gt_N <- maf_hg18_19[gt_N, on = .(hg18_snpID)]
final <- k[gt_N_hg19 %in% gt_N$hg19_snpID, ] 
fwrite(final,
       sprintf("C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_%s/outcome/OQN_MixFinngenPval_7_EUR_bon.txt",r2_threshold),
       row.names = F, col.names = T, sep = "\t") 


cat(sprintf("pass eQTL bon %d snp, retain %d in OQN_MixFinngenPval_7_EUR_bon.txt\n", snp_number, nrow(final)))


# 重算 FDR, bon, qval ----

final[!is.na(gt_N_hg38_finngen_pval),
  `:=`(
    gt_N_hg38_finngen_FDR =
      p.adjust(gt_N_hg38_finngen_pval, method = "BH"),
    gt_N_hg38_finngen_Bonfi =
      fifelse(gt_N_hg38_finngen_pval < 0.05 / .N, 1, 0)
  )
]
final[!is.na(mostLD_or_finngen_pval),
  `:=`(
    mostLD_or_finngen_FDR =
      p.adjust(mostLD_or_finngen_pval, method = "BH"),
    mostLD_or_finngen_Bonfi =
      fifelse(mostLD_or_finngen_pval < 0.05 / .N, 1, 0)
  )
  
]
final[!is.na(MOST_snp_nearest_pvalue),
  `:=`(
    MOST_snp_hg38_finngen_FDR =
      p.adjust(MOST_snp_nearest_pvalue, method = "BH"),
    MOST_snp_hg38_finngen_Bonfi =
      fifelse(MOST_snp_nearest_pvalue < 0.05 / .N, 1, 0)
  )
  
]

fwrite(final,
       sprintf("C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_%s/outcome/OQN_MixFinngenPval_7_EUR_bon.txt",r2_threshold),
       row.names = F, col.names = T, sep = "\t") 





## FIN ----


### Merge ----


  ld_list <- vector("list", 22)

  for (i in 1:22) {
    file_path <- sprintf("C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_%s/trash/no_intersect_LD/FIN/FIN_LDchr%d.ld", r2_threshold, i)

    # --- 檢查檔案是否存在，不存在則跳過 ---
    if (!file.exists(file_path)) {
      message(sprintf(" 第 %d 號染色體檔案不存在，已跳過。", i))
      next
    }

    a <- fread(file_path, header = T)

    # 把 1:123_2 轉成 1:123
    a[, SNP_B_noRepeat := str_extract(SNP_B, "^[^_]+")]
    # 挑出SNP_A 是 not in finngen，SNP_B 是 in finngen 的組合 (以確認所有SNP_A 都是 not in finngen)
    a <- a[which(SNP_B_noRepeat %in% finngen$hg19_snpID), ]
    a[, SNP_B_noRepeat := NULL]

    # 按snp_a, R2從大排到小，用row index 的方式選sub data，省記憶體
    # 挑前380 個LD高的，因為有371 snp 附近有146個 LD=1 的 snp
    setorder(a, SNP_A, BP_A, -R2)
    ld_list[[i]] <- a[, head(.SD, 380), by = SNP_A]
  }

df <- rbindlist(ld_list, use.names = TRUE)


# plink1.9 cmd 不會計算自己跟自己的 LD
df[, dist := abs(BP_A- BP_B)]
setorder(df, CHR_A, BP_A, -R2, dist)

# seq_len(.N) 產生 1: nrow(df)，grp 是一個「每 180 row 一組」的 group index，grp= 0,...,0,1,...,1,...
# 每個 grp 中，比較是否等於第一個 pval
df_1 <- df[
  , grp := (seq_len(.N) - 1) %/% 380
][
  , keep := R2 == R2[1], by = grp
][
  keep == TRUE
]

# check 前面選 380合理，最多附近有 371個snp 相同 LD (不一定是LD=1)
if(df_1[,.N, by=SNP_A][, max(N)]>380){
  stop("head(.SD, 380) 不該只選前380 個")
}


# 跟finngen 合併 ####

# 把 16:2045121_2 改成 16:2045121，才能在 finngen 找到 info
setnames(df_1, old = "SNP_B", new = "hg19_snpID")
df_1[, SNP_A := str_extract(SNP_A, "^\\d+:\\d+")]
df_1[, hg19_snpID := str_extract(hg19_snpID, "^\\d+:\\d+")]
df_2 <- finngen[df_1, on= .(hg19_snpID)]


df_2[,"mostLD_snp_hg19" := paste0(hg19_snpID,":",ref,":",alt)]


setnames(df_2,
         old = c("snp_hg38","CHR_A","SNP_A","BP_A","R2","pval"),
         new = c("mostLD_snp_hg38", "CHR","snp_gt_N","start","LD","mostLD_snp_hg19_finngen_pval"))

setorder(df_2, CHR, start, -LD, dist, mostLD_snp_hg19_finngen_pval)

df_2[,"mostLD_snp_nearest_hg19" := mostLD_snp_hg19]
df_2[,c("ref","alt","CHR_B","BP_B","dist","keep","grp","hg19_snpID") := NULL]


# 每個gt_N 取最近的高LD snp，再把同樣高LD snp 但距離遠的併入 #### 

# 每個gt_N 取最近的高LD snp
df_most <- df_2[, .SD[1], by=snp_gt_N]
# 其他同樣高LD snp
df_others <- setdiff(df_2,df_most)

# df_others 中將相同 snp_gt_N 和 pval 的 MOST_snp_hg38 用逗號串接起來
df_others_collapsed <- df_others[, .(
  others_combined = paste(mostLD_snp_hg19, collapse = ", ")
), by = .(snp_gt_N, LD)]

# 將 df_most 的內容併入 df_others_collapsed
df_most[df_others_collapsed, 
        on = .(snp_gt_N, LD), 
        mostLD_snp_hg19 := paste0(mostLD_snp_hg19, ", ", i.others_combined)]



rm(a,i,df_others_collapsed,df_others,ld_list,df)



### FDR ----



most_pval_info <- sprintf("C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_%s/outcome/OQN_maf_gt_N_pvalue_FDR_MixFinngenPval_6.txt",r2_threshold) %>% fread()

# 紀錄不在 finngen snp
not_in_finngen <- most_pval_info[most_pval_info$alt=="",c("gt_N_hg19","gt_N_hg38")]

# df_most$snp_gt_N 是 hg19
setnames(df_most,
         old = c("snp_gt_N","mostLD_snp_hg38"),
         new = c("gt_N_hg19","mostLD_snp_nearest_hg38"))

# 這32 snp 不在finngen，且附近1000 kb 沒有 finngen 的 snp (前面算R2就連R2=0都會算出)
setdiff(not_in_finngen$gt_N_hg19, df_most$gt_N_hg19 )





fwrite(data.table(k = setdiff(not_in_finngen$gt_N_hg19, df_most$gt_N_hg19 )),
       sprintf("C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_%s/outcome/not_in_finngen_1000gFIN_hg19SNP.txt",r2_threshold),
       row.names = F, col.names = F, sep = "\n") 



# 跟 maf_gt_N_pvalue_FDR_MixFinngenPval_6.txt 紀錄的gt_N, most_pval 結果合併 ####
most_pval_info <- sprintf("C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_%s/outcome/OQN_maf_gt_N_pvalue_FDR_MixFinngenPval_6.txt",r2_threshold) %>% fread()
most_pval_info[,c("CHR","start"):= NULL]
a <- df_most[most_pval_info, on=.(gt_N_hg19)]


# 確認了，只有80 snp 不在finngen，且附近1000 kb 沒有 finngen 的 snp，才會 mostLD_snp_hg19_finngen_pval, gt_N_finngen_pval 兩個都是 NA，其他snp 都只有一個是 NA
nofin_LD <-  setdiff(not_in_finngen$gt_N_hg19, df_most$gt_N_hg19 )
test <- intersect(which(is.na(a$mostLD_snp_hg19_finngen_pval)), which(is.na(a$gt_N_finngen_pval))) %>%
  a[.,gt_N_hg19]
all(test==nofin_LD)





# 把 gt_N_finngen_pval , mostLD_snp_hg19_finngen_pval 合併算 FDR ####
most_LD <- df_most[,.SD[1],by=mostLD_snp_nearest_hg19] %>% 
  select(gt_N_hg19,mostLD_snp_hg19_finngen_pval,mostLD_snp_nearest_hg19)

setnames(most_LD,old = "mostLD_snp_hg19_finngen_pval", new = 'mostLD_or_finngen_pval')

most_pval_info_unique <- most_pval_info[,.SD[1],by=gt_N_hg19] %>%
  select(gt_N_hg19,gt_N_finngen_pval)
setnames(most_pval_info_unique,old = 'gt_N_finngen_pval', new = 'mostLD_or_finngen_pval')
most_pval_info_unique[,mostLD_snp_nearest_hg19 := NA_character_ ]


# 扣掉 4095 snp 不在finngen, pval =NA 的
most_pval_info_unique <- most_pval_info_unique[which(!is.na(mostLD_or_finngen_pval)),]


# 合併
record_in_and_LD_pval <- rbind(most_pval_info_unique,most_LD)


# all col NA row number
colSums(is.na(record_in_and_LD_pval))
# all col 空值 row number
colSums((record_in_and_LD_pval==""))


record_in_and_LD_pval[,mostLD_or_finngen_FDR :=p.adjust(mostLD_or_finngen_pval,
                                                        method = "BH") %>% 
                        format(digits = 10,scientific = T) %>% 
                        as.numeric()]

record_in_and_LD_pval[,mostLD_or_finngen_qvalue := qvalue(mostLD_or_finngen_pval)$qvalues]

record_in_and_LD_pval[,mostLD_or_finngen_Bonfi := 
                        ifelse(mostLD_or_finngen_pval<(0.05/nrow(record_in_and_LD_pval)),
                               1,0)]


# record 紀錄 excel 數量，因為前面把 LDsnp, gt_N snp 取唯一後相加，所以不會有 snp 重複出現的問題
record_in_and_LD_pval[mostLD_or_finngen_FDR<0.05,] %>%  nrow()
record_in_and_LD_pval[mostLD_or_finngen_qvalue<0.05,] %>%  nrow()
record_in_and_LD_pval[mostLD_or_finngen_Bonfi==1,] %>%  nrow()





# 不同 gt_N snp 對到同個 LD snp，併入 ####

target_cols <- c("mostLD_or_finngen_pval", 
                 "mostLD_or_finngen_FDR", 
                 "mostLD_or_finngen_qvalue", 
                 "mostLD_or_finngen_Bonfi")

# record_in_and_LD_pval 併入 df_most 
df_result <- record_in_and_LD_pval[, c("mostLD_snp_nearest_hg19", target_cols), with = FALSE][
  df_most, on = "mostLD_snp_nearest_hg19"
]

# 選出紀錄在 finngen的 8599 snp  FDR, qvalue 數值，跟 不在finngen 但有 LD snp 的 3366 snp FDR 數值合併
test <- rbind(record_in_and_LD_pval[which(is.na(mostLD_snp_nearest_hg19))],
              df_result %>% 
                select(names(record_in_and_LD_pval)))
# 更新資料
a <- test[a,on=.(gt_N_hg19,mostLD_snp_nearest_hg19)]



a[,start:= NULL]
setnames(a,
         old=c("end","gt_N_finngen_pval", "snp_finngen_FDR",
               "snp_finngen_qvalue","snp_finngen_Bonfi"),
         new=c("pos","gt_N_hg38_finngen_pval", "gt_N_hg38_finngen_FDR",
               "gt_N_hg38_finngen_qvalue","gt_N_hg38_finngen_Bonfi"))

setcolorder(a,
            c("CHR","pos","gt_N_hg38","gt_N_hg19","alt", "ref",
              "gt_N_hg38_finngen_pval", "gt_N_hg38_finngen_FDR",
              "gt_N_hg38_finngen_qvalue","gt_N_hg38_finngen_Bonfi","MOST_snp_hg38",
              "MOST_snp_nearest_hg38","MOST_snp_nearest_pvalue","MOST_snp_hg38_finngen_FDR",
              "MOST_snp_hg38_finngen_Bonfi","LD","mostLD_snp_hg19","mostLD_snp_nearest_hg19",
              "mostLD_snp_nearest_hg38","mostLD_snp_hg19_finngen_pval",
              "mostLD_or_finngen_pval","mostLD_or_finngen_FDR",
              "mostLD_or_finngen_qvalue","mostLD_or_finngen_Bonfi")
)
setkey(a,CHR,pos)



a[,CHR := str_extract(gt_N_hg38, ".*(?=:)")]


# 新增 MOST_snp_nearest_hg19 欄位 ####
most_snp <- a %>% 
  select(MOST_snp_nearest_hg38) %>% 
  unique() %>%
  separate(MOST_snp_nearest_hg38,
           into = c("chr", "pos", "ref", "alt"),
           sep = ":", remove = FALSE) %>% 
  as.data.table()

most_snp[,chr := as.numeric(chr)]
most_snp[,pos := as.numeric(pos)]


setorder(most_snp,chr,pos)
test <- data.table(chr=paste0("chr",most_snp$chr),
                   start=most_snp$pos, 
                   end=most_snp$pos ) 

test_eur <- sprintf("C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_%s/trash/most_snp.txt",r2_threshold) %>% fread()
names(test_eur) <- names(test)

if(!isTRUE(all.equal(test, test_eur))){
  stop("most pval snp in EUR, FIN are different")
}


# 用 EUR liftover 轉換檔案
snp_hg19 <- sprintf("C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_%s/trash/most_SNPhg38TOhg19.bed",r2_threshold) %>% 
  fread(header=F)
snp_hg19[,chr := str_remove(V1, "chr")]

unmapped_path <- sprintf(
  "C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_%s/trash/hg38TOhg19_unmapped.bed",
  r2_threshold
)

if (!file.exists(unmapped_path) || file.info(unmapped_path)$size == 0) {
  a_unmapped <- data.table(V1 = character(), V2 = integer(), V3 = integer())
} else {
  a_unmapped <- fread(
    unmapped_path,
    header = FALSE,
    comment.char = "#"
  )
}

trans_success <- paste0(snp_hg19$chr, ":", snp_hg19$V2)
trans_fail <- paste0(a_unmapped$V1, ":", a_unmapped$V2, "-", a_unmapped$V2) %>%
  str_extract("(\\d+:\\d+)(?=\\-)")

most_snp[, MOST_snp_nearest_hg19 := NA_character_]
index_hg19_suc <- !paste0(most_snp$chr, ":", most_snp$pos) %in% trans_fail
most_snp[index_hg19_suc, MOST_snp_nearest_hg19 := trans_success]

most_snp <- as.data.table(most_snp)  
most_snp[,c("chr","pos","ref","alt" ) := NULL]
a <- most_snp[a, on= .(MOST_snp_nearest_hg38)]


fwrite(a,
       sprintf("C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_%s/outcome/OQN_MixFinngenPval_7_FIN.txt",r2_threshold),
       row.names = F, col.names = T, sep = "\t")


names(most_snp)




### Final FDR ----

k <- sprintf("C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_%s/outcome/OQN_MixFinngenPval_7_FIN.txt",r2_threshold) %>% fread()

# 挑出 LD snp
ld_snp <- k[mostLD_snp_nearest_hg19!="",
            c("mostLD_snp_nearest_hg38","mostLD_snp_hg19_finngen_pval")] %>% 
  unique()

# 挑出 gt_N record in finngen GWAS snp
gt_N <- k[mostLD_snp_nearest_hg19=="",
          c("gt_N_hg38","gt_N_hg38_finngen_pval")] %>% 
  unique()

# 6:31267725 等185 個 snp 重複出現在 finngen, 挑最顯著的
setorder(gt_N,gt_N_hg38_finngen_pval)
gt_N <- gt_N[, .SD[1], by="gt_N_hg38"]

# 刪掉 32 個沒有 LD snp 且沒紀錄在 GWAS 的
gt_N <- gt_N[!is.na(gt_N_hg38_finngen_pval),]


# LD snp 中有 731 snp 出現在 gt_N，刪掉
combine_ldsnp_gtN <- merge(ld_snp,gt_N,
                           by.x = c("mostLD_snp_nearest_hg38","mostLD_snp_hg19_finngen_pval"),
                           by.y = c("gt_N_hg38","gt_N_hg38_finngen_pval"))

ld_snp <- ld_snp[!mostLD_snp_nearest_hg38 %in% combine_ldsnp_gtN$mostLD_snp_nearest_hg38,]



# 合併不重複的 203 LD snp, 9421 原本就在 finngen 的 snp
setnames(ld_snp, old = c("mostLD_snp_nearest_hg38","mostLD_snp_hg19_finngen_pval"),
         new = c("snp","pval"))
setnames(gt_N, old = c("gt_N_hg38","gt_N_hg38_finngen_pval"),
         new = c("snp","pval"))
final <- rbind(gt_N, ld_snp)

final[,mostLD_or_finngen_FDR :=p.adjust(pval,method = "BH") %>% 
        format(digits = 10,scientific = T) %>% 
        as.numeric()]

final[,mostLD_or_finngen_qvalue := qvalue(pval)$qvalues]
final[,mostLD_or_finngen_Bonfi := 
        ifelse(pval<(0.05/nrow(final)),
               1,0)]


cols_to_fill <- c("mostLD_or_finngen_FDR", "mostLD_or_finngen_qvalue", "mostLD_or_finngen_Bonfi")
k[, (cols_to_fill) := NULL] 


# 對 final$snp 對應 k$gt_N_hg38 的 row，final的 cols_to_fill 變數填入 k
k[final, on = .(gt_N_hg38 = snp), 
  (cols_to_fill) := .(i.mostLD_or_finngen_FDR, 
                      i.mostLD_or_finngen_qvalue, 
                      i.mostLD_or_finngen_Bonfi)]


# 第二次填入，對 final$snp 對應 k$mostLD_snp_nearest_hg38 的 row，final的 cols_to_fill 變數填入 k
k[final, on = .(mostLD_snp_nearest_hg38 = snp), 
  (cols_to_fill) := .(i.mostLD_or_finngen_FDR, 
                      i.mostLD_or_finngen_qvalue, 
                      i.mostLD_or_finngen_Bonfi)]




# 加 R2 info
cis_snp <- fread("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/cis_hg18_19_snp_MAF.txt",header=T)
cis_snp[!cis_snp$REF %in% c("A", "T", "C", "G") | !cis_snp$ALT %in% c("A", "T", "C", "G"),
        rsID := NA_character_]

cis_snp <- cis_snp[,c("hg19_snpID","rsID","R2","ER2","impute_type")]
setnames(cis_snp, old = "hg19_snpID", new = "gt_N_hg19")
k <- cis_snp[k, on= .(gt_N_hg19)]

setcolorder(k,c("CHR",	"pos",	"gt_N_hg38",	"gt_N_hg19","rsID","R2","ER2","impute_type",	"alt",	"ref",
                "gt_N_hg38_finngen_pval","gt_N_hg38_finngen_FDR",	"gt_N_hg38_finngen_qvalue",
                "gt_N_hg38_finngen_Bonfi","MOST_snp_hg38",	"MOST_snp_nearest_hg38",	"MOST_snp_nearest_hg19",
                "MOST_snp_nearest_pvalue","MOST_snp_hg38_finngen_FDR",	"MOST_snp_hg38_finngen_Bonfi",
                "LD",	"mostLD_snp_hg19",	"mostLD_snp_nearest_hg19",	"mostLD_snp_nearest_hg38",
                "mostLD_snp_hg19_finngen_pval",	"mostLD_or_finngen_pval",	"mostLD_or_finngen_FDR",
                "mostLD_or_finngen_qvalue",	"mostLD_or_finngen_Bonfi"))
fwrite(k,
       sprintf("C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_%s/outcome/OQN_MixFinngenPval_7_FIN.txt",r2_threshold),
       row.names = F, col.names = T, sep = "\t") 


# 造過 eQTL bon 門檻的 3 手法結果 ----
gt_N <- sprintf(
  "C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_%s/outcome/OQN_N_bon_R2_%s.txt",
  r2_threshold, r2_threshold
) %>% fread()
setnames(gt_N, old="SNP", new="hg18_snpID")
snp_number <- nrow(gt_N)

# 用 maf_hg18_19 對照，抓 hg19 id
maf_hg18_19 <- fread("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/cis_hg18_19_snp_MAF.txt", header = T)[
  ,
  .(hg18_snpID, hg19_snpID)
]

gt_N <- maf_hg18_19[gt_N, on = .(hg18_snpID)]
final <- k[gt_N_hg19 %in% gt_N$hg19_snpID, ] 
fwrite(final,
       sprintf("C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_%s/outcome/OQN_MixFinngenPval_7_FIN_bon.txt",r2_threshold),
       row.names = F, col.names = T, sep = "\t") 



cat(sprintf("pass eQTL bon %d snp, retain %d in OQN_MixFinngenPval_7_FIN_bon.txt\n", snp_number, nrow(final)))

# 重算 FDR, bon, qval ----
# 自己估 pi_hat
# pi0_hat <- min(max(
#     sum(x$`p-value` > 0.7) / ((1-0.7) * length(x$`p-value`)),
#     0), 1)
# x[,qvalue := qvalue(`p-value`,pi0=pi0_hat)$qvalues]


final[!is.na(gt_N_hg38_finngen_pval),
  `:=`(
    gt_N_hg38_finngen_FDR =
      p.adjust(gt_N_hg38_finngen_pval, method = "BH"),
    gt_N_hg38_finngen_Bonfi =
      fifelse(gt_N_hg38_finngen_pval < 0.05 / .N, 1, 0)
  )
]
final[!is.na(mostLD_or_finngen_pval),
  `:=`(
    mostLD_or_finngen_FDR =
      p.adjust(mostLD_or_finngen_pval, method = "BH"),
    mostLD_or_finngen_Bonfi =
      fifelse(mostLD_or_finngen_pval < 0.05 / .N, 1, 0)
  )
  
]
final[!is.na(MOST_snp_nearest_pvalue),
  `:=`(
    MOST_snp_hg38_finngen_FDR =
      p.adjust(MOST_snp_nearest_pvalue, method = "BH"),
    MOST_snp_hg38_finngen_Bonfi =
      fifelse(MOST_snp_nearest_pvalue < 0.05 / .N, 1, 0)
  )
  
]

fwrite(final,
  sprintf("C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_%s/outcome/OQN_MixFinngenPval_7_FIN_bon.txt", r2_threshold),
  row.names = F, col.names = T, sep = "\t"
)







