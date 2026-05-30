# ============================================================
# TITLE: Imputation Result Dealing
# SUBTITLE: Handle strand flips, allele switches, and PLINK post-imputation cleanup.
# SEARCH TAGS: imputation, result, PLINK, strand flip, allele switch, cleanup
# NOTE: Code below is unchanged; only this navigation header was added.
# ============================================================

# OUTLINE: Load required packages ----
# package ----
library(data.table)
library(dplyr)
library(tidyverse)

# OUTLINE: Define command runner ----
# function ----
## plink code ----

# 這樣plink 執行出錯，就會停下來

run_cmd <- function(cmd) {
  
  exit_code <- system(cmd)
  
  # 發生錯誤，system 回傳非 0 的數字
  if (exit_code != 0) {
    stop(paste("cmd 執行失敗"))
  }
}




# flip, switch處理
# OUTLINE: Classify strand flip and allele switch SNPs ----
## distinguish & flip ----

# 1. distinguish different failed type: Allele switch, Strand flip and both
# 2. deal Strand flip
a1 <- fread("D:/oral_cancer/deal_strand_flip/snps-excluded.txt",header = T)

# .*: 代表匹配任意字元直到冒號 : 為止
# ([A-Z]) 代表空格後的第一個大寫字母。我們之後用 \\1 取出
# /([A-Z]) 代表/ 後的另一個大寫字母，捕捉群組 \\2
a1$REF_reference <-  sub(".*: ([A-Z])/([A-Z]).*", "\\1", a1$INFO)
a1$ALT_reference <-  sub(".*: ([A-Z])/([A-Z]).*", "\\2", a1$INFO) 

switch <- (str_detect(a1$INFO,pattern = "Allele switch.") & 
             !str_detect(a1$INFO,pattern = "Strand flip and Allele switch.")) %>% 
  a1[.,]

strand <- (str_detect(a1$INFO,pattern = "Strand flip.") & 
             !str_detect(a1$INFO,pattern = "Strand flip and Allele switch.")) %>% 
  a1[.,]

both <- str_detect(a1$INFO,pattern = "Strand flip and Allele switch.") %>% 
  a1[.,]

# 輸出flip_SNPlist.txt 紀錄Strand flip snp id
fwrite(c(strand$ID,
         both$ID) %>% 
         as.data.table(),
       col.names = F,
       "D:/oral_cancer/deal_strand_flip/flip_SNPlist.txt")


setwd("D:/oral_cancer/deal_strand_flip") 
# 把 flip_SNPlist.txt 檔案中的 snp 進行allele flip，像是從A/G 變成 T/C
system("plink --bfile D:/oral_cancer/imputation_preprocessing/ForImputation --flip flip_SNPlist.txt --recode vcf --out FinalData")


# OUTLINE: Apply allele switch correction ----
## switch ----
setwd("D:/oral_cancer/deal_strand_flip")

# 輸出switch_SNPlist.txt 紀錄switch snp id, ref panel ALT
fwrite(data.table(id = c(both$ID, switch$ID),
                  REF = c(both$REF_reference, switch$REF_reference)),
       col.names = F,
       sep = "\t",
       "D:/oral_cancer/deal_strand_flip/switch_SNPlist.txt")

system("plink --vcf FinalData.vcf --a2-allele switch_SNPlist.txt --recode vcf --out FinalData_switch_done")


# OUTLINE: Split corrected VCF by chromosome ----
## split chr ----
setwd("D:/oral_cancer/deal_strand_flip")

# --keep-allele-order 保持REF, ALT順序。PLINK 在每次子集合(例如 --chr)，預設會把 A1 設成該子集合中的minor allele
for (chr in 1:22) {
  sprintf("plink --vcf FinalData_switch_done.vcf --keep-allele-order --chr %d --recode vcf --out ForImputation_chr%d", chr, chr) %>% 
    system()
  
}

# 釋放記憶體
rm(list=ls())
gc()


# imputation result 整理 ----
# OUTLINE: Remove repeated imputed SNP positions ----
## delete repeat snp ----

# 1. delete repeat snp
# 2. save as chr22_t_select.tped
# record(snp重複個數, 資料被刪個數, 最終沒重複的snp個數)
record <- matrix(0, nrow = 1, ncol = 3)


for (i in 1:22) {
  
  snp_imputed <- sprintf("D:/oral_cancer/imputation_result/chr_%d/chr%d.info", i, i) %>%
    fread() 

  # .*AF= 代表從AF= 之後的字元讀取
  # ([^;]+) 小括弧代表，保存匹配到的內容，稍後可以用 \\1 取出
  # [^;] 代表保存不是分號的字元，後面多個+ 代表1 次或多次。所以合起來就是，AF= 後面直到分號為止的所有字元
  snp_imputed[, AF     := as.numeric(sub(".*AF=([^;]+).*", "\\1", INFO))]
  snp_imputed[, MAF    := as.numeric(sub(".*MAF=([^;]+).*", "\\1", INFO))]
  snp_imputed[, AVG_CS := as.numeric(sub(".*AVG_CS=([^;]+).*", "\\1", INFO))]
  snp_imputed[, R2     := as.numeric(sub(".*R2=([^;]+).*", "\\1", INFO))]
  snp_imputed[, row_num := 1:nrow(snp_imputed)]

  # .vcf.gz 轉成 chr22_t.map, chr22_t.tfam, chr22_t.tped
  sprintf("plink --vcf D:/oral_cancer/imputation_result/chr_%d/chr%d.dose.vcf.gz --recode transpose --out D:/oral_cancer/imputation_result/chr_%d/chr%d_t",
          i, i, i, i) %>% 
    system()

  tped <- sprintf("D:/oral_cancer/imputation_result/chr_%d/chr%d_t.tped",
                  i, i) %>% 
    fread()

  
  # 按照R2 從大到小排列，再刪掉POS 重複的row，如果R2 一樣，保留最上面的row
  imputed_unique <- snp_imputed[order(-R2)][, .SD[1], by = POS] 
  
  
  # 確認 .tped 跟.info snp 順序一致
  if(all(tped$V4==snp_imputed$POS)){
    
    # 確認 imputed_unique 沒有重複的snp
    if(nrow(imputed_unique) == unique(tped$V4) %>% 
       length()){
      
      # 幾個snp有重複
      snp_repeat <- snp_imputed[, .N, by = POS][N > 1] %>%
        nrow()
      
      # 幾比資料被刪掉
      data_delete <- nrow(snp_imputed)-nrow(imputed_unique)
  
      # 更新紀錄
      record <- rbind(record,
                      c(snp_repeat, data_delete,
                        nrow(imputed_unique)))
      
      # tped_unique 是snp沒有重複，且選出R2最大的snp 的tped 檔案
      tped_unique <- imputed_unique$row_num %>%
        sort() %>% 
        tped[.,]
      
      # 存檔
      fwrite(tped_unique,
             sprintf("D:/oral_cancer/imputation_result/chr_%d/chr%d_t_select.tped",i, i),
             row.names = F, col.names = F, sep = "\t")
             
      
    }else{
      cat("imputed_unique snp also repeated！")
      
    }
  }else{
    cat(".tped and .info snp order are different！")
    
  }  
} 









# <!-- "C:/Users/user/Desktop/暫存/0723.pdf" -->
# OUTLINE: Combine chromosome binaries ----
## Combine chromesome ----


setwd("D:/oral_cancer")

# 用 .tped 跟 .tfam 檔案轉成3個binary file
for (i in 1:22) {
  sprintf("plink --tped imputation_result/chr_%d/chr%d_t_select.tped --tfam imputation_result/chr_%d/chr%d_t.tfam --make-bed --out imputation_result_selectSNP/chr%d_final",i,i,i,i,i) %>% 
    system()
  
}

# 避免部分snp使用 rs id，統一改成 chr:pos
for (i in 1:22) {
  test <- sprintf("D:/oral_cancer/imputation_result_selectSNP/chr%d_final.bim",i) %>% 
    fread(, header = F)
  
  names(test) <- c("chr", "id","cM","pos","A1","A2")
  test[,id := paste0(chr, ":", pos)]
  
  fwrite(test, sprintf("D:/oral_cancer/imputation_result_selectSNP/chr%d_final.bim",i),
         sep = "\t", col.names = F)
  
  rm(test)
  gc()
}



# 把chr2~22 的 binary file name 放進txt 檔案裡
mergelist <- sprintf("chr%d_final.%s",
                     rep(2:22, each = 3),
                     c("bed", "bim", "fam")) %>% 
  matrix(ncol = 3, byrow=T)

mergelist <- data.frame(mergelist)
fwrite(mergelist, "D:/oral_cancer/imputation_result_selectSNP/mergelist.txt",
       row.names = F, col.names = F, sep = "\t")


# 合併不同chr
setwd("D:/oral_cancer/imputation_result_selectSNP")
system("plink --bfile chr1_final --merge-list mergelist.txt --make-bed --out chr1-22_imputation")

## record ----

# record matrix 紀錄每個chr 的 snp重複個數, 資料被刪個數, 最終沒重複的snp個數 ----


record <- record[-1,]
record <- record %>% 
  as.data.table()

names(record) <- c("repeat_snp_number","delete_row_number","unique_snp_number")

fwrite(record, "D:/oral_cancer/imputation_result_selectSNP/snp_number_record.txt",
       row.names = T, col.names = T, sep = "\t")



# 釋放記憶體
rm(list=ls())
gc()


# OUTLINE: Clean duplicated SNPs in VCF files ----
# Deal vcf ----
# 處理vcf 檔案出現重複snp 的問題
# 把有 DS的 vcf 檔案，重複snp 挑出R2大的、snp id 改成 chr:pos，存檔後再手動加入前面id說明，跑 parallel


for (i in 1:22) {
  
  
# 更新vcf，刪掉POS 重複的snp ----
  
  sprintf("\"C:/Program Files/7-Zip/7z.exe\" e D:/oral_cancer/imputation_result/chr_%d/chr%d.dose.vcf.gz -oD:/oral_cancer/imputation_result/chr_%d",i,i,i ) %>%
    run_cmd()
  
  
  snp_imputed <- sprintf("D:/oral_cancer/imputation_result/chr_%d/chr%d.dose.vcf", i, i) %>%
    fread() 
  
  snp_imputed[, ID := paste0(`#CHROM`,":", POS)]
  snp_imputed[, R2 := as.numeric(sub(".*R2=([^;]+).*", "\\1", INFO))]
  
  # 按照R2 從大到小排列，再刪掉POS 重複的row，如果R2 一樣，保留最上面的row
  snp_imputed <- snp_imputed[order(-R2)][, .SD[1], by = POS] 
  snp_imputed[,c("R2") := NULL]
  
  # 調換前兩個變數順序，維持 vcf 格式
  all_vars <- names(snp_imputed)
  new_order_v2 <- all_vars
  new_order_v2[1:2] <- all_vars[2:1]
  setcolorder(snp_imputed, new_order_v2)
  
  # 按照pos 大小排序，維持 vcf 格式
  setorder(snp_imputed, POS)
  
  
  # 更新vcf
  fwrite(snp_imputed, sprintf("D:/oral_cancer/imputation_result/chr_%d/chr%d.dose.vcf", i, i),
         row.names = F, col.names = T, sep = "\t")
  
  
  
  
# 更新vcf，加上前面id說明####
  
  # 讀取原始壓縮檔的 header
  # gzfile 可以直接讀寫壓縮檔的內容，readLines 讀取檔案，只讀前100 row
  header_lines <- sprintf("D:/oral_cancer/imputation_result/chr_%d/chr%d.dose.vcf.gz", i,i) %>% 
    gzfile() %>% 
    readLines(,n=100)
  
  # 抓vcf header(開頭 ## 的部分)
  header_lines <- header_lines[grepl("^##", header_lines)]
  
  # 建立檔案 chr22_with_header.vcf，寫入 header_lines
  writeLines(header_lines,
             sprintf("D:/oral_cancer/imputation_result/chr_%d/chr%d_snp_unique.dose.vcf", i, i) )
  
  # 複製chr%d.dose.vcf ，放到 chr%d_snp_unique.dose.vcf 檔案內容的後面
  file.append(sprintf("D:/oral_cancer/imputation_result/chr_%d/chr%d_snp_unique.dose.vcf", i, i),
              sprintf("D:/oral_cancer/imputation_result/chr_%d/chr%d.dose.vcf", i, i) )

  # 刪掉沒有vcf 標頭說明的檔案
  file.remove(sprintf("D:/oral_cancer/imputation_result/chr_%d/chr%d.dose.vcf", i, i))
}




# OUTLINE: Convert cleaned VCF files to PLINK binary ----
## vcf -> binary ----

for(i in 1:22){
  sprintf("plink --vcf D:/oral_cancer/imputation_result/chr_%d/chr%d_snp_unique.dose.vcf --keep-allele-order --make-bed --out D:/oral_cancer/imputation_result/trash/chr%d_final",i,i,i) %>% 
    system()
  
}
  
mergelist <- sprintf("chr%d_final.%s",
                     rep(2:22, each = 3),
                     c("bed", "bim", "fam")) %>% matrix(ncol = 3, byrow=T)

mergelist <- data.frame(mergelist)
fwrite(mergelist, "D:/oral_cancer/imputation_result/trash/mergelist.txt",
       row.names = F, col.names = F, sep = "\t")


# 合併不同chr
setwd("D:/oral_cancer/imputation_result/trash")
system("plink --bfile chr1_final --merge-list mergelist.txt --make-bed --out chr1-22_imputation")





# supplementary ----

## .N ----

# code from <https://www.rdocumentation.org/packages/data.table/versions/1.10.0/topics/special-symbols>
DT = data.table(x=rep(c("b","a","c"),each=3), v=c(1,1,1,2,2,1,1,2,2), y=c(1,3,6), a=1:9, b=9:1)
X = data.table(x=c("c","b"), v=8:7, foo=c(4,2))

# 欄位x 各種元素出現次數
DT[, .N, by=x]                         

# 欄位x 各種元素第i次出現的row
DT[, .SD[i], by=x]                     

# 欄位x 各種元素在其他欄位的 sum，欄位N 代表欄位x 各種元素出現次數
DT[, c(.N, lapply(.SD, sum)), by=x]    # get rows *and* sum columns 'v' and 'y' by group

# 欄位x 各種元素在 DT 出現第i次的 row index
DT[, .I[i], by=x]    
