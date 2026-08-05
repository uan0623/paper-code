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

finngen <- fread("D:/oral_cancer/TWAS/summary_stat/FinnGen/C3_ORALCAVITY_EXALLC")
gt_N <- fread("D:/oral_cancer/snp_repeat_Finngen/outcome/maf_gt_N_pvalue_FDR_info.txt")


# 刪掉REF, ALT 不是 a,t,c,g 的 rsID
gt_N[!gt_N$REF %in% c("A", "T", "C", "G") | !gt_N$ALT %in% c("A", "T", "C", "G"),
  rsID := NA_character_]

fwrite(data.table(a=paste0("chr",gt_N$CHR),
                  b=str_extract(gt_N$SNP, "(?<=\\:)\\d+"), 
                  c=str_extract(gt_N$SNP, "(?<=\\:)\\d+")) ,
       "D:/oral_cancer/snp_repeat_Finngen/process_trash/maf_gt_N_pvalue_FDR_snp.txt",
         row.names = F, col.names = F, sep = "\t")


# 檔案 maf_gt_N_pvalue_FDR_snp.txt 放進leftover轉換成 maf_gt_N_pvalue_FDR_snp_success.bed
a <- fread("D:/oral_cancer/snp_repeat_Finngen/process_trash/maf_gt_N_pvalue_FDR_snp_success.bed")

# "chr22:123-124" -> "22:123"
trans_success <-  paste0(a$V1,":",a$V2,"-",a$V2) %>% 
  str_extract( "(\\d+:\\d+)(?=\\-)")

# 轉失敗的是 c("1:148546405", "1:148567630", "19:59684885")，把除了這些以外的 row 依序放入轉換成功的snp in hg38
gt_N[,snp_hg38 := NA_character_]
index_hg38_suc <- ! gt_N$SNP %in% c("1:148546405", "1:148567630", "19:59684885")
gt_N[index_hg38_suc, snp_hg38 := trans_success]


finngen[, snp:= paste0(`#chrom`,":",pos)]

# finngen 裡有一個snp 不同 alt,pval 的情況，對重複snp 取 pval最小的 
setorder(finngen, pval)
finngen_unique <- finngen[, .SD[1], by = snp]

repeat_snp <- finngen[, .N, by = snp][N > 1]$snp
repeat_snp <- finngen[snp %in% repeat_snp,]
setkey(repeat_snp,snp)



# 把想要的 "`#chrom`","pos","snp" 移前面
setcolorder(repeat_snp,
            c("#chrom","pos","snp", 
              setdiff(names(repeat_snp),
                      c("#chrom","pos","snp"))
              )
            )

fwrite(repeat_snp,
       "D:/oral_cancer/snp_repeat_Finngen/outcome/maf_gt_N_pvalue_FDR_repeatSNP.txt",
         row.names = F, col.names = T, sep = "\t")




# 找出 gt_N$snp_hg38  ==  finngen_unique$snp 的 row，把這些row 的 finngen_unique$pval 寫進 gt_N 新變數 finngen_pval，沒對上的 finngen_pval= NA，使用前確定 finngen_unique$snp 元素唯一出現
gt_N[
  finngen_unique,
  on = c(snp_hg38 = "snp"),
  finngen_pval := i.pval
]

# 3個snp 沒轉成功，snp number: 12024 -> 12021，其中 8599 個有出現在 finngen OC GWAS
uniqueN(gt_N$SNP, na.rm = T)
uniqueN(gt_N$snp_hg38, na.rm = T)
(unique(gt_N$snp_hg38) %in% finngen_unique$snp) %>% which() %>% length()



# 確認 3點 "1:148546405", "1:148567630", "19:59684885" snp_hg38=na not null，再存檔
fwrite(gt_N,
       "D:/oral_cancer/snp_repeat_Finngen/outcome/maf_gt_N_pvalue_FDR_MixFinngenPval.txt",
         row.names = F, col.names = T, sep = "\t")


# 共12021 snp 轉換成 hg38成功，NA不算
uniqueN(gt_N$snp_hg38, na.rm = T)

# 共8599 snp 有對到 finngen_pval，NA不算
(gt_N[!is.na(finngen_pval),snp_hg38]) %>% unique() %>% length()



gt_N <- fread("D:/oral_cancer/snp_repeat_Finngen/outcome/maf_gt_N_pvalue_FDR_MixFinngenPval.txt")


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

# 對每個 finngen_unique 資料，找他落在哪些gt_N$snp_hg38 的區間裡，nomatch = 0L 代表不輸出匹配失敗的
overlap <- foverlaps(finngen_unique, gt_N, nomatch = 0L)
# overlap[, c("i.start", "i.end") := NULL]







overlap[, CHR := as.character(CHR)]

# 取得所有染色體
chr_list <- unique(overlap$CHR)

# 逐染色體處理
for(chr_i in chr_list){
  
  # 只取當前染色體的資料
  finngen_chr <- overlap[CHR == chr_i]
  
  # 如果有 start/end 欄位，用 foverlaps 前先刪掉 NA
  # finngen_chr <- finngen_chr[!is.na(start) & !is.na(end)]
  
  # 存檔
  fwrite(finngen_chr,
         sprintf("D:/oral_cancer/snp_repeat_Finngen/outcome/gt_N_ALLcisSNP/finngen_cisSNP_gt_N_chr%s.txt", chr_i),
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
  overlaps <- sprintf("D:/oral_cancer/snp_repeat_Finngen/outcome/gt_N_ALLcisSNP/finngen_cisSNP_gt_N_chr%d.txt", i) %>%
    fread() 
  overlaps[, start:= str_extract(snp_gt_N , "(?<=\\:)\\d+") %>% as.numeric()]
  overlaps[, end:= start]
  setkey(overlaps, CHR, start)
  
  fwrite(overlaps,
         sprintf("D:/oral_cancer/snp_repeat_Finngen/outcome/gt_N_ALLcisSNP/finngen_cisSNP_gt_N_chr%d.txt", i),
         sep = "\t",
         row.names = FALSE,
         col.names = TRUE)
}





# 每個 snp +-1000KB 選個pval 最小的 snp ####
data_list <- vector("list", 22)

for (chr_i in 1:22) {
  
  
  a <- sprintf("D:/oral_cancer/snp_repeat_Finngen/outcome/gt_N_ALLcisSNP/finngen_cisSNP_gt_N_chr%d.txt", chr_i) %>% 
    fread()
  
  setkey(a, snp_gt_N, pval)
  data_list[[chr_i]] <- a[,  .SD[1], by = snp_gt_N] 
}

df <- rbindlist(data_list, use.names = TRUE)

fwrite(df,
       "D:/oral_cancer/snp_repeat_Finngen/outcome/gt_N_MOSTsnp/finngen_sigSNP_chr1-22.txt",
        row.names = F, col.names = T, sep = "\t")

data_list <- vector("list", 22)
for (chr_i in 1:22) {
a <- sprintf("D:/oral_cancer/snp_repeat_Finngen/outcome/gt_N_ALLcisSNP/finngen_cisSNP_gt_N_chr%d.txt", chr_i) %>%
fread()
setkey(a, snp_gt_N, pval)
data_list[[chr_i]] <- a[,  head(.SD, 10), by = snp_gt_N]
}
df <- rbindlist(data_list, use.names = TRUE)
df[, start:= str_extract(snp_gt_N , "(?<=\\:)\\d+") %>% as.numeric()]
df[, end:= start]
df[, dist := abs(start- i.start)]

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




# 最多出現3 個一樣都是最顯著的pval
table(df_1[,.N, by=snp_gt_N]$N)
df_1[,"MOST_snp_hg38" := paste0(snp,":",ref,":",alt)]
df_1[,"MOST_snp_nearest_hg38" := MOST_snp_hg38]
df_1[,c("i.start" ,"i.end","snp","ref","alt","mlogp","beta","sebeta","dist","keep","grp") := NULL]


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
       "D:/oral_cancer/snp_repeat_Finngen/outcome/maf_gt_N_pvalue_FDR_MixFinngenPval_5.txt",
        row.names = F, col.names = T, sep = "\t")

df_most <- fread("D:/oral_cancer/snp_repeat_Finngen/outcome/maf_gt_N_pvalue_FDR_MixFinngenPval_5.txt")

finngen <- fread(
  "D:/oral_cancer/snp_repeat_Finngen/outcome/C3_ORALCAVITY_EXALLC_2.txt",
  select = c("hg38_snpID","hg19_snpID", "pval","alt","ref")
)
setnames(finngen,old = "hg38_snpID",new = "snp_hg38")
finngen[,"snp_record_finngen" := paste0(snp_hg38,":",ref,":",alt)]

# 將 finngen 的內容併入 df_most
df_most <- finngen[df_most, on = .(snp_hg38)]


setnames(df_most,
         old=c("pval", "snp_record_finngen"),
         new=c("snp_hg38_finngen_pval", "snp_hg38_finngen"))

df_most_1 <- df_most[which(!is.na(snp_hg38_finngen_pval)),]
setkey(df_most_1,CHR, start)



# 算 gt_N 8599 snp FDR ####
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


MOSTpval_info <- MOSTpval_info[df_most, on = .(MOST_snp_nearest_hg38)]



fwrite(data.table(a=paste0("chr",MOSTpval_info$CHR),
                  b=str_extract(MOSTpval_info$snp_hg38, "(?<=\\:)\\d+"), 
                  c=str_extract(MOSTpval_info$snp_hg38, "(?<=\\:)\\d+")) ,
       "D:/oral_cancer/snp_repeat_Finngen/process_trash/MixFinngenPval_5_SNPhg38.txt",
         row.names = F, col.names = F, sep = "\t")


# 把 MOSTpval_info$snp_hg38 轉成 hg19 ，檔案 MixFinngenPval_5_SNPhg38.txt 放進leftover轉換成 MixFinngenPval_5_SNPhg38_success.bed
a <- fread("D:/oral_cancer/snp_repeat_Finngen/process_trash/MixFinngenPval_5_SNPhg38_success.bed")

names(a) <- c("chr","start","end","hg38","unknown")
a[,chr_nochr := str_extract(chr, "\\d+")]
MOSTpval_info_hg19 <- paste0(a$chr_nochr,":",a$start)

MOSTpval_info[,hg19_snpID := MOSTpval_info_hg19]




# 刪變數 ####

# MOST_snp_nearest_hg38 比較近
MOSTpval_info[,c( "MOST_snp_hg38_finngen_pval","snp_hg38_finngen") := NULL]
setnames(MOSTpval_info,
         old=c("snp_hg38", "hg19_snpID", "snp_hg38_finngen_pval"),
         new=c("gt_N_hg38", "gt_N_hg19", "gt_N_finngen_pval"))

a <- c("gt_N_hg38","gt_N_hg19","gt_N_finngen_pval","CHR","start",
       "end","alt","ref","MOST_snp_hg38","MOST_snp_nearest_hg38")

# 把想要的 "`#chrom`","pos","snp" 移前面
setcolorder(MOSTpval_info,
            c(a,
              setdiff(names(MOSTpval_info),a)
              )
            )

fwrite(MOSTpval_info,
       "D:/oral_cancer/snp_repeat_Finngen/outcome/gt_N_MOSTsnp/maf_gt_N_pvalue_FDR_MixFinngenPval_6.txt",
        row.names = F, col.names = T, sep = "\t")



most_pval_info <- fread("D:/oral_cancer/snp_repeat_Finngen/outcome/gt_N_MOSTsnp/maf_gt_N_pvalue_FDR_MixFinngenPval_6.txt")

# 紀錄不在 finngen snp
not_in_finngen <- most_pval_info[most_pval_info$alt=="",c("gt_N_hg19","gt_N_hg38")]


fwrite(data.table(k = not_in_finngen$gt_N_hg19),
       "D:/oral_cancer/snp_repeat_Finngen/outcome/not_in_finngen/not_in_finngen_hg19SNP.txt",
       row.names = F, col.names = F, sep = "\n") 

gt_N_MOSTsnp <- fread("D:/oral_cancer/snp_repeat_Finngen/outcome/gt_N_MOSTsnp/finngen_sigSNP_chr1-22.txt")

# 併進 gt_N，gt_N_MOSTsnp$snp_gt_N 對應到 gt_N$snp_hg38
gt_N <- fread("D:/oral_cancer/snp_repeat_Finngen/outcome/maf_gt_N_pvalue_FDR_MixFinngenPval.txt")

gt_N_MOSTsnp <- gt_N_MOSTsnp %>%
  select(snp_gt_N,snp)

names(gt_N_MOSTsnp) <- c("snp_hg38","1000kb_MOST_snp")
a <- merge(gt_N, gt_N_MOSTsnp, by = "snp_hg38",  all.x = TRUE)
fwrite(a,
       "D:/oral_cancer/snp_repeat_Finngen/outcome/maf_gt_N_pvalue_FDR_MixFinngenPval_2.txt",
        row.names = F, col.names = T, sep = "\t")





# 釋放記憶體
rm(list=ls())
gc()


run_cmd <- function(cmd) {
  
  exit_code <- system(cmd)
  
  # 發生錯誤，system 回傳非 0 的數字
  if (exit_code != 0) {
    stop(paste("cmd 執行失敗"))
  }
}



for (i in 1:22) {
  
  ## vcf to bed
  sprintf("plink1.9 --vcf D:/oral_cancer/PCA_1000G_20130502/1000g_DATA_forPCA/ALL.chr%d.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz --make-bed --out C:/Peter/PCA_1000G_20130502/decompress_binary/chr%d", i, i ) %>% 
    run_cmd()
}


race <- read.table("C:/Peter/PCA_1000G_20130502/integrated_call_samples_v3.20130502.ALL.panel.txt",header = T)
race <- subset(race, super_pop == "EUR" ) 

# 不能用integrated_call_samples_v3.20250704.ALL.ped 裡的family id，同個樣本 decompress_binary/chr%d 裡的family id，跟 integrated_call_samples_v3.20250704.ALL.ped 裡的family id 會不同
fwrite( data.table(family_id = race$sample,
                   indivudual_id = race$sample), "C:/Peter/PCA_1000G_20130502/sample_EUR.txt",
         row.names = F, col.names = F, sep = "\t")

for (i in 1:22) {

  ## select race
  sprintf("plink1.9 --bfile C:/Peter/PCA_1000G_20130502/decompress_binary/chr%d --keep C:/Peter/PCA_1000G_20130502/sample_EUR.txt --make-bed --out C:/Peter/PCA_1000G_20130502/trash/chr_%d", i, i ) %>%
    run_cmd()
}

finngen <- fread("D:/oral_cancer/TWAS/summary_stat/FinnGen/C3_ORALCAVITY_EXALLC")
finngen[,hg38_snpID :=paste0(finngen$`#chrom`, ":",finngen$pos) ]


# 製作 leftover需要的格式，like chr1 123 123  ####
fwrite(data.table(a= paste0("chr",finngen$`#chrom`)[1:10000000],
                  b = finngen$pos[1:10000000],
                  c = finngen$pos[1:10000000]),
       "D:/oral_cancer/snp_repeat_Finngen/process_trash/finngen_for_leftover/finngen_snphg38_1.txt",
       row.names = F, col.names = F, sep = "\t")

fwrite(data.table(a= paste0("chr",finngen$`#chrom`)[10000001:15000000],
                  b = finngen$pos[10000001:15000000],
                  c = finngen$pos[10000001:15000000]),
       "D:/oral_cancer/snp_repeat_Finngen/process_trash/finngen_for_leftover/finngen_snphg38_2.txt",
       row.names = F, col.names = F, sep = "\t")

fwrite(data.table(a= paste0("chr",finngen$`#chrom`)[15000001:nrow(finngen)],
                  b = finngen$pos[15000001:nrow(finngen)],
                  c = finngen$pos[15000001:nrow(finngen)]),
       "D:/oral_cancer/snp_repeat_Finngen/process_trash/finngen_for_leftover/finngen_snphg38_3.txt",
       row.names = F, col.names = F, sep = "\t")



# 確認刪掉轉失敗的點，都會轉成功，code就不留了 ####


# 新增欄位 hg19_snpID ####

# 合併切成3份轉換，失敗的點
failure_snp <- c()
for (i in 1:3) {
  a <-   sprintf("D:/oral_cancer/snp_repeat_Finngen/process_trash/finngen_for_leftover/trans_failure_%d.txt", i) %>%
      read.table()
  
  names(a) <- c("chr","pos","pos_b")
  id <- str_extract(a$chr, "(?<=chr)\\d+") %>%
    paste0(":",a$pos)
  
  failure_snp <- append(failure_snp, id)
}


# 合併切成3份轉換，成功的點(轉換後)
file_paths <- sprintf("D:/oral_cancer/snp_repeat_Finngen/process_trash/finngen_for_leftover/trans_success_%d.bed", 1:3)

suc_snp_hg19 <- lapply(file_paths, fread, header = F) %>% 
  rbindlist()


# 刪掉轉成 chrX, chrY, chr8_gl000196_random, chrUn_gl000247  的 ####
names(suc_snp_hg19) <- c("chr","pos","pos_a","before_trans_id","score")

# 轉失敗的hg19_snpID=NA，成功的賦予值0
finngen[,  hg19_snpID := NA_character_]
finngen[ !hg38_snpID %in% failure_snp,
         hg19_snpID := paste0(suc_snp_hg19$chr,":",suc_snp_hg19$pos) ]




# update ####
# 把 hg19_snpID 欄位從 chr1:123 -> 1:123
# 且令非chr1-22 的為 NA，像是 X,Y,或 chr{染色體編號}_{序列編號}_random, chrUn_{序列編號} 
a <- str_extract(finngen$hg19_snpID, "(?<=^chr)\\d+(?=:)")
not1_22_index <- which(is.na(a))

finngen[,hg19_snpID_nochr:= NA_character_]
finngen[setdiff(1:nrow(finngen), not1_22_index),
        hg19_snpID_nochr:= str_remove(hg19_snpID, "chr")]
finngen[, hg19_snpID:= NULL]
setnames(finngen, old="hg19_snpID_nochr", new="hg19_snpID")

fwrite(finngen,
       "D:/oral_cancer/snp_repeat_Finngen/outcome/C3_ORALCAVITY_EXALLC_2.txt",
       row.names = F, col.names = T, sep = "\t")

finngen <- fread(
    "D:/oral_cancer/snp_repeat_Finngen/outcome/C3_ORALCAVITY_EXALLC_2.txt"
)

# finngen 裡有一個snp 不同 alt,pval 的情況，對重複snp 取 pval最小的 
setorder(finngen, pval)
finngen_unique <- finngen[, .SD[1], by = hg19_snpID ]
fwrite(finngen_unique,
       "D:/oral_cancer/snp_repeat_Finngen/outcome/C3_ORALCAVITY_EXALLC_2_hg19SNPunique.txt",
       row.names = F, col.names = T, sep = "\t")

finngen <- fread("D:/oral_cancer/snp_repeat_Finngen/outcome/C3_ORALCAVITY_EXALLC_1.txt")

# 計算 na, chrX, chrY, chr8_gl000196_random, chrUn_gl000247, ... 各多少 row ####
finngen_hg19 <- str_extract(finngen$hg19_snpID, "^[^:]+")

# 非chr1-22的 row index，包含 NA, chrX, chrY, chr8_gl000196_random, chrUn_gl000247 
CHRnot_1to22 <- which(!finngen_hg19 %in% paste0("chr", c(1:22)))

# 去掉NA
a <- setdiff(CHRnot_1to22,
                        which(is.na(finngen$hg19_snpID)))

str_extract(finngen[a,hg19_snpID], "^[^:]+") %>%
  table()



hg19_chr1_22 <- finngen[-CHRnot_1to22,hg19_snpID] %>%
  unique()


# 把 finngen$hg19_snpID chr1-22 的 id 存下來
df <- data.table( id = hg19_chr1_22,
                 pos = str_extract(hg19_chr1_22 , "(?<=\\:)\\d+") %>% as.numeric(),
                 chr = str_extract(hg19_chr1_22, "(?<=^chr)\\d+(?=:)") %>% as.numeric() )

setkey(df,chr, pos)


fwrite(data.table(k = paste0(df$chr, ":", df$pos) ),
       "D:/oral_cancer/snp_repeat_Finngen/outcome/intersect_finngen_1000G/finngen_hg19SNP.txt",
       row.names = F, col.names = F, sep = "\n") 


maf_hg1819 <- fread("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/cis_hg18_19_snp_MAF.txt")

gt_N <- fread("D:/oral_cancer/snp_repeat_Finngen/outcome/maf_gt_N_pvalue_FDR_MixFinngenPval_2.txt")

# 確認所有 gt_N snp(hg18) in maf_hg1819$hg18_snpID
all(gt_N$SNP %in% maf_hg1819$hg18_snpID)

a <- merge(gt_N, maf_hg1819 %>%
        select("hg18_snpID", "hg19_snpID"),
        by.x = "SNP",by.y ="hg18_snpID",   all.x = TRUE)


fwrite(data.table(kk = a$hg19_snpID %>% 
                   unique()),
       "D:/oral_cancer/snp_repeat_Finngen/outcome/intersect_finngen_1000G/eQTL_SNPhg19.txt",
       row.names = F, col.names = T, sep = "\t")

# IN_DIR="/mnt/c/Peter/PCA_1000G_20130502/trash"

# parallel -j 3 '
#     echo "Processing chromosome {}..."
#     # 使用 awk 處理重複位點，並賦予唯一 ID
#     awk '\''{
#         id = $1":"$4;
#         count[id]++;
#         if (count[id] == 1) {
#             $2 = id;
#         } else {
#             $2 = id "_" count[id];
#         }
#         print $0;
#     }'\'' '"$IN_DIR"'/chr_{}.bim > '"$IN_DIR"'/deal_repeatSNP/chr_{}_rename.bim
    
#     cp '"$IN_DIR"'/chr_{}.bed '"$IN_DIR"'/deal_repeatSNP/chr_{}_rename.bed
#     cp '"$IN_DIR"'/chr_{}.fam '"$IN_DIR"'/deal_repeatSNP/chr_{}_rename.fam
# ' ::: {1..22}

# export EQTL_LIST="/mnt/d/oral_cancer/snp_repeat_Finngen/outcome/not_in_finngen/not_in_finngen_hg19SNP.txt"
# export BFILE_DIR="/mnt/c/Peter/PCA_1000G_20130502/trash/deal_repeatSNP"

# for chr in {1..22}; do
#     echo "正在處理第 ${chr} 號染色體..."
    
#     plink1.9 \
#         --bfile "${BFILE_DIR}/chr_${chr}_rename" \
#         --ld-snp-list "${EQTL_LIST}" \
#         --r2 \
#         --ld-window-kb 1000 \
#         --ld-window 999999 \
#         --ld-window-r2 0 \
#         --list-all \
#         --out "/mnt/c/Peter/PCA_1000G_20130502/no_intersect_LD/LDchr${chr}" &
    
#     if [[ $(($chr % 3)) -eq 0 ]]; then
#         wait
#     fi
# done
# wait
# echo "全部處理完成！"

ld_list <- vector("list", 22)
finngen <- fread("D:/oral_cancer/snp_repeat_Finngen/outcome/C3_ORALCAVITY_EXALLC_2.txt")

for (i in 1:22) {
  a <-  sprintf("C:/Peter/PCA_1000G_20130502/no_intersect_LD/LDchr%d.ld", i) %>%
    fread(header = T)
  
  # 把 1:123_2 轉成 1:123
  a[,SNP_B_noRepeat:= str_extract(SNP_B, "^[^_]+")]
  # 挑出SNP_A 是 not in finngen，SNP_B 是 in finngen 的組合 (以確認所有SNP_A 都是 not in finngen)
  a <- a[which(SNP_B_noRepeat %in% finngen$hg19_snpID),]
  a[,SNP_B_noRepeat:= NULL]
  
  # 按snp_a, R2從大排到小，用row index 的方式選sub data，省記憶體
  # 挑前151個LD高的，因為有67 snp 附近有146個 LD=1 的 snp
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

# check 前面選 181合理，最多附近有 148個snp 相同 LD (不一定是LD=1)
table(df_1[,.N, by=SNP_A]$N)



# 跟finngen 合併 ####
finngen <- finngen %>% 
  select("hg38_snpID","hg19_snpID","ref","alt", "pval")


# 把 16:2045121_2 改成 16:2045121，才能在 finngen 找到 info
setnames(df_1, old = "SNP_B", new = "hg19_snpID")
df_1[, SNP_A := str_extract(SNP_A, "^\\d+:\\d+")]
df_1[, hg19_snpID := str_extract(hg19_snpID, "^\\d+:\\d+")]
df_2 <- finngen[df_1, on= .(hg19_snpID)]


df_2[,"mostLD_snp_hg19" := paste0(hg19_snpID,":",ref,":",alt)]
setnames(df_2,
         old = c("hg38_snpID","CHR_A","SNP_A","BP_A","R2","pval"),
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



rm(finngen,a,i,df_others_collapsed,df_others,ld_list,df)

# check don't have any NA，正常不會出現 NA ####
# all col NA row number
colSums(is.na(df_most))
# all col 空值 row number
colSums((df_most==""))



most_pval_info <- fread("D:/oral_cancer/snp_repeat_Finngen/outcome/gt_N_MOSTsnp/maf_gt_N_pvalue_FDR_MixFinngenPval_6.txt")

# 紀錄不在 finngen snp
not_in_finngen <- most_pval_info[most_pval_info$alt=="",c("gt_N_hg19","gt_N_hg38")]

# df_most$snp_gt_N 是 hg19
setnames(df_most,
         old = c("snp_gt_N","mostLD_snp_hg38"),
         new = c("gt_N_hg19","mostLD_snp_nearest_hg38"))


# 這56 snp 不在finngen，且附近1000 kb 沒有 finngen 的 snp (前面算R2就連R2=0都會算出)
setdiff(not_in_finngen$gt_N_hg19, df_most$gt_N_hg19 )

fwrite(data.table(k = setdiff(not_in_finngen$gt_N_hg19, df_most$gt_N_hg19 )),
       "D:/oral_cancer/snp_repeat_Finngen/outcome/not_in_finngen/not_in_finngen_1000gEUR_hg19SNP.txt",
       row.names = F, col.names = F, sep = "\n") 

# 跟 maf_gt_N_pvalue_FDR_MixFinngenPval_6.txt 紀錄的gt_N, most_pval 結果合併 ####
most_pval_info <- fread("D:/oral_cancer/snp_repeat_Finngen/outcome/gt_N_MOSTsnp/maf_gt_N_pvalue_FDR_MixFinngenPval_6.txt")
most_pval_info[,c("CHR","start"):= NULL]
a <- df_most[most_pval_info, on=.(gt_N_hg19)]


# 確認了，只有56 snp 不在finngen，且附近1000 kb 沒有 finngen 的 snp，才會 mostLD_snp_hg19_finngen_pval, gt_N_finngen_pval 兩個都是 NA，其他snp 都只有一個是 NA
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


# 扣掉 3422 snp 不在finngen, pval =NA 的
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




# check mostLD_or_finngen_pval na 數量是 56, mostLD_snp_hg19 na 數量是 8833+56 (8833 record in finngen, 56 not record in 且 maf=0 or no LD snp), gt_N_hg38_finngen_pval na 數量是 3422 (not record in finngen) ####
# all col NA row number
colSums(is.na(a))
# all col 空值 row number
colSums((a==""))

fwrite(a,
       "D:/oral_cancer/snp_repeat_Finngen/outcome/MixFinngenPval_7_EUR.txt",
        row.names = F, col.names = T, sep = "\t")

race <- read.table("C:/Peter/PCA_1000G_20130502/integrated_call_samples_v3.20130502.ALL.panel.txt",header = T)
race <- as.data.table(race)
race <- subset(race, pop == "FIN" ) 

# 不能用integrated_call_samples_v3.20250704.ALL.ped 裡的family id，同個樣本 decompress_binary/chr%d 裡的family id，跟 integrated_call_samples_v3.20250704.ALL.ped 裡的family id 會不同
fwrite( data.table(family_id = race$sample,
                   indivudual_id = race$sample), "C:/Peter/PCA_1000G_20130502/FIN_sample/sample_FIN.txt",
        row.names = F, col.names = F, sep = "\t")

for (i in 1:22) {
    
    ## select race
    sprintf("plink1.9 --bfile C:/Peter/PCA_1000G_20130502/decompress_binary/chr%d --keep C:/Peter/PCA_1000G_20130502/FIN_sample/sample_FIN.txt --make-bed --out C:/Peter/PCA_1000G_20130502/FIN_sample/trash/chr_%d", i, i ) %>%
        run_cmd()
}

IN_DIR="/mnt/c/Peter/PCA_1000G_20130502/FIN_sample/trash"

parallel -j 3 '
    echo "Processing chromosome {}..."
    # 使用 awk 處理重複位點，並賦予唯一 ID
    awk '\''{
        id = $1":"$4;
        count[id]++;
        if (count[id] == 1) {
            $2 = id;
        } else {
            $2 = id "_" count[id];
        }
        print $0;
    }'\'' '"$IN_DIR"'/chr_{}.bim > '"$IN_DIR"'/deal_repeatSNP/chr_{}_rename.bim

cp '"$IN_DIR"'/chr_{}.bed '"$IN_DIR"'/deal_repeatSNP/chr_{}_rename.bed
cp '"$IN_DIR"'/chr_{}.fam '"$IN_DIR"'/deal_repeatSNP/chr_{}_rename.fam
' ::: {1..22}

export EQTL_LIST="/mnt/d/oral_cancer/snp_repeat_Finngen/outcome/not_in_finngen/not_in_finngen_hg19SNP.txt"
export BFILE_DIR="/mnt/c/Peter/PCA_1000G_20130502/FIN_sample/trash/deal_repeatSNP"

for chr in {1..22}; do
    echo "正在處理第 ${chr} 號染色體..."
    
    plink1.9 \
        --bfile "${BFILE_DIR}/chr_${chr}_rename" \
        --ld-snp-list "${EQTL_LIST}" \
        --r2 \
        --ld-window-kb 1000 \
        --ld-window 999999 \
        --ld-window-r2 0 \
        --list-all \
        --out "/mnt/c/Peter/PCA_1000G_20130502/FIN_sample/no_intersect_LD/LDchr${chr}" &
    
    if [[ $(($chr % 3)) -eq 0 ]]; then
        wait
    fi
done
wait
echo "全部處理完成！"

finngen <- fread("D:/oral_cancer/snp_repeat_Finngen/outcome/C3_ORALCAVITY_EXALLC_2.txt")
ld_list <- vector("list", 22)

for (i in 1:22) {
  a <-  sprintf("C:/Peter/PCA_1000G_20130502/FIN_sample/no_intersect_LD/LDchr%d.ld", i) %>%
    fread(header = T)
  
  # 把 1:123_2 轉成 1:123
  a[,SNP_B_noRepeat:= str_extract(SNP_B, "^[^_]+")]
  # 挑出SNP_A 是 not in finngen，SNP_B 是 in finngen 的組合 (以確認所有SNP_A 都是 not in finngen)
  a <- a[which(SNP_B_noRepeat %in% finngen$hg19_snpID),]
  a[,SNP_B_noRepeat:= NULL]
  
  # 按snp_a, R2從大排到小，用row index 的方式選sub data，省記憶體
  # 挑前380 個LD高的，因為有371 snp 附近有146個 LD=1 的 snp
  setorder(a, SNP_A, BP_A, -R2)
  ld_list[[i]] <- a[, head(.SD, 380), by = SNP_A]
}

df <- rbindlist(ld_list, use.names = TRUE)
nrow(df)%%380==0

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
table(df_1[,.N, by=SNP_A]$N)



# 跟finngen 合併 ####
finngen <- fread(
  "D:/oral_cancer/snp_repeat_Finngen/outcome/C3_ORALCAVITY_EXALLC_2.txt",
  select = c("hg38_snpID","hg19_snpID","ref","alt", "pval")
)


# 把 16:2045121_2 改成 16:2045121，才能在 finngen 找到 info
setnames(df_1, old = "SNP_B", new = "hg19_snpID")
df_1[, SNP_A := str_extract(SNP_A, "^\\d+:\\d+")]
df_1[, hg19_snpID := str_extract(hg19_snpID, "^\\d+:\\d+")]
df_2 <- finngen[df_1, on= .(hg19_snpID)]


df_2[,"mostLD_snp_hg19" := paste0(hg19_snpID,":",ref,":",alt)]
setnames(df_2,
         old = c("hg38_snpID","CHR_A","SNP_A","BP_A","R2","pval"),
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



rm(finngen,a,i,df_others_collapsed,df_others,ld_list,df)

# check don't have any NA，正常不會出現 NA ####
# all col NA row number
colSums(is.na(df_most))
# all col 空值 row number
colSums((df_most==""))



most_pval_info <- fread("D:/oral_cancer/snp_repeat_Finngen/outcome/gt_N_MOSTsnp/maf_gt_N_pvalue_FDR_MixFinngenPval_6.txt")

# 紀錄不在 finngen snp
not_in_finngen <- most_pval_info[most_pval_info$alt=="",c("gt_N_hg19","gt_N_hg38")]

# df_most$snp_gt_N 是 hg19
setnames(df_most,
         old = c("snp_gt_N","mostLD_snp_hg38"),
         new = c("gt_N_hg19","mostLD_snp_nearest_hg38"))

# 這80 snp 不在finngen，且附近1000 kb 沒有 finngen 的 snp (前面算R2就連R2=0都會算出)
setdiff(not_in_finngen$gt_N_hg19, df_most$gt_N_hg19 )

fwrite(data.table(k = setdiff(not_in_finngen$gt_N_hg19, df_most$gt_N_hg19 )),
       "D:/oral_cancer/snp_repeat_Finngen/outcome/not_in_finngen/not_in_finngen_1000gFIN_hg19SNP.txt",
       row.names = F, col.names = F, sep = "\n") 

# 跟 maf_gt_N_pvalue_FDR_MixFinngenPval_6.txt 紀錄的gt_N, most_pval 結果合併 ####
most_pval_info <- fread("D:/oral_cancer/snp_repeat_Finngen/outcome/gt_N_MOSTsnp/maf_gt_N_pvalue_FDR_MixFinngenPval_6.txt")
most_pval_info[,c("CHR","start"):= NULL]
a <- df_most[most_pval_info, on=.(gt_N_hg19)]


# 確認了，只有80 snp 不在finngen，且附近1000 kb 沒有 finngen 的 snp，才會 mostLD_snp_hg19_finngen_pval, gt_N_finngen_pval 兩個都是 NA，其他snp 都只有一個是 NA
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


# 扣掉 3422 snp 不在finngen, pval =NA 的
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


# check mostLD_or_finngen_pval na 數量是 80, mostLD_snp_hg19 na 數量是 8833+80 (8833 record in finngen, 80 not record in 且 maf=0 or no LD snp), gt_N_hg38_finngen_pval na 數量是 3422 (not record in finngen) ####
# all col NA row number
colSums(is.na(a))
# all col 空值 row number
colSums((a==""))


fwrite(a,
       "D:/oral_cancer/snp_repeat_Finngen/outcome/MixFinngenPval_7_FIN.txt",
        row.names = F, col.names = T, sep = "\t")





# 加上 chr
a <-fread("D:/oral_cancer/snp_repeat_Finngen/outcome/MixFinngenPval_7_EUR.txt")
all(str_extract(a$gt_N_hg38, ".*(?=:)")==str_extract(a$gt_N_hg19, ".*(?=:)"))
a[,CHR:= str_extract(gt_N_hg38, ".*(?=:)") %>% as.numeric()]

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
       "D:/oral_cancer/snp_repeat_Finngen/outcome/MixFinngenPval_7_EUR.txt",
        row.names = F, col.names = T, sep = "\t")

k <- fread("D:/oral_cancer/snp_repeat_Finngen/outcome/MixFinngenPval_7_EUR.txt")
nrow(k[gt_N_hg38_finngen_FDR<0.05,])
nrow(k[gt_N_hg38_finngen_qvalue<0.05,])
nrow(k[gt_N_hg38_finngen_Bonfi==1,])

kk <- k[,.SD[1],by="MOST_snp_nearest_hg38"]
nrow(kk[MOST_snp_hg38_finngen_FDR<0.05,])
nrow(kk[MOST_snp_hg38_finngen_Bonfi==1,])

nrow(k[mostLD_or_finngen_FDR<0.05,])
nrow(k[mostLD_or_finngen_qvalue<0.05,])
nrow(k[mostLD_or_finngen_Bonfi==1,])


ld_list <- vector("list", 22)

for (i in 1:22) {
  a <-  sprintf("c:/Peter/PCA_1000G_20130502/LD/deal_repeatSNP/LDchr%d.ld", i) %>%
    fread(header = T)
  
  # 按snp_a, R2從大排到小，用row index 的方式選sub data，省記憶體
  # 挑前k個，因為只挑一個，可能是自己跟自己的 LD
  setorder(a, SNP_A, BP_A, -R2)
  a <- a[!which(SNP_A==SNP_B),]
  a <- a[which(R2==1),]
  ld_list[[i]] <- a
}

df <- rbindlist(ld_list, use.names = TRUE)

# 排除跟自己算LD
df <- df[!which(SNP_A==SNP_B),]
setorder(df, CHR_A, BP_A, -R2)

LD1_freq <- df[,.N, by=SNP_A]
hist(LD1_freq$N)


# get finngen p-value ####
finngen <- fread(
  "D:/oral_cancer/snp_repeat_Finngen/outcome/C3_ORALCAVITY_EXALLC_2.txt",
  select = c("hg19_snpID", "pval")
)

# finngen 裡有一個snp 不同 alt,pval 的情況，對重複snp 取 pval最小的 
finngen <- finngen[, .SD[1], by = hg19_snpID] 


# 跟merge(df, finngen,, by.x = "SNP_B",by.y ="hg19_snpID") 一樣，但比較快，適合用大資料
# nomatch = NULL 可以剔除掉沒匹配到的row
setkey(finngen, hg19_snpID,pval)
setkey(df, SNP_B)
df_1 <- finngen[df]   


# 只保留需要欄位，並重命名
df_1 <- df_1[, .(
  gt_N_SNP_hg19 = SNP_A,
  intersect_finngen_1000gEUR_1000kbLDmax_SNP_hg19 = hg19_snpID,
  R2,
  pval
)]

df_1 <- df_1[!which(is.na(pval)),]
df_1 <- df_1[!which((pval)==""),]
LD1_freq_df1 <- df_1[,.N, by=gt_N_SNP_hg19]
hist(LD1_freq_df1$N, breaks = 40,
     main = expression("Freq Plot"),
     xlab = expression("Number of SNP with LD=1 in" %+-% "1000KB"), ylab = "Freq")


# for 1-22 loop ####
# 在 cmd 執行
FOR /L %i IN (1,1,22) DO ( plink1.9 --bfile "C:/Peter/PCA_1000G_20130502/FIN_sample/trash/deal_repeatSNP/chr_%i_rename" --extract "D:/oral_cancer/snp_repeat_Finngen/outcome/not_in_finngen/not_in_finngen_1000gFIN_hg19SNP.txt" --freq --out "D:/oral_cancer/snp_repeat_Finngen/outcome/not_in_finngen/frq/FIN_chr_%i")


FOR /L %i IN (1,1,22) DO ( plink1.9 --bfile "C:/Peter/PCA_1000G_20130502/trash/deal_repeatSNP/chr_%i_rename" --extract "D:/oral_cancer/snp_repeat_Finngen/outcome/not_in_finngen/not_in_finngen_1000gEUR_hg19SNP.txt" --freq --out "D:/oral_cancer/snp_repeat_Finngen/outcome/not_in_finngen/frq/EUR_chr_%i")






# combine ####
file_paths <- list.files(path = "D:/oral_cancer/snp_repeat_Finngen/outcome/not_in_finngen/frq/", 
                         pattern = "^EUR_chr_\\d+\\.frq$", 
                         full.names = TRUE)

all_maf <- lapply(file_paths, fread, header = T) %>% 
    rbindlist()

fwrite(all_maf, "D:/oral_cancer/snp_repeat_Finngen/outcome/not_in_finngen/frq/all_EUR_combined.txt", sep = "\t")


# combine ####
file_paths <- list.files(path = "D:/oral_cancer/snp_repeat_Finngen/outcome/not_in_finngen/frq/", 
                         pattern = "^FIN_chr_\\d+\\.frq$", 
                         full.names = TRUE)

all_maf <- lapply(file_paths, fread, header = T) %>% 
    rbindlist()

fwrite(all_maf, "D:/oral_cancer/snp_repeat_Finngen/outcome/not_in_finngen/frq/all_FIN_combined.txt", sep = "\t")


a <- fread("D:/oral_cancer/snp_repeat_Finngen/outcome/not_in_finngen/frq/all_FIN_combined.txt")
k <- fread("D:/oral_cancer/snp_repeat_Finngen/outcome/not_in_finngen/not_in_finngen_1000gFIN_hg19SNP.txt",header=F)
setdiff(a$SNP,k$V1)
setdiff(k$V1,a$SNP)

gt_dt <- fread("D:/oral_cancer/snp_repeat_Finngen/outcome/not_in_finngen/not_in_finngen_1000gFIN_hg19SNP.txt",header=F)
names(gt_dt) <- c("snp")

gt_dt[,pos:= str_extract(snp, "(?<=\\:)\\d+") %>% as.numeric()]
gt_dt[,chr:= str_extract(snp, ".*(?=:)") %>% as.numeric()]

ld_dt <- fread("D:/oral_cancer/snp_repeat_Finngen/outcome/MixFinngenPval_7_FIN.txt")
ld_dt <- ld_dt %>% select(mostLD_snp_nearest_hg19) 

ld_dt[,pos:= str_extract(mostLD_snp_nearest_hg19, "(?<=:)[^:]+") %>% as.numeric()]
ld_dt[,chr:= str_extract(mostLD_snp_nearest_hg19, "^[^:]+") %>% as.numeric()]
ld_dt <- unique(ld_dt)

names(ld_dt) <- names(gt_dt)
ld_dt[,type:="LD"]
gt_dt[,type:="gt_N"]

# 2. 定義 1MB (1,000,000 bp) 的邊界
gt_dt[, `:=`(start_window = pos - 1e6, end_window = pos + 1e6)]

# 3. 找在同一條染色體，且 LD 的 pos 落在 gt_N 的 window 內的資料
# nomatch = NULL 只保留有匹配到的
# x 是 LD, i 是 gt_N
matches <- ld_dt[gt_dt,
                 on = .(chr, pos >= start_window, pos <= end_window),  
                 nomatch = NULL, 
                 .(x.snp, i.snp)] 

# 取得有匹配到的 gt_N SNP 清單
hit_snps <- unique(matches$i.snp)

gt_dt[, has_LD_within_1MB := snp %in% hit_snps]

# 有 LD 的 74 snp，其中"2:98080989" "6:26722543" 沒紀錄在 1000G ，從code setdiff(k$V1,a$SNP) 得知
 gt_dt[has_LD_within_1MB == TRUE]
# 沒有 LD 的 6 snp
gt_dt[has_LD_within_1MB == FALSE]

fwrite(gt_dt,
       "D:/oral_cancer/snp_repeat_Finngen/outcome/not_in_finngen/1000gEUR_LD_ORnot.txt",
        row.names = F, col.names = T, sep = "\t")







# 跟 maf_gt_N_pvalue_FDR_MixFinngenPval_6.txt 紀錄的gt_N, most_pval 結果合併 ####
most_pval_info <- fread("D:/oral_cancer/snp_repeat_Finngen/outcome/gt_N_MOSTsnp/maf_gt_N_pvalue_FDR_MixFinngenPval_6.txt")
most_pval_info[,c("CHR","start"):= NULL]
a <- df_most[most_pval_info, on=.(gt_N_hg19)]


# 確認了，只有80 snp 不在finngen，且附近1000 kb 沒有 finngen 的 snp，才會 mostLD_snp_hg19_finngen_pval, gt_N_finngen_pval 兩個都是 NA，其他snp 都只有一個是 NA
nofin_LD <-  setdiff(not_in_finngen$gt_N_hg19, df_most$gt_N_hg19 )
test <- intersect(which(is.na(a$mostLD_snp_hg19_finngen_pval)), which(is.na(a$gt_N_finngen_pval))) %>%
  a[.,gt_N_hg19]
all(test==nofin_LD)






# 把gt_N_finngen_pval , mostLD_snp_hg19_finngen_pval 合併算 FDR ####
most_LD <- df_most %>% 
  select(gt_N_hg19,mostLD_snp_hg19_finngen_pval)
setnames(most_LD,old = "mostLD_snp_hg19_finngen_pval", new = 'mostLD_or_finngen_pval')

most_pval_info_unique <- most_pval_info[,.SD[1],by=gt_N_hg19] %>%
  select(gt_N_hg19,gt_N_finngen_pval)
setnames(most_pval_info_unique,old = 'gt_N_finngen_pval', new = 'mostLD_or_finngen_pval')

# 扣掉 3422 snp 不在finngen, pval =NA 的
most_pval_info_unique <- most_pval_info_unique[which(!is.na(mostLD_or_finngen_pval)),]





# 合併
record_in_and_LD_pval <- rbind(most_pval_info_unique,most_LD)

# 確認沒 NA 或 空
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


# 將 record_in_and_LD_pval 的內容併入 a
a <- record_in_and_LD_pval[a, on=.(gt_N_hg19)]





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


# check mostLD_or_finngen_pval na 數量是 80, mostLD_snp_hg19 na 數量是 8833+80 (8833 record in finngen, 80 not record in 且 maf=0 or no LD snp), gt_N_hg38_finngen_pval na 數量是 3422 (not record in finngen) ####
# all col NA row number
colSums(is.na(a))
# all col 空值 row number
colSums((a==""))


fwrite(a,
       "D:/oral_cancer/snp_repeat_Finngen/outcome/MixFinngenPval_7_FIN.txt",
        row.names = F, col.names = T, sep = "\t")

k <- a$gt_N_hg38_finngen_pval
# 幾種不同的值
uniqueN(k, na.rm=T)
# NA 多少row
length(which(is.na(k)))
# 空值 多少row
length(which(k==""))



k <- a
# all col NA row number
colSums(is.na(k))
# all col 空值 row number
colSums((k==""))
