# ============================================================
# TITLE: LD Clumping QN
# SUBTITLE: Combine QN eQTL p-values with FinnGen p-values for LD clumping.
# SEARCH TAGS: LD clumping, QN, eQTL, FinnGen, 1000G FIN, R2
# NOTE: Code below is unchanged; only this navigation header was added.
# ============================================================
# 簡介：
# 用 1000G FIN_sample/40 OC sample 算LD，保留跟 finngen minor= ref, major=alt 或是對調的 snp，用 finngen GWAS pval/eQTL pval 做 LD clumping




# 以下是 code 做了什麼 (LD clumping 需要 genotype data, summary stat)

# 1. 40 OC sample 的 genotype (hg18) 算 LD。把 finngen GWAS & eQTL summary stat 改成 hg18，做 LD clumping，算 Bonferroni

# 2. 1000G FIN 的 genotype (hg19) 算 LD，把 finngen GWAS summary stat 改成 hg19，做 LD clumping，算 Bonferroni

# 3. 看 12024 eQTL snp +- 1000kb 範圍，挑出在 finngen GWAS 的 snp，不考慮LD，因為真正的 causal effect 可能在 eQTL snp 附近。對這些 snp 算 Bonferroni。以及 用1000G FIN 的 genotype (hg19) + GWAS summary stat 做 LD clumping，r^2=03，算 Bonferroni

# 4. 看 12024 eQTL snp +- 500kb 範圍，挑出在 finngen GWAS 的 snp，不考慮LD，因為真正的 causal effect 可能在 eQTL snp 附近。對這些 snp 算 Bonferroni。以及 用1000G FIN 的 genotype (hg19) + GWAS summary stat 做 LD clumping，r^2=03，算 Bonferroni

# (4. 檔案步驟)
# 找附近 500kb snp 產生 500kb_finngenPval_bon_noClumping.txt (.1)，
# 把 finngen GWAS allele 跟 1000G FIN genotype allele 留下 minor= ref, major=alt 或是對調的 500kb_finngenPval_bon.txt (.2)，
# 從各個染色體的 binary file 取出上一步這些 snp 的 genotype data，產生 gt_N_500kb_in_finngen_hg19SNP_chr1-22 (.3)，
# 做 LD clumping 產生 ld_clumping_500kb_finngenPval_1000GfinLD.clumped，計算 FDR, bon, qvalue 產生 ld_clumping_500kb_finngenPval_1000GfinLD.txt  (.4)


# 6. 補上有在finngen 顯著的snp 5:1333830 (hg18) eQTL info



# OUTLINE: Load required packages ----
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
library(matrixStats)  # rowMins
library(qvalue)


# OUTLINE: Method 1 OC genotype QN clumping ----
## 1.  ----

### Get genotype data ----
# - genearate gt_N eQTL 出現在 finngen GWAS 的 8599 snp binary file

gt_N <- fread("C:/Peter/QN_before_eQTL/outcome/QN_MixFinngenPval_7_EUR.txt")

# 找出gt_N 紀錄在 finngne 的 8599 SNP (8833 rsID) in hg19
gt_N[ref!="",gt_N_hg19] %>% 
  length()

hg18_19_tab <- fread("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/cis_hg18_19_snp.txt")
gt_N_in_finngen_hg18SNP <- hg18_19_tab[ hg19_snpID %in% gt_N[ref!="",gt_N_hg19] , hg18_snpID]

# 轉成 hg18 存檔
fwrite(data.table(k = gt_N_in_finngen_hg18SNP),
       "C:/Peter/LD_clumping/QN_before_eQTL/trash/gt_N_in_finngen_hg18SNP.txt",
       row.names = F, col.names = F, sep = "\n") 

# 從原始 binary file 取出這些snp
system("plink --bfile D:/oral_cancer/expression/trash/chr1-22_imputation_hg18 --extract C:/Peter/LD_clumping/QN_before_eQTL/trash/gt_N_in_finngen_hg18SNP.txt --make-bed --out C:/Peter/LD_clumping/QN_before_eQTL/trash/gt_N_in_finngen_hg18SNP")



# 找出gt_N 紀錄再 finngne 的 8599 SNP (8833 rsID) in hg38
gt_N_in_finngen_SNP <- gt_N[ref!="",c("gt_N_hg38","gt_N_hg19")] 
names(hg18_19_tab) <- c("hg18_snpID", "gt_N_hg19")

# 把 hg18_19_tab 併入 gt_N_in_finngen_SNP
gt_N_in_finngen_SNP <- hg18_19_tab[gt_N_in_finngen_SNP, on= .(gt_N_hg19)]






### Get finngen summary stat ----

# 使用 finngen pval 當成 summary stat
# - 把 GWAS 改成 hg18
# - 挑出有出現在 eQTL 的8599 SNP (8833 rsID)，只保留 minor= ref, major=alt 或是對調的，剩下 8527 rsID

finngen <- fread("D:/oral_cancer/TWAS/summary_stat/FinnGen/C3_ORALCAVITY_EXALLC")
finngen[,chr_pos := paste0(`#chrom`,":",pos)]

# 保留有在 gt_N 的 GWAS
finngen <- finngen[chr_pos %in% gt_N_in_finngen_SNP$gt_N_hg38,]
finngen[,c("rsids","nearest_genes") := NULL]
setnames(finngen, c("#chrom","chr_pos","pval"),c("chr","gt_N_hg38","P"))

# 合併獲取 hg18 ID
gt_N_in_finngen_SNP <- gt_N_in_finngen_SNP[finngen, on= .(gt_N_hg38)]
# 去掉 hg19, hg38 等沒用資訊
gt_N_in_finngen_SNP[,c("gt_N_hg19","gt_N_hg38","af_alt","af_alt_cases","af_alt_controls","chr","pos"):= NULL]


# 保留 minor= ref, major=alt 或是對調的
bim <- fread("C:/Peter/LD_clumping/QN_before_eQTL/trash/gt_N_in_finngen_hg18SNP.bim")
names(bim) <- c("chr","SNP","genetic_dist","pos","minor","major")
setnames(gt_N_in_finngen_SNP,"hg18_snpID","SNP")

gt_N_in_finngen_SNP <- merge(bim,gt_N_in_finngen_SNP,by="SNP")

final_data <- gt_N_in_finngen_SNP[
  (gt_N_in_finngen_SNP$minor == gt_N_in_finngen_SNP$ref & gt_N_in_finngen_SNP$major == gt_N_in_finngen_SNP$alt) | 
  (gt_N_in_finngen_SNP$minor == gt_N_in_finngen_SNP$alt & gt_N_in_finngen_SNP$major == gt_N_in_finngen_SNP$ref), 
]

final_data <- unique(final_data)

# 檢查重複的 snp，同時出現 minor= ref, major=alt 跟對調的，beta 是否只差負號
test <- final_data[, .N, by=SNP][N>1]$SNP

# 對重複snp 取 pval最小的
setorder(final_data, P)
final_data <- final_data[, .SD[1], by = SNP]


# 保留這些 snp GWAS in hg18 存檔
fwrite(final_data,
       "C:/Peter/LD_clumping/QN_before_eQTL/trash/gt_N_in_finngen_summarySTAT_hg18.txt",
       row.names = F, col.names = T, sep = "\t") 

system("plink --bfile C:/Peter/LD_clumping/QN_before_eQTL/trash/gt_N_in_finngen_hg18SNP --clump C:/Peter/LD_clumping/QN_before_eQTL/trash/gt_N_in_finngen_summarySTAT_hg18.txt --clump-p1 0.99999 --clump-r2 0.3 --clump-kb 1000 --out C:/Peter/LD_clumping/QN_before_eQTL/outcome/ld_clumping_finngenPval_40OC_LD")


a <- fread("C:/Peter/LD_clumping/QN_before_eQTL/outcome/ld_clumping_finngenPval_40OC_LD.clumped")
a[, Bonfi := ifelse(P<(0.05/nrow(a)), 1,0)]
a[,qvalue := qvalue(P)$qvalues]
a[,FDR :=p.adjust(P, method = "BH") %>% 
      format(digits = 10,scientific = T) %>% 
      as.numeric()]

fwrite(a,
       "C:/Peter/LD_clumping/QN_before_eQTL/outcome/ld_clumping_finngenPval_40OC_LD.txt",
       row.names = F, col.names = T, sep = "\t") 



### Get eQTL summary stat ----

# 用 eQTL pval 當成 summary stat
gt_N_in_finngen_hg18SNP <- fread("C:/Peter/LD_clumping/QN_before_eQTL/trash/gt_N_in_finngen_hg18SNP.txt",header=F)
gt_N <- fread("C:/Peter/QN_before_eQTL/outcome/QN_maf_gt_N_pvalue_FDR_info.txt")
gt_N <- gt_N[SNP %in% gt_N_in_finngen_hg18SNP$V1,]

# 有些 snp 對到不同的 expression，對重複snp 取 pval最小的
setorder(gt_N, `p-value`)
gt_N <- gt_N[, .SD[1], by = SNP] %>% 
  select(c("CHR","Probe","Gene","SNP","p-value"))
setnames(gt_N,"p-value", "P")

# 保留這些 snp GWAS  存檔
fwrite(gt_N,
       "C:/Peter/LD_clumping/QN_before_eQTL/trash/gt_N_in_finngen_summarySTAT_eQTL_hg18.txt",
       row.names = F, col.names = T, sep = "\t") 


system("plink --bfile C:/Peter/LD_clumping/QN_before_eQTL/trash/gt_N_in_finngen_hg18SNP --clump C:/Peter/LD_clumping/QN_before_eQTL/trash/gt_N_in_finngen_summarySTAT_eQTL_hg18.txt --clump-p1 0.99999 --clump-r2 0.3 --clump-kb 1000 --out C:/Peter/LD_clumping/QN_before_eQTL/outcome/ld_clumping_eQTLPval_40OC_LD")

a <- fread("C:/Peter/LD_clumping/QN_before_eQTL/outcome/ld_clumping_eQTLPval_40OC_LD.clumped")
a[, Bonfi := ifelse(P<(0.05/nrow(a)), 1,0)]
a[,qvalue := qvalue(P)$qvalues]
a[,FDR :=p.adjust(P, method = "BH") %>% 
      format(digits = 10,scientific = T) %>% 
      as.numeric()]

fwrite(a,
       "C:/Peter/LD_clumping/QN_before_eQTL/outcome/ld_clumping_eQTLPval_40OC_LD.txt",
       row.names = F, col.names = T, sep = "\t") 

# OUTLINE: Method 2 1000G FIN QN clumping ----
## 2. ----

# 1000G use hg19 snp，官網有寫
# <https://apexbtic.icgeb.res.in/dbtapex/resources/genome/category/grch37/index.html>

# data from 1000G
# <https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/>

### Get genotype data ----

# C:/Peter/PCA_1000G_20130502/FIN_sample/trash/deal_repeatSNP/chr1-22_rename_FIN 檔案是 260117_repeatSNP_FinnGen r.notebook 做的檔案。做了

# 1. 1000 G data 挑出 FIN sample 
# 2. 因為1000g data 同一個位點出現多次，則自動加上 _2, _3 等後綴



gt_N <- fread("C:/Peter/QN_before_eQTL/outcome/QN_MixFinngenPval_7_EUR.txt")

# 找出gt_N 紀錄在 finngne 的 8599 SNP (8833 rsID) in hg19
hg18_19_tab <- fread("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/cis_hg18_19_snp.txt")
gt_N_in_finngen_hg19SNP <- hg18_19_tab[ hg19_snpID %in% gt_N[ref!="",gt_N_hg19] , hg19_snpID]

# 轉成 hg19 存檔
fwrite(data.table(k = gt_N_in_finngen_hg19SNP),
       "D:/oral_cancer/snp_repeat_Finngen/process_trash/QN_gt_N_in_finngen_hg19SNP.txt",
       row.names = F, col.names = F, sep = "\n") 

# 從原始 binary file 取出這些snp
shell('for /L %i in (1,1,22) do plink --bfile "C:/Peter/PCA_1000G_20130502/FIN_sample/trash/deal_repeatSNP/chr_%i_rename" --extract "D:/oral_cancer/snp_repeat_Finngen/process_trash/QN_gt_N_in_finngen_hg19SNP.txt" --make-bed --out "C:/Peter/LD_clumping/QN_before_eQTL/trash/gt_N_in_finngen_hg19SNP_chr%i"')



### merge ----

# 把chr2~22 的 binary file name 放進txt 檔案裡
mergelist <- sprintf("C:/Peter/LD_clumping/QN_before_eQTL/trash/gt_N_in_finngen_hg19SNP_chr%d.%s",
                     rep(2:22, each = 3),
                     c("bed", "bim", "fam")) %>% 
  matrix(ncol = 3, byrow=T) 

mergelist <- data.frame(mergelist)
fwrite(mergelist, "C:/Peter/LD_clumping/QN_before_eQTL/trash/mergelist.txt",
       row.names = F, col.names = F, sep = "\t")


# 合併不同chr
system("plink --bfile C:/Peter/LD_clumping/QN_before_eQTL/trash/gt_N_in_finngen_hg19SNP_chr1 --merge-list C:/Peter/LD_clumping/QN_before_eQTL/trash/mergelist.txt --make-bed --out C:/Peter/LD_clumping/QN_before_eQTL/trash/gt_N_in_finngen_hg19SNP_chr1-22")



### Get finngen summary stat ----

# 使用 finngen pval 當成 summary stat

# - 挑出有出現在 eQTL 的8599 SNP (8833 rsID)，只保留 minor= ref, major=alt 或是對調的，剩下 8527 rsID

# 保留有在 gt_N 的 8599 SNP (8833 rsID) GWAS stat
finngen <- fread("D:/oral_cancer/snp_repeat_Finngen/outcome/C3_ORALCAVITY_EXALLC_2.txt")[
  hg19_snpID %in% gt_N_in_finngen_SNP$gt_N_hg19
][ ,.(chr = `#chrom`, hg19_snpID, P = pval, ref, alt, beta,mlogp,sebeta)]


gt_N <- fread("C:/Peter/QN_before_eQTL/outcome/QN_MixFinngenPval_7_FIN.txt")
# 找出gt_N 紀錄再 finngne 的 8599 SNP (8833 rsID) in hg38
gt_N_in_finngen_SNP <- gt_N[ref!="",c("gt_N_hg38","gt_N_hg19")]


# 合併與 bim 對齊
bim <- fread("C:/Peter/LD_clumping/QN_before_eQTL/trash/gt_N_in_finngen_hg19SNP_chr1-22.bim",
             col.names = c("chr", "hg19_snpID", "genetic_dist", "pos", "minor", "major"))

# bim 併入 final_data
final_data <- bim[finngen, on= .(chr,hg19_snpID)]


# 過濾 Allele 、針對每個 SNP 取 P-value 最小的那一筆
final_data <- final_data[
  (minor == ref & major == alt) | (minor == alt & major == ref)
][ order(P), .SD[1], by = hg19_snpID ]

setnames(final_data,"hg19_snpID","SNP")


# 保留這些 snp GWAS in hg19 存檔
fwrite(final_data, "C:/Peter/LD_clumping/QN_before_eQTL/trash/gt_N_in_finngen_summarySTAT_hg19.txt", 
       row.names = F, col.names = T, sep = "\t") 


### clumping ----

system("plink --bfile C:/Peter/LD_clumping/QN_before_eQTL/trash/gt_N_in_finngen_hg19SNP_chr1-22 --clump C:/Peter/LD_clumping/QN_before_eQTL/trash/gt_N_in_finngen_summarySTAT_hg19.txt --clump-p1 0.99999 --clump-r2 0.3 --clump-kb 1000 --out C:/Peter/LD_clumping/QN_before_eQTL/outcome/ld_clumping_finngenPval_1000GfinLD")


a <- fread("C:/Peter/LD_clumping/QN_before_eQTL/outcome/ld_clumping_finngenPval_1000GfinLD.clumped")
a[, Bonfi := ifelse(P<(0.05/nrow(a)), 1,0)]
a[,qvalue := qvalue(P)$qvalues]
a[,FDR :=p.adjust(P, method = "BH") %>% 
      format(digits = 10,scientific = T) %>% 
      as.numeric()]

fwrite(a,
       "C:/Peter/LD_clumping/QN_before_eQTL/outcome/ld_clumping_finngenPval_1000GfinLD.txt",
       row.names = F, col.names = T, sep = "\t") 



### record ----

a <- fread("C:/Peter/LD_clumping/QN_before_eQTL/outcome/ld_clumping_finngenPval_40OC_LD.clumped")
b <- fread("C:/Peter/LD_clumping/QN_before_eQTL/trash/gt_N_in_finngen_summarySTAT_hg18.txt")
gt_N <- fread("C:/Peter/QN_before_eQTL/outcome/QN_MixFinngenPval_7_EUR.txt")
c <- fread("C:/Peter/LD_clumping/QN_before_eQTL/outcome/ld_clumping_finngenPval_40OC_LD.txt")
# 法1,2,3 snp beginning 數量
uniqueN(gt_N[ref!="",gt_N_hg19])
# 法1 allele篩選後 snp 數量
uniqueN(b$SNP)
# 法1 clumping 後 snp 數量
uniqueN(a$SNP)
# 法1 過 bon 門檻的 snp 數量
length(c[Bonfi==1,SNP])


a <- fread("C:/Peter/LD_clumping/QN_before_eQTL/outcome/ld_clumping_eQTLPval_40OC_LD.clumped")
c <- fread("C:/Peter/LD_clumping/QN_before_eQTL/outcome/ld_clumping_eQTLPval_40OC_LD.txt")
# 法2 allele篩選後 snp 數量
uniqueN(gt_N[ref!="",gt_N_hg19])
# 法2 clumping 後 snp 數量
uniqueN(a$SNP)
# 法2 過 bon 門檻的 snp 數量
length(c[Bonfi==1,SNP])



a <- fread("C:/Peter/LD_clumping/QN_before_eQTL/outcome/ld_clumping_finngenPval_1000GfinLD.clumped")
b <- fread("C:/Peter/LD_clumping/QN_before_eQTL/trash/gt_N_in_finngen_summarySTAT_hg19.txt")
c <- fread("C:/Peter/LD_clumping/QN_before_eQTL/outcome/ld_clumping_finngenPval_1000GfinLD.txt")
# 法3 allele篩選後 snp 數量
uniqueN(b$SNP)
# 法3 clumping 後 snp 數量
uniqueN(a$SNP)
# 法3 過 bon 門檻的 snp 數量
length(c[Bonfi==1,SNP])
# 法3 過 FDR 門檻的 snp 數量
length(c[FDR<0.05,SNP])



# OUTLINE: Method 3 1000kb nearby SNP analysis ----
## 3.  ----


# 找附近 500kb snp 產生 1000kb_finngenPval_bon_noClumping.txt (.1)，
# 把 finngen GWAS allele 跟 1000G FIN genotype allele 留下 minor= ref, major=alt 或是對調的 1000kb_finngenPval_bon.txt (.2)，
# 從各個染色體的 binary file 取出上一步這些 snp 的 genotype data，產生 gt_N_1000kb_in_finngen_hg19SNP_chr1-22 (.3)，
# 做 LD clumping 產生 ld_clumping_1000kb_finngenPval_1000GfinLD.clumped，計算 FDR, bon, qvalue 產生 ld_clumping_1000kb_finngenPval_1000GfinLD.txt  (.4)



finngen_hg19 <- fread("D:/oral_cancer/snp_repeat_Finngen/outcome/C3_ORALCAVITY_EXALLC_2.txt")
gt_N_hg19 <- fread("C:/Peter/QN_before_eQTL/outcome/QN_intersect_finngen_1000G/QN_eQTL_SNPhg19.txt")

finngen_hg19 <- finngen_hg19 %>% 
  select(`#chrom`,hg19_snpID,pos,pval,ref,alt)

setnames(finngen_hg19, old = "#chrom", new = "CHR")
finngen_hg19[,start := pos ]
finngen_hg19[,end := pos ]

# 調整範圍正負1000KB (10^6 bp)
names(gt_N_hg19) <- c("hg19_snpID")
gt_N_hg19[, pos := str_extract(hg19_snpID, "(?<=\\:)\\d+") %>% as.numeric()]
gt_N_hg19[, start := ifelse(pos-1e6>0, pos-1e6, 0)]
gt_N_hg19[, end := pos+1e6]
gt_N_hg19[, CHR := str_extract(hg19_snpID, ".*(?=:)") %>% as.numeric()]


# 要先排序
setkey(finngen_hg19, CHR, start, end)
setkey(gt_N_hg19, CHR, start, end)
names(gt_N_hg19) <- c("gt_N_hg19","pos","start","end","CHR")

finngen_hg19[,c("pos"):= NULL]
# 語法解釋：找 gt_N 中，CHR 一樣，且 pos 在 finngen 的 lower 與 upper 之間的位點
overlap <- finngen_hg19[gt_N_hg19, 
                     on = .(CHR, start >= start, start <= end), 
                     nomatch = 0L,
                     allow.cartesian = TRUE]



gt_N_surround_snp <- unique(overlap,by="hg19_snpID") 
setnames(gt_N_surround_snp,old="hg19_snpID", new="fin_GWAS_hg19")
gt_N_surround_snp[, c("pos","start","end","start.1","gt_N_hg19") := NULL ]
gt_N_surround_snp <- na.omit(gt_N_surround_snp)
gt_N_surround_snp <- gt_N_surround_snp[rowSums(gt_N_surround_snp == "") == 0]


fwrite(gt_N_surround_snp,
       "C:/Peter/LD_clumping/QN_before_eQTL/outcome/1000kb_finngenPval_bon_noClumping.txt",
       row.names = F, col.names = T, sep = "\t") 




### 留下跟 binary allele 符合的 ----

# 檢查重複的 snp，同時出現 minor= ref, major=alt 跟對調的，beta 是否只差負號


# 迴圈抓
data_list <- vector("list", 22)

for (i in 1:22) {
  bim <- sprintf("C:/Peter/PCA_1000G_20130502/FIN_sample/trash/deal_repeatSNP/chr_%d_rename.bim",i) %>% 
    fread(.,col.names = c("CHR", "fin_GWAS_hg19", "genetic_dist", "pos", "minor", "major"))

  bim <- bim[gt_N_surround_snp, on= .(CHR,fin_GWAS_hg19), nomatch = NULL]
  
  data_list[[i]] <- bim[
    (minor == ref & major == alt) | (minor == alt & major == ref)
    ][ order(pval), .SD[1], by = fin_GWAS_hg19 ]
  
}
final_data <- rbindlist(data_list, use.names = TRUE)

fwrite(final_data,
       "C:/Peter/LD_clumping/QN_before_eQTL/outcome/1000kb_finngenPval_bon.txt",
       row.names = F, col.names = T, sep = "\t") 

print(nrow(final_data))



# 1000G data (hg19)+ finngen GWAS (hg38)，轉成 hg19

### Get genotype   ----
sur_snp <- fread("C:/Peter/LD_clumping/QN_before_eQTL/outcome/1000kb_finngenPval_bon.txt")

# 取出 附近 1000kb snp hg19 存檔
fwrite(data.table(k = sur_snp$fin_GWAS_hg19),
       "C:/Peter/LD_clumping/QN_before_eQTL/trash/gt_N_1000kb_in_finngen_hg19SNP.txt",
       row.names = F, col.names = F, sep = "\n") 

# 從原始 binary file 取出這些snp
shell('for /L %i in (1,1,22) do plink --bfile "C:/Peter/PCA_1000G_20130502/FIN_sample/trash/deal_repeatSNP/chr_%i_rename" --extract "C:/Peter/LD_clumping/QN_before_eQTL/trash/gt_N_1000kb_in_finngen_hg19SNP.txt" --make-bed --out "C:/Peter/LD_clumping/QN_before_eQTL/trash/gt_N_1000kb_in_finngen_hg19SNP_chr%i"')

# 把chr2~22 的 binary file name 放進txt 檔案裡
mergelist <- sprintf("C:/Peter/LD_clumping/QN_before_eQTL/trash/gt_N_1000kb_in_finngen_hg19SNP_chr%d.%s",
                     rep(2:22, each = 3),
                     c("bed", "bim", "fam")) %>% 
    matrix(ncol = 3, byrow=T) 

mergelist <- data.frame(mergelist)
fwrite(mergelist, "C:/Peter/LD_clumping/QN_before_eQTL/trash/1000kb_mergelist.txt",
       row.names = F, col.names = F, sep = "\t")


# 合併不同chr
system("plink --bfile C:/Peter/LD_clumping/QN_before_eQTL/trash/gt_N_1000kb_in_finngen_hg19SNP_chr1 --merge-list C:/Peter/LD_clumping/QN_before_eQTL/trash/1000kb_mergelist.txt --make-bed --out C:/Peter/LD_clumping/QN_before_eQTL/trash/gt_N_1000kb_in_finngen_hg19SNP_chr1-22")

### clumping ----

# --clump-snp-field, --clump-field 分別指定 snp, pval 在資料的哪個欄位
system("plink --bfile C:/Peter/LD_clumping/QN_before_eQTL/trash/gt_N_1000kb_in_finngen_hg19SNP_chr1-22 --clump C:/Peter/LD_clumping/QN_before_eQTL/outcome/1000kb_finngenPval_bon.txt --clump-p1 0.99999 --clump-r2 0.3 --clump-kb 1000 --clump-snp-field fin_GWAS_hg19 --clump-field pval --out C:/Peter/LD_clumping/QN_before_eQTL/outcome/LD_0.3/ld_clumping_1000kb_finngenPval_1000GfinLD")

af_clumping <- fread("C:/Peter/LD_clumping/QN_before_eQTL/outcome/LD_0.3/ld_clumping_1000kb_finngenPval_1000GfinLD.clumped")

af_clumping[,Bonfi := ifelse(P<(0.05/nrow(af_clumping)),
                             1,0)]
af_clumping[,qvalue := qvalue(P)$qvalues]
af_clumping[,FDR :=p.adjust(P, method = "BH") %>% 
      format(digits = 10,scientific = T) %>% 
      as.numeric()]

fwrite(af_clumping,
       "C:/Peter/LD_clumping/QN_before_eQTL/outcome/LD_0.3/ld_clumping_1000kb_finngenPval_1000GfinLD.txt",
       row.names = F, col.names = T, sep = "\t") 


### record ----

# 紀錄不同方式snp 數量變化

# [1] 3843420
# [1] 2725011
# [1] 799764
# [1] 0
# [1] 0
# [1] 0

c <- fread("C:/Peter/LD_clumping/QN_before_eQTL/outcome/1000kb_finngenPval_bon_noClumping.txt")
# 法4,5 snp beginning 數量
uniqueN(c$fin_GWAS_hg19)


b <-  fread("C:/Peter/LD_clumping/QN_before_eQTL/outcome/1000kb_finngenPval_bon.txt")
# 法5 500kb allele篩選後 snp 數量
uniqueN(b$fin_GWAS_hg19)

a <-  fread("C:/Peter/LD_clumping/QN_before_eQTL/outcome/LD_0.3/ld_clumping_1000kb_finngenPval_1000GfinLD.txt")
# 法5 clumping 後 snp 數量
nrow(a)


# 法5 過 bon 門檻的 snp 數量
length(a[Bonfi==1,SNP])
# 法5 過 FDR 門檻的 snp 數量
length(a[FDR<0.05,SNP])
# 法5 過 qvalue 門檻的 snp 數量
length(a[qvalue<0.05,SNP])

# 釋放記憶體
rm(list=ls())
gc()




# OUTLINE: Method 4 500kb nearby SNP analysis ----
## 4. ----
# 用意：真正的 causal effect 可能在 eQTL snp 附近


# 找附近 500kb snp 產生 500kb_finngenPval_bon_noClumping.txt (.1)，
# 把 finngen GWAS allele 跟 1000G FIN genotype allele 留下 minor= ref, major=alt 或是對調的 500kb_finngenPval_bon.txt (.2)，
# 從各個染色體的 binary file 取出上一步這些 snp 的 genotype data，產生 gt_N_500kb_in_finngen_hg19SNP_chr1-22 (.3)，
# 做 LD clumping 產生 ld_clumping_500kb_finngenPval_1000GfinLD.clumped，計算 FDR, bon, qvalue 產生 ld_clumping_500kb_finngenPval_1000GfinLD.txt  (.4)

### (.1) ----
finngen_hg19 <- fread("D:/oral_cancer/snp_repeat_Finngen/outcome/C3_ORALCAVITY_EXALLC_2.txt")
gt_N_hg19 <- fread("C:/Peter/QN_before_eQTL/outcome/QN_intersect_finngen_1000G/QN_eQTL_SNPhg19.txt")


finngen_hg19 <- finngen_hg19 %>% 
  select(`#chrom`,hg19_snpID,pos,pval,ref,alt)

setnames(finngen_hg19, old = "#chrom", new = "CHR")
finngen_hg19[,start := pos ]
finngen_hg19[,end := pos ]

# 調整範圍正負500KB (5e5 bp)
names(gt_N_hg19) <- c("hg19_snpID")
gt_N_hg19[, pos := str_extract(hg19_snpID, "(?<=\\:)\\d+") %>% as.numeric()]
gt_N_hg19[, start := ifelse(pos-5e5>0, pos-5e5, 0)]
gt_N_hg19[, end := pos+5e5]
gt_N_hg19[, CHR := str_extract(hg19_snpID, ".*(?=:)") %>% as.numeric()]


# 要先排序
setkey(finngen_hg19, CHR, start, end)
setkey(gt_N_hg19, CHR, start, end)
names(gt_N_hg19) <- c("gt_N_hg19","pos","start","end","CHR")

# 語法解釋：找 gt_N 中，CHR 一樣，且 pos 在 finngen 的 lower 與 upper 之間的位點
overlap <- finngen_hg19[gt_N_hg19, 
                     on = .(CHR, pos >= start, pos <= end), 
                     nomatch = 0L,
                     allow.cartesian = TRUE]



gt_N_surround_snp <- unique(overlap,by="hg19_snpID") 
setnames(gt_N_surround_snp,old="hg19_snpID", new="fin_GWAS_hg19")
# 把 pos 刪掉，下面會得到 bim$pos
gt_N_surround_snp[, c("pos","start","end","pos.1","gt_N_hg19","i.pos") := NULL ]

gt_N_surround_snp <- na.omit(gt_N_surround_snp)
gt_N_surround_snp <- gt_N_surround_snp[rowSums(gt_N_surround_snp == "") == 0]


fwrite(gt_N_surround_snp,
       "C:/Peter/LD_clumping/QN_before_eQTL/outcome/500kb_finngenPval_bon_noClumping.txt",
       row.names = F, col.names = T, sep = "\t") 




### (.2) ----

# 迴圈抓
data_list <- vector("list", 22)

for (i in 1:22) {
  bim <- sprintf("C:/Peter/PCA_1000G_20130502/FIN_sample/trash/deal_repeatSNP/chr_%d_rename.bim",i) %>% 
    fread(.,col.names = c("CHR", "fin_GWAS_hg19", "genetic_dist", "pos", "minor", "major"))

  bim <- bim[gt_N_surround_snp, on= .(CHR,fin_GWAS_hg19), nomatch = NULL]
  
  data_list[[i]] <- bim[
    (minor == ref & major == alt) | (minor == alt & major == ref)
    ][ order(pval), .SD[1], by = fin_GWAS_hg19 ]
  
}
final_data <- rbindlist(data_list, use.names = TRUE)


fwrite(final_data,
       "C:/Peter/LD_clumping/QN_before_eQTL/outcome/500kb_finngenPval_bon.txt",
       row.names = F, col.names = T, sep = "\t") 

print(nrow(final_data))




### (.3) ----
# 1000G data (hg19)+ finngen GWAS (hg38)，轉成 hg19

sur_snp <- fread("C:/Peter/LD_clumping/QN_before_eQTL/outcome/500kb_finngenPval_bon.txt")

# 取出 附近 500kb snp hg19 存檔
fwrite(data.table(k = sur_snp$fin_GWAS_hg19),
       "C:/Peter/LD_clumping/QN_before_eQTL/trash/gt_N_500kb_in_finngen_hg19SNP.txt",
       row.names = F, col.names = F, sep = "\n") 

# 從原始 binary file 取出這些 snp
shell('for /L %i in (1,1,22) do plink --bfile "C:/Peter/PCA_1000G_20130502/FIN_sample/trash/deal_repeatSNP/chr_%i_rename" --extract "C:/Peter/LD_clumping/QN_before_eQTL/trash/gt_N_500kb_in_finngen_hg19SNP.txt" --make-bed --out "C:/Peter/LD_clumping/QN_before_eQTL/trash/gt_N_500kb_in_finngen_hg19SNP_chr%i"')

# 把chr2~22 的 binary file name 放進txt 檔案裡
mergelist <- sprintf("C:/Peter/LD_clumping/QN_before_eQTL/trash/gt_N_500kb_in_finngen_hg19SNP_chr%d.%s",
                     rep(2:22, each = 3),
                     c("bed", "bim", "fam")) %>% 
    matrix(ncol = 3, byrow=T) 

mergelist <- data.frame(mergelist)
fwrite(mergelist, "C:/Peter/LD_clumping/QN_before_eQTL/trash/500kb_mergelist.txt",
       row.names = F, col.names = F, sep = "\t")


# 合併不同chr
system("plink --bfile C:/Peter/LD_clumping/QN_before_eQTL/trash/gt_N_500kb_in_finngen_hg19SNP_chr1 --merge-list C:/Peter/LD_clumping/QN_before_eQTL/trash/500kb_mergelist.txt --make-bed --out C:/Peter/LD_clumping/QN_before_eQTL/trash/gt_N_500kb_in_finngen_hg19SNP_chr1-22")

### (.4) ----


# --clump-snp-field, --clump-field 分別指定 snp, pval 在資料的哪個欄位
system("plink --bfile C:/Peter/LD_clumping/QN_before_eQTL/trash/gt_N_500kb_in_finngen_hg19SNP_chr1-22 --clump C:/Peter/LD_clumping/QN_before_eQTL/outcome/500kb_finngenPval_bon.txt --clump-p1 0.99999 --clump-r2 0.3 --clump-kb 500 --clump-snp-field fin_GWAS_hg19 --clump-field pval --out C:/Peter/LD_clumping/QN_before_eQTL/outcome/LD_0.3/ld_clumping_500kb_finngenPval_1000GfinLD")

af_clumping <- fread("C:/Peter/LD_clumping/QN_before_eQTL/outcome/LD_0.3/ld_clumping_500kb_finngenPval_1000GfinLD.clumped")

af_clumping[,Bonfi := ifelse(P<(0.05/nrow(af_clumping)),
                             1,0)]
af_clumping[,qvalue := qvalue(P)$qvalues]
af_clumping[,FDR :=p.adjust(P, method = "BH") %>% 
      format(digits = 10,scientific = T) %>% 
      as.numeric()]

fwrite(af_clumping,
       "C:/Peter/LD_clumping/QN_before_eQTL/outcome/LD_0.3/ld_clumping_500kb_finngenPval_1000GfinLD.txt",
       row.names = F, col.names = T, sep = "\t") 

# 釋放記憶體
rm(list=ls())
gc()

### record ----

# 紀錄不同方式snp 數量變化


c <- fread("C:/Peter/LD_clumping/QN_before_eQTL/outcome/500kb_finngenPval_bon_noClumping.txt")
# 法4,5 snp beginning 數量
uniqueN(c$fin_GWAS_hg19)


b <-  fread("C:/Peter/LD_clumping/QN_before_eQTL/outcome/500kb_finngenPval_bon.txt")
# 法5 500kb allele篩選後 snp 數量
uniqueN(b$fin_GWAS_hg19)

a <-  fread("C:/Peter/LD_clumping/QN_before_eQTL/outcome/LD_0.3/ld_clumping_500kb_finngenPval_1000GfinLD.txt")
# 法5 clumping 後 snp 數量
nrow(a)


# 法5 過 bon 門檻的 snp 數量
length(a[Bonfi==1,SNP])
# 法5 過 FDR 門檻的 snp 數量
length(a[FDR<0.05,SNP])
# 法5 過 qvalue 門檻的 snp 數量
length(a[qvalue<0.05,SNP])

# OUTLINE: Sensitivity clumping by LD threshold ----
## 5. ----

### 1000kb r2 0.2  ----

# --clump-snp-field, --clump-field 分別指定 snp, pval 在資料的哪個欄位
system("plink --bfile C:/Peter/LD_clumping/QN_before_eQTL/trash/gt_N_1000kb_in_finngen_hg19SNP_chr1-22 --clump C:/Peter/LD_clumping/QN_before_eQTL/outcome/1000kb_finngenPval_bon.txt --clump-p1 0.99999 --clump-r2 0.2 --clump-kb 1000 --clump-snp-field fin_GWAS_hg19 --clump-field pval --out C:/Peter/LD_clumping/QN_before_eQTL/outcome/LD_0.2/ld_clumping_1000kb_finngenPval_1000GfinLD")

af_clumping <- fread("C:/Peter/LD_clumping/QN_before_eQTL/outcome/LD_0.2/ld_clumping_1000kb_finngenPval_1000GfinLD.clumped")

af_clumping[,Bonfi := ifelse(P<(0.05/nrow(af_clumping)),
                             1,0)]
af_clumping[,qvalue := qvalue(P)$qvalues]
af_clumping[,FDR :=p.adjust(P, method = "BH") %>% 
      format(digits = 10,scientific = T) %>% 
      as.numeric()]

fwrite(af_clumping,
       "C:/Peter/LD_clumping/QN_before_eQTL/outcome/LD_0.2/ld_clumping_1000kb_finngenPval_1000GfinLD.txt",
       row.names = F, col.names = T, sep = "\t") 


# record ----


a <-  fread("C:/Peter/LD_clumping/QN_before_eQTL/outcome/LD_0.2/ld_clumping_1000kb_finngenPval_1000GfinLD.txt")
# 法5 clumping 後 snp 數量
nrow(a)


# 法5 過 bon 門檻的 snp 數量
length(a[Bonfi==1,SNP])
# 法5 過 FDR 門檻的 snp 數量
length(a[FDR<0.05,SNP])
# 法5 過 qvalue 門檻的 snp 數量
length(a[qvalue<0.05,SNP])


### 500kb r2 0.2  ----

# 把 r2 改成 0.2

# clumping ----

# --clump-snp-field, --clump-field 分別指定 snp, pval 在資料的哪個欄位
system("plink --bfile C:/Peter/LD_clumping/QN_before_eQTL/trash/gt_N_500kb_in_finngen_hg19SNP_chr1-22 --clump C:/Peter/LD_clumping/QN_before_eQTL/outcome/500kb_finngenPval_bon.txt --clump-p1 0.99999 --clump-r2 0.2 --clump-kb 500 --clump-snp-field fin_GWAS_hg19 --clump-field pval --out C:/Peter/LD_clumping/QN_before_eQTL/outcome/LD_0.2/ld_clumping_500kb_finngenPval_1000GfinLD")

af_clumping <- fread("C:/Peter/LD_clumping/QN_before_eQTL/outcome/LD_0.2/ld_clumping_500kb_finngenPval_1000GfinLD.clumped")

af_clumping[,Bonfi := ifelse(P<(0.05/nrow(af_clumping)),
                             1,0)]
af_clumping[,qvalue := qvalue(P)$qvalues]
af_clumping[,FDR :=p.adjust(P, method = "BH") %>% 
      format(digits = 10,scientific = T) %>% 
      as.numeric()]

fwrite(af_clumping,
       "C:/Peter/LD_clumping/QN_before_eQTL/outcome/LD_0.2/ld_clumping_500kb_finngenPval_1000GfinLD.txt",
       row.names = F, col.names = T, sep = "\t") 


#  record ----

a <-  fread("C:/Peter/LD_clumping/QN_before_eQTL/outcome/LD_0.2/ld_clumping_500kb_finngenPval_1000GfinLD.txt")
# 法5 clumping 後 snp 數量
nrow(a)


# 法5 過 bon 門檻的 snp 數量
length(a[Bonfi==1,SNP])
# 法5 過 FDR 門檻的 snp 數量
length(a[FDR<0.05,SNP])
# 法5 過 qvalue 門檻的 snp 數量
length(a[qvalue<0.05,SNP])


# 釋放記憶體
rm(list=ls())
gc()

### 1000kb r2 0.1 ----
# 把 r2 改成 0.1

# clumping ----

# --clump-snp-field, --clump-field 分別指定 snp, pval 在資料的哪個欄位
system("plink --bfile C:/Peter/LD_clumping/QN_before_eQTL/trash/gt_N_1000kb_in_finngen_hg19SNP_chr1-22 --clump C:/Peter/LD_clumping/QN_before_eQTL/outcome/1000kb_finngenPval_bon.txt --clump-p1 0.99999 --clump-r2 0.1 --clump-kb 1000 --clump-snp-field fin_GWAS_hg19 --clump-field pval --out C:/Peter/LD_clumping/QN_before_eQTL/outcome/LD_0.1/ld_clumping_1000kb_finngenPval_1000GfinLD")

af_clumping <- fread("C:/Peter/LD_clumping/QN_before_eQTL/outcome/LD_0.1/ld_clumping_1000kb_finngenPval_1000GfinLD.clumped")

af_clumping[,Bonfi := ifelse(P<(0.05/nrow(af_clumping)),
                             1,0)]
af_clumping[,qvalue := qvalue(P)$qvalues]
af_clumping[,FDR :=p.adjust(P, method = "BH") %>% 
      format(digits = 10,scientific = T) %>% 
      as.numeric()]

fwrite(af_clumping,
       "C:/Peter/LD_clumping/QN_before_eQTL/outcome/LD_0.1/ld_clumping_1000kb_finngenPval_1000GfinLD.txt",
       row.names = F, col.names = T, sep = "\t") 


# record ----


a <-  fread("C:/Peter/LD_clumping/QN_before_eQTL/outcome/LD_0.1/ld_clumping_1000kb_finngenPval_1000GfinLD.txt")
# 法5 clumping 後 snp 數量
nrow(a)


# 法5 過 bon 門檻的 snp 數量
length(a[Bonfi==1,SNP])
# 法5 過 FDR 門檻的 snp 數量
length(a[FDR<0.05,SNP])
# 法5 過 qvalue 門檻的 snp 數量
length(a[qvalue<0.05,SNP])


### 500kb r2 0.1   ----

# clumping ----

# --clump-snp-field, --clump-field 分別指定 snp, pval 在資料的哪個欄位
system("plink --bfile C:/Peter/LD_clumping/QN_before_eQTL/trash/gt_N_500kb_in_finngen_hg19SNP_chr1-22 --clump C:/Peter/LD_clumping/QN_before_eQTL/outcome/500kb_finngenPval_bon.txt --clump-p1 0.99999 --clump-r2 0.1 --clump-kb 500 --clump-snp-field fin_GWAS_hg19 --clump-field pval --out C:/Peter/LD_clumping/QN_before_eQTL/outcome/LD_0.1/ld_clumping_500kb_finngenPval_1000GfinLD")

af_clumping <- fread("C:/Peter/LD_clumping/QN_before_eQTL/outcome/LD_0.1/ld_clumping_500kb_finngenPval_1000GfinLD.clumped")

af_clumping[,Bonfi := ifelse(P<(0.05/nrow(af_clumping)),
                             1,0)]
af_clumping[,qvalue := qvalue(P)$qvalues]
af_clumping[,FDR :=p.adjust(P, method = "BH") %>% 
      format(digits = 10,scientific = T) %>% 
      as.numeric()]

fwrite(af_clumping,
       "C:/Peter/LD_clumping/QN_before_eQTL/outcome/LD_0.1/ld_clumping_500kb_finngenPval_1000GfinLD.txt",
       row.names = F, col.names = T, sep = "\t") 


#  record ----



a <-  fread("C:/Peter/LD_clumping/QN_before_eQTL/outcome/LD_0.1/ld_clumping_500kb_finngenPval_1000GfinLD.txt")
# 法5 clumping 後 snp 數量
nrow(a)


# 法5 過 bon 門檻的 snp 數量
length(a[Bonfi==1,SNP])
# 法5 過 FDR 門檻的 snp 數量
length(a[FDR<0.05,SNP])
# 法5 過 qvalue 門檻的 snp 數量
length(a[qvalue<0.05,SNP])



# OUTLINE: Add chr5 rs4975538 eQTL information ----
## 6. ----

# 補上有在finngen 顯著的snp 5:1333830 (hg18) eQTL info

qn_n <- fread("C:/Peter/result/QN_maf_gt_N_cis_chr5.txt")
qn_t <- fread("C:/Peter/result/QN_maf_gt_T_cis_chr5.txt")
gt_n <- fread("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/Normal/maf_gt_N_cis_chr5.txt")
gt_t <- fread("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/Tumor/maf_gt_T_cis_chr5.txt")

qn_n <- qn_n[SNP=="5:1333830",]
qn_t <- qn_t[SNP=="5:1333830",]
gt_n <- gt_n[SNP=="5:1333830",]
gt_t <- gt_t[SNP=="5:1333830",]

qn_n[,gene := gsub("^([^_]+_[^_]+).*", "\\1",gene)]
qn_n <- unique(qn_n)
qn_t[,gene := gsub("^([^_]+_[^_]+).*", "\\1",gene)]
qn_t <- unique(qn_t)
gt_n[,gene := gsub("^([^_]+_[^_]+).*", "\\1",gene)]
gt_n <- unique(gt_n)
gt_t[,gene := gsub("^([^_]+_[^_]+).*", "\\1",gene)]
gt_t <- unique(gt_t)


fwrite(qn_n,
       "C:/Peter/LD_clumping/QN_before_eQTL/outcome/fin_mostSNP_eQTL/qn_n_5.txt",
       row.names = F, col.names = T, sep = "\t") 

fwrite(qn_t,
       "C:/Peter/LD_clumping/QN_before_eQTL/outcome/fin_mostSNP_eQTL/qn_t_5.txt",
       row.names = F, col.names = T, sep = "\t") 

fwrite(gt_n,
       "C:/Peter/LD_clumping/QN_before_eQTL/outcome/fin_mostSNP_eQTL/NOqn_n_5.txt",
       row.names = F, col.names = T, sep = "\t") 

fwrite(gt_t,
       "C:/Peter/LD_clumping/QN_before_eQTL/outcome/fin_mostSNP_eQTL/NOqn_t_5.txt",
       row.names = F, col.names = T, sep = "\t") 



 
### function  ----

mix_info <- function(inputname, outputname){
  
  if (!is.character(inputname) || length(inputname) != 1 || !is.character(outputname) || length(outputname) != 1) {
    stop("inputname, outputname 必須是一個字串")
  }
  
  fdr <- fread(inputname, header=T)
  
  #merge, due to some snp MAF<0.05, 
  fdr <- merge(fdr, pvalue_oqn, by.x = "gene", by.y = "PROBE_ID")
  fdr <- merge(fdr, probe_pos, by.x = "gene", by.y = "Gene")
  fdr <- merge(fdr, maf, by.x = "SNP", by.y = "hg18_snpID")
  fdr <- merge(fdr, probe_info, by.x = "gene", by.y = "PROBE_ID")
  
  
  
  # 調整
  fdr[,c("start","end","INFO","ER2","chr") := NULL]
  setnames(fdr, old = c("gene", "CHROMOSOME"), new = c("Probe", "CHR"))
  
 
  setcolorder(fdr,c("CHR","ProbeID","Probe","Gene","PROBE_COORDINATES",
                    "SNP","rsID","R2","impute_type","MAF","REF","ALT","beta","t-stat",
                    "p-value","FDR",
                    "pval_QN","t_QN","pval_QN_BH","sig_QN","sig_QN_Bonfi"))
  
  fdr <- fdr[order(FDR),]
  
  # 調整數值
  fdr[,MAF := round(MAF, digits = 4)]
  fdr[,beta := round(beta, digits = 4)]
  fdr[,`t-stat` := round(`t-stat`, digits = 4)]
  fdr[,`p-value` := format(`p-value`, digits = 4, scientific = T)]
  fdr[,FDR := format(FDR, digits = 4, scientific = T)]
  fdr[,pval_QN := format(pval_QN, digits = 4, scientific = T)]
  fdr[,t_QN := format(t_QN, digits = 4, scientific = T)]
  fdr[,R2 := round(R2, digits = 4)]
  
  
  fwrite(fdr, outputname,
         row.names = F, col.names = T, sep = "\t")
}








pvalue_oqn <- fread("C:/Peter/QN_before_eQTL/trash/QN_exp_different.txt",header=T)

# 把各檔案cis snp number 資料刪掉，info 檔案用不到
pvalue_oqn <- pvalue_oqn[,c("QN_maf_gt_N_pvalue_FDR_sigCis-SNP_number","QN_maf_gt_T_pvalue_FDR_sigCis-SNP_number",
                            "QN_maf_gt_N_MOSTpvalue_FDR_sigCis-SNP_number","QN_maf_gt_T_MOSTpvalue_FDR_sigCis-SNP_number") := NULL]
  
  
probe_pos <- fread("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/probe_pos.txt",header=T)
maf <- fread("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/cis_hg18_19_snp_MAF.txt", header=T)
maf[,hg19_snpID := NULL]

probe_info <- fread("D:/oral_cancer/expression/expression_data/probe24526infoB.txt", header=T)
probe_info[,c("TargetID","CHROMOSOME") := NULL]


file_name <- c("C:/Peter/LD_clumping/QN_before_eQTL/outcome/fin_mostSNP_eQTL/NOqn_n_5",
               "C:/Peter/LD_clumping/QN_before_eQTL/outcome/fin_mostSNP_eQTL/NOqn_t_5",
               "C:/Peter/LD_clumping/QN_before_eQTL/outcome/fin_mostSNP_eQTL/qn_n_5",
               "C:/Peter/LD_clumping/QN_before_eQTL/outcome/fin_mostSNP_eQTL/qn_t_5")


sig_number <- fread("C:/Peter/QN_before_eQTL/trash/QN_exp_different.txt",header=T)
sig_number <- sig_number %>% 
  select(c("PROBE_ID",
           "QN_maf_gt_N_pvalue_FDR_sigCis-SNP_number",
           "QN_maf_gt_T_pvalue_FDR_sigCis-SNP_number"))
  
NOqn_sig_number <- fread("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/exp_different.txt",header=T)
NOqn_sig_number <- NOqn_sig_number %>% 
  select(c("PROBE_ID",
           "maf_gt_N_pvalue_FDR_sigCis-SNP_number",
           "maf_gt_T_pvalue_FDR_sigCis-SNP_number"))

sig_number_final <- merge(NOqn_sig_number,sig_number, by="PROBE_ID")

# 對這些檔案生成info 檔
for (idx in 1:length(file_name)) {
  i <- file_name[idx]
  
  mix_info(paste0(i,".txt"),
           paste0(i,"_info.txt"))
  
  # 新增 cis sig snp
  a <- paste0(i,"_info.txt") %>% 
    fread()
  
  # 挑出1, idx+1 col
  b <- sig_number_final[, .SD, .SDcols = c(1, idx+1)]
  names(b) <- c("PROBE_ID", "sig_Cis-SNP_number")
  a <- merge(a,b,by.x="Probe",by.y="PROBE_ID")
  
  
  
  #setcolorder(a,c("ProbeID","Probe","Gene","CHR","PROBE_COORDINATES","sig_Cis-SNP_number","SNP","rsID","impute_type","R2","MAF","REF","ALT","beta","t-stat","p-value","FDR","t_QN","pval_QN","pval_QN_BH","sig_QN","sig_QN_Bonfi"))
  
  fwrite(a, paste0(i,"_info.txt"),
         row.names = F, col.names = T, sep = "\t")
}






# OUTLINE: Compute additional FDR summaries ----
### 補算 FDR ----
# 原本流程只保留 FDR<0.05 snp，重用該染色體的結果跑 FDR

compute_FDR <- function(inputname,outputname, qvalue_type){
  
  if (!is.character(inputname) || length(inputname) != 1 || !is.character(outputname) || length(outputname) != 1) {
    stop("inputname, outputname 必須是一個字串")
  }
  
  x <- fread(inputname,header = T)
  
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
  
  
  x <- x[order(FDR),]
  x[, gene2 := gsub("^([^_]+_[^_]+).*", "\\1",gene)]
  
  
  fwrite(x, outputname,row.names = F, col.names = T, sep = "\t")
}






file_name_all <- c("QN_maf_gt_N_pvalue",
               "QN_maf_gt_T_pvalue")



for (i in file_name_all) {
  compute_FDR(paste0("C:/Peter/QN_before_eQTL/trash/",i,".txt"),
           paste0("C:/Peter/QN_before_eQTL/trash/",i,"_FDR2.txt"),2)
}



# 只用chr5 的全部結果算FDR ----

fdr_n <- fread("C:/Peter/QN_before_eQTL/trash/QN_maf_gt_N_pvalue_FDR2.txt")
fdr_n <- fdr_n[SNP=="5:1333830",]
fdr_n <- fdr_n %>% select("gene2","FDR","qvalue")
qn_n <- fread("C:/Peter/LD_clumping/QN_before_eQTL/outcome/fin_mostSNP_eQTL/qn_n_5_info.txt")
qn_n[,FDR:= NULL]

fdr_n <- merge(qn_n, fdr_n,by.x="Probe",by.y="gene2")

fwrite(fdr_n,
       "C:/Peter/LD_clumping/QN_before_eQTL/outcome/fin_mostSNP_eQTL/qn_n_5_FDR.txt",
       row.names = F, col.names = T, sep = "\t") 




fdr_t <- fread("C:/Peter/QN_before_eQTL/trash/QN_maf_gt_T_pvalue_FDR2.txt")
fdr_t <- fdr_t[SNP=="5:1333830",]
fdr_t <- fdr_t %>% select("gene2","FDR","qvalue")
qn_t <- fread("C:/Peter/LD_clumping/QN_before_eQTL/outcome/fin_mostSNP_eQTL/qn_t_5_info.txt")
qn_t[,FDR:= NULL]

fdr_t <- merge(qn_t, fdr_t,by.x="Probe",by.y="gene2")

fwrite(fdr_t,
       "C:/Peter/LD_clumping/QN_before_eQTL/outcome/fin_mostSNP_eQTL/qn_t_5_FDR.txt",
       row.names = F, col.names = T, sep = "\t") 






# 只用 5:1333830 結果算FDR ----
a <- fread("C:/Peter/LD_clumping/QN_before_eQTL/outcome/fin_mostSNP_eQTL/qn_t_5_FDR.txt")
a[,c("qvalue","FDR") := NULL]
a[,FDR :=p.adjust(`p-value`, method = "BH") %>% 
      format(digits = 10,scientific = T) %>% 
      as.numeric()]

fwrite(a,
       "C:/Peter/LD_clumping/QN_before_eQTL/outcome/fin_mostSNP_eQTL/qn_t_5_FDR.txt",
       row.names = F, col.names = T, sep = "\t") 


a <- fread("C:/Peter/LD_clumping/QN_before_eQTL/outcome/fin_mostSNP_eQTL/qn_n_5_FDR.txt")
a[,c("qvalue","FDR") := NULL]
a[,FDR :=p.adjust(`p-value`, method = "BH") %>% 
      format(digits = 10,scientific = T) %>% 
      as.numeric()]

fwrite(a,
       "C:/Peter/LD_clumping/QN_before_eQTL/outcome/fin_mostSNP_eQTL/qn_n_5_FDR.txt",
       row.names = F, col.names = T, sep = "\t") 









