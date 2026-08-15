# ============================================================
# TITLE: Imputation Pre-processing
# SUBTITLE: Build PLINK MAP/PED inputs from oral cancer genotype calls.
# SEARCH TAGS: imputation, preprocessing, PLINK, MAP, PED, genotype
# NOTE: Code below is unchanged; only this navigation header was added.
# ============================================================

# OUTLINE: Load required packages ----
## package ----
library(data.table)
library(dplyr)
library(stringr)




# OUTLINE: Build initial PLINK MAP and PED ----
## Make .map and .ped ----

setwd("D:/oral_cancer/imputation_preprocessing") 
OCdata <- fread("GT9810_FinalReport_clean.txt",header = T)
colnames(OCdata) <- c("snp","sample","A1","A2")


OCdata$allele <-  paste(OCdata$A1,OCdata$A2)

# 用矩陣byrow 依序放allele
OC_ped <- matrix(OCdata$allele %>% 
                   unlist(),
                 ncol = 345111,
                 byrow = TRUE)

OC_ped <- as.data.table(OC_ped)
colnames(OC_ped) <- paste0("snp",1:345111)


OC_ped$FID <- unique(OCdata$sample)
OC_ped$IID <- unique(OCdata$sample)
OC_ped$PaID <- 0
OC_ped$MaID <- 0
OC_ped$sex <- 1
OC_ped$pheno <- 2

# 調整變數順序
setcolorder(OC_ped, c("FID", "IID","PaID","MaID","sex","pheno",paste0("snp",1:345111))) 

fwrite(OC_ped, "OC_case.ped", sep = "\t", col.names = FALSE)
rm(OC_ped)


setwd("D:/oral_cancer/imputation_preprocessing") 

snp_data <- fread("SNP Table.txt")
snp_data$GD <- 0
setcolorder(snp_data, c("Chr","Name","GD","Position"))
names(snp_data)
fwrite(snp_data, "OC_case.map", sep = "\t", col.names = FALSE)
rm(snp_data)


# OUTLINE: Run initial PLINK QC filters ----
## QC(1) ----

setwd("D:/oral_cancer/imputation_preprocessing")

system("plink --file D:/oral_cancer/imputation_preprocessing/OC_case --make-bed --out D:/oral_cancer/imputation_preprocessing/OC_case")

# Exclude duplicate SNPs, no duplicate SNPs
system("plink --bfile OC_case --missing --out output1")

# Remove duplicate samples
system("plink --bfile OC_case --chr 1-22 --indep-pairwise 50000 5000 0.2 --out output2")

# no PI_HAT between samples >0.8
system("plink --bfile OC_case --extract output2.prune.in --genome --min 0.8 --out output3")
system("plink --bfile OC_case --extract output2.prune.in --genome --out aa")

# Exclude SNPs not on chr 1-24
system("plink --bfile OC_case --chr 1-24 --make-bed --out OC_case3")

# Exclude Indels or multi-base alleles
snp_allele <- fread("OC_case3.bim")
colnames(snp_allele) <- c("chr", "id","IDoNotKnow","pos","A1","A2")

# Exclude "-"
snp_not_measure <- which(snp_allele$A1=="-" | snp_allele$A2=="-") %>% 
  unique()

# no snp is multi-base (like AGAA)
nchar(snp_allele$A1) %>% 
  table()
nchar(snp_allele$A2) %>% 
  table()

# save SNPs is "-"
fwrite(snp_allele$id[snp_not_measure] %>% 
         as.data.frame(),
       "snp_not_measure.txt", sep = "\n", col.names = FALSE)

# Exclude Indels or multi-base alleles
system("plink --bfile OC_case3 --exclude snp_not_measure.txt --make-bed --out OC_case4")
rm(snp_not_measure,snp_allele)

# OUTLINE: Run sample and SNP QC filters ----
### QC(2) ----

setwd("D:/oral_cancer/imputation_preprocessing")

# Exclude SNPs with call rates < 95%
system("plink --bfile OC_case4 --geno 0.05 --make-bed --out OC_case5")

# Exclude SNPs with HWE p-values<$10^{-4}$
system("plink --bfile OC_case5 --mind 0.05 --hwe 0.0001 --make-bed --out output5")
system("plink --bfile OC_case5 --extract output5.bim --make-bed --out OC_case6")

# Remove samples with sex discrepancy
system("plink --bfile OC_case6 --chr 23,24 --maf 0.05 --make-bed --out output6")

sex_chr_maf <- fread("D:/oral_cancer/imputation_preprocessing/output6.bim")

# save SNPs in 23,24 chr is maf >= 0.05
fwrite(sex_chr_maf$V2 %>% 
         as.data.frame(),
       "SNPmaf.txt", sep = "\n", col.names = FALSE)

system("plink --bfile OC_case6 --extract SNPmaf.txt --indep-pairwise 50000 5000 0.2 --out output7")
system("plink --bfile OC_case6 --extract output7.prune.in --check-sex ycount --out output8")
rm(sex_chr_maf)

# Remove samples with call rates < 95%
system("plink --bfile OC_case6 --chr 1-22 --maf 0.05 --make-bed --out output9")
system("plink --bfile OC_case6 --extract output9.bim --mind 0.05 --make-bed --out output10")
system("plink --bfile OC_case6 --keep output10.fam --make-bed --out OC_case8")


# Remove samples with outlying heterozygosity (3/6 SD)
system("plink --bfile OC_case8 --chr 1-22 --maf 0.05 --het --out output12")
system("plink --bfile OC_case8 --chr 1-22 --maf 0.05 --make-bed --out output13")
system("plink --bfile output13 --missing --out output14")

# no sample heterozygosity over 3 times standard deviation
heterozygosity <- fread("D:/oral_cancer/imputation_preprocessing/output12.het")
heterozygosity$he <- 1-heterozygosity$`O(HOM)`/heterozygosity$`N(NM)`
upper <- mean(heterozygosity$he)+3*sd(heterozygosity$he)
lower <- mean(heterozygosity$he)-3*sd(heterozygosity$he)

which(heterozygosity$he<lower | heterozygosity$he > upper) %>% 
  length()

rm(heterozygosity,upper,lower)


system("plink --bfile OC_case8  --make-bed --out afGenotypeQC")

# QC 好的 snp ----
afGenotypeQC <- fread("D:/oral_cancer/imputation_preprocessing/afGenotypeQC.bim")

# 把QC後的snp改成輸入格式為，chr1:8207084-8207084
v36_transformed <- paste0("chr",
                          afGenotypeQC$V1,
                          ":",
                          afGenotypeQC$V4,
                          "-",
                          afGenotypeQC$V4)


fwrite(v36_transformed %>% 
         as.data.frame(),
       "afQC_hg18.txt", sep = "\n", col.names = FALSE)




# OUTLINE: Remove liftover failed SNPs ----
## Delete snp ----

# 刪掉轉換失敗的snp，版本轉換網址 <https://genome.ucsc.edu/cgi-bin/hgLiftOver>
setwd("D:/oral_cancer/imputation_preprocessing") 


# 從網址下載hglft_genome.err.txt，紀錄轉換失敗的snp
failed <- fread("D:/oral_cancer/imputation_preprocessing/hglft_genome.err.txt",header=F)
v36_transformed <- fread("D:/oral_cancer/imputation_preprocessing/afQC_hg18.txt",header=F)

# 篩選開頭是 chr 的元素
a <- grep("^chr",
          failed$V1,
          value = TRUE)

# 挑出轉換成功的snp，重新轉換一次，確保沒有轉換失敗的
success <- setdiff(v36_transformed$V1,a)
fwrite(success %>% 
         as.data.frame(),
       "afQC_hg18.txt", sep = "\n", col.names = FALSE)

v36_transformed$order <- c(1:nrow(v36_transformed))


# 從網站下載，轉換成功的snp file，命名 afQC_hg19.txt



# OUTLINE: Generate hg19 MAP file ----
## .map generating ----

# 生出轉換成功的 .map
setwd("D:/oral_cancer/imputation_preprocessing") 


afGenotypeQC <- fread("D:/oral_cancer/imputation_preprocessing/afGenotypeQC.bim",header = F)
# GD is genetic distance
names(afGenotypeQC) <- c("chr","snp","GD","pos","A1","A2")

# afQC_hg18 是所有成功轉換的snp，pos version is hg18
# afQC_hg19 是所有成功轉換的snp，pos version is hg19
afQC_hg18 <- fread("afQC_hg18.txt",header = F)
afQC_hg19 <- fread("afQC_hg19.txt",header = F)

# 挑出 chr7:105746325-105746325 中 chr, pos
afQC_hg19 <- str_match(afQC_hg19$V1, "chr(\\d+):(\\d+)-") %>% 
  as.data.table()
afQC_hg19 <- afQC_hg19[,-1]
names(afQC_hg19) <- c("chr","pos")

afQC_hg18 <- str_match(afQC_hg18$V1, "chr(\\d+):(\\d+)-") %>% 
  as.data.table()
afQC_hg18 <- afQC_hg18[,-1]
names(afQC_hg18) <- c("chr","pos")
afQC_hg18$success <- 1


afGenotypeQC$pos <- afGenotypeQC$pos %>% 
  as.character()
afGenotypeQC$chr <- afGenotypeQC$chr %>% 
  as.character()

#  用chr, pos合併且保持snp 順序
afQC_hg18_bim <- left_join(afGenotypeQC,afQC_hg18,by = c("chr","pos"))
afQC_hg18_bim <- afQC_hg18_bim[!is.na(success)]

# pos 轉換版本
afQC_hg18_bim$pos <- afQC_hg19$pos

# build map
afQC_hg19_map <- data.table(Chr = afQC_hg18_bim$chr,
                           Name = afQC_hg18_bim$snp,
                           GD = 0,
                           Position = afQC_hg18_bim$pos)
fwrite(afQC_hg19_map, "afQC_hg19.map", sep = "\t", col.names = FALSE)


# OUTLINE: Generate hg19 PED file ----
## .ped generating ----

# 生出轉換成功的 .ped
setwd("D:/oral_cancer/imputation_preprocessing") 

# build ped
system("plink --bfile afGenotypeQC --recode --out afGenotypeQC")


# ped 讀取會把allele1, allele2 分開，把allele 合併起來。afGenotypeQC_ped_1 紀錄QC後還沒轉換的所有樣本的allele
afGenotypeQC_ped <- fread("D:/oral_cancer/imputation_preprocessing/afGenotypeQC.ped",header = F)
afGenotypeQC_ped_1 <- sapply(seq(7, ncol(afGenotypeQC_ped), by=2),
                             function(i){
                               paste(afGenotypeQC_ped[[i]],
                                     afGenotypeQC_ped[[i+1]],
                                     sep = " ")} )
afGenotypeQC_ped_1 <- afGenotypeQC_ped_1 %>% 
  as.data.table()


# 用afGenotypeQC.bim 和afQC_hg19_map，找出轉換成功的snp id。因為afGenotypeQC.ped 的snp順序跟afGenotypeQC.bim 一樣，再回到afGenotypeQC.ped 刪掉轉換失敗的
afQC_hg19_ped <- (afGenotypeQC$snp %in% afQC_hg19_map$Name) %>% 
  which() %>% 
  afGenotypeQC_ped_1[,.]

names(afQC_hg19_ped) <- paste0("snp",1:332614)
afQC_hg19_ped$FID <- unique(afGenotypeQC_ped$V1)
afQC_hg19_ped$IID <- unique(afGenotypeQC_ped$V1)
afQC_hg19_ped$PaID <- 0
afQC_hg19_ped$MaID <- 0
afQC_hg19_ped$sex <- 1
afQC_hg19_ped$pheno <- 2

# 調整變數順序
setcolorder(afQC_hg19_ped, c("FID", "IID", "PaID", "MaID", "sex", "pheno",
                             paste0("snp",1:332614))) 

fwrite(afQC_hg19_ped, "afQC_hg19.ped", sep = "\t", col.names = FALSE)
rm(aa,afGenotypeQC,afGenotypeQC_ped,afGenotypeQC_ped_1,afQC_hg18,afQC_hg19, success_snp_order)

# 把 map, ped 轉成 bed, bim, fam
system("plink --file afQC_hg19 --make-bed --out afQC_hg19")

# OUTLINE: Split final data by chromosome ----
## split chr ----

setwd("D:/oral_cancer/imputation_preprocessing") 

system("plink --bfile afQC_hg19 --chr 1-22 --maf 0.05 --make-bed --out ForImputation")

for (chr in 1:22) {
  cmd <- sprintf("plink --bfile ForImputation --chr %d --recode vcf --out ForImputation_chr%d", chr, chr)
  system(cmd)
}





