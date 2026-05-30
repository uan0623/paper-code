# ============================================================
# TITLE: Compute LD for GWAS Catalog
# SUBTITLE: Calculate LD around GWAS Catalog variants and summarize nearby signals.
# SEARCH TAGS: LD, GWAS Catalog, SNP, R2, linkage disequilibrium
# NOTE: Code below is unchanged; only this navigation header was added.
# ============================================================

library(data.table)
library(dplyr)
library(tidyverse)
library(stringr) 
library(ggplot2)
library(ggrepel)
library(matrixStats)  # rowMins
library(qvalue)


# OUTLINE: Add ancestry/race labels ----
## add race ----
# all snp position version are hg19(GRCH38)

loci_record <- fread("D:/oral_cancer/expression/expression_data/EFO_0005570_associations_export.tsv")

study_mapping <- c("Diversity and longitudinal records: Genetic architecture of disease associations and polygenic risk in the Taiwanese Han population." = "Taiwanese_Han",
                   "A genome-wide association study identifies two novel susceptible regions for squamous cell carcinoma of the head and neck."  = "non-Hispanic_white",
                   "Host Genetic Associations with Salivary Microbiome in Oral Cancer." = "Taiwanese",
                   "Diversity and scale: Genetic architecture of 2068 traits in the VA Million Veteran Program." = "AFR_AMR_EAS_EUR",
                   "Genome-wide association analyses identify new susceptibility loci for oral cavity and pharyngeal cancer." = "EUR_and_America",
                   "Germline determinants of humoral immune response to HPV-16 protect against oropharyngeal cancer." = "EUR_and_America")

loci_record[,race:= study_mapping[STUDY]]


# OUTLINE: Map GWAS Catalog variants to GRCh37 ----
## get GRCh37 ID ----


# 去除重複
rs_list <- unique(loci_record$SNPS)

get_grch37_coord <- function(rsid) {
  
  numeric_id <- gsub("rs", "", rsid)
  url <- paste0("https://api.ncbi.nlm.nih.gov/variation/v0/beta/refsnp/", numeric_id)
  
  resp <- GET(url)
  if (resp$status_code != 200) {
    return(data.frame(rsid=rsid, chr=NA, pos=NA))
  }
  
  json <- tryCatch(
    fromJSON(content(resp, "text", encoding="UTF-8")),
    error=function(e) return(NULL)
  )
  if (is.null(json)) {
    return(data.frame(rsid=rsid, chr=NA, pos=NA))
  }
  
  placements <- json$primary_snapshot_data$placements_with_allele
  if (length(placements) == 0) {
    return(data.frame(rsid=rsid, chr=NA, pos=NA))
  }
  
  # ----------- 找 GRCh37 使用 seq_id_traits_by_assembly -------
  
  grch37_index <- NULL
  
  for (i in 1:nrow(placements)) {
    pa <- placements$placement_annot$seq_id_traits_by_assembly[[i]]
    
    if (length(pa)!=0) {
      assembly_names <- pa$assembly_name
      if (assembly_names == "GRCh37.p13") {
        grch37_index <- i
        break
      }
    }
  }
  
  if (is.null(grch37_index)) {
    return(data.frame(rsid=rsid, chr=NA, pos=NA))
  }
  
  # ----------- 抽取 GRCh37 座標 -------
  
  target <- placements[grch37_index,]
  
  # seq_id 例如 "NC_000020.10" → chr 20
  chr_raw <- target$seq_id
  chr <- gsub("NC_0+|\\..*$", "", chr_raw)  # 最穩定抽法
  
  # 取第一個 allele 的位置（GRCh37）
  allele <- target$alleles[[1]]
  pos <- allele$allele$spdi$position + 1   # NCBI 是 0-based，要轉為 1-based
  
  data.frame(
    rsid = rsid,
    chr = chr,
    pos = pos
  )
}

results <- list()
for (i in rs_list) {
  results[[i]] <- get_grch37_coord(i)
  Sys.sleep(0.5)  # 每次呼叫間隔 0.5 秒 → 每秒最多 2 次
}

final_df <- bind_rows(results) %>% 
  unique()
rownames(final_df) <- NULL


# 把hg19/GRCH37 id 併進資料
final_df <- as.data.table(final_df)
final_df[,SNP_GRCH37 := paste0(chr,":",pos)]
final_df[,c("chr","pos") := NULL] 
loci_record <- merge(loci_record,final_df,by.x="SNPS", by.y="rsid")




fwrite(loci_record,
       "D:/oral_cancer/expression/expression_data/EFO_0005570_associations_export.tsv",
         row.names = F, col.names = T, sep = "\t")

# OUTLINE: Compute LD for overlapping variants ----
## compute LD ----
# EFO_0005570_associations_export.tsv 調整欄位順序，且新增變數，紀錄GWAScatalog 已知的口腔癌位點+ 我們eQTL發現的位點 in GRCh37 ID，變成 EFO_0005570_associations_export.txt

### 找 GWAS catalog 跟我們研究重複的snp ----

file_name <- c("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/maf_gt_N_pvalue_FDR_info.txt",
               "D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/maf_gt_T_pvalue_FDR_info.txt")

loci_record <- fread("D:/oral_cancer/expression/expression_data/EFO_0005570_associations_export.tsv")


snp_list <- c()
for (i in 1:length(file_name)) {
  a <- fread(file_name[i])
  snp_list <- c(snp_list,a$rsID)
}


loci_record[,in_our_eQTL := loci_record$SNPS %in% snp_list]

# adjust col order
setcolorder(loci_record, 
            c("SNPS","SNP_GRCH37","in_our_eQTL","race","MAPPED_GENE","MAPPED_TRAIT",
              "DISEASE/TRAIT","SNP_GENE_IDS",
              setdiff(names(loci_record), c("SNPS","SNP_GRCH37","in_our_eQTL","race",
                                            "MAPPED_GENE","MAPPED_TRAIT","DISEASE/TRAIT","SNP_GENE_IDS")))
            )

fwrite(loci_record,
       "D:/oral_cancer/expression/expression_data/EFO_0005570_associations_export.txt",
         row.names = F, col.names = T, sep = "\t")



# OUTLINE: Run LD calculation ----
### compute LD ----

#### 取出 hg18 的snp id ----
# 使用 rsID, hg18, hg19 cis snp 對照表，找出 GWAS catalog + 我們 eQTL FDR顯著結果的交集 snp in hg18

rm(list = ls())


file_name <- c("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/maf_gt_N_pvalue_FDR_info.txt",
               "D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/maf_gt_T_pvalue_FDR_info.txt")
loci_record <- fread("D:/oral_cancer/expression/expression_data/EFO_0005570_associations_export.tsv")
snp_list <- c()
for (i in 1:length(file_name)) {
  a <- fread(file_name[i])
  snp_list <- c(snp_list,a$rsID)
}

snp_list <- snp_list[!is.na(snp_list) & snp_list != ""] 

maf_snp <- fread("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/cis_hg18_19_snp_MAF.txt")
GWAScatalog_ourResult_rsID <- maf_snp[rsID %in% snp_list]

# 取出 hg18 的snp id
fwrite(data.table(snp = GWAScatalog_ourResult_rsID$hg18_snpID),
       "D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/loci_record_GWAScatalog.txt",
         row.names = F, col.names = F, sep = "\n")


#### 算 LD ----
# 使用轉換成 hg18 的 binary file
# D:/oral_cancer/expression/trash/chr1-22_imputation_hg18



system("plink --bfile D:/oral_cancer/expression/trash/chr1-22_imputation_hg18 --extract D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/loci_record_GWAScatalog.txt --make-bed --out D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/loci_record_GWAScatalog")

# 計算任兩個1000 KB snp LD
system("plink --bfile D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/loci_record_GWAScatalog --r2 --ld-window-kb 1000 --ld-window 999999 --ld-window-r2 0 --out D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/loci_record_GWAScatalog_LD_all")

# 做 ld pruning, 刪掉50000 snp 內 ld>0.2 的 snp
system("plink --bfile D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/loci_record_GWAScatalog --indep-pairwise 50000 5000 0.2 --out D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/loci_record_GWAScatalog_ld")



## add  ----
# - info 檔案加新變數，判斷哪些 FDR<0.05 eQTL snp 跟 GWAScatalog 位點 LD<0.2
# - 刪掉REF, ALT 不是 a,t,c,g 的 rsID

a <- fread("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/maf_gt_N_pvalue_FDR_info.txt")
prunein <- fread("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/loci_record_GWAScatalog_ld.prune.in",header=F)

a[,GWAScatalog_ld := a$SNP %in% prunein$V1]

# 刪掉REF, ALT 不是 a,t,c,g 的 rsID
a[!a$REF %in% c("A", "T", "C", "G") | !a$ALT %in% c("A", "T", "C", "G"),
  rsID := NA_character_]

fwrite(a,
       "D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/maf_gt_N_pvalue_FDR_info_ld.txt",
         row.names = F, col.names = T, sep = "\t")


# T ----
a <- fread("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/maf_gt_T_pvalue_FDR_info.txt")
prunein <- fread("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/loci_record_GWAScatalog_ld.prune.in",header=F)

a[,GWAScatalog_ld := a$SNP %in% prunein$V1]

# 刪掉REF, ALT 不是 a,t,c,g 的 rsID
a[!a$REF %in% c("A", "T", "C", "G") | !a$ALT %in% c("A", "T", "C", "G"),
  rsID := NA_character_]
fwrite(a,
       "D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/maf_gt_T_pvalue_FDR_info_ld.txt",
         row.names = F, col.names = T, sep = "\t")


# OUTLINE: Link GWAS Catalog variants to local results ----
## Link to our result ----

# OUTLINE: Define FDR calculation ----
### compute_FDR 函數 ----
# - 用 BH 估計 FDR
# - 先只跑 normal part

# 計算並挑出FDR<0.05 的row
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
  
  
  # x <- x %>% filter(FDR< 0.05)
  x <- x[order(FDR),]
  x[, gene2 := gsub("^([^_]+_[^_]+).*", "\\1",gene)]
  
  
  fwrite(x, outputname,row.names = F, col.names = T, sep = "\t")
}


#### 計算 FDR ----

# 對這些檔案生成有FDR 檔
compute_FDR("D:/oral_cancer/expression/trash/maf_gt_N_pvalue.txt",
           "D:/oral_cancer/expression/trash/maf_gt_N_pvalue_noFDRfilter.txt",2)




# OUTLINE: Define ordered quantile normalization ----
### OQN 函數 ----
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






### normal_quantile 函數 ----


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


# OUTLINE: Merge annotation and result information ----
### mix_info function  ----
rm(list=ls())
gc()

# 找這75位點，在我們eQTL 的結果
loci_record <- fread("D:/oral_cancer/expression/expression_data/EFO_0005570_associations_export.txt")
maf_snp <- fread("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/cis_hg18_19_snp_MAF.txt")
GWAScatalog_maf <- maf_snp[rsID %in% loci_record$SNPS]


mix_info <- function(inputname, outputname){
  
  if (!is.character(inputname) || length(inputname) != 1 || !is.character(outputname) || length(outputname) != 1) {
    stop("inputname, outputname 必須是一個字串")
  }
  
  fdr <- fread(inputname, header=T)
  
  # 只要這75位點有在我們eQTL 的
  fdr <- fdr[SNP %in% GWAScatalog_maf$hg18_snpID]
  
  # 對於多區間的probe 名稱"ILMN_2215025_2" 取出 "ILMN_2215025"
  fdr[, gene2 := gsub("^([^_]+_[^_]+).*", "\\1",gene)]
  
  #merge, due to some snp MAF<0.05, 
  fdr <- merge(fdr, pvalue_oqn, by.x = "gene2", by.y = "PROBE_ID")
  fdr <- merge(fdr, probe_pos, by.x = "gene", by.y = "Gene")
  fdr <- merge(fdr, maf, by.x = "SNP", by.y = "hg18_snpID")
  fdr <- merge(fdr, probe_info, by.x = "gene2", by.y = "PROBE_ID")
  
  
  
  # 調整
  fdr[,PROBE_COORDINATES := paste0(start,"-",end)]
  fdr[,c("start","end","gene","INFO","ER2","chr") := NULL]
  setnames(fdr, old = c("gene2", "CHROMOSOME"), new = c("Probe", "CHR"))
  
 
  setcolorder(fdr,c("CHR","ProbeID","Probe","Gene","PROBE_COORDINATES",
                    "SNP","rsID","R2","impute_type","MAF","REF","ALT","beta","t-stat",
                    "p-value","FDR","qvalue","sig_pval_Bonfi","pval_OQN","t_OQN",
                    "pval_OQN_BH","sig_OQN","sig_OQN_Bonfi",
                    "pval_QN","t_QN","pval_QN_BH","sig_QN","sig_QN_Bonfi"))
  
  fdr <- fdr[order(FDR),]
  
  # 調整數值
  fdr[,MAF := round(MAF, digits = 4)]
  fdr[,beta := round(beta, digits = 4)]
  fdr[,`t-stat` := round(`t-stat`, digits = 4)]
  fdr[,`p-value` := format(`p-value`, digits = 4, scientific = T)]
  fdr[,FDR := format(FDR, digits = 4, scientific = T)]
  fdr[,pval_OQN := format(pval_OQN, digits = 4, scientific = T)]
  fdr[,pval_QN := format(pval_QN, digits = 4, scientific = T)]
  fdr[,t_OQN := format(t_OQN, digits = 4, scientific = T)]
  fdr[,t_QN := format(t_QN, digits = 4, scientific = T)]
  fdr[,R2 := round(R2, digits = 4)]
  fdr[,qvalue := format(qvalue, digits = 4, scientific = T)]
  
  
  fwrite(fdr, outputname,
         row.names = F, col.names = T, sep = "\t")
}




#### 增加變數 ----
# 資料準備 ----
pvalue_oqn <- fread("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/exp_different.txt",header=T)

# 把各檔案cis snp number 資料刪掉，info 檔案用不到
pvalue_oqn <- pvalue_oqn[,c(paste0("maf_gt_", c("N","T"), "_pvalue_FDR_sigCis-SNP_number"),
                            paste0("maf_gt_", c("N","T"), "_MOSTpvalue_FDR_sigCis-SNP_number"),
                            paste0("maf_ds_", c("N","T"), "_pvalue_FDR_sigCis-SNP_number"),
                            paste0("maf_ds_", c("N","T"), "_MOSTpvalue_FDR_sigCis-SNP_number")) := NULL]
  
  
probe_pos <- fread("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/probe_pos.txt",header=T)
maf <- fread("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/cis_hg18_19_snp_MAF.txt", header=T)
maf[,hg19_snpID := NULL]

probe_info <- fread("D:/oral_cancer/expression/expression_data/probe24526infoB.txt", header=T)
probe_info[,c("TargetID","CHROMOSOME","PROBE_COORDINATES") := NULL]






# 跑 #####

file_name <- c("D:/oral_cancer/expression/trash/maf_gt_N_pvalue_noFDRfilter")

file_name_new <- c("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/maf_gt_N_pvalue_noFDRfilter")

sig_number <- fread("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/exp_different.txt",header=T)
sig_number <- sig_number %>% 
  select(c("PROBE_ID",
           "maf_gt_N_pvalue_FDR_sigCis-SNP_number",
           "maf_gt_T_pvalue_FDR_sigCis-SNP_number",
           "maf_gt_N_MOSTpvalue_FDR_sigCis-SNP_number",
           "maf_gt_T_MOSTpvalue_FDR_sigCis-SNP_number",
           "maf_ds_N_pvalue_FDR_sigCis-SNP_number",
           "maf_ds_T_pvalue_FDR_sigCis-SNP_number",
           "maf_ds_N_MOSTpvalue_FDR_sigCis-SNP_number",
           "maf_ds_T_MOSTpvalue_FDR_sigCis-SNP_number"))
  

# 對這些檔案生成info 檔
for (idx in 1:length(file_name)) {
  i <- file_name[idx]
  j <- file_name_new[idx]
  
  mix_info(paste0(i,".txt"),
           paste0(j,"_info.txt"))
  
  # 新增 cis sig snp
  a <- paste0(j,"_info.txt") %>% 
    fread()
  
  # 挑出1, idx+1 col
  b <- sig_number[, .SD, .SDcols = c(1, idx+1)]
  names(b) <- c("PROBE_ID", "sig_Cis-SNP_number")
  a <- merge(a,b,by.x="Probe",by.y="PROBE_ID")
  
  
  
  setcolorder(a,c("ProbeID","Probe","Gene","CHR","PROBE_COORDINATES","sig_Cis-SNP_number","SNP","rsID","impute_type","R2","MAF","REF","ALT","beta","t-stat","p-value","FDR","qvalue","sig_pval_Bonfi","t_OQN","pval_OQN","pval_OQN_BH","sig_OQN","sig_OQN_Bonfi","t_QN","pval_QN","pval_QN_BH","sig_QN","sig_QN_Bonfi"))
  
  fwrite(a, paste0(j,"_info.txt"),
         row.names = F, col.names = T, sep = "\t")
}




