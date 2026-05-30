# 目的
# 生成 chr1-22 的 eQTL result OQN_maf_gt_N_pvalue.txt, differentially express 的 OQN_exp_different.txt


## package ----
rm(list = ls())
gc()

library(data.table)
library(dplyr)
library(tidyverse)
library(stringr) 
library(ggplot2)
library(ggrepel)
library(matrixStats)  # rowMins
library(qvalue)



  # differentially express ----
  library(bestNormalize)
  exp <- fread("D:/oral_cancer/expression/expression_data/ExpRes11sv20120601Gene.txt", header = T)

  # N,T 共80 份 一起 OQN
  exp_OQN <- (apply(as.matrix(exp[, 3:82]), 2, function(x) {
    predict(orderNorm(x))
  })) %>%
    as.data.table()

  k <- exp_OQN[, 1:40]
  names(k) <- c(sprintf("0%dB", 1:9), sprintf("%dB", 10:40))
  k <- cbind(exp[, c(1:2)], k)
  fwrite(k, "C:/Peter/OQN_FIXpeople_before_eQTL/trash/OQN_ExpRes11sv20120601Gene_T.txt",
    row.names = F, col.names = T, sep = "\t"
  )

  k <- exp_OQN[, 41:80]
  names(k) <- c(sprintf("0%dB", 1:9), sprintf("%dB", 10:40))
  k <- cbind(exp[, c(1:2)], k)
  fwrite(k, "C:/Peter/OQN_FIXpeople_before_eQTL/trash/OQN_ExpRes11sv20120601Gene_N.txt",
    row.names = F, col.names = T, sep = "\t"
  )

  
  OQN_exp_N <- fread("C:/Peter/OQN_FIXpeople_before_eQTL/trash/OQN_ExpRes11sv20120601Gene_N.txt")
  OQN_exp_T <- fread("C:/Peter/OQN_FIXpeople_before_eQTL/trash/OQN_ExpRes11sv20120601Gene_T.txt")

  exp_diff_afOQN <- OQN_exp_N[, -c(1, 2)] - OQN_exp_T[, -c(1, 2)]
  n <- ncol(exp_diff_afOQN)
  dbar <- rowMeans(exp_diff_afOQN)
  sd_d <- apply(exp_diff_afOQN, 1, sd)
  tval <- dbar / (sd_d / sqrt(n))
  pval <- 2 * pt(-abs(tval), df = n - 1)

  pvalue_OQN <- cbind(OQN_exp_N[, 1:2], pval, tval)
  names(pvalue_OQN) <- c("Gene", "PROBE_ID", "pval_OQN", "t_OQN")
  probe <- fread("D:/oral_cancer/expression/expression_data/probe24526infoB.txt",
    header = T
  )[, .(chr = CHROMOSOME, PROBE_ID, PROBE_COORDINATES)]
  
  pvalue_OQN <- probe[pvalue_OQN, on=.(PROBE_ID)]
  pvalue_OQN <- pvalue_OQN[chr %in% (c(1:22) %>% as.character())]


  pvalue_OQN[, pval_OQN_BH := p.adjust(pval_OQN, method = "BH") %>%
    format(digits = 4, scientific = T) %>%
    as.numeric()]
  pvalue_OQN[, sig_OQN := ifelse(pval_OQN_BH < 0.05, 1, 0)]
  pvalue_OQN[, sig_OQN_Bonfi := ifelse(pval_OQN < 0.05 / nrow(pvalue_OQN), 1, 0)]

  # 調整數值
  pvalue_OQN[, pval_OQN := format(pval_OQN, digits = 4, scientific = T)]


  # differentially express pval 排序
  pvalue_OQN[, rank_pval := frank(pval_OQN, ties.method = "dense")]
  fwrite(pvalue_OQN, "C:/Peter/OQN_FIXpeople_before_eQTL/outcome/OQN_exp_different.txt",
    row.names = F, col.names = T, sep = "\t"
  )



### merge eQTL result ----

# build OQN_maf_gt_T_pvalue.txt
# GT 

# Normal part
file_paths <- sprintf("C:/Peter/OQN_FIXpeople_before_eQTL/outcome/OQN_maf_gt_N_cis_chr%d.txt", 1:22)

gt <- lapply(file_paths, fread, header = T) %>% 
    rbindlist(use.names = TRUE)



fwrite(gt, "C:/Peter/OQN_FIXpeople_before_eQTL/trash/OQN_maf_gt_N_pvalue.txt",
       row.names = F, col.names = T, sep = "\t")




# 同樣pvalue 大小，篩選距離最近的 
gene_pos <- fread("D:/oral_cancer/expression/expression_data/probe_pos_mulInterval.txt", header = T)


# 按基因、pvalue 排序
setkey(gt, gene, `p-value`)
gt <- merge(gt, gene_pos,by.x="gene",by.y="Gene")
gt[, snp_pos := str_extract(SNP , "(?<=\\:)\\d+") %>% as.numeric()]

# 算距離
gt[, dist := pmin(abs(snp_pos - start), abs(snp_pos - end))]
setkey(gt, gene, `p-value`, dist)
gt <- gt[, .SD[1L], by = "gene"]

# 對於 gene ILMN_2353642_1, ILMN_2353642_2，只留下一個 snp
gt[, gene2 := gsub("^([^_]+_[^_]+).*", "\\1",gene)]
setkey(gt, gene2, `p-value`, dist)
gt <- gt[, .SD[1L], by = "gene2"]
gt[,c("CHROMOSOME","start","end","snp_pos","dist","gene2") := NULL]

# 這挑出的最近 snp，可能是錯的，因為沒考慮落在區間內的，也許有 snp 不在 probe 範圍內，卻是最近的
fwrite(gt, "C:/Peter/OQN_FIXpeople_before_eQTL/trash/OQN_maf_gt_N_MOSTpvalue.txt",
       row.names = F, col.names = T, sep = "\t")
  



  
  
# Tumor part
  
file_paths <- sprintf("C:/Peter/OQN_FIXpeople_before_eQTL/outcome/OQN_maf_gt_T_cis_chr%d.txt", 1:22)

gt <- lapply(file_paths, fread, header = T) %>% 
    rbindlist(use.names = TRUE)


  
fwrite(gt, "C:/Peter/OQN_FIXpeople_before_eQTL/trash/OQN_maf_gt_T_pvalue.txt",
       row.names = F, col.names = T, sep = "\t")
  

  
# 同樣pvalue 大小，篩選距離最近的 
gene_pos <- fread("D:/oral_cancer/expression/expression_data/probe_pos_mulInterval.txt", header = T)

# 按基因、pvalue 排序
setkey(gt, gene, `p-value`)
gt <- merge(gt, gene_pos,by.x="gene",by.y="Gene")
gt[, snp_pos := str_extract(SNP , "(?<=\\:)\\d+") %>% as.numeric()]

# 算距離
gt[, dist := pmin(abs(snp_pos - start), abs(snp_pos - end))]
setkey(gt, gene, `p-value`, dist)
gt <- gt[, .SD[1L], by = "gene"]

# 對於 gene ILMN_2353642_1, ILMN_2353642_2，只留下一個 snp
gt[, gene2 := gsub("^([^_]+_[^_]+).*", "\\1",gene)]
setkey(gt, gene2, `p-value`, dist)
gt <- gt[, .SD[1L], by = "gene2"]
gt[,c("CHROMOSOME","start","end","snp_pos","dist","gene2") := NULL]

# 這挑出的最近 snp，可能是錯的，因為沒考慮落在區間內的，也許有 snp 不在 probe 範圍內，卻是最近的
fwrite(gt, "C:/Peter/OQN_FIXpeople_before_eQTL/trash/OQN_maf_gt_T_MOSTpvalue.txt",
       row.names = F, col.names = T, sep = "\t")





rm(gt)
gc()




## R2 filter cisSNP ----

# # *find R2*
# maf_snp <- fread("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/cis_hg18_19_snp_MAF.txt")

# # 用 list 儲存 22 個結果
# all_merged <- vector("list", 22)

# for (i in 1:22) {
  
#   snp_imputed <- sprintf("D:/oral_cancer/imputation_result/chr_%d/chr%d_snp_unique.dose.vcf", i, i) %>%
#     fread() 
  
#   setkey(snp_imputed, ID)
#   setkey(maf_snp, hg19_snpID)
#   names(snp_imputed)[ which(names(snp_imputed)== "#CHROM")] <-  "chr"
  
#   # nomatch = 0 代表只保留匹配的行 (即 Inner Join)
#   # 保留 snp_imputed 的所有欄位，list 不特別指出。要留的 maf_snp 欄位使用 i. 前綴
#   merged_data <- maf_snp[snp_imputed, nomatch = 0,on = c(hg19_snpID = "ID"),
#                               j = list(
#                                  chr = i.chr, 
#                                  hg18_snpID, hg19_snpID, rsID, MAF, REF, ALT,
#                                  INFO = i.INFO
#                              )]
  
  
#   merged_data[, R2 := fifelse(
#     grepl("R2=", INFO),
#     as.numeric(sub(".*?R2=([^;]+).*", "\\1", INFO)),
#     NA
#   )]
  
#   merged_data[, ER2 := fifelse(
#     grepl("ER2=", INFO),
#     as.numeric(sub(".*?ER2=([^;]+).*", "\\1", INFO)),
#     NA
#   )]
  
#   merged_data[, impute_type := 
#                 fifelse( grepl("TYPED", INFO), "TYPED",
#                 fifelse( grepl("IMPUTED", INFO), "IMPUTED",NA)
#   )]
  
#   setkey(merged_data, hg19_snpID)
#   all_merged[[i]] <- merged_data

# }

# final_merged <- rbindlist(all_merged, use.names = TRUE)

# fwrite(final_merged,
#        "D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/cis_hg18_19_snp_MAF.txt",
#        row.names = F, col.names = T, sep = "\t")


# # *挑出R2>0.8 snp*
# # 挑出原本有的snp或是impute r2>0.8 的
 
# cis_snp <- fread("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/cis_hg18_19_snp_MAF.txt",header=T)

# cis_snp_filter <- cis_snp %>% filter(R2>0.6)
# setkey(cis_snp_filter,hg19_snpID)
# fwrite(cis_snp_filter,
#        "D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/cis_snp_R2_0.6.txt",
#        row.names = F, col.names = T, sep = "\t")



# for (i in c(0.7,0.8,0.9)) {
  
#   cis_snp_filter <- cis_snp_filter %>% filter(R2>i)
#   setkey(cis_snp_filter,hg19_snpID)
#   fwrite(cis_snp_filter,
#          sprintf("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/cis_snp_R2_%s.txt",
#                  i),
#          row.names = F, col.names = T, sep = "\t")


# }








## VEP ----
# 只對 MAF>0.05 snp 跑rsID, REF, ALT

# *select SNP*




# maf_hg18_19 <- fread("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/cis_hg18_19_snp_MAF.txt",header=T)

# maf_hg18_19 <- maf_hg18_19 %>% 
#   filter(MAF>0.05)

# fwrite(data.table(chr = str_extract(maf_hg18_19$hg19_snpID, ".*(?=:)"),
#                   pos_start = str_extract(maf_hg18_19$hg19_snpID, "(?<=\\:)(\\d+)"),
#                   pos_end = str_extract(maf_hg18_19$hg19_snpID, "(?<=\\:)(\\d+)")
#                   ),
#        "D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/maf_cis_snpID_hg19.txt",
#        row.names = F, col.names = F, sep = "\t")




# ### 轉成rsID 
# # chr1-22 要將近12 hr
# # (不生成summary file 可能更快)


# # maf_cis_snpID_hg19.txt 檔案放 22   16050115   16050115，挑出 maf_cis_snp
# cmd <- paste0(
#   "wsl bash -c 'set -e; for i in {1..22}; do ",
#   "echo Processing chr${i}...; ",
#   "tabix -R /mnt/d/oral_cancer/expression/outcome/multiple_nucleotide_variant/maf_cis_snpID_hg19.txt /mnt/d/oral_cancer/expression/trash/chr${i}_cis_snp.vcf.gz > /mnt/d/oral_cancer/expression/trash/maf_chr${i}.vcf; ",
#   "done'"
# )
# system(cmd)


# # 對maf_cis_snp 轉成 rsID
# cmd <- paste0(
#   "wsl bash -c 'set -e; for i in {1..22}; do ",
#   "echo Running VEP on chr${i}...; ",
#   "/root/ensembl-vep/vep --cache --offline ",
#   "--dir_cache ~/.vep ",
#   "--assembly GRCh37 ",
#   "--fasta ~/.vep/homo_sapiens/115_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa ",
#   "--input_file /mnt/d/oral_cancer/expression/trash/maf_chr${i}.vcf ",
#   "--output_file /mnt/d/oral_cancer/expression/outcome/multiple_nucleotide_variant/chr${i}_vep.vcf ",
#   "--vcf --force_overwrite --check_existing ",
#   "--custom ~/.vep/homo_sapiens-chr${i}.vcf.gz,dbSNP,vcf,exact,0,ID; ",
#   "echo Done chr${i}; ",
#   "done'"
# )
# system(cmd)







# 挑出 rsID, chr:pos, REF, ALT

# 挑出rsID, chr:pos, REF, ALT in ubuntu cmd

# # set -e
# # for i in {1..22}; do
# #   echo "Running VEP on chr${i}..."
# #   awk -v OFS="\t" '!/^#/ {rs="NA"; if (match($8, /rs[0-9]+/, a)) rs=a[0]; print rs, $1 ":" $2, $4, $5}' \
# #   /mnt/d/oral_cancer/expression/outcome/multiple_nucleotide_variant/chr${i}_vep.vcf > /mnt/d/oral_cancer/expression/trash/chr${i}_vep.txt
# #   echo "Done chr${i}"
# # done



# # *pos rsID 對照表*

# # 合併
# file_paths <- sprintf("D:/oral_cancer/expression/trash/chr%d_vep.txt", 1:22)

# ds <- lapply(file_paths, fread, header = F) %>% 
#     rbindlist(use.names = F)

# names(ds) <- c("rsID","pos","REF","ALT")

# fwrite(ds, "D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/chr1-22_vep.txt",
#        row.names = F, col.names = T, sep = "\t")

# snp_pos <- fread("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/cis_hg18_19_snp_MAF.txt",header=T)
# rsID <- fread("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/chr1-22_vep.txt",header=T) 
# mix_rs_pos <- merge(snp_pos,rsID, by.x="hg19_snpID", by.y="pos", all.x = T)

# mix_rs_pos <- mix_rs_pos %>% 
#   unique()

# setcolorder(mix_rs_pos,c("hg18_snpID","hg19_snpID","rsID","MAF","REF","ALT"))
# fwrite(mix_rs_pos, "D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/cis_hg18_19_snp_MAF.txt",
#        row.names = F, col.names = T, sep = "\t")




# VEP check ----

# df <- fread("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/cis_hg18_19_snp_MAF.txt")


# # 去除重複
# rs_list <- unique(df$rsID)

# get_grch37_coord <- function(rsid) {
  
#   numeric_id <- gsub("rs", "", rsid)
#   url <- paste0("https://api.ncbi.nlm.nih.gov/variation/v0/beta/refsnp/", numeric_id)
  
#   resp <- GET(url)
#   if (resp$status_code != 200) {
#     return(data.frame(rsid=rsid, chr=NA, pos=NA))
#   }
  
#   json <- tryCatch(
#     fromJSON(content(resp, "text", encoding="UTF-8")),
#     error=function(e) return(NULL)
#   )
#   if (is.null(json)) {
#     return(data.frame(rsid=rsid, chr=NA, pos=NA))
#   }
  
#   placements <- json$primary_snapshot_data$placements_with_allele
#   if (length(placements) == 0) {
#     return(data.frame(rsid=rsid, chr=NA, pos=NA))
#   }
  
#   # ----------- 找 GRCh37 使用 seq_id_traits_by_assembly 
  
#   grch37_index <- NULL
  
#   for (i in 1:nrow(placements)) {
#     pa <- placements$placement_annot$seq_id_traits_by_assembly[[i]]
    
#     if (length(pa)!=0) {
#       assembly_names <- pa$assembly_name
#       if (assembly_names == "GRCh37.p13") {
#         grch37_index <- i
#         break
#       }
#     }
#   }
  
#   if (is.null(grch37_index)) {
#     return(data.frame(rsid=rsid, chr=NA, pos=NA))
#   }
  
#   # ----------- 抽取 GRCh37 座標 
  
#   target <- placements[grch37_index,]
  
#   # seq_id 例如 "NC_000020.10" → chr 20
#   chr_raw <- target$seq_id
#   chr <- gsub("NC_0+|\\..*$", "", chr_raw)  # 最穩定抽法
  
#   # 取第一個 allele 的位置（GRCh37）
#   allele <- target$alleles[[1]]
#   pos <- allele$allele$spdi$position + 1   # NCBI 是 0-based，要轉為 1-based
  
#   data.frame(
#     rsid = rsid,
#     chr = chr,
#     pos = pos
#   )
# }

# results <- list()
# for (i in rs_list) {
#   results[[i]] <- get_grch37_coord(i)
#   Sys.sleep(0.5)  # 每次呼叫間隔 0.5 秒 → 每秒最多 2 次
# }

# final_df <- bind_rows(results) %>% 
#   unique()
# rownames(final_df) <- NULL


# # 把hg19/GRCH37 id 併進資料
# final_df <- as.data.table(final_df)
# final_df[,SNP_GRCH37 := paste0(chr,":",pos)]
# final_df[,c("chr","pos") := NULL] 


# # 刪不是 ATCG
# maf <- fread("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/cis_hg18_19_snp_MAF.txt", header=T)
# maf[!maf$REF %in% c("A", "T", "C", "G") | !maf$ALT %in% c("A", "T", "C", "G"),
#      rsID := NA_character_]

# fwrite(maf,"D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/cis_hg18_19_snp_MAF.txt",
#          row.names = F, col.names = T, sep = "\t")







