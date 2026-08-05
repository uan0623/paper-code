## Setup ----

options(error = stop)

library(data.table)
library(dplyr)
library(tidyverse)
library(stringr) 
library(ggplot2)
library(ggrepel)
library(matrixStats)  # rowMins
library(qvalue)




  # 定義資料夾名稱 ####
  base_dir <- sprintf(
    "C:/Peter/repeatSNP_clumping_FIXpeople/r2_filter_%s/eQTL_%s",
    r2_threshold, eqtl_threshold
  )

  # 刪舊資料夾
  if (dir.exists(base_dir)) {
  message("Deleting: ", base_dir)
  unlink(base_dir, recursive = TRUE, force = TRUE)
} else {
  message("Directory does not exist: ", base_dir)
}


  clump_path <- function(pop, section, ...) {
    if (!pop %in% c("EUR", "FIN")) {
      stop("pop must be 'EUR' or 'FIN'")
    }
    if (!section %in% c("outcome", "trash")) {
      stop("section must be 'outcome' or 'trash'")
    }

    file.path(base_dir, pop, section, ...)
  }
  pops <- c("EUR", "FIN")
  LD <- c(0.1, 0.2, 0.3)
  kb_ranges <- c(500, 1000)
  all_paths <- c()


  for (pop in pops) {
    trash_path <- file.path(pop, "trash/tmp_chr")
    all_paths <- c(all_paths, trash_path)

    for (para in LD) {
      for (range_kb in kb_ranges) {
        path <- file.path(pop, "outcome", sprintf("%dkb_LD%g", range_kb, para))
        all_paths <- c(all_paths, path)
      }
    }
  }

  all_paths <- file.path(base_dir, all_paths)


# 沒有資料夾，就自動新增 ####
for (i in all_paths) {
  if(!dir.exists(i)) dir.create(i, recursive = TRUE)
  
}




# 1000G EUR/FIN

for (pop in pops) {


message("==========================================")
message(sprintf("1000G %s", pop))
message("==========================================")



## SNP Select ----
if (eqtl_threshold=="bon") {
   dt <- fread(sprintf(
     "C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_%s/outcome/OQN_MixFinngenPval_7_%s_bon.txt",
     r2_threshold, pop
   ))

} else if (eqtl_threshold=="FDR") {
   dt <- fread(sprintf(
     "C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_%s/outcome/OQN_MixFinngenPval_7_%s.txt",
     r2_threshold, pop
   ))
} 

dt[, MOST_snp_nearest_alt := tstrsplit(MOST_snp_nearest_hg38, ":", keep = 4)[[1]] ]
dt[, MOST_snp_nearest_ref := tstrsplit(MOST_snp_nearest_hg38, ":", keep = 3)[[1]] ]
dt[, MOST_snp_nearest_chr := tstrsplit(MOST_snp_nearest_hg38, ":", keep = 1)[[1]] ]


# eQTL 出現在fin gwas的 snp，過 FDR, qval 的 snp info ###

gt_N_bon <- dt[gt_N_hg38_finngen_Bonfi  ==1, .(CHR, SNP_hg19=gt_N_hg19, alt, ref)]
gt_N_fdr <- dt[gt_N_hg38_finngen_FDR < 0.05, .(CHR, SNP_hg19=gt_N_hg19, alt, ref)]
gt_N_qval <- dt[gt_N_hg38_finngen_qvalue < 0.05, .(CHR, SNP_hg19=gt_N_hg19, alt, ref)]

if(nrow(gt_N_bon)!=0){
  fwrite(gt_N_bon,
         clump_path(pop, "outcome", "1_bon.txt"),
         row.names = F, col.names = T, sep = "\t") 
  
} else(message("gt_N_bon no snp"))

if(nrow(gt_N_fdr)!=0){
  fwrite(gt_N_fdr,
         clump_path(pop, "outcome", "1_fdr.txt"),
         row.names = F, col.names = T, sep = "\t") 
  
} else(message("gt_N_fdr no snp"))

if(nrow(gt_N_qval)!=0){
  fwrite(gt_N_qval,
         clump_path(pop, "outcome", "1_qval.txt"),
         row.names = F, col.names = T, sep = "\t") 
  
} else(message("gt_N_qval no snp"))










# eQTL 沒出現在fin gwas的，用 LD snp 替代 
ld_fdr <- dt[mostLD_or_finngen_FDR < 0.05, .(CHR, gt_N_hg19, alt, ref,mostLD_snp_nearest_hg19)]
if(nrow(ld_fdr)!=0){
  ld_fdr[mostLD_snp_nearest_hg19=="",
         mostLD_snp_nearest_hg19:= paste0(gt_N_hg19,":",ref,":",alt)]
  ld_fdr[,LDSNP_replace := ifelse(alt=="",1,0)]
  
  
  # SNP_hg19 表示 LD SNP 或是 gt_N snp
  setnames(ld_fdr,old = "mostLD_snp_nearest_hg19",new = "SNP_hg19")
  ld_fdr[, `:=`(
    SNP_alt = alt,
    SNP_ref = ref,
    SNP_hg19_new = gt_N_hg19
  )]
  ld_fdr[!is.na(SNP_hg19), `:=`(
    SNP_alt = tstrsplit(SNP_hg19, ":", keep = 4)[[1]],
    SNP_ref = tstrsplit(SNP_hg19, ":", keep = 3)[[1]],
    SNP_hg19_new = sub(":[^:]+:[^:]+$", "", SNP_hg19)
  )]
  ld_fdr[, SNP_hg19 := SNP_hg19_new]
  ld_fdr[, SNP_hg19_new := NULL]
  
  # 統一欄位名稱，並存檔
  fwrite(  ld_fdr[, .(CHR, SNP_hg19, alt=SNP_alt, ref=SNP_ref)],
           clump_path(pop, "outcome", "2_fdr.txt"),
           row.names = F, col.names = T, sep = "\t") 
}else(message("ld_fdr no snp"))


ld_qval <- dt[mostLD_or_finngen_qvalue < 0.05, .(CHR, gt_N_hg19, alt, ref,mostLD_snp_nearest_hg19)]
if(nrow(ld_qval)!=0){
  ld_qval[mostLD_snp_nearest_hg19=="",
          mostLD_snp_nearest_hg19:= paste0(gt_N_hg19,":",ref,":",alt)]
  ld_qval[,LDSNP_replace := ifelse(alt=="",1,0)]
  
  setnames(ld_qval,old = "mostLD_snp_nearest_hg19",new = "SNP_hg19")
  ld_qval[, `:=`(
    SNP_alt = alt,
    SNP_ref = ref,
    SNP_hg19_new = gt_N_hg19
  )]
  ld_qval[!is.na(SNP_hg19), `:=`(
    SNP_alt = tstrsplit(SNP_hg19, ":", keep = 4)[[1]],
    SNP_ref = tstrsplit(SNP_hg19, ":", keep = 3)[[1]],
    SNP_hg19_new = sub(":[^:]+:[^:]+$", "", SNP_hg19)
  )]
  ld_qval[, SNP_hg19 := SNP_hg19_new]
  ld_qval[, SNP_hg19_new := NULL]
  
  fwrite(  ld_qval[, .(CHR, SNP_hg19, alt=SNP_alt, ref=SNP_ref)],
           clump_path(pop, "outcome", "2_qval.txt"),
           row.names = F, col.names = T, sep = "\t") 
  
}else(message("ld_qval no snp"))


ld_bon <- dt[mostLD_or_finngen_Bonfi ==1, .(CHR, gt_N_hg19, alt, ref,mostLD_snp_nearest_hg19)] 
if(nrow(ld_bon)!=0){
  ld_bon[mostLD_snp_nearest_hg19=="",
         mostLD_snp_nearest_hg19:= paste0(gt_N_hg19,":",ref,":",alt)]
  ld_bon[,LDSNP_replace := ifelse(alt=="",1,0)]
  
  setnames(ld_bon,old = "mostLD_snp_nearest_hg19",new = "SNP_hg19")
  ld_bon[, `:=`(
    SNP_alt = alt,
    SNP_ref = ref,
    SNP_hg19_new = gt_N_hg19
  )]
  ld_bon[!is.na(SNP_hg19), `:=`(
    SNP_alt = tstrsplit(SNP_hg19, ":", keep = 4)[[1]],
    SNP_ref = tstrsplit(SNP_hg19, ":", keep = 3)[[1]],
    SNP_hg19_new = sub(":[^:]+:[^:]+$", "", SNP_hg19)
  )]
  ld_bon[, SNP_hg19 := SNP_hg19_new]
  ld_bon[, SNP_hg19_new := NULL]

  fwrite(  ld_bon[, .(CHR, SNP_hg19, alt=SNP_alt, ref=SNP_ref)],
           clump_path(pop, "outcome", "2_bon.txt"),
           row.names = F, col.names = T, sep = "\t") 
  
}else(message("ld_bon no snp"))







# 13516 eQTL snp 找 finngen pval 最顯著的 snp info 

most_p <- dt[,.SD[1],by="MOST_snp_nearest_hg19"]

# 統一欄位名稱，並存檔
fwrite( most_p[MOST_snp_hg38_finngen_FDR<0.05,
               .(CHR = MOST_snp_nearest_chr,
                 SNP_hg19 = MOST_snp_nearest_hg19,
                 alt = MOST_snp_nearest_alt,
                 ref = MOST_snp_nearest_ref)],
        clump_path(pop, "outcome", "3_fdr.txt"),
        row.names = F, col.names = T, sep = "\t") 

fwrite( most_p[MOST_snp_hg38_finngen_Bonfi==1,
               .(CHR = MOST_snp_nearest_chr,
                 SNP_hg19 = MOST_snp_nearest_hg19,
                 alt = MOST_snp_nearest_alt,
                 ref = MOST_snp_nearest_ref)],
        clump_path(pop, "outcome", "3_bon.txt"),
        row.names = F, col.names = T, sep = "\t") 






## SNP Save ----

prefix <- paste0(clump_path(pop, "outcome"), "/")
files <- c("3_bon.txt", "3_fdr.txt","2_fdr.txt", "2_qval.txt", "1_qval.txt", "1_fdr.txt") 
file_name <- paste0(prefix, files)

files <- c("3_bon_SNP.txt", "3_fdr_SNP.txt","2_fdr_SNP.txt", "2_qval_SNP.txt", "1_qval_SNP.txt", "1_fdr_SNP.txt") 
prefix <- paste0(clump_path(pop, "trash"), "/")
file_output <- paste0(prefix, files)

# 判斷 3 手法檔案是否存在，不在就別跑，避免報錯
existing_idx <- file.exists(file_name)
if (any(!existing_idx)) {
  message("Skip missing SNP list files: \n",
          paste(file_name[!existing_idx], collapse = "\n"))
}
file_output <- file_output[existing_idx]
file_name <- file_name[existing_idx]



for (i in 1:length(file_name)) {
  a <- fread(file_name[i])
  
  fwrite( data.table(kk=a$SNP_hg19 %>% unique()),
          file_output[i],
          row.names = F, col.names = F, sep = "\n") 
  
  
}


## PLINK Run ----

prefix <- paste0(clump_path(pop, "trash"), "/")
files <- c("3_bon_SNP.txt", "3_fdr_SNP.txt","2_fdr_SNP.txt", "2_qval_SNP.txt", "1_qval_SNP.txt", "1_fdr_SNP.txt") 
file_output <- paste0(prefix, files)




# 定義目錄（建議統一管理）
tmp_dir    <- clump_path(pop, "trash", "tmp_chr")    # 存放 1-22 染色體的小檔案

final_dir  <- clump_path(pop, "outcome") # 存放合併後的結果
bfile_path <- if (pop == "EUR") {
  "C:/Peter/PCA_1000G_20130502/trash/deal_repeatSNP"
} else {
  "C:/Peter/PCA_1000G_20130502/FIN_sample/trash/deal_repeatSNP"
} # 原始數據位置
# 取得乾淨的名稱 (例如: 1_qval)
new_names <- sub("_SNP\\.txt$", "", basename(file_output))


# 判斷 3 手法檔案是否存在，不在就別跑，避免報錯
existing_idx <- file.exists(file_output)
if (any(!existing_idx)) {
  message("Skip missing SNP list files: \n",
          paste(file_output[!existing_idx], collapse = "\n"))
}
file_output <- file_output[existing_idx]
new_names <- new_names[existing_idx]


# 對6 檔案做
for (j in seq_along(file_output)) {
  
  current_prefix <- new_names[j]
  current_list   <- file_output[j]
  
  message("==========================================")
  message(sprintf("🚀 開始處理清單 [%d/%d]: %s", j, length(file_output), current_prefix))
  message("==========================================")
  
  
  extract_snp <- fread(current_list,header=F)
  names(extract_snp) <- "snp"
  extract_snp[, chr:= str_extract(snp, ".*(?=:)") %>% as.numeric()]
  extract_chr <- extract_snp$chr %>% 
    unique() %>% 
    sort()
  
  
  
  # 從 22 染色體 依序取出 snp
  for (i in extract_chr) {
    # 構建輸入與輸出路徑
    input_bfile <- file.path(bfile_path, sprintf("chr_%d_rename", i))
    output_tmp  <- file.path(tmp_dir, sprintf("%s_chr%d", current_prefix, i))
    
    # 構建 PLINK 指令
    cmd_extract <- sprintf('plink --bfile "%s" --extract "%s" --make-bed --out "%s"',
                           input_bfile, current_list, output_tmp)
    
    # 執行（這裡使用 shell() 或 system() 皆可，Windows 建議 shell()）
    message(sprintf("  -> 正在提取 Chr %d...", i))
    # intern = TRUE 過程不顯示，存成變數
    plink_extract <- suppressWarnings(shell(paste(cmd_extract, "2>&1"), intern = TRUE))
    plink_status <- attr(plink_extract, "status")

    if (
      any(grepl("No variants remaining after --extract\\.", plink_extract)) ||
        identical(plink_status, 12L)
    ) {
      message("No variants remaining after --extract.")
      next
    }
    if (!is.null(plink_status) && plink_status != 0) {
      stop("PLINK --extract failed with status ", plink_status, ":\n",
           paste(plink_extract, collapse = "\n"))
    }
  }
  
  
  
  # 1. 搜尋 tmp_dir 下所有符合該 prefix 且後綴為 .bed 的檔案
  # pattern 範例: ^3_bon_chr[0-9]+\.bed$
  all_bed_files <- list.files(path = tmp_dir, 
                              pattern = sprintf("^%s_chr[0-9]+\\.bed$", current_prefix), 
                              full.names = TRUE)
  
  # 2. 將 .bed 副檔名去掉，因為 PLINK --merge-list 需要的是檔案前綴
  all_prefixes <- sub("\\.bed$", "", all_bed_files)
  
  if (length(all_prefixes) < 1) {
    message(sprintf("⚠️ 找不到 %s 的暫存檔，跳過合併。", current_prefix))
    next
  }
  
  # 3. 準備合併：取出第一個檔案作為主檔案，其餘放進 merge-list
  first_chr <- all_prefixes[1]
  remaining_chrs <- all_prefixes[-1]
  
  merge_list_path <- file.path(tmp_dir, paste0(current_prefix, "_mergelist.txt"))
  writeLines(remaining_chrs, merge_list_path)
  
  # 4. 執行合併 (如果只有一個檔案則不需要 --merge-list)
  final_out <- file.path(final_dir, sprintf("%s_chr1-22", current_prefix))
  message(sprintf(" 📦 偵測到 %d 個檔案，正在合併 %s ...", length(all_prefixes), current_prefix))
  
  if (length(all_prefixes) > 1) {
    cmd_merge <- sprintf('plink --bfile "%s" --merge-list "%s" --make-bed --out "%s"',
                         first_chr, merge_list_path, final_out)
  } else {
    # 只有一個檔案時，直接用 --make-bed 複製過去即可
    cmd_merge <- sprintf('plink --bfile "%s" --make-bed --out "%s"',
                         first_chr, final_out)
  }
  
  shell(cmd_merge,intern = TRUE)
  message(sprintf("✅ 清單 %s 處理完成！", current_prefix))
}
  









## Summary Stat ----




prefix <- paste0(clump_path(pop, "outcome"), "/")
files <- c("3_bon.txt", "3_fdr.txt","2_fdr.txt", "2_qval.txt", "1_qval.txt", "1_fdr.txt") 
file_name <- paste0(prefix, files)

# 判斷 3 手法檔案是否存在，不在就別跑，避免報錯
existing_idx <- file.exists(file_name)
if (any(!existing_idx)) {
  message("Skip summary stat for missing SNP list files: \n",
          paste(file_name[!existing_idx], collapse = "\n"))
}
file_name <- file_name[existing_idx]



for (i in file_name) {
  
  base <- basename(i)
  prefix <- sub("\\.txt$", "", base) # 取得 1_qval, 1_fdr 等前綴
  message(paste("--- 正在處理檔案:", base, "---"))
  
  # --- 步驟 1: 根據檔名設定變動參數 ---
  
  # 預設欄位名稱 
  col_map <- list(chr = "CHR", snp = "SNP_hg19", ref = "ref", alt = "alt")
  bim_suffix <- "_chr1-22.bim"
  
  
  
  # get("SNP_ref") 就是取出變數 SNP_ref 的值，
  if (!file.exists(i)) {
    message("Skip missing input file: ", i)
    next
  }
  snp <- fread(i)
  snp <- snp[, .(chr = get(col_map$chr), hg19_snpID = get(col_map$snp), 
                 ref = get(col_map$ref), alt = get(col_map$alt))]
  
  # 讀取 BIM
  bim_path <- clump_path(pop, "outcome", paste0(prefix, bim_suffix))
  if (!file.exists(i)) {
    message("Skip missing input file: ", i)
    next
  }
  if (!file.exists(bim_path)) {
    message("Skip missing BIM file: ", bim_path)
    next
  }
  bim <- fread(bim_path, col.names = c("chr", "hg19_snpID", "genetic_dist", "pos", "minor", "major"))
  
  # 合併與過濾 (假設 finngen 已經存在於環境中)
  final_data <- bim[finngen, on = .(chr, hg19_snpID), nomatch = 0]
  
  # 邏輯過濾：Allele 匹配 & 取 P-value 最小
  final_data <- final_data[
    (minor == ref & major == alt) | (minor == alt & major == ref)
  ][order(P), .SD[1], by = hg19_snpID]
  
  setnames(final_data, "hg19_snpID", "SNP")
  setorder(final_data, chr, pos) 
  
  # --- 步驟 3: 存檔 ---
  output_path <- clump_path(pop, "outcome", paste0(prefix, "_summarySTAT.txt"))
  fwrite(final_data, output_path, row.names = FALSE, sep = "\t")
  
  message(paste("✅ 完成存檔:", output_path))
}



## Clumping ----

prefix <- paste0(clump_path(pop, "outcome"), "/")
files <- c("3_bon.txt", "3_fdr.txt","2_fdr.txt", "2_qval.txt", "1_qval.txt", "1_fdr.txt") 
file_name <- paste0(prefix, files)
trash_prefix <- paste0(clump_path(pop, "trash"), "/")

# 判斷 3 手法檔案是否存在，不在就別跑，避免報錯
existing_idx <- file.exists(file_name)
if (any(!existing_idx)) {
  message("Skip summary stat for missing SNP list files: \n",
          paste(file_name[!existing_idx], collapse = "\n"))
}
file_name <- file_name[existing_idx]



for (i in file_name) {
  base <- basename(i)
  prefix <- sub("\\.txt$", "", base) # 取得 1_qval, 1_fdr 等前綴
  message(paste("--- 正在處理檔案:", base, "---"))
  
  
  # 預設欄位名稱
  col_map <- list(chr = "CHR", snp = "SNP_hg19", ref = "ref", alt = "alt")
  bim_suffix <- "_chr1-22" 
  
  
  
  # 1000, 500 kb LD=0.1, 0.2, 0.3 clumping
  for (para in c(0.1, 0.2, 0.3)) {
    # 設 binary file, summary stat 路徑
    binary_path <- clump_path(pop, "outcome", paste0(prefix, bim_suffix))
    summary_path <- clump_path(pop, "outcome", paste0(prefix, "_summarySTAT.txt"))
    if (!file.exists(paste0(binary_path, ".bed")) ||
        !file.exists(paste0(binary_path, ".bim")) ||
        !file.exists(paste0(binary_path, ".fam"))) {
      message("Skip missing PLINK binary set: ", binary_path)
      next
    }
    if (!file.exists(summary_path)) {
      message("Skip missing summary stat file: ", summary_path)
      next
    }
    
    
    
    # 1000kb 
    final_out <- clump_path(pop, "outcome", sprintf("1000kb_LD%s", para), prefix)
    cmd_extract <- sprintf('plink --bfile "%s" --clump "%s" --clump-p1 0.99999 --clump-r2 %s --clump-kb 1000 --clump-snp-field SNP --clump-field P --out "%s"',
                           binary_path, summary_path, para, final_out)
    shell(cmd_extract,intern = TRUE)
    
    
    
    # 500kb 
    final_out <- clump_path(pop, "outcome", sprintf("500kb_LD%s", para), prefix)
    cmd_extract <- sprintf('plink --bfile "%s" --clump "%s" --clump-p1 0.99999 --clump-r2 %s --clump-kb 500 --clump-snp-field SNP --clump-field P --out "%s"',
                           binary_path, summary_path, para, final_out)
    shell(cmd_extract,intern = TRUE)
  }
  
  message(sprintf("  ->  %s... 1000kb, 500kb clumping 結束", base))
}


}



