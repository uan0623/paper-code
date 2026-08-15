# TITLE: TWAS Workflow
# SUBTITLE: Prepare genotype and summary statistics for FUSION TWAS analysis.
# SEARCH TAGS: TWAS, FUSION, GWAS, genotype, summary statistics
# NOTE: Code below is unchanged; only this navigation header was added.


# - fusion 官網 <http://gusevlab.org/projects/fusion/>
# - 教學 <https://cloufield.github.io/GWASTutorial/21_twas/#installation>

# - 使用的 individual data： 40 OC patient genotype data
# "D:/oral_cancer/imputation_result/trash/chr1-22_imputation" 複製到 "C:/TWAS/chr1-22_imputation"

# - 使用的 summary stat： finngen OC GWAS
# "D:/oral_cancer/TWAS/summary_stat/FinnGen/C3_ORALCAVITY_EXALLC"


# - top1,lasso,enet 3種模型算權重檔案：
# "/root/fusion_project/Normal_part/weight_directory"


# - TWAS 執行結果檔案：
# /root/fusion_project/Normal_part/twas/C3_ORAL_TWAS_Result_chr${chr}.dat
# /root/fusion_project/Normal_part/twas/C3_ORAL_TWAS_Result_chr1-22.dat

# 欄位意思：
# <http://gusevlab.org/projects/fusion/#compute-your-own-predictive-models>
# 的 Output: Gene-disease association

# package ----
library(data.table)
library(dplyr)
library(stringr)
library(qvalue)
install.packages(c("devtools", "optparse", "glmnet", "Rcpp", "Matrix", "methods"))

# # 複製別人github 的 plink2R (under windows cmd)
# git clone https://github.com/gabraham/plink2R.git

# # 從下載位址安裝
# devtools::install_local("D:/fusion_TWAS/plink2R-master/plink2R")
# library(devtools)
# install.packages("cli", type = "source")


# # funcion ----
# run_cmd <- function(cmd) {

#   exit_code <- system(cmd)

#   # 發生錯誤，system 回傳非 0 的數字
#   if (exit_code != 0) {
#     stop(paste("cmd 執行失敗"))
#   }
# }


# ## 安裝linux plink, gcta64 ----

# # 1. 安裝linux plink, gcta64
# # 2. 執行

# # 以下code 在 ubuntu cmd 輸入 ----

# # 建立資料夾 my_bin
# mkdir -p /root/my_bin

# # 目錄切換到 my_bin
# cd /root/my_bin

# # https://github.com/jianyangqt/gcta/releases/download/v1.94.1/gcta-1.94.1-linux-x86_64-static 下載 linux 版本的 gcta64 至 my_bin
# wget https://github.com/jianyangqt/gcta/releases/download/v1.94.1/gcta-1.94.1-linux-x86_64-static

# # gcta-1.94.1-linux-x86_64-static 改名成 gcta64
# mv gcta-1.94.1-linux-x86_64-static gcta64

# # https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20250819.zip 下載 linux 版本的 plink 至 my_bin
# wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20250819.zip

# # 解壓縮
# unzip plink_linux_x86_64_20250819.zip

# # 只留plink2，其他手動刪掉

# # 讓檔案變成可執行
# chmod +x /root/my_bin/gcta64
# chmod +x /root/my_bin/plink

# # 加進 path
# echo 'export PATH="/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/root/my_bin:$PATH"' >> /root/.bashrc
# source /root/.bashrc


# # 建立一個放程式的資料夾和放資料的資料夾
# mkdir /root/fusion_project/Normal_part

# # 進入資料夾, 之後的所有動作都在這裡面進行
# cd /root/fusion_project/Normal_part

# # 移動檔案
# mv /root/my_bin/fusion_twas-master/FUSION.assoc_test.R /root/fusion_project
# mv /root/my_bin/fusion_twas-master/FUSION.compute_weights.R /root/fusion_project
















# 算 weight 前置作業
# 算 weight 需要準備每個 probe 的 cis snp binary file, expression file，gcta 會根據cis snp binary file 對每個 snp 估計權重

# OUTLINE: Build probe position file ----
## 弄 probe pos file ----

# probe pos 都轉成一個區間，cis snp 範圍1000KB 之下，我們的probe 單區間跟多區間結果一樣
probe <- fread("D:/oral_cancer/expression/expression_data/probe24526infoB.txt", header = T)

# 判斷"^\\d+$"是否被包含在 CHROMOSOME 的每個元素中，"^\\d+$" 代表從字串開頭到結尾，只能有數字
probe <- probe[grepl("^\\d+$", CHROMOSOME)]

probe$CHROMOSOME <- as.numeric(probe$CHROMOSOME)


# by = list(TargetID, CHROMOSOME) 代表以TargetID, CHROMOSOME 的值將所有row 做分組
# 新增range 變數，把多重區間拆成一個個字串，像是 1-7:10-12 -> "1-7" "10-12"
probe_bounds <- probe[, list(
  range = str_extract_all(PROBE_COORDINATES, "\\d+-\\d+")
), by = list(TargetID, CHROMOSOME)]


# 把list "1-7" "10-12" 變成字串
probe_bounds <- probe_bounds[, list(
  start = as.integer(sub("-.*", "", unlist(range))),
  end   = as.integer(sub(".*-", "", unlist(range)))
), by = list(TargetID, CHROMOSOME)]


# 新增PROBE_ID, ProbeID. str_count 找出每個row "-" 出現次數，代表區間拆分後，這 ProbeID 要出現幾次
probe_bounds[, PROBE_ID := rep(
  probe$PROBE_ID,
  str_count(probe$PROBE_COORDINATES, "-")
)]
probe_bounds[]
probe_bounds[duplicated(probe_bounds, by = c("TargetID", "PROBE_ID"))]

# 調整多區間->一個大區間，因為多區間之間間隔不超過1000KB
setkey(probe_bounds, CHROMOSOME, start)
a <- duplicated(probe_bounds, by = c("TargetID", "PROBE_ID")) %>% which()
probe_bounds[a - 1, "end"] <- probe_bounds[a, "end"]

probe_bounds <- probe_bounds[!duplicated(probe_bounds, by = c("TargetID", "PROBE_ID"))]

fwrite(probe_bounds %>% select(PROBE_ID, CHROMOSOME, start, end),
  "D:/oral_cancer/expression/expression_data/probe24526infoB_OneInterval.txt",
  row.names = F, col.names = T, sep = "\t"
)













# OUTLINE: Build and split probe SNP table ----
## 造probe snp table 且拆分 ----
# "D:/oral_cancer/TWAS/trash/cis_hg19_snp_R2filter.txt" 檔案紀錄那些 cis snp，typed 或是 imputed 且 R2>0.8 的 snp_hg19，來自檔案
# "D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/cis_hg18_19_snp_R2filter.txt"

gene <- fread("D:/oral_cancer/expression/expression_data/probe24526infoB_OneInterval.txt")
snp <- fread("D:/oral_cancer/TWAS/trash/cis_hg19_snp_R2filter.txt")
snp[, start := hg19_pos]
snp[, end := hg19_pos]

# 基因位置延伸 1000 kb
gene[, gene_start := start - 1e6]
gene[, gene_end := end + 1e6]

# 若有負數位置，拉回 1
gene[gene_start < 1, gene_start := 1]

# 設定 key（必須做，否則不會快）
setkey(snp, chr, start, end)
setkey(gene, CHROMOSOME, gene_start, gene_end)

# 用 data.table 的 foverlaps 快速找區間重疊

overlap <- foverlaps(
  snp,
  gene,
  by.x = c("chr", "start", "end"),
  by.y = c("CHROMOSOME", "gene_start", "gene_end"),
  type = "within",
  nomatch = 0L
)

# 輸出 GENE_ID 與 SNP_ID

fwrite(overlap[, .(PROBE_ID, hg19_snpID)],
  "D:/oral_cancer/TWAS/trash/probe_to_snp.txt",
  row.names = F, col.names = T, sep = "\t"
)


#  拆成一個個probe 檔案，每個probe 檔案紀錄對應的cis snp ----
cis <- fread("D:/oral_cancer/TWAS/trash/probe_to_snp.txt")

# 為每個 gene 建 snplist 檔
dir.create("C:/TWAS/trash/snplists")

# 生每個probe cis snp 檔案
cis[, fwrite(.SD,
  file = paste0("C:/TWAS/trash/snplists/cisSNP_hg19_R2filter_", PROBE_ID, ".txt"),
  row.names = F, col.names = F, sep = "\n"
),
by = PROBE_ID
]













# OUTLINE: Export expression by probe ----
## Collect exp by probe ----
# 造一個個probe 檔案對應的expression ----
exp <- fread("D:/oral_cancer/expression/expression_data/ExpRes11sv20120601Gene_N.txt", header = TRUE)

# 取得所有樣本 ID 作為 pheno 檔案的FID/IID
sample_ids <- names(exp)[-c(1:2)]

for (i in 1:nrow(exp)) {
  pheno_df <- data.table(
    FID = sample_ids, # 樣本家族 ID
    IID = sample_ids, # 樣本個體 ID
    phe = exp[i, -c(1, 2)] %>% as.vector() # 表現量數值
  )

  fwrite(pheno_df,
    paste0("C:/TWAS/trash/expression_by_probe/expression_for_", exp$PROBE_ID[i], ".txt"),
    row.names = FALSE,
    col.names = TRUE,
    sep = "\t"
  )
}


# 造probe list file ----
fwrite(data.table(a = exp$PROBE_ID),
  paste0("C:/TWAS/trash/probe_list.txt"),
  row.names = F,
  col.names = F,
  sep = "\n"
)

# 移至 ubuntu (code under ubuntu)
cp - r / mnt / c / TWAS / trash / expression_by_probe / root / fusion_project / Normal_part /
  cp / mnt / c / TWAS / trash / probe_list.txt / root / fusion_project / Normal_part /












  # OUTLINE: Export cis SNPs by probe ----
  ## Collect cis-snp by probe ----
  # 用plink extract 從all cis snp binary file 抽取一個個probe 的 cis snp，腳本執行花費1 hr
  # 造一個個probe 檔案對應的cis-snp ----
  # 　D:/oral_cancer/imputation_result/trash/chr1-22_imputation 先複製到c槽，執行較快

  cis_snp_rsq <- fread("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/cis_hg18_19_snp_R2filter.txt", header = TRUE)
setkey(cis_snp_rsq, hg18_snpID)
cis_snp_rsq[, hg19_pos := str_extract(hg19_snpID, "(?<=\\:)\\d+")]

fwrite(cis_snp_rsq %>% select(hg19_snpID, chr, hg19_pos), "D:/oral_cancer/TWAS/trash/cis_hg19_snp_R2filter.txt",
  row.names = F, col.names = T, sep = "\t"
)

system("plink --bfile C:/TWAS/chr1-22_imputation --extract D:/oral_cancer/TWAS/trash/cis_hg19_snp_R2filter.txt --make-bed --out C:/TWAS/cis_hg19_snp_R2filter")















# # code under ubuntu
# # 複製檔案至ubuntu
# cp /mnt/c/TWAS/cis_hg19_snp_R2filter.* /root/fusion_project/
# cp -r /mnt/c/TWAS/trash/snplists /root/fusion_project/Normal_part/

# # 安裝平行運算
# sudo apt update
# sudo apt install parallel

# # 建立檔案
# nano /root/fusion_project/run_parallel_extract.sh

## 腳本 ----

# SNPLIST_DIR="/root/fusion_project/Normal_part/snplists"
# OUT_DIR="/root/fusion_project/Normal_part/cisSNP_binary_file"
# BFILE="/root/fusion_project/cis_hg19_snp_R2filter"

# mkdir -p "$OUT_DIR"

# # /cisSNP_hg19_R2filter_ILMN_*.txt 當作單一元素放入陣列，total 代表陣列長度
# file_list=("$SNPLIST_DIR"/cisSNP_hg19_R2filter_ILMN_*.txt)
# total=${#file_list[@]}
# count=0

# export OUT_DIR BFILE total

# # 定義每個 probe 要跑的 function
# run_one() {
#   # 接收函數的第一個參數，也就是完整路徑，像是/root/.../cisSNP_hg19_R2filter_ILMN_1651228.txt
#     file="$1"

#     # 去掉路徑,只留檔名
#     fname=$(basename "$file")
#     # 去掉 cisSNP_hg19_R2filter_ 與尾巴 .txt
#     id=$(echo "$fname" | sed 's/cisSNP_hg19_R2filter_//; s/\.txt//')
#     out_prefix="$OUT_DIR/Cissnp_in_${id}"

#     # 執行 PLINK
#     plink \
#       --bfile "$BFILE" \
#       --extract "$file" \
#       --make-bed \
#       --out "$out_prefix"

#     # 更新進度
#     (
#         flock -x 200
#         count_file="$OUT_DIR/.progress_count"
#         if [[ ! -f "$count_file" ]]; then echo 0 > "$count_file"; fi
#         c=$(cat "$count_file")
#         c=$((c+1))
#         echo "$c" > "$count_file"

#         percent=$((100 * c / total))
#         echo "[進度] $c / $total genes (${percent}%)"
#     ) 200>/tmp/plink_parallel.lock
# }

# export -f run_one

# # 用 GNU Parallel 平行執行
# parallel -j 12 run_one ::: "${file_list[@]}"


# # 執行
# # 按以下按鍵，存檔後退出
# Ctrl+O, Enter, Ctrl+X

# # 賦予執行權限並執行
# chmod +x /root/fusion_project/run_parallel_extract.sh
# /root/fusion_project/run_parallel_extract.sh


# OUTLINE: Record cis SNP counts ----
## 紀錄 cis snp number ----
# 挑出原本有的snp或是impute r2>0.8 的，snplists 裡的檔案寫著每個probe cis snp ，先複製到 C:/TWAS/snplists

# 取得資料夾內所有檔案完整路徑
files <- list.files("C:/TWAS/snplists", full.names = TRUE)


# 計算每個檔案的 row 數
result <- rbindlist(lapply(files, function(f) {
  n <- fread(f, header = FALSE, nThread = 4)[, .N]
  data.table(
    file = basename(f),
    rows = n
  )
}))


names(result) <- c("file_name", "cisSNP_number")
setkey(result, file_name)
result[, probe := sub(".*(ILMN_[0-9]+).*", "\\1", file_name)]


fwrite(result,
  "\\\\wsl.localhost\\Ubuntu-22.04\\root\\fusion_project\\Normal_part\\snplists\\cisSNP_number_record.txt",
  row.names = F, col.names = T, sep = "\t"
)


# #  算 weight 的腳本
# cd /root/fusion_project
# nano run_fusion_compute_weights.sh


# # 腳本內容
# GENE_LIST="/root/fusion_project/Normal_part/probe_list.txt"
# FUSION_SCRIPT="/root/fusion_project/FUSION.compute_weights.R"
# GCTA_PATH="/root/my_bin/gcta64"
# PLINK_PATH="/root/my_bin/plink"
# OUTPUT_DIR="/root/fusion_project/Normal_part/weight_directory"
# TMP_BASE="/root/fusion_project/Normal_part/tmp_fusion"
# mkdir -p "$OUTPUT_DIR"
# mkdir -p "$TMP_BASE"

# # ================================
# # 檢查基因列表檔案
# # ================================
# if [[ ! -f "$GENE_LIST" ]]; then
#     echo "錯誤: 找不到基因列表檔案 $GENE_LIST"
#     exit 1
# fi

# # 修正 Windows CRLF
# sed -i 's/\r$//' "$GENE_LIST"

# # ================================
# # 計算總基因數，用於進度
# # ================================
# TOTAL_GENES=$(grep -cve '^\s*$' "$GENE_LIST")
# echo "總共基因數量：$TOTAL_GENES"

# # 計數器
# count=0

# # ================================
# # 逐行處理每個基因
# # ================================
# while read -r GENE; do
#     # 空行跳過
#     if [[ -z "$GENE" ]]; then
#         continue
#     fi

#     count=$((count+1))
#     percent=$((100 * count / TOTAL_GENES))
#     echo "------------------------------------------------------"
#     echo "[$percent%] 處理基因 $count / $TOTAL_GENES : $GENE"

#     # 檔案路徑
#     BFILE_PATH="/root/fusion_project/Normal_part/cisSNP_binary_file/Cissnp_in_${GENE}"
#     PHENO_PATH="/root/fusion_project/Normal_part/expression_by_probe/expression_for_${GENE}.txt"
#     OUT_PATH="${OUTPUT_DIR}/${GENE}"
#     TMP_PATH="${TMP_BASE}/${GENE}"
#     mkdir -p "$TMP_PATH"

#     # 檔案檢查
#     if [[ ! -f "${BFILE_PATH}.bed" ]]; then
#         echo "警告: 缺少 ${BFILE_PATH}.bed，跳過 ${GENE}"
#         continue
#     fi
#     if [[ ! -f "$PHENO_PATH" ]]; then
#         echo "警告: 缺少 ${PHENO_PATH}，跳過 ${GENE}"
#         continue
#     fi

#     # ================================
#     # 執行 FUSION 權重計算 (安靜模式)
#     # ================================
#     Rscript "$FUSION_SCRIPT" \
#         --bfile "$BFILE_PATH" \
#         --pheno "$PHENO_PATH" \
#         --models top1,lasso,enet \
#         --PATH_gcta "$GCTA_PATH" \
#         --PATH_plink "$PLINK_PATH" \
#         --out "$OUT_PATH" \
#         --tmp "$TMP_PATH" \
#         --hsq_p 1 \
#         --verbose 0

#     exit_code=$?

#     # 成功或失敗處理
#     if [[ $exit_code -eq 0 ]]; then
#         # 只刪除存在的資料夾
#         if [[ -d "$TMP_PATH" ]]; then
#             rm -rf "$TMP_PATH"
#         fi
#         echo "成功：$GENE"
#     else
#         echo "失敗：$GENE"
#     fi

# done < "$GENE_LIST"

# echo "------------------------------------------------------"
# echo "所有基因處理完畢。"

# # 存檔
# # 按以下按鍵，存檔後退出
# Ctrl+O, Enter, Ctrl+X

# # 賦予執行權限並執行
# chmod +x /root/fusion_project/run_fusion_compute_weights.sh
# /root/fusion_project/run_fusion_compute_weights.sh















# # 整理 weight result 的腳本 (跑)
# # 建立腳本
# nano /home/jcc623/fusion_project/summarize_weights.R


# # 腳本內容
# pkgs <- c("parallel")

# for (p in pkgs) {
#     if (!requireNamespace(p, quietly = TRUE))
#         install.packages(p)
#     suppressMessages(library(p, character.only = TRUE))
# }

# # 設定資料夾
# wgt_dir <- "/home/jcc623/fusion_project/Normal_part/weight_directory"
# output_csv <- "/home/jcc623/fusion_project/Normal_part/fusion_weights_summary.csv"

# # 取得所有 .wgt.RDat 檔案
# files <- list.files(wgt_dir, pattern = "\\.wgt\\.RDat$", full.names = TRUE)
# cat("開始平行處理", length(files), "個權重檔案...\n")

# # 設定 CPU 使用核心數（可調整）
# NCORE <- max(1, detectCores() - 4)
# cat("使用核心數:", NCORE, "\n")

# # 定義處理 function（每個檔案跑一次）
# process_one <- function(f) {
#   # 建立獨立環境載入 RDat，避免變數殘留
#   env <- new.env()

#   tryCatch({
#     load(f, envir = env) # 將資料載入到 env 環境中

#     # 提取 Gene_ID (移除前綴與後綴)
#     # 例如：TCGA-HNSC.TUMOR.ABCA8_10351.wgt.RDat -> ABCA8_10351
#     fname <- basename(f)
#     gene_id <- gsub("\\.wgt\\.RDat$", "", fname)

#     # 1. 提取遺傳力 (Heritability)
#     hsq_val <- if (!is.null(env$hsq)) env$hsq[1] else NA
#     hsq_pval <- if (!is.null(env$hsq.pv)) env$hsq.pv[1] else NA

#     # 2. 提取模型表現
#     # 初始化所有模型欄位為 NA
#     top1_r2 <- top1_p <- lasso_r2 <- lasso_p <- enet_r2 <- enet_p <- NA
#     best_model <- NA
#     best_r2 <- NA
#     best_p <- NA
#     n_snps <- if (!is.null(env$wgt.matrix)) nrow(env$wgt.matrix) else NA

#     if (!is.null(env$cv.performance)) {
#       cv <- env$cv.performance

#       # 提取個別模型數值 (根據 row name 提取，較安全)
#       if ("top1" %in% colnames(cv)) {
#         top1_r2 <- cv[1, "top1"]
#         top1_p  <- cv[2, "top1"]
#       }
#       if ("lasso" %in% colnames(cv)) {
#         lasso_r2 <- cv[1 ,"lasso"]
#         lasso_p  <- cv[2, "lasso"]
#       }
#       if ("enet" %in% colnames(cv)) {
#         enet_r2 <- cv[1, "enet"]
#         enet_p  <- cv[2, "enet"]
#       }

#       # 判斷最佳模型 (R2 最大者)
#       best_idx <- which.max(cv[1 ,])
#       best_model <-  colnames(cv)[best_idx]
#       best_r2 <- cv[1, best_idx]
#       best_p <- cv[2, best_idx]
#     }

#     # 3. 回傳完整 Data Frame
#     return(data.frame(
#       Gene_ID    = gene_id,
#       N_SNPs     = n_snps,
#       Hsq        = hsq_val,
#       Hsq_Pval   = hsq_pval,
#       Top1_R2    = top1_r2,
#       Top1_Pval  = top1_p,
#       Lasso_R2   = lasso_r2,
#       Lasso_Pval = lasso_p,
#       Enet_R2    = enet_r2,
#       Enet_Pval  = enet_p,
#       Best_Model = best_model,
#       Best_R2    = best_r2,
#       Best_Pval  = best_p,
#       stringsAsFactors = FALSE
#     ))

#   }, error = function(e) {
#     cat("錯誤: 無法處理", f, "\n")
#     return(NULL)
#   })
# }

# # 平行執行
# results <- mclapply(files, process_one, mc.cores = NCORE)

# # 合併所有結果
# summary_table <- do.call(rbind, results)

# # 輸出 CSV
# write.csv(summary_table, file = output_csv, row.names = FALSE, quote = FALSE)

# cat("------------------------------------------------------\n")
# cat("平行處理完成！\n")
# cat("共整理了", nrow(summary_table), "個基因。\n")
# cat("結果已儲存至:", output_csv, "\n")

# # 存檔

# # 按以下按鍵，存檔後退出
# Ctrl+O, Enter, Ctrl+X

# # 執行
# Rscript /home/jcc623/fusion_project/summarize_weights.R













# # 製造紀錄 probe pos 檔案
# nano /root/fusion_project/make_pos_file.R

# ## 腳本 ----
# pkgs <- c("data.table", "dplyr")

# for (p in pkgs) {
#     if (!requireNamespace(p, quietly = TRUE))
#         install.packages(p)
#     suppressMessages(library(p, character.only = TRUE))
# }


# MASTER_ANNOTATION_FILE <- "/mnt/d/oral_cancer/expression/expression_data/probe24526infoB_OneInterval.txt"
# WEIGHTS_DIR <- "/root/fusion_project/Normal_part/weight_directory"
# OUTPUT_POS_FILE <- paste0(WEIGHTS_DIR, "/weights.pos")

# if (!file.exists(MASTER_ANNOTATION_FILE)) {
#   stop(paste("錯誤：找不到主要註釋檔案:", MASTER_ANNOTATION_FILE,
#              "\n請先確認路徑是否正確，或自行建立此檔案。"))
# }
# if (!dir.exists(WEIGHTS_DIR)) {
#   stop(paste("錯誤：找不到權重目錄:", WEIGHTS_DIR))
# }

# master_df <- fread(MASTER_ANNOTATION_FILE,header=T)
# names(master_df) <- c("ID", "CHR", "P0", "P1")

# wgt_files <- list.files(WEIGHTS_DIR, pattern = "\\.wgt\\.RDat$", full.names = FALSE)

# if (length(wgt_files) == 0) {
#   stop("錯誤：權重目錄中找不到任何 *.wgt.RDat 檔案。請先執行 FUSION.compute_weights.R")
# }

# # 從檔名中提取基因 ID (這是唯一的連接鍵)
# successful_genes <- data.frame(
#   WGT = wgt_files,
#   ID = gsub("\\.wgt\\.RDat$", "", wgt_files),
#   stringsAsFactors = FALSE
# )

# final_pos_df <- merge(
#   successful_genes,
#   master_df,
#   by = "ID",
#   all.x = FALSE,
#   all.y = FALSE
# )

# setcolorder(final_pos_df,c("WGT", "ID", "CHR", "P0", "P1"))

# fwrite(final_pos_df,
#        OUTPUT_POS_FILE,
#        row.names = F, col.names = T, sep = "\t")

# cat("--------------------------------------------------\n")
# cat("檔案路徑:", OUTPUT_POS_FILE, "\n")
# cat("總共包含", nrow(final_pos_df), "個基因。\n")


# # 按以下按鍵，存檔後退出
# Ctrl+O, Enter, Ctrl+X

# # 執行
# Rscript /root/fusion_project/make_pos_file.R


















# TWAS
# summary file, ref panel file, weight file 當中snp 命名方式、版本要相同，這裡一致用 hg19 的 chr:pos 命名

# OUTLINE: Convert rsID to chromosome position ----
## rsID -> chr:pos ----
### for summary data  ----
sumstat <- fread("D:/oral_cancer/TWAS/summary_stat/FinnGen/C3_ORALCAVITY_EXALLC")
names(sumstat) <- c(
  "chr", "pos", "A2", "A1", "SNP", "nearest_genes", "pval", "mlogp", "beta",
  "sebeta", "af_alt", "af_alt_cases", "af_alt_controls"
)
sumstat[, SNP := paste0(chr, ":", pos)]
sumstat[, Z := beta / sebeta]
sumstat <- sumstat %>%
  select("SNP", "A1", "A2", "Z")

sumstat <- sumstat[apply(sumstat, 1, function(x) all(!is.na(x) & x != "")), ]


fwrite(sumstat,
  "D:/oral_cancer/TWAS/summary_stat/FinnGen/C3_ORALCAVITY_EXALLC.txt",
  row.names = F, col.names = T, sep = "\t"
)


# 移至 ubuntu (code under ubuntu)
# mv /mnt/d/oral_cancer/TWAS/summary_stat/FinnGen/C3_ORALCAVITY_EXALLC.txt /root/fusion_project/C3_ORALCAVITY_EXALLC.txt


### for ref data  ----

for (i in 1:22) {
  chr1 <- sprintf("D:/oral_cancer/TWAS/summary_stat/LDREF_hg19/1000G.EUR.%d.bim", i) %>%
    fread(header = FALSE)

  names(chr1) <- c("chr", "id", "maf", "pos", "A1", "A2")
  chr1[, id := paste0(chr, ":", pos)]
  fwrite(chr1,
    sprintf("D:/oral_cancer/TWAS/summary_stat/LDREF_hg19/1000G.EUR.%d.bim", i),
    row.names = F, col.names = F, sep = "\t"
  )
}


# # 移至 ubuntu (code under ubuntu)
# mv /mnt/d/oral_cancer/TWAS/summary_stat/LDREF_hg19 /root/fusion_project/LDREF_hg19


# ## 腳本 ----
# nano /root/fusion_project/run_fusion_assoc_test.sh














# ## 平行話版本 ----


# # 腳本內容
# FUSION_R="/root/fusion_project/FUSION.assoc_test.R"
# WEIGHTS_DIR="/root/fusion_project/Normal_part/weight_directory/"
# WEIGHTS_LIST="/root/fusion_project/Normal_part/weight_directory/weights.pos"
# SUMSTATS="/root/fusion_project/C3_ORALCAVITY_EXALLC.txt"
# OUT_PREFIX="/root/fusion_project/Normal_part/twas/C3_ORAL_TWAS_Result"

# mkdir -p "$OUT_PREFIX"

# # 檢查權重清單是否為空
# if [ ! -s "$WEIGHTS_LIST" ]; then
#     echo "錯誤: 權重清單為空！"
#     exit 1
# fi

# # 進度檔
# PROGRESS_FILE="${OUT_PREFIX}/.progress"
# echo 0 > "$PROGRESS_FILE"
# TOTAL=22

# # ================================
# # 定義單染色體執行函數
# # ================================
# run_chr() {
#     chr="$1"
#     OUT_NAME="${OUT_PREFIX}_chr${chr}.dat"

#     Rscript "$FUSION_R" \
#         --sumstats "$SUMSTATS" \
#         --weights "$WEIGHTS_LIST" \
#         --weights_dir "$WEIGHTS_DIR" \
#         --ref_ld_chr "/root/fusion_project/LDREF_hg19/1000G.EUR." \
#         --chr "$chr" \
#         --out "$OUT_NAME"

#     # 更新進度
#     (
#         flock -x 200
#         c=$(cat "$PROGRESS_FILE")
#         c=$((c+1))
#         echo "$c" > "$PROGRESS_FILE"
#         percent=$((100 * c / TOTAL))
#         echo "[進度] $c / $TOTAL 染色體完成 (${percent}%)"
#     ) 200>/tmp/twas_parallel.lock
# }

# export -f run_chr
# export FUSION_R WEIGHTS_LIST WEIGHTS_DIR SUMSTATS OUT_PREFIX PROGRESS_FILE TOTAL

# # ================================
# # 使用 GNU Parallel 平行執行
# # ================================
# parallel -j 4 run_chr ::: {1..22}

# # ================================
# # 合併結果
# # ================================
# echo "========================================"
# echo "正在合併 22 條染色體的結果..."
# echo "========================================"

# FINAL_FILE="${OUT_PREFIX}_FINAL.dat"

# # 抓 chr1 header
# head -n 1 "${OUT_PREFIX}_chr1.dat" > "$FINAL_FILE"

# # append 其他染色體資料
# for chr in {1..22}; do
#     if [ -f "${      }_chr${chr}.dat" ]; then
#         tail -n +2 "${OUT_PREFIX}_chr${chr}.dat" >> "$FINAL_FILE"
#     fi
# done

# echo "TWAS完成！"














# # 合併TWAS file
# nano /root/fusion_project/combine_TWAS_result.sh

# # 腳本內容

# INPUT_PATTERN="/root/fusion_project/Normal_part/twas/*"
# # 合併後的輸出檔案名稱
# OUTPUT_FILE="/root/fusion_project/Normal_part/twas/C3_ORAL_TWAS_Result_chr1-22.dat"

# # 當一個萬用字元模式沒有匹配到任何檔案時，防止將/root/fusion_project/Normal_part/twas/* 當成檔案名稱處理，將 /root/fusion_project/Normal_part/twas/* 替換為空字串
# shopt -s nullglob

# files=( $INPUT_PATTERN )

# if [ ${#files[@]} -eq 0 ]; then
#     echo "錯誤：在當前目錄中找不到任何符合 '$INPUT_PATTERN' 模式的檔案。"
#     echo "請確保腳本與 .dat 檔案位於同一目錄或提供正確的路徑。"
#     exit 1
# fi

# echo "開始合併以下檔案："
# # 2. 迴圈處理所有找到的 .dat 檔案
# for file in $INPUT_PATTERN; do
#     # 檢查是否為當前的輸出檔案 (避免重複寫入自己)
#     if [ "$file" != "$OUTPUT_FILE" ]; then
#         echo "  - 合併：$file"
#         # 使用 cat 將檔案內容追加 (append >>) 到輸出檔案
#         cat "$file" >> "$OUTPUT_FILE"
#     fi
# done

# # 3. 輸出結果
# echo "---"
# echo "✅ 合併完成！"
# echo "所有資料已成功寫入到檔案：'$OUTPUT_FILE'"
# echo "檔案大小：$(du -h "$OUTPUT_FILE" | awk '{print $1}')"

# # 恢復 nullglob 預設狀態
# shopt -u nullglob


# # 存檔
# # 按以下按鍵，存檔後退出
# Ctrl+O, Enter, Ctrl+X

# # 賦予權限並執行
# chmod +x /root/fusion_project/combine_TWAS_result.sh
# /root/fusion_project/combine_TWAS_result.sh














# OUTLINE: Read TWAS output ----
## read ----

# 讀權重檔案 "C:\Users\user\Desktop\TCGA-HNSC.TUMOR\TCGA-HNSC.TUMOR.TRIT1_54802.wgt.RDat"
target_env <- new.env()

# 2. 將檔案載入到這個特定環境中
load("D:/root/fusion_project/Normal_part/weight_directory/ILMN_1770811.wgt.RDat", envir = target_env)

# 3. 從小房間裡把 cv.performance 拿出來，存成你的新變數
my_perf <- target_env$cv.performance

# 4. 查看結果
print(my_perf)

# lasso model, 幾個 snp 係數不為0
which(target_env$wgt.matrix[, "lasso"] != 0)
