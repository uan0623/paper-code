## package ----
library(ggplot2)
library(data.table)
library(dplyr)
library(stringr)
library(readxl) # read_xlsx ft




# 讀權重檔案 ----
# target_env <- new.env()
# # 2. 將檔案載入到這個特定環境中


# load("C:/Peter/gene_enrichment/code_project/data/GTEx_Heart_Left/GTExv8.ALL.Heart_Left_Ventricle/ENSG00000000419.12.wgt.RDat", envir = target_env)
# # 3. 從小房間裡把 cv.performance 拿出來，存成你的新變數
# my_perf <- target_env$cv.performance
# # 4. 查看結果
# print(my_perf)
# # lasso model, 幾個 snp 係數不為0
# which(target_env$wgt.matrix[,"lasso"]!=0)
# target_env$hsq.pv
#  head(target_env$wgt.matrix)


## 選跟 eQTL sig. 同數量的，跟 HNSC 交集數量 hist ----


run_tcga_eqtl_enrichment <- function(
  tissue_name,
  weight_path,
  output_dir
) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  cat(
    "After remove R2<0.8 snp, run eQTL genes:",
    uniqueN(all_eQTL_gene),
    "\n"
  )

  tcga_weight <- fread(weight_path)
  tcga_weight[, Gene := gsub("_.*$", "", Gene_ID)]
  cat(
    tissue_name, "genes:",
    length(tcga_weight$Gene),
    "\n"
  )


  tcga_weight <- tcga_weight[!is.na(Gene) & Gene != ""]
  tcga_weight <- tcga_weight[!is.na(Top1_Pval) & Top1_Pval != ""]
  cat(
    tissue_name, "After remove NULL and NA GeneID/Top1_Pval, genes:",
    uniqueN(tcga_weight$Gene),
    "\n"
  )

  gt_N <- fread("C:/Peter/rawData_eQTL/r2_filter_0.8/outcome/raw_N_FDR_R2_0.8.txt")

  eQTL_gene <- unique(gt_N$Gene)

  tcga_sigGene <- tcga_weight[Top1_Pval < 0.01, Gene] %>%
    unique()
  cat(
    tissue_name, "sig. genes:",
    uniqueN(tcga_sigGene),
    "\n"
  )
  cat(
    "eQTL sig. genes:",
    uniqueN(eQTL_gene),
    "\n"
  )


  cat(
    tissue_name, "significant genes recorded in eQTL:",
    tcga_sigGene[tcga_sigGene %in% all_eQTL_gene] %>% uniqueN(), "\n"
  )


  # 交集數量

  common_gene <- all_eQTL_gene
  sample_gene_number <- tcga_sigGene[tcga_sigGene %in% all_eQTL_gene] %>% uniqueN()
  observed_intersect <- intersect(eQTL_gene, tcga_sigGene) %>% uniqueN()

  cat(
    "intersect of eQTL,", tissue_name, "sig gene:",
    observed_intersect, "\n"
  )


  # 從共同基因隨機挑，看多少 在 oral 顯著 ----
  library(parallel)


  oral_sig_sub <- eQTL_gene
  is_sig <- common_gene %in% oral_sig_sub
  random_times <- 1e9

  n_core <- max(1, detectCores() - 2)

  chunk_times <- rep(random_times %/% n_core, n_core)
  remainder <- random_times %% n_core

  if (remainder > 0) {
    chunk_times[seq_len(remainder)] <- chunk_times[seq_len(remainder)] + 1
  }

  cl <- makeCluster(n_core)
  # 設定平行亂數 seed，不同核心會自動用不同 seed
  clusterSetRNGStream(cl, 123)

  clusterExport(
    cl,
    varlist = c("is_sig", "sample_gene_number", "observed_intersect"),
    envir = environment()
  )
  # parLapply 類似 lapply，chunk_times list 數字分別傳入 n_sim
  core_result <- parLapply(cl, chunk_times, function(n_sim) {
    out <- integer(n_sim)

    for (i in seq_len(n_sim)) {
      # is_sig 判斷 common_gene 是否在 oral_sig_sub，一串的 T,F，sum 會把T當1, F當0
      out[i] <- sum(sample(is_sig, sample_gene_number))
    }

    out
  })

  stopCluster(cl)

  repeat_gene_number <- unlist(core_result, use.names = FALSE)

  total_random_times <- length(repeat_gene_number)
  total_sum <- sum(repeat_gene_number)
  min_repeat <- min(repeat_gene_number)
  max_repeat <- max(repeat_gene_number)
  ge_observed_count <- sum(repeat_gene_number >= observed_intersect)

  mean_repeat <- total_sum / total_random_times
  empirical_p <- (ge_observed_count + 1) / (total_random_times + 1)

  cat(
    "sample", sample_gene_number, "genes, mean sig in", tissue_name, "genes number:",
    mean_repeat, "\n"
  )

  cat(
    "sample", sample_gene_number, "genes, sig in", tissue_name, "genes number range:",
    min_repeat, max_repeat, "\n"
  )

  cat(
    "sample", sample_gene_number, "genes, sig in", tissue_name, "genes number >=",
    observed_intersect, "times:",
    ge_observed_count, "\n"
  )

  cat(
    "empirical p-value:",
    empirical_p, "\n"
  )

  cat(rep("\n", 2))
  png(file.path(output_dir, sprintf("random_%sgene.png", sample_gene_number)),
    width = 8, height = 6, units = "in", res = 300
  )
  hist(repeat_gene_number,
    main = sprintf("%s times Intersection in %s Sig.", random_times, tissue_name),
    xlab = "Intersection Size each Times",
    col = "skyblue",
    breaks = seq(min_repeat - 0.5, max_repeat + 0.5, by = 1)
  )
  abline(v = observed_intersect, col = "red", lty = 1)
  dev.off()
}





# 執行重抽 ----
all_eQTL <- fread("C:/Peter/rawData_eQTL/trash/raw_maf_gt_N_pvalue.txt", header = T)[
  ,
  .(pval = `p-value`, SNP, probe = gene)
]
R2_filter <- fread("C:/Peter/oral_cancer/expression/outcome/multiple_nucleotide_variant/cis_snp_R2_0.8.txt", header = T)

setkey(all_eQTL, pval)
all_eQTL <- all_eQTL[SNP %in% R2_filter$hg18_snpID, ]
all_eQTL <- all_eQTL[, .SD[1], by = probe]
all_eQTL[, probe := gsub("^([^_]+_[^_]+).*", "\\1", probe)]


probe_info <- fread("C:/Peter/rawData_eQTL/r2_filter_0.8/outcome/raw_exp_different_r2_0.8.txt")[
  ,
  .(probe = PROBE_ID, Gene)
]
all_eQTL <- probe_info[all_eQTL, on = .(probe)]
all_eQTL_gene <- unique(all_eQTL$Gene)
all_eQTL_probe <- unique(all_eQTL$probe)




tissue_type <- c("TCGA-HNSC", "TCGA-LUAD", "TCGA-LUSC")

plot_output_dir <- sprintf(
  "C:/Peter/gene_enrichment/code_project/outcome/rawData_eQTL_N/for_gene/%s",
  tissue_type
)


for (i in seq_along(tissue_type)) {
  run_tcga_eqtl_enrichment(
    tissue_name = tissue_type[i],
    weight_path = sprintf(
      "C:/Peter/gene_enrichment/code_project/data/%s.TUMOR/fusion_project/fusion_weights_summary.csv",
      tissue_type[i]
    ),
    output_dir = plot_output_dir[i]
  )
}


# 有m白球，n黑球，抽k球，抽後不放回，白球數量>q的機率。lower.tail=T 則是，白球數量<=q的機率
phyper(
  # 兩邊都顯著數量-1
  q = 44 - 1,
  # oral sig number
  m = 177,
  # background - oral sig number
  n = 15992 - 177,
  # TCGA sig number
  k = 1540,
  lower.tail = FALSE
) %>% sprintf("%.2e", .)


# 1e9 次抽取結果，4 hr per tissue ----


# After remove R2<0.8 snp, run eQTL genes: 15992
# TCGA-HNSC genes: 2767
# TCGA-HNSC After remove NULL and NA GeneID/Top1_Pval, genes: 2763
# TCGA-HNSC sig. genes: 2069
# eQTL sig. genes: 177
# TCGA-HNSC significant genes recorded in eQTL: 1540
# intersect of eQTL, TCGA-HNSC sig gene: 44
# sample 1540 genes, mean sig in TCGA-HNSC genes number: 17.04477
# sample 1540 genes, sig in TCGA-HNSC genes number range: 0 46
# sample 1540 genes, sig in TCGA-HNSC genes number >= 44 times: 3
# empirical p-value: 4e-09
# Hypergeometric dist. pval: 2.67e-09


# After remove R2<0.8 snp, run eQTL genes: 15992
# TCGA-LUAD genes: 2978
# TCGA-LUAD After remove NULL and NA GeneID/Top1_Pval, genes: 2974
# TCGA-LUAD sig. genes: 2248
# eQTL sig. genes: 177
# TCGA-LUAD significant genes recorded in eQTL: 1644
# intersect of eQTL, TCGA-LUAD sig gene: 47
# sample 1644 genes, mean sig in TCGA-LUAD genes number: 18.1958
# sample 1644 genes, sig in TCGA-LUAD genes number range: 0 45
# sample 1644 genes, sig in TCGA-LUAD genes number >= 47 times: 0
# empirical p-value: 1e-09
# Hypergeometric dist. pval: 6.27e-10


# After remove R2<0.8 snp, run eQTL genes: 15992
# TCGA-LUSC genes: 2548
# TCGA-LUSC After remove NULL and NA GeneID/Top1_Pval, genes: 2545
# TCGA-LUSC sig. genes: 1862
# eQTL sig. genes: 177
# TCGA-LUSC significant genes recorded in eQTL: 1366
# intersect of eQTL, TCGA-LUSC sig gene: 47
# sample 1366 genes, mean sig in TCGA-LUSC genes number: 15.11891
# sample 1366 genes, sig in TCGA-LUSC genes number range: 0 42
# sample 1366 genes, sig in TCGA-LUSC genes number >= 47 times: 0
# empirical p-value: 1e-09
# Hypergeometric dist. pval: 1.03e-12