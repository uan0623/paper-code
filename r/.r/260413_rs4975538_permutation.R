# TITLE: rs4975538 Permutation
# SUBTITLE: Run permutation analysis for rs4975538 and related association signals.
# SEARCH TAGS: rs4975538, permutation, association, SNP, eQTL
# NOTE: Code below is unchanged; only this navigation header was added.



# package ----
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
library(bestNormalize)



  # finngen GWAS most sig. rs4975538  eQTL result ----
  
  maf_hg18_19 <- fread("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/cis_hg18_19_snp_MAF.txt",
    header = T
  )
  exp_diff <- fread("C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_0.8/outcome/OQN_exp_different_r2_0.8.txt", header = T)
  maf_sub <- maf_hg18_19[hg18_snpID == "5:1333830", ][
    ,
    c("hg18_snpID", "hg19_snpID", "rsID", "MAF", "REF", "ALT", "R2", "impute_type")
  ]

 
  
  for (part in c("N", "T")) {
     pval_all <- sprintf(
       "C:/Peter/OQN_FIXpeople_before_eQTL/trash/OQN_maf_gt_%s_pvalue.txt",
       part
     ) %>% fread()
     
      af_clump_fdr <- pval_all[SNP %in% "5:1333830", ]

      setnames(af_clump_fdr, old = c("gene", "SNP"), c("PROBE_ID", "hg18_snpID"))
      af_clump_fdr <- exp_diff[af_clump_fdr, on = .(PROBE_ID)]
      af_clump_fdr <- maf_sub[af_clump_fdr, on = .(hg18_snpID)]

      af_clump_fdr[, FDR := NULL]
      af_clump_fdr[, FDR_20asso := p.adjust(`p-value`, method = "BH") %>%
        format(digits = 10, scientific = T) %>%
        as.numeric()]

      setkey(af_clump_fdr, `p-value`)

    a <- sprintf(
       "C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_0.8/outcome/OQN_%s_FDR_R2_0.8.txt",
       part
     ) %>% fread()
     
    setkey(a, `p-value`)

    # 找 pval 至少多少能過 bon 門檻
    pval_thres <- a[sum(a$sig_pval_Bonfi) + 1, `p-value`]
    if (all(af_clump_fdr$`p-value` > pval_thres)) {

      # 加上 sig_pval_Bonfi
      af_clump_fdr[, sig_pval_Bonfi := 0]
    } else {
      stop("5:1333830 association pass bon threshold ! ")
    }

    fwrite(af_clump_fdr,
      sprintf("C:/Peter/rs4975538_permutation/outcome/OQN_FIXpeople_before_eQTL/rs4975538_%s_eQTL.txt",
       part),
      row.names = F, col.names = T, sep = "\t"
    )
  }




  rm(list = ls())
  gc()





  # probe cor ----
  exp_N_OQN <- fread("C:/Peter/OQN_FIXpeople_before_eQTL/OQN_Exp_mulInterval_N.txt", header = T)
  probe_rs4975538 <- fread("C:/Peter/rs4975538_permutation/outcome/OQN_FIXpeople_before_eQTL/rs4975538_N_eQTL.txt")

  probe_N <- unique(probe_rs4975538$PROBE_ID)
  exp_N_OQN <- exp_N_OQN[PROBE_ID %in% probe_N, ]
  exp_N_OQN <- t(exp_N_OQN) %>% as.data.table()

  names(exp_N_OQN) <- as.character(exp_N_OQN[1, ])
  exp_N_OQN <- exp_N_OQN[-c(1), ]
  exp_N_OQN[, (names(exp_N_OQN)) := lapply(.SD, as.numeric)]
  setcolorder(exp_N_OQN, probe_N)

  # before OQN exp
  exp <- fread("D:/oral_cancer/expression/expression_data/ExpRes11sv20120601Gene.txt", header = T)
  exp_N <- exp[, c(2, 43:82)]
  names(exp_N) <- c(
    "PROBE_ID",
    c(sprintf("0%dB", 1:9), sprintf("%dB", 10:40))
  )

  probe_N <- unique(probe_rs4975538$PROBE_ID)
  exp_N <- exp_N[PROBE_ID %in% probe_N, ]
  exp_N <- t(exp_N) %>% as.data.table()

  names(exp_N) <- as.character(exp_N[1, ])
  exp_N <- exp_N[-c(1), ]
  exp_N[, (names(exp_N)) := lapply(.SD, as.numeric)]
  setcolorder(exp_N, probe_N)


  plot_abs_correlation <- function(exp_data, title) {
    cor_mat <- cor(exp_data, use = "pairwise.complete.obs")
    diag(cor_mat) <- 1
    cor_dt <- as.data.table(as.table(cor_mat))
    setnames(cor_dt, c("probe_x", "probe_y", "correlation"))
    cor_dt[, abs_cor_group := cut(
      abs(correlation),
      breaks = c(-Inf, 0.3, 0.5, Inf),
      labels = c("|r| < 0.3", "0.3 <= |r| < 0.5", "|r| >= 0.5"),
      right = FALSE
    )]
    cor_dt[, probe_x := factor(probe_x, levels = colnames(cor_mat))]
    cor_dt[, probe_y := factor(probe_y, levels = colnames(cor_mat))]

    ggplot(cor_dt, aes(x = probe_x, y = probe_y, fill = abs_cor_group)) +
      geom_tile(color = "white") +
      geom_text(aes(label = sprintf("%.2f", correlation), color = abs_cor_group), size = 3.8) +
      scale_fill_manual(
        values = c(
          "|r| < 0.3" = "#f7fbff",
          "0.3 <= |r| < 0.5" = "#6baed6",
          "|r| >= 0.5" = "#08306b"
        ),
        labels = c(
          expression(abs(r) < 0.3),
          expression(0.3 <= abs(r) * " < " * 0.5),
          expression(abs(r) >= 0.5)
        ),
        drop = FALSE,
        name = "absolute cor."
      ) +
      scale_color_manual(
        values = c(
          "|r| < 0.3" = "black",
          "0.3 <= |r| < 0.5" = "black",
          "|r| >= 0.5" = "white"
        ),
        guide = "none"
      ) +
      coord_fixed() +
      labs(title = title, x = NULL, y = NULL) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 11),
        axis.text.y = element_text(size = 11),
        panel.grid = element_blank()
      )
  }

  png("C:/Peter/rs4975538_permutation/outcome/OQN_FIXpeople_before_eQTL/rs4975538_asso_probe_correlation.png",
    width = 14, height = 12, units = "in", res = 300
  )
  
  plot_abs_correlation(exp_N_OQN, title = "rs4975538 Normal Part eQTL Gene Cor. After OQN") 
  # 關閉檔案
  dev.off()


  png("C:/Peter/rs4975538_permutation/outcome/OQN_FIXpeople_before_eQTL/rs4975538_asso_probe_correlation_raw.png",
    width = 14, height = 12, units = "in", res = 300
  )
  plot_abs_correlation(exp_N, title = "rs4975538 Normal Part eQTL Gene Cor.")

  # 關閉檔案
  dev.off()






  # eQTL permut ----
  # 1. 挑出 5:1333830 (rs4975538) 的 pos, snp binary file
  # 2. exp data 多區間合併

  all_paths <- c(
    "C:/Peter/rs4975538_permutation/outcome/OQN_FIXpeople_before_eQTL/",
    "C:/Peter/rs4975538_permutation/outcome/OQN_FIXpeople_before_eQTL/"
  )

  # 沒有資料夾，就自動新增
  for (i in all_paths) {
    if (!dir.exists(i)) dir.create(i, recursive = TRUE)
  }






# snp gt data ----
# a <- fread("C:/Peter/OQN_FIXpeople_before_eQTL/cis_snp_gt_maf_chr5.txt")
# a <- a[ID %in% "5:1333830",]
# fwrite(a,
#        "C:/Peter/rs4975538_permutation/outcome/cisSNP.txt",
#         row.names = F, col.names = T, sep = "\t")






# exp data ----
probe_rs4975538 <- fread("C:/Peter/rs4975538_permutation/outcome/OQN_FIXpeople_before_eQTL/rs4975538_N_eQTL.txt")
probe <- unique(probe_rs4975538$PROBE_ID)
permu_times <- 10000

for (j in c("N","T")) {
  
  exp <- sprintf("C:/Peter/OQN_FIXpeople_before_eQTL/OQN_Exp_mulInterval_%s.txt",j) %>% 
    fread()
  # 只留下 rs4975538 附近的 gene，N,T 都是
  exp <- exp[PROBE_ID %in% probe,]
  
  res_list <- vector("list", length = permu_times + 1)
  res_list[[1]] <- exp 

  for (permu_th in 1:permu_times) {
  
  # 打亂，但保持 expression 之間的關係
  set.seed(permu_th)
  exp_permu <- cbind(exp[,c(1:2)],
                 exp[, sample(.SD), .SDcols = names(exp)[-c(1:2)]] )
  names(exp_permu) <- c("PROBE_ID",c(sprintf("0%dB",1:9), sprintf("%dB", 10:40)) )
  exp_permu[, PROBE_ID := paste0(PROBE_ID,"_permu",permu_th)]
  res_list[[permu_th + 1]] <- exp_permu
  
  }
  
  # 一次性合併，比迴圈 rbind 快
  final <- rbindlist(res_list, use.names = TRUE)

  fwrite(final,
         sprintf("C:/Peter/rs4975538_permutation/outcome/OQN_FIXpeople_before_eQTL/OQN_Exp_%s_permu.txt",j),
         row.names = F, col.names = T, sep = "\t")
  
}


  library(MatrixEQTL)
  probe_rs4975538 <- fread("C:/Peter/rs4975538_permutation/outcome/OQN_FIXpeople_before_eQTL/rs4975538_N_eQTL.txt")
  probe <- unique(probe_rs4975538$PROBE_ID)

  full_snpspos = fread("C:/Peter/OQN_FIXpeople_before_eQTL/snp_pos.txt", header = TRUE)
  full_genepos = fread("D:/oral_cancer/expression/expression_data/probe_pos_mulInterval.txt", header = TRUE)
  current_snpspos = full_snpspos[ID %in% "5:1333830", ]
  current_genepos = full_genepos[Gene %in% probe, ]


  permu_all <- rbindlist(lapply(1:permu_times, function(i) {
    # 複製原始資料
    dt_tmp <- copy(current_genepos)

    # 修改 PROBE_ID，加上 _permu 加上序號
    dt_tmp[, Gene := paste0(Gene, "_permu", i)]
    return(dt_tmp)
  }))


  snps = SlicedData$new()
  snps$fileDelimiter = "\t"
  snps$fileOmitCharacters = "NA"
  snps$fileSkipRows = 1 # 跳過表頭
  snps$fileSkipColumns = 1 # 第一欄是SNP ID
  snps$fileSliceSize = 50000 # 每次讀入多少行
  snps_path = "C:/Peter/rs4975538_permutation/outcome/cisSNP.txt"
  snps$LoadFile(snps_path)


  gene_T = SlicedData$new()
  gene_T$fileDelimiter = "\t"
  gene_T$fileOmitCharacters = "NA"
  gene_T$fileSkipRows = 1
  gene_T$fileSkipColumns = 1
  gene_T$fileSliceSize = 20000
  gene_T$LoadFile("C:/Peter/rs4975538_permutation/outcome/OQN_FIXpeople_before_eQTL/OQN_Exp_T_permu.txt")

  gene_N = SlicedData$new()
  gene_N$fileDelimiter = "\t"
  gene_N$fileOmitCharacters = "NA"
  gene_N$fileSkipRows = 1
  gene_N$fileSkipColumns = 1
  gene_N$fileSliceSize = 20000
  gene_N$LoadFile("C:/Peter/rs4975538_permutation/outcome/OQN_FIXpeople_before_eQTL/OQN_Exp_N_permu.txt")

  # 執行 Normal
  test <- Matrix_eQTL_main(
    snps = snps,
    gene = gene_N,
    snpspos = current_snpspos, # 使用過濾後的座標
    genepos = permu_all, # 使用過濾後的座標
    cvrt = SlicedData$new(), # 沒有 covariates
    output_file_name = NULL,
    output_file_name.cis = "C:/Peter/rs4975538_permutation/outcome/OQN_FIXpeople_before_eQTL/gt_N_permu.txt",
    useModel = modelLINEAR,
    errorCovariance = numeric(),
    pvOutputThreshold.cis = 1, #  # cis snp pvalue <= 1
    pvOutputThreshold = 0,
    cisDist = 1e6,
    verbose = T,
    pvalue.hist = F # 關閉繪圖可微幅加速
  )

  # 執行 Tumor
  test_1 <- Matrix_eQTL_main(
    snps = snps,
    gene = gene_T,
    snpspos = current_snpspos, # 使用過濾後的座標
    genepos = permu_all, # 使用過濾後的座標
    cvrt = SlicedData$new(), # 沒有 covariates
    output_file_name = NULL,
    output_file_name.cis = "C:/Peter/rs4975538_permutation/outcome/OQN_FIXpeople_before_eQTL/gt_T_permu.txt",
    useModel = modelLINEAR,
    errorCovariance = numeric(),
    pvOutputThreshold.cis = 1, #  # cis snp pvalue <= 1
    pvOutputThreshold = 0,
    cisDist = 1e6,
    verbose = T,
    pvalue.hist = F # 關閉繪圖可微幅加速
  )


  # 釋放記憶體
  rm(list = ls())
  gc()

  # OUTLINE: Summarize permutation results ----
  # 整理結果 ----
  # pval_permu 加進資料
  probe_rs4975538 <- fread("C:/Peter/rs4975538_permutation/outcome/OQN_FIXpeople_before_eQTL/rs4975538_N_eQTL.txt")
  probe <- unique(probe_rs4975538$PROBE_ID)


  for (part in c("N", "T")) {
    df_raw <- sprintf("C:/Peter/rs4975538_permutation/outcome/OQN_FIXpeople_before_eQTL/rs4975538_%s_eQTL.txt", part) %>%
      fread()
    setkey(df_raw, `p-value`)

    gt <- sprintf("C:/Peter/rs4975538_permutation/outcome/OQN_FIXpeople_before_eQTL/gt_%s_permu.txt", part) %>%
      fread()
    names(gt) <- c("SNP", "probe", "beta", "t_stat", "pval", "FDR")
    # "ILMN_1794303_permu264" -> 264
    gt[, permu_idx := str_extract(probe, "(?<=permu)\\d+") %>%
      as.numeric()]


    # 3. 在每一組 permu 中，按 pval 由小到大排序，並給予 rank (1~20)
    setorder(gt, permu_idx, pval)
    gt[, rank := 1:.N, by = permu_idx]

    pval_list <- split(gt$pval, gt$rank)

    df_raw <- df_raw[order(`p-value`)]
    raw_adj <- numeric(20)
    for (j in 1:20) {
      pval_all <- pval_list[[j]]
      raw_adj[j] <- mean(pval_all <= df_raw$`p-value`[j])
    }

    adj_stepdown <- cummax(raw_adj)

    # 加 pval_permu 進原始資料
    df_raw[, pval_permu := adj_stepdown]




    
    if (part == "T") {
      sig_cols <- c(
        sprintf("0.8_%s_sigCis-SNP_number", part),
        sprintf("0.8_%s_sigCis-SNP_number (bon)", part)
      )
      df_raw[, c(
        "hg19_snpID",
        sprintf("0.8_N_sigCis-SNP_number%s", c("", " (bon)"))
      ) := NULL]
    } else {
       sig_cols <- c(
        sprintf("0.8_%s_sigCis-SNP_number", part),
        sprintf("0.8_%s_sigCis-SNP_number (bon)", part)
      )
      df_raw[, c(
        "hg19_snpID",
        sprintf("0.8_T_sigCis-SNP_number%s", c("", " (bon)"))
      ) := NULL]
    }


    setcolorder(
      df_raw,
      c(
        "chr", "PROBE_ID", "Gene", "PROBE_COORDINATES", sig_cols, "hg18_snpID", "rsID", "impute_type",
        "R2", "MAF", "REF", "ALT", "beta", "t-stat", "p-value", "pval_permu", "FDR_20asso",
        "sig_pval_Bonfi", "t_OQN", "pval_OQN", "pval_OQN_BH", "sig_OQN", "sig_OQN_Bonfi"
      ) 
    )

    fwrite(df_raw,
      sprintf("C:/Peter/rs4975538_permutation/outcome/OQN_FIXpeople_before_eQTL/rs4975538_%s_eQTL.txt", part),
      row.names = F, col.names = T, sep = "\t"
    )
  }




# 問題 ----

# 基因有可能跨不同 chr?
#    TargetID ProbeID     PROBE_ID CHROMOSOME                   PROBE_COORDINATES
#      <char>   <int>       <char>      <num>                              <char>
# 1:   GABPB2 1170706 ILMN_2161007          1                 149358958-149359007
# 2:   GABPB2 1690215 ILMN_2331701         15                   48365611-48365660
# 3:   GABPB2 7200431 ILMN_1761147         15                   48369074-48369123
# 4:     PRG2 1580195 ILMN_1722290         19                       764090-764139
# 5:     PRG2 3840072 ILMN_1729314         11 56911880-56911914:56912626-56912640





