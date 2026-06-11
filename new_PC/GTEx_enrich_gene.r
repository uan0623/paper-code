## package ----
library(data.table)
library(dplyr)
library(tidyverse)
library(stringr)
library(ggplot2)
library(ggrepel)
library(matrixStats) # rowMins
library(qvalue)
library(bestNormalize)
library(rtracklayer)


# dir_dt <- expand.grid(
#   id_type = c("for_ENSG", "for_gene"),
#   dataset = c(
#     "lungTWAS",
#     "HNSC",
#     "TCGA-LUAD",
#     "TCGA-LUSC",
#     "GTEx-salivary",
#     "GTEx-esophagus",
#     "GTEx-thyroid",
#     "GTEx-lung"
#   ),
#   stringsAsFactors = FALSE
# )

# plot_output_dir <- file.path(
#   "D:/Peter/gene_enrichment/code_project/outcome/rawData_eQTL_N",
#   dir_dt$id_type,
#   dir_dt$dataset
# )


# # 刪舊檔案，刪掉資料夾內所有檔案
# sapply(plot_output_dir, function(x) {
#   if (dir.exists(x)) {
#     message("Deleting: ", x)
#     unlink(x, recursive = TRUE, force = TRUE)
#   } else {
#     message("Directory does not exist sfiles: ", x)
#   }
# })


# 使用 TWAS data ----
# 引用 gene_enrichment 檔案的 建立腳本 summarize_weights.R


# 建立腳本
# nano /mnt/d/Peter/gene_enrichment/code_project/data/GTEx_thyroid/fusion_project/summarize_weights.R


# # 整理 GTEx data
# pkgs <- c("parallel")

# for (p in pkgs) {
#   if (!requireNamespace(p, quietly = TRUE)) {
#     install.packages(p)
#   }
#   suppressMessages(library(p, character.only = TRUE))
# }


# # 設定資料夾
# wgt_dir <- "/mnt/d/Peter/gene_enrichment/code_project/data/GTEx_thyroid/GTExv8.ALL.Thyroid"
# output_csv <- "/mnt/d/Peter/gene_enrichment/code_project/data/GTEx_thyroid/fusion_project/fusion_weights_summary.csv"

# # 取得所有 .wgt.RDat 檔案
# files <- list.files(wgt_dir, pattern = "\\.wgt\\.RDat$", full.names = TRUE)
# cat("開始平行處理", length(files), "個權重檔案...\n")

# # 設定 CPU 使用核心數（可調整）
# NCORE <- max(1, detectCores() - 6)
# cat("使用核心數:", NCORE, "\n")

# # 定義處理 function（每個檔案跑一次）
# process_one <- function(f) {
#   # 建立獨立環境載入 RDat，避免變數殘留
#   env <- new.env()

#   tryCatch({
#     load(f, envir = env) # 將資料載入到 env 環境中

#     # 提取 Gene_ID (移除前綴與後綴)
#     # 例如：GTEx-lung.ABCA8_10351.wgt.RDat -> ABCA8_10351
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


# 按以下按鍵，存檔後退出
# Ctrl+O, Enter, Ctrl+X














# 下載 GTEx v8 用的 annotation ----
# 到 <https://www.gencodegenes.org/human/release_26.html> 下載，
# gtf_file <- "D:/Peter/gene_enrichment/code_project/data/gencode.v26.annotation.gtf.gz"
# gtf <- rtracklayer::import(gtf_file)
# gene_anno <- as.data.table(gtf)
# gene_anno <- gene_anno[type == "gene"]
# gene_anno <- unique(gene_anno[, .(
#   ensg_version = gene_id,
#   query = sub("\\..*$", "", gene_id),
#   symbol = gene_name,
#   gene_type = gene_type,
#   chr = as.character(seqnames),
#   start = start,
#   end = end
# )])

# gene_anno <- gene_anno[chr %in% sprintf("chr%s", 1:22), ]
# fwrite(
#   data.table(
#     a = gene_anno$chr,
#     b = gene_anno$start,
#     c = gene_anno$end
#   ),
#   "D:/Peter/gene_enrichment/code_project/data/gene_hg38.txt",
#   row.names = F, col.names = F, sep = "\t"
# )


# # liftOver 轉換 pos
# cmd <- paste(
#   "/mnt/d/Peter/oral_cancer/liftover/liftOver",
#   "/mnt/d/Peter/gene_enrichment/code_project/data/gene_hg38.txt",
#   "/mnt/d/Peter/oral_cancer/liftover/hg38ToHg18.over.chain",
#   "/mnt/d/Peter/gene_enrichment/code_project/data/gene_hg38TOhg18.bed",
#   "/mnt/d/Peter/gene_enrichment/code_project/data//hg38TOhg18_unmapped.bed"
# )

# system2("wsl", cmd)

# a <- fread("D:/Peter/gene_enrichment/code_project/data/gene_hg38TOhg18.bed", header = F)
# # 忽略轉失敗 .bed 檔案的 #Deleted in new 該行

# unmapped_path <- "D:/Peter/gene_enrichment/code_project/data//hg38TOhg18_unmapped.bed"

# if (!file.exists(unmapped_path) || file.info(unmapped_path)$size == 0) {
#   a_unmapped <- data.table(V1 = character(), V2 = integer(), V3 = integer())
# } else {
#   a_unmapped <- fread(
#     unmapped_path,
#     header = FALSE,
#     comment.char = "#"
#   )
# }

# col_order <- names(gene_anno)
# names(a_unmapped) <- c("chr", "start", "end")
# names(a) <- c("chr_hg18", "start_hg18", "end_hg18")

# # 新增轉換後 pos
# a_unmapped[, c("chr_hg18", "start_hg18", "end_hg18") := "no"]
# final <- merge(gene_anno, a_unmapped, all.x = T)

# success_idx <- is.na(final$chr_hg18)
# final[
#   success_idx,
#   c("chr_hg18", "start_hg18", "end_hg18") := a[, .(chr_hg18, start_hg18, end_hg18)]
# ]

# # 轉換紀錄
# # fwrite(final,
# #   "D:/Peter/gene_enrichment/code_project/data/hg38_hg18.txt",
# #   row.names = F, col.names = T, sep = "\t"
# # )

# final[, c("chr", "start", "end") := NULL]
# setnames(final, old = c("chr_hg18", "start_hg18", "end_hg18"), new = c("chr", "start", "end"))
# setcolorder(final, col_order)
# final[which(gene_anno$start == "no"), c("chr", "start", "end") := NA_character_]
# final[, start := as.numeric(start)]
# final[, end := as.numeric(end)]


# fwrite(final,
#   "D:/Peter/gene_enrichment/code_project/data/gencode_v26_gene_anno_hg18.txt",
#   row.names = F, col.names = T, sep = "\t"
# )








## ft ----
find_nearest_interval <- function(dt, unit, query_start, query_end, for_gene = TRUE) {
  if (for_gene) {
    dt_sub <- copy(dt[gene == unit])
  } else {
    dt_sub <- copy(dt[ENSEMBL == unit])
  }



  dt_sub[, start := as.numeric(start)]
  dt_sub[, end := as.numeric(end)]

  # 避免 gtex data 位置有 na
  dt_sub <- dt_sub[!is.na(start) & !is.na(end)]
  if (nrow(dt_sub) == 0) {
    return(data.table())
  }

  query_dt <- data.table(
    query_start = query_start,
    query_end = query_end
  )

  query_dt[, query_start := as.numeric(query_start)]
  query_dt[, query_end := as.numeric(query_end)]

  # 我們資料的位置 query_start(因為基因 pos 會有很多組) 跟 gtex data 位置 tmp 比距離
  result <- rbindlist(lapply(seq_len(nrow(query_dt)), function(i) {
    tmp <- copy(dt_sub)
    tmp <- tmp[!is.na(start)]

    if (nrow(tmp) == 0) {
      return(data.table())
    }

    qs <- query_dt$query_start[i]
    qe <- query_dt$query_end[i]

    # 區間距離：有重疊 distance = 0
    tmp[, distance := fifelse(
      end < qs,
      qs - end,
      fifelse(
        start > qe,
        start - qe,
        0
      )
    )]

    # 重疊長度：沒有重疊 overlap_len = 0
    tmp[, overlap_start := pmax(start, qs)]
    tmp[, overlap_end := pmin(end, qe)]
    tmp[, overlap_len := pmax(0, overlap_end - overlap_start + 1)]

    # 中心點距離
    tmp[, query_mid := (qs + qe) / 2]
    tmp[, interval_mid := (start + end) / 2]
    tmp[, mid_distance := abs(interval_mid - query_mid)]

    tmp[, `:=`(
      query_start = qs,
      query_end = qe,
      query_index = i
    )]

    # 排序規則：
    # 1. distance 最小
    # 2. overlap_len 最大
    # 3. mid_distance 最小
    tmp[order(distance, -overlap_len, mid_distance)][1]
  }))

  if (nrow(result) == 0) {
    return(data.table())
  }
  # 多組 query interval 最後再挑一次最佳結果
  result[order(distance, -overlap_len, mid_distance)][1]
}



density_plot <- function(df, value_col, main, x_range = c(0, 1)) {
  all_types <- unique(df$type)
  random_types <- all_types[all_types != "eQTL"]

  color_values <- c(
    "eQTL" = "red",
    setNames(rep("gray70", length(random_types)), random_types)
  )


  ggplot(df, aes(x = (.data[[value_col]]), color = type, fill = type)) +
    geom_density(
      kernel = "gaussian", adjust = 0.5, alpha = 0.1, n = 10^6
    ) +
    scale_color_manual(values = color_values) +
    scale_fill_manual(values = color_values) +
    labs(
      title = main,
      x = "P-value",
      y = "Density",
      color = "Group",
      fill = "Group"
    ) +
    theme_minimal() +
    coord_cartesian(xlim = x_range) +
    theme(legend.position = "none")
}



# fusion_weights ENSG -> gene ----
# install.packages("BiocManager")
# BiocManager::install("mygene")
library(mygene)

# ft
gene_enrichment <- function(
  tissue_name,
  weight_path,
  output_dir
) {
  gtex_weight <- fread(weight_path)
  gtex_weight[, query := sub("\\..*", "", Gene_ID)]

  gene_anno <- fread("D:/Peter/gene_enrichment/code_project/data/gencode_v26_gene_anno_hg18.txt")

  # gene 名稱併入fusion weight
  gtex_weight <- gene_anno[gtex_weight, on = .(query), nomatch = 0]

  # 加上 ensg 出現次數
  query_times <- gtex_weight[, .(
    n_ensg = uniqueN(query)
  ), by = symbol]
  gtex_weight <- gtex_weight[query_times, on = .(symbol), nomatch = 0]


  fwrite(gtex_weight,
    file.path(output_dir, "manyENSG_to_oneGENE.txt"),
    row.names = F, col.names = T, sep = "\t"
  )
  setnames(gtex_weight, old = "symbol", new = "gene")

  gtex_weight <- gtex_weight[!is.na(gene) & !is.na(query), ]
  cat(
    tissue_name, "raw data genes and ENSG genes (remove NA):",
    uniqueN(gtex_weight$gene),
    uniqueN(gtex_weight$query),
    "\n"
  )

  gtex_weight <- gtex_weight[!is.na(gene) & !is.na(Top1_Pval), ]
  cat(
    tissue_name, "After remove NA Top1_Pval, genes and ENSg gene:",
    uniqueN(gtex_weight$gene), uniqueN(gtex_weight$query),
    "\n"
  )


  # 我們資料會有 一個 gene 多個 probe, 區間的情況，把 gtex_weight 跟我們資料每組 gene 區間算距離後，取出最近的 gtex_weight，
  # 再從中選出最近的。 例如，我們資料一個 gene 有 3 區間，gtex_weight gene 重複2次，算3*2組距離，找出離我們資料任一區間，
  # 最近的 gtex_weight gene
  probe_pos <- fread("D:/Peter/oral_cancer/expression/expression_data/probe_pos_mulInterval.txt")
  gene_probe <- fread("D:/Peter/oral_cancer/expression/expression_data/probe24526infoB.txt")[
    ,
    .(TargetID, PROBE_ID)
  ]
  probe_pos <- gene_probe[probe_pos, on = .(PROBE_ID = Gene), nomatch = 0]
  repeat_ensg <- gtex_weight[n_ensg > 1 & gene %in% probe_pos$TargetID, gene] %>%
    unique()


  repeat_row <- list()
  for (i in seq_along(repeat_ensg)) {
    repeat_row[[i]] <- find_nearest_interval(
      dt = gtex_weight,
      unit = repeat_ensg[i],
      query_start = probe_pos[TargetID == repeat_ensg[i], start],
      query_end = probe_pos[TargetID == repeat_ensg[i], end],
      for_gene = T
    )
  }

  repeat_row_result <- rbindlist(repeat_row) %>% as.data.table()
  col_order <- names(gtex_weight)
  repeat_row_result <- repeat_row_result[, ..col_order]

  # repeat gene choose one
  # 處理出現在 probe_pos 的 repeat ensg 就好，其他的不重要，因為後續跟我們資料 取交集，會被刪掉
  repeat_row_result[, choose := 1]
  col_order <- names(gtex_weight)
  gtex_weight <- repeat_row_result[gtex_weight, on = col_order]

  gtex_weight[is.na(choose), choose := 2]
  setkey(gtex_weight, choose)
  gtex_weight <- gtex_weight[, .SD[1], by = gene]

  cat(
    "After choose most nearest gene for repeat gene, genes and ENSG genes:",
    uniqueN(gtex_weight$gene), uniqueN(gtex_weight$query),
    "\n"
  )

  gt_N <- fread("D:/Peter/rawData_eQTL/r2_filter_0.8/outcome/raw_N_FDR_R2_0.8.txt")
  
  all_eQTL_gene <- unique(all_eQTL$Gene)
  eQTL_gene <- unique(all_eQTL[Gene %in% gt_N$Gene, Gene])

  gtex_sigGene <- gtex_weight[Top1_Pval < 0.01, gene] %>%
    unique()

  cat(
    tissue_name, "sig. genes:",
    uniqueN(gtex_sigGene),
    "\n"
  )
  cat(
    "eQTL sig. genes:",
    uniqueN(eQTL_gene),
    "\n"
  )

  common_gene <- intersect(all_eQTL_gene, gtex_weight$gene)
  sample_gene_number <- uniqueN(eQTL_gene[eQTL_gene %in% gtex_weight$gene])
  observed_intersect <- intersect(eQTL_gene, gtex_sigGene) %>% uniqueN()

  cat(
    "intersect of eQTL,", tissue_name, "total genes:",
    intersect(all_eQTL_gene, gtex_weight$gene) %>% uniqueN(), "\n"
  )
  cat(
    tissue_name, "significant genes recorded in eQTL:",
    gtex_sigGene[gtex_sigGene %in% all_eQTL_gene] %>% uniqueN(), "\n"
  )
  cat(
    "eQTL significant genes recorded in", tissue_name, ":",
    uniqueN(eQTL_gene[eQTL_gene %in% gtex_weight$gene]), "\n"
  )
  cat(
    "intersect of eQTL,", tissue_name, "sig genes:",
    observed_intersect, "\n"
  )

  repeat_gene_number <- c()
  random_times <- 1e4
  for (i in 1:random_times) {
    set.seed(i)
    repeat_gene_number[i] <- (sample(
      common_gene,
      sample_gene_number
    ) %in% gtex_sigGene[gtex_sigGene %in% all_eQTL_gene]) %>%
      which() %>%
      length()
  }

  cat(
    "sample", sample_gene_number, "genes, mean sig in", tissue_name, "genes number:",
    mean(repeat_gene_number), "\n"
  )
  cat(
    "sample", sample_gene_number, "genes, sig in", tissue_name, "genes number range:",
    range(repeat_gene_number), "\n"
  )
  cat(
    "sample", sample_gene_number, "genes, sig in", tissue_name, "genes number >=",
    observed_intersect, "times:",
    length(which(repeat_gene_number >= observed_intersect)), "\n"
  )

  png(file.path(output_dir, sprintf("random_%sgene.png", sample_gene_number)),
    width = 8, height = 6, units = "in", res = 300
  )
  hist(repeat_gene_number,
    main = sprintf("%s times Gene Intersection Size between %s Sig. and eQTL Sig.",random_times, tissue_name),
    xlab = "Intersection Size each Times",
    col = "skyblue",
    breaks = c(min(repeat_gene_number):max(repeat_gene_number))
  )
  abline(v = observed_intersect, col = "red", lty = 1)
  dev.off()

  a <- gtex_weight[gene %in% eQTL_gene, ][, .(gene_name = gene, pval = Top1_Pval)]
  setkey(a, pval)
  a <- a[, .SD[1], by = gene_name]
  a[, type := "eQTL"]




  # 畫 5次隨機 dist
  random_list <- list()
  random_times <- 5
  for (i in 1:random_times) {
    set.seed(i)
    b_tmp <- gtex_weight[gene %in% sample(common_gene, sample_gene_number), ][, .(gene_name = gene, pval = Top1_Pval)]
    b_tmp[, type := paste0("random_", i)]
    random_list[[i]] <- b_tmp
  }

  all_random_b <- rbindlist(random_list)
  combine_all <- rbind(a[, .(gene_name, pval, type)], all_random_b)
  all_pval <- combine_all[, pval]
  x_auto <- quantile(all_pval, probs = 0.8, na.rm = TRUE)
  x_auto <- max(x_auto, 1e-4)
  p_start <- ceiling(log10(x_auto))
  pow10_large <- 10^seq(0, p_start, by = -1)

  for (x_cutoff in c(pow10_large, x_auto)) {
    png(file.path(output_dir, sprintf("%sgene_Top1_PvalDIST_%s.png", sample_gene_number, format(x_cutoff, scientific = TRUE))),
      width = 8, height = 6, units = "in", res = 300
    )
    print(
      density_plot(
        combine_all, "pval",
        sprintf("Density Plot for %s Gene Top1_Pval vs %s Random", sample_gene_number, random_times),
        x_range = c(0, x_cutoff)
      )
    )
    dev.off()
  }




  # 畫 1e6 次隨機 dist的平均
  random_times <- 1e4
  random_list <- vector("list", random_times)
  for (i in 1:random_times) {
    set.seed(i)
    b_tmp <- gtex_weight[gene %in% sample(common_gene, sample_gene_number), ][, .(gene_name = gene, pval = Top1_Pval)]
    b_tmp[, type := paste0("random_", i)]
    random_list[[i]] <- b_tmp
  }

  all_random_b <- rbindlist(random_list)
  combine_all <- rbind(a[, .(gene_name, pval, type)], all_random_b)
  all_pval <- combine_all[, pval]
  x_auto <- quantile(all_pval, probs = 0.8, na.rm = TRUE)
  x_auto <- max(x_auto, 1e-4)
  p_start <- ceiling(log10(x_auto))
  pow10_large <- 10^seq(0, p_start, by = -1)


  plot_x_range <- c(0, 1)
  random_density <- combine_all[type != "eQTL",
    {
      density_fit <- density(pval, from = plot_x_range[1], to = plot_x_range[2], n = 2^15)
      .(x = density_fit$x, density = density_fit$y)
    },
    by = type
  ]

  random_density_summary <- random_density[, .(
    mean_density = mean(density),
    ci_lower = mean(density) - qt(0.975, .N - 1) * sd(density) / sqrt(.N),
    ci_upper = mean(density) + qt(0.975, .N - 1) * sd(density) / sqrt(.N),
    q025 = quantile(density, 0.025),
    q975 = quantile(density, 0.975)
  ), by = x]

  eqtl_density_fit <- density(
    combine_all[type == "eQTL", pval],
    from = plot_x_range[1],
    to = plot_x_range[2],
    n = 2^15
  )
  eqtl_density <- data.table(x = eqtl_density_fit$x, density = eqtl_density_fit$y)

  all_pval <- combine_all[, pval]
  x_auto <- quantile(all_pval, probs = 0.9, na.rm = TRUE)
  x_auto <- max(x_auto, 1e-4)
  p_start <- ceiling(log10(x_auto))
  pow10_large <- 10^seq(0, p_start, by = -1)

  for (x_cutoff in c(pow10_large, x_auto)) {
    png(file.path(output_dir,
     sprintf("%sgene_confi_xlim_%s.png",
     sample_gene_number,
      format(x_cutoff, scientific = TRUE))),
      width = 8, height = 6, units = "in", res = 300
    )
    print(
      ggplot() +
        geom_ribbon(
          data = random_density_summary,
          aes(x = x, ymin = ci_lower, ymax = ci_upper),
          fill = "blue",
          alpha = 0.35
        ) +
        geom_ribbon(
          data = random_density_summary,
          aes(x = x, ymin = q025, ymax = q975),
          fill = "gray70",
          alpha = 0.35
        ) +
        geom_line(
          data = random_density_summary,
          aes(x = x, y = mean_density, color = "Random mean"),
          linewidth = 0.9
        ) +
        geom_line(
          data = eqtl_density,
          aes(x = x, y = density, color = "eQTL"),
          linewidth = 0.9
        ) +
        scale_color_manual(values = c("eQTL" = "red", "Random mean" = "gray40")) +
        coord_cartesian(xlim = c(0, x_cutoff)) +
        labs(
          title = sprintf("Mean Density Plot for %s Gene Top1_Pval vs %s Random", sample_gene_number, random_times),
          x = "p-value",
          y = "Density",
          color = "Group"
        ) +
        theme_minimal()
    )
    dev.off()
  }
}








# 跑多種組織需要的 csv ----
#  Rscript /mnt/d/Peter/gene_enrichment/code_project/data/summarize_weights_LOTStissue.R



base_dir <- "D:/Peter/gene_enrichment/code_project/data"
out_dir <- "D:/Peter/gene_enrichment/code_project/outcome/rawData_eQTL_N/for_gene"



gtex <- data.frame(
  tissue_dir = c("GTEx_salivary",
   "GTEx_esophagus", "GTEx_thyroid", "GTEx_lung",
    "GTEx_Adipose_Subcutaneous",
    "GTEx_Brain_Cerebellum",
    "GTEx_Esophagus_Gastroesophageal",
    "GTEx_Esophagus_Muscularis",
    "GTEx_Heart_Left",
    "GTEx_Muscle_Skeletal"
  )
)

output_csv <- file.path(
  base_dir,
  gtex$tissue_dir,
  "fusion_project",
  "fusion_weights_summary.csv"
)

plot_output_dir <- file.path(
  out_dir,
  gtex$tissue_dir
)


# 沒有則新建資料夾
sapply(plot_output_dir, function(x) {
  if (!dir.exists(x)) dir.create(x, recursive = TRUE)
})





# 執行 ----
all_eQTL <- fread("D:/Peter/rawData_eQTL/trash/raw_maf_gt_N_pvalue.txt", header = T)[
  ,
  .(pval = `p-value`, SNP, probe = gene)
]
R2_filter <- fread("D:/Peter/oral_cancer/expression/outcome/multiple_nucleotide_variant/cis_snp_R2_0.8.txt", header = T)

setkey(all_eQTL, pval)
all_eQTL <- all_eQTL[SNP %in% R2_filter$hg18_snpID, ]
all_eQTL <- all_eQTL[, .SD[1], by = probe]
all_eQTL[, probe := gsub("^([^_]+_[^_]+).*", "\\1", probe)]

probe_info <- fread("D:/Peter/rawData_eQTL/r2_filter_0.8/outcome/raw_exp_different_r2_0.8.txt")[
  ,
  .(probe = PROBE_ID, Gene)
]
all_eQTL <- probe_info[all_eQTL, on = .(probe), nomatch = 0]
all_eQTL_gene <- unique(all_eQTL$Gene)
all_eQTL_probe <- unique(all_eQTL$probe)





for (i in seq_along(gtex$tissue_dir)) {
  gene_enrichment(
    tissue_name = gtex$tissue_dir[i],
    weight_path = output_csv[i],
    output_dir = plot_output_dir[i]
  )
}


# GTEx-salivary raw data genes and ENSG genes (remove NA): 23814 23847
# GTEx-salivary After remove NA Top1_Pval, genes and ENSg gene: 23814 23847
# After choose most nearest gene for repeat gene, genes and ENSG genes: 23814 23814
# GTEx-salivary sig. genes: 4282
# eQTL sig. genes: 177
# intersect of eQTL, GTEx-salivary total genes: 12059
# GTEx-salivary significant genes recorded in eQTL: 2058
# eQTL significant genes recorded in GTEx-salivary : 136
# intersect of eQTL, GTEx-salivary sig genes: 42
# sample 136 genes, mean sig in GTEx-salivary genes number: 23.2556
# sample 136 genes, sig in GTEx-salivary genes number range: 9 42
# sample 136 genes, sig in GTEx-salivary genes number >= 42 times: 1

# GTEx-esophagus raw data genes and ENSG genes (remove NA): 22588 22607
# GTEx-esophagus After remove NA Top1_Pval, genes and ENSg gene: 22588 22607
# After choose most nearest gene for repeat gene, genes and ENSG genes: 22588 22588
# GTEx-esophagus sig. genes: 10347
# eQTL sig. genes: 177
# intersect of eQTL, GTEx-esophagus total genes: 12014
# GTEx-esophagus significant genes recorded in eQTL: 5903
# eQTL significant genes recorded in GTEx-esophagus : 140
# intersect of eQTL, GTEx-esophagus sig genes: 94
# sample 140 genes, mean sig in GTEx-esophagus genes number: 68.8856
# sample 140 genes, sig in GTEx-esophagus genes number range: 44 94
# sample 140 genes, sig in GTEx-esophagus genes number >= 94 times: 1

# GTEx-thyroid raw data genes and ENSG genes (remove NA): 24722 24767
# GTEx-thyroid After remove NA Top1_Pval, genes and ENSg gene: 24722 24767
# After choose most nearest gene for repeat gene, genes and ENSG genes: 24722 24722
# GTEx-thyroid sig. genes: 12640
# eQTL sig. genes: 177
# intersect of eQTL, GTEx-thyroid total genes: 12272
# GTEx-thyroid significant genes recorded in eQTL: 6610
# eQTL significant genes recorded in GTEx-thyroid : 142
# intersect of eQTL, GTEx-thyroid sig genes: 95
# sample 142 genes, mean sig in GTEx-thyroid genes number: 76.5144
# sample 142 genes, sig in GTEx-thyroid genes number range: 55 98
# sample 142 genes, sig in GTEx-thyroid genes number >= 95 times: 11

# GTEx-lung raw data genes and ENSG genes (remove NA): 24645 24687
# GTEx-lung After remove NA Top1_Pval, genes and ENSg gene: 24645 24687
# After choose most nearest gene for repeat gene, genes and ENSG genes: 24645 24645
# GTEx-lung sig. genes: 10378
# eQTL sig. genes: 177
# intersect of eQTL, GTEx-lung total genes: 12282
# GTEx-lung significant genes recorded in eQTL: 5473
# eQTL significant genes recorded in GTEx-lung : 142
# intersect of eQTL, GTEx-lung sig genes: 82
# sample 142 genes, mean sig in GTEx-lung genes number: 63.2441
# sample 142 genes, sig in GTEx-lung genes number range: 42 87
# sample 142 genes, sig in GTEx-lung genes number >= 82 times: 5

# GTEx_Adipose_Subcutaneous raw data genes and ENSG genes (remove NA): 23361 23395
# GTEx_Adipose_Subcutaneous After remove NA Top1_Pval, genes and ENSg gene: 23361 23395
# After choose most nearest gene for repeat gene, genes and ENSG genes: 23361 23361
# GTEx_Adipose_Subcutaneous sig. genes: 10993
# eQTL sig. genes: 177
# intersect of eQTL, GTEx_Adipose_Subcutaneous total genes: 12075
# GTEx_Adipose_Subcutaneous significant genes recorded in eQTL: 5965
# eQTL significant genes recorded in GTEx_Adipose_Subcutaneous : 135
# intersect of eQTL, GTEx_Adipose_Subcutaneous sig genes: 95
# sample 135 genes, mean sig in GTEx_Adipose_Subcutaneous genes number: 66.709
# sample 135 genes, sig in GTEx_Adipose_Subcutaneous genes number range: 46 89
# sample 135 genes, sig in GTEx_Adipose_Subcutaneous genes number >= 95 times: 0

# GTEx_Brain_Cerebellum raw data genes and ENSG genes (remove NA): 23892 23922
# GTEx_Brain_Cerebellum After remove NA Top1_Pval, genes and ENSg gene: 23892 23922
# After choose most nearest gene for repeat gene, genes and ENSG genes: 23892 23892
# GTEx_Brain_Cerebellum sig. genes: 8303
# eQTL sig. genes: 177
# intersect of eQTL, GTEx_Brain_Cerebellum total genes: 11911
# GTEx_Brain_Cerebellum significant genes recorded in eQTL: 4355
# eQTL significant genes recorded in GTEx_Brain_Cerebellum : 133
# intersect of eQTL, GTEx_Brain_Cerebellum sig genes: 74
# sample 133 genes, mean sig in GTEx_Brain_Cerebellum genes number: 48.6346
# sample 133 genes, sig in GTEx_Brain_Cerebellum genes number range: 28 73
# sample 133 genes, sig in GTEx_Brain_Cerebellum genes number >= 74 times: 0

# GTEx_Esophagus_Gastroesophageal raw data genes and ENSG genes (remove NA): 22679 22703
# GTEx_Esophagus_Gastroesophageal After remove NA Top1_Pval, genes and ENSg gene: 22679 22703
# After choose most nearest gene for repeat gene, genes and ENSG genes: 22679 22679
# GTEx_Esophagus_Gastroesophageal sig. genes: 8062
# eQTL sig. genes: 177
# intersect of eQTL, GTEx_Esophagus_Gastroesophageal total genes: 11896
# GTEx_Esophagus_Gastroesophageal significant genes recorded in eQTL: 4328
# eQTL significant genes recorded in GTEx_Esophagus_Gastroesophageal : 136
# intersect of eQTL, GTEx_Esophagus_Gastroesophageal sig genes: 79
# sample 136 genes, mean sig in GTEx_Esophagus_Gastroesophageal genes number: 49.4927
# sample 136 genes, sig in GTEx_Esophagus_Gastroesophageal genes number range: 28 74
# sample 136 genes, sig in GTEx_Esophagus_Gastroesophageal genes number >= 79 times: 0

# GTEx_Esophagus_Muscularis raw data genes and ENSG genes (remove NA): 22481 22503
# GTEx_Esophagus_Muscularis After remove NA Top1_Pval, genes and ENSg gene: 22481 22503
# After choose most nearest gene for repeat gene, genes and ENSG genes: 22481 22481
# GTEx_Esophagus_Muscularis sig. genes: 10263
# eQTL sig. genes: 177
# intersect of eQTL, GTEx_Esophagus_Muscularis total genes: 11919
# GTEx_Esophagus_Muscularis significant genes recorded in eQTL: 5649
# eQTL significant genes recorded in GTEx_Esophagus_Muscularis : 136
# intersect of eQTL, GTEx_Esophagus_Muscularis sig genes: 87
# sample 136 genes, mean sig in GTEx_Esophagus_Muscularis genes number: 64.4645
# sample 136 genes, sig in GTEx_Esophagus_Muscularis genes number range: 43 84
# sample 136 genes, sig in GTEx_Esophagus_Muscularis genes number >= 87 times: 0

# GTEx_Heart_Left raw data genes and ENSG genes (remove NA): 19870 19884
# GTEx_Heart_Left After remove NA Top1_Pval, genes and ENSg gene: 19870 19884
# After choose most nearest gene for repeat gene, genes and ENSG genes: 19870 19870
# GTEx_Heart_Left sig. genes: 7374
# eQTL sig. genes: 177
# intersect of eQTL, GTEx_Heart_Left total genes: 11126
# GTEx_Heart_Left significant genes recorded in eQTL: 4269
# eQTL significant genes recorded in GTEx_Heart_Left : 131
# intersect of eQTL, GTEx_Heart_Left sig genes: 75
# sample 131 genes, mean sig in GTEx_Heart_Left genes number: 50.2148
# sample 131 genes, sig in GTEx_Heart_Left genes number range: 30 71
# sample 131 genes, sig in GTEx_Heart_Left genes number >= 75 times: 0

# GTEx_Muscle_Skeletal raw data genes and ENSG genes (remove NA): 19963 19978
# GTEx_Muscle_Skeletal After remove NA Top1_Pval, genes and ENSg gene: 19963 19978
# After choose most nearest gene for repeat gene, genes and ENSG genes: 19963 19963
# GTEx_Muscle_Skeletal sig. genes: 9673
# eQTL sig. genes: 177
# intersect of eQTL, GTEx_Muscle_Skeletal total genes: 11429
# GTEx_Muscle_Skeletal significant genes recorded in eQTL: 5725
# eQTL significant genes recorded in GTEx_Muscle_Skeletal : 133
# intersect of eQTL, GTEx_Muscle_Skeletal sig genes: 89
# sample 133 genes, mean sig in GTEx_Muscle_Skeletal genes number: 66.6538
# sample 133 genes, sig in GTEx_Muscle_Skeletal genes number range: 41 88
# sample 133 genes, sig in GTEx_Muscle_Skeletal genes number >= 89 times: 0
