rm(list = ls())
gc()


## package ----
library(ggplot2)
library(data.table)
library(dplyr)
library(stringr)
library(readxl) # read_xlsx ft


## ft ----
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


plot_output_dir <- sprintf(
  "D:/Peter/gene_enrichment/code_project/outcome/rawData_eQTL_N/for_gene/%s",
  c("lungTWAS", "HNSC", "TCGA-LUAD", "TCGA-LUSC")
)

# 刪舊檔案，刪掉資料夾內所有檔案
# sapply(plot_output_dir, function(x) {
#   if (dir.exists(x)) {
#     message("Deleting: ", x)
#     unlink(x, recursive = TRUE, force = TRUE)
#   } else {
#     message("Directory does not exist sfiles: ", x)
#   }
# })

sapply(plot_output_dir, function(x) {
  if (!dir.exists(x)) dir.create(x, recursive = TRUE)
})

# 取出刪掉 R2<0.8 snp 後，仍有跑出 eQTL 的 probe
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
all_eQTL <- probe_info[all_eQTL, on = .(probe)]
all_eQTL_gene <- unique(all_eQTL$Gene)
all_eQTL_probe <- unique(all_eQTL$probe)


gt_N <- fread("D:/Peter/rawData_eQTL/r2_filter_0.8/outcome/raw_N_FDR_R2_0.8.txt")


lung_twas <- fread(
  "D:/Peter/gene_enrichment/code_project/data/lung_國衛院/20260322_115samples_24216probes_ModelPerformance.txt",
  select = c("ProbeID", "GeneSymbol", "top1_p")
)
cat("lung probe: ", uniqueN(lung_twas$ProbeID), "\n")
cat("lung gene: ", uniqueN(lung_twas$GeneSymbol), "\n")

lung_twas <- lung_twas[!is.na(top1_p), ]
cat("lung probe with pval (not NA): ", uniqueN(lung_twas$ProbeID), "\n")
cat("lung gene with pval (not NA): ", uniqueN(lung_twas$GeneSymbol), "\n")

lung_sig <- lung_twas[top1_p < 0.01]
cat("lung sig probe: ", uniqueN(lung_sig$ProbeID), "\n")
cat("lung sig gene: ", uniqueN(lung_sig$GeneSymbol), "\n")

eQTL_probe <- unique(gt_N$Probe)
eQTL_gene <- unique(gt_N$Gene)
cat("eQTL sig probe: ", uniqueN(eQTL_probe), "\n")
cat("eQTL sig gene: ", uniqueN(eQTL_gene), "\n")


# 挑出我們資料有紀錄、lung 紀錄的 probe/gene
common_probe <- intersect(all_eQTL_probe, lung_twas$ProbeID)
common_gene <- intersect(all_eQTL_gene, lung_twas$GeneSymbol)

cat("intersect of eQTL, lung total probe:", intersect(all_eQTL_probe, lung_twas$ProbeID) %>% uniqueN(), "\n")
cat("intersect of eQTL, lung total gene:", intersect(all_eQTL_gene, lung_twas$GeneSymbol) %>% uniqueN(), "\n")

# 在一邊顯著，同時有被紀錄在另一邊 數量
cat("eQTL significant probes recorded in lung TWAS:", uniqueN(gt_N$Probe[gt_N$Probe %in% lung_twas$ProbeID]), "\n")
cat("eQTL significant genes recorded in lung TWAS:", uniqueN(gt_N$Gene[gt_N$Gene %in% lung_twas$GeneSymbol]), "\n")
cat("lung significant probes recorded in eQTL:", lung_sig$ProbeID[lung_sig$ProbeID %in% all_eQTL_probe] %>% uniqueN(), "\n")
cat("lung significant genes recorded in eQTL:", lung_sig$GeneSymbol[lung_sig$GeneSymbol %in% all_eQTL_gene] %>% uniqueN(), "\n")


# 交集數量
observed_intersect_probe <- intersect(eQTL_probe, lung_twas[top1_p < 0.01, ProbeID]) %>% uniqueN()
observed_intersect_gene <- intersect(eQTL_gene, lung_twas[top1_p < 0.01, GeneSymbol]) %>% uniqueN()
cat("intersect of eQTL, lung sig probe:", observed_intersect_probe, "\n")
cat("intersect of eQTL, lung sig gene:", observed_intersect_gene, "\n")

sample_probe_number <- uniqueN(gt_N$Probe[gt_N$Probe %in% lung_twas$ProbeID])
sample_gene_number <- uniqueN(gt_N$Gene[gt_N$Gene %in% lung_twas$GeneSymbol])







## 選跟 eQTL sig. 同數量的，跟 lung_twas 交集數量 hist ----
# 隨機選 105 probe, 91 gene，多少顯著 in lung?
common_idx_probe <- match(common_probe, lung_sig$ProbeID[lung_sig$ProbeID %in% all_eQTL_probe])
common_idx_gene <- match(common_gene, lung_sig$GeneSymbol[lung_sig$GeneSymbol %in% all_eQTL_gene])

random_times <- 1e6
repeat_probe_number <- c()
repeat_gene_number <- c()
for (i in 1:random_times) {
  set.seed(i)
  idx <- sample(common_idx_probe, sample_probe_number)
  repeat_probe_number[i] <- idx[!is.na(idx)] %>% length()

  idx <- sample(common_idx_gene, sample_gene_number)
  repeat_gene_number[i] <- idx[!is.na(idx)] %>% length()
}

cat("抽", sample_probe_number, "個，平均顯著 in lung probe number:", mean(repeat_probe_number), "\n")
cat("抽", sample_gene_number, "個，平均顯著 in lung gene number:", mean(repeat_gene_number), "\n")
cat("抽", sample_probe_number, "個，顯著 in lung probe number 範圍:", range(repeat_probe_number), "\n")
cat("抽", sample_gene_number, "個，顯著 in lung gene number 範圍:", range(repeat_gene_number), "\n")

cat(
  "sample", sample_probe_number, "genes, sig in lung probe number >=",
  observed_intersect_probe, "times:",
  length(which(repeat_probe_number >= observed_intersect_probe)), "\n"
)
cat(
  "sample", sample_gene_number, "genes, sig in lung gene number >=",
  observed_intersect_gene, "times:",
  length(which(repeat_gene_number >= observed_intersect_gene)), "\n"
)

png(file.path(plot_output_dir[1], sprintf("random_%sprobe.png", sample_probe_number)),
  width = 8, height = 6, units = "in", res = 300
)
hist(repeat_probe_number,
  main = c("1e4 times Probe Intersection Size between Lung Sig. and eQTL Sig."),
  xlab = "Intersection Size each Times",
  col = "skyblue"
)
dev.off()

png(file.path(plot_output_dir[1], sprintf("random_%sgene.png", sample_gene_number)),
  width = 8, height = 6, units = "in", res = 300
)
hist(repeat_gene_number,
  main = c("1e4 times Gene Intersection Size between Lung Sig. and eQTL Sig."),
  xlab = "Intersection Size each Times",
  col = "skyblue"
)
dev.off()



# 對 probe ----
# 挑出 eQTL_probe 在 lung_twas 的 pval
setorder(gt_N, `p-value`)
a <- lung_twas[ProbeID %in% eQTL_probe, ][, .(ProbeID, pval = top1_p)]
a[, type := "eQTL"]
random_times <- 5
random_list <- list()




for (i in 1:random_times) {
  set.seed(i)
  # 抽樣並選取欄位
  b_tmp <- lung_twas[ProbeID %in% sample(common_probe, sample_probe_number), ][, .(ProbeID, pval = top1_p)]

  # 給予獨特的標籤，例如 "random_1", "random_2" ...
  b_tmp[, type := paste0("random_", i)]

  random_list[[i]] <- b_tmp
}

all_random_b <- rbindlist(random_list)
combine_all <- rbind(a[, .(ProbeID, pval, type)], all_random_b)



png(file.path(plot_output_dir[1], sprintf("%sprobe_pval_frequency.png", sample_probe_number)),
  width = 8, height = 6, units = "in", res = 300
)
print(
  density_plot(
    combine_all, "pval",
    sprintf("Density Plot for %s Probe pval vs 1e4 Random", sample_probe_number),
    x_range = c(0, 1)
  )
)
dev.off()



# probe: 取平均畫 dist，加上 quantile band  ----
lung_probe_small <- lung_twas[
  ProbeID %in% common_probe & !is.na(top1_p),
  .(ProbeID, pval = top1_p)
]
setorder(lung_probe_small, pval)
lung_probe_small <- lung_probe_small[, .SD[1], by = ProbeID]
common_idx <- match(common_probe, lung_probe_small$ProbeID)
common_idx <- common_idx[!is.na(common_idx)]
random_times <- 1e6


plot_x_range <- c(0, 1)
# 不同 type 的 row分組，每組的 pdf 從 from 到 to 切成n份，用 預設 gaussian kernel 估計 pdf
density_n <- 2^12
quantile_sample_times <- min(random_times, 10000L)
density_x <- NULL
mean_density <- NULL
m2_density <- NULL
density_quantile_sample <- matrix(NA_real_, nrow = density_n, ncol = quantile_sample_times)

for (i in seq_len(random_times)) {
  set.seed(i)
  idx <- sample(common_idx, sample_probe_number)
  pval_tmp <- lung_probe_small[idx, pval]

  density_fit <- density(
    pval_tmp,
    from = plot_x_range[1],
    to = plot_x_range[2],
    n = density_n
  )
  y <- density_fit$y

  if (is.null(mean_density)) {
    density_x <- density_fit$x
    mean_density <- y
    m2_density <- numeric(length(y))
  } else {
    delta <- y - mean_density
    mean_density <- mean_density + delta / i
    m2_density <- m2_density + delta * (y - mean_density)
  }

  if (i <= quantile_sample_times) {
    density_quantile_sample[, i] <- y
  } else {
    replace_idx <- sample.int(i, 1)
    if (replace_idx <= quantile_sample_times) {
      density_quantile_sample[, replace_idx] <- y
    }
  }
}

density_quantile <- t(apply(
  density_quantile_sample,
  1,
  quantile,
  probs = c(0.025, 0.975),
  na.rm = TRUE
))
sd_density <- sqrt(m2_density / (random_times - 1))
se_density <- sd_density / sqrt(random_times)
random_density_summary <- data.table(
  x = density_x,
  mean_density = mean_density,
  ci_lower = mean_density - qt(0.975, random_times - 1) * se_density,
  ci_upper = mean_density + qt(0.975, random_times - 1) * se_density,
  q025 = density_quantile[, 1],
  q975 = density_quantile[, 2]
)



eqtl_density_fit <- density(
  a[, pval],
  from = plot_x_range[1],
  to = plot_x_range[2],
  n = density_n
)
eqtl_density <- data.table(x = eqtl_density_fit$x, density = eqtl_density_fit$y)

all_pval <- c(a[, pval], lung_probe_small[common_idx, pval])
x_auto <- quantile(all_pval, probs = 0.9, na.rm = TRUE)
x_auto <- max(x_auto, 1e-4)
p_start <- ceiling(log10(x_auto))
pow10_large <- 10^seq(0, p_start, by = -1)

for (x_cutoff in c(pow10_large, x_auto)) {
  png(file.path(plot_output_dir[1], sprintf("%sprobe_confi_xlim_%s.png", sample_probe_number, format(x_cutoff, digits = 2, scientific = TRUE))),
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
        aes(x = x, y = density, color = "TWAS"),
        linewidth = 0.9
      ) +
      scale_color_manual(values = c("TWAS" = "red", "Random mean" = "gray40")) +
      coord_cartesian(xlim = c(0, x_cutoff)) +
      labs(
        title = sprintf("Mean Density Plot for %s probe pval vs %s Random", sample_probe_number, random_times),
        x = "P-value",
        y = "Density",
        color = "Group"
      ) +
      theme_minimal()
  )
  dev.off()
}


# 對 gene ----
setorder(gt_N, `p-value`)
a <- lung_twas[GeneSymbol %in% eQTL_gene, ][, .(GeneSymbol, pval = top1_p)]
a[, type := "eQTL"]

# gene 重複出現，取最顯著的
setkey(a, pval)
a <- a[, .SD[1], by = GeneSymbol]
random_times <- 5
random_list <- list()

for (i in 1:random_times) {
  set.seed(i)
  b_tmp <- lung_twas[GeneSymbol %in% sample(common_gene, sample_gene_number), ][, .(GeneSymbol, pval = top1_p)]
  b_tmp[, type := paste0("random_", i)]
  random_list[[i]] <- b_tmp
}

all_random_b <- rbindlist(random_list)
combine_all <- rbind(a[, .(GeneSymbol, pval, type)], all_random_b)


png(file.path(plot_output_dir[1], sprintf("%sgene_pval_frequency.png", sample_gene_number)),
  width = 8, height = 6, units = "in", res = 300
)
print(
  density_plot(combine_all, "pval",
    sprintf("Density Plot for %s Gene pval vs 1e4 Random", sample_gene_number),
    x_range = c(0, 1)
  )
)
dev.off()




# gene: 取平均畫 dist，加上 quantile band  ----
lung_gene_small <- lung_twas[
  GeneSymbol %in% common_gene & !is.na(top1_p),
  .(GeneSymbol, pval = top1_p)
]
setorder(lung_gene_small, pval)
lung_gene_small <- lung_gene_small[, .SD[1], by = GeneSymbol]
common_idx <- match(common_gene, lung_gene_small$GeneSymbol)
common_idx <- common_idx[!is.na(common_idx)]
random_times <- 1e6


plot_x_range <- c(0, 1)
# 不同 type 的 row分組，每組的 pdf 從 from 到 to 切成n份，用 預設 gaussian kernel 估計 pdf
density_n <- 2^12
quantile_sample_times <- min(random_times, 10000L)
density_x <- NULL
mean_density <- NULL
m2_density <- NULL
density_quantile_sample <- matrix(NA_real_, nrow = density_n, ncol = quantile_sample_times)

for (i in seq_len(random_times)) {
  set.seed(i)
  idx <- sample(common_idx, sample_gene_number)
  pval_tmp <- lung_gene_small[idx, pval]

  density_fit <- density(
    pval_tmp,
    from = plot_x_range[1],
    to = plot_x_range[2],
    n = density_n
  )
  y <- density_fit$y

  if (is.null(mean_density)) {
    density_x <- density_fit$x
    mean_density <- y
    m2_density <- numeric(length(y))
  } else {
    delta <- y - mean_density
    mean_density <- mean_density + delta / i
    m2_density <- m2_density + delta * (y - mean_density)
  }

  if (i <= quantile_sample_times) {
    density_quantile_sample[, i] <- y
  } else {
    replace_idx <- sample.int(i, 1)
    if (replace_idx <= quantile_sample_times) {
      density_quantile_sample[, replace_idx] <- y
    }
  }
}

density_quantile <- t(apply(
  density_quantile_sample,
  1,
  quantile,
  probs = c(0.025, 0.975),
  na.rm = TRUE
))
sd_density <- sqrt(m2_density / (random_times - 1))
se_density <- sd_density / sqrt(random_times)
random_density_summary <- data.table(
  x = density_x,
  mean_density = mean_density,
  ci_lower = mean_density - qt(0.975, random_times - 1) * se_density,
  ci_upper = mean_density + qt(0.975, random_times - 1) * se_density,
  q025 = density_quantile[, 1],
  q975 = density_quantile[, 2]
)



eqtl_density_fit <- density(
  a[, pval],
  from = plot_x_range[1],
  to = plot_x_range[2],
  n = density_n
)
eqtl_density <- data.table(x = eqtl_density_fit$x, density = eqtl_density_fit$y)

all_pval <- c(a[, pval], lung_gene_small[common_idx, pval])
x_auto <- quantile(all_pval, probs = 0.9, na.rm = TRUE)
x_auto <- max(x_auto, 1e-4)
p_start <- ceiling(log10(x_auto))
pow10_large <- 10^seq(0, p_start, by = -1)

for (x_cutoff in c(pow10_large, x_auto)) {
  png(file.path(plot_output_dir[1], sprintf("%sgene_confi_xlim_%s.png", sample_gene_number, format(x_cutoff, digits = 2, scientific = TRUE))),
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
        aes(x = x, y = density, color = "TWAS"),
        linewidth = 0.9
      ) +
      scale_color_manual(values = c("TWAS" = "red", "Random mean" = "gray40")) +
      coord_cartesian(xlim = c(0, x_cutoff)) +
      labs(
        title = sprintf("Mean Density Plot for %s Gene pval vs %s Random", sample_gene_number, random_times),
        x = "P-value",
        y = "Density",
        color = "Group"
      ) +
      theme_minimal()
  )
  dev.off()
}



## 整理 HNSC data ----
# 建立腳本
# nano /mnt/d/Peter/gene_enrichment/code_project/data/TCGA-HNSC.TUMOR/fusion_project/summarize_weights.R

# pkgs <- c("parallel")

# for (p in pkgs) {
#     if (!requireNamespace(p, quietly = TRUE))
#         install.packages(p)
#     suppressMessages(library(p, character.only = TRUE))
# }

# # 設定資料夾
# wgt_dir <- "/mnt/d/Peter/gene_enrichment/code_project/data/TCGA-HNSC.TUMOR/TCGA-HNSC.TUMOR"
# output_csv <- "/mnt/d/Peter/gene_enrichment/code_project/data/TCGA-HNSC.TUMOR/fusion_project/fusion_weights_summary.csv"

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
#     # 例如：TCGA-HNSC.TUMOR.ABCA8_10351.wgt.RDat -> ABCA8_10351
#     fname <- basename(f)
#     gene_id <- gsub("^TCGA-HNSC\\.TUMOR\\.|\\.wgt\\.RDat$", "", fname)

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

# 執行
# Rscript  /mnt/d/Peter/gene_enrichment/code_project/data/TCGA-HNSC.TUMOR/fusion_project/summarize_weights.R


# 建立腳本 ----
# nano /mnt/d/Peter/gene_enrichment/code_project/data/TCGA-LUSC.TUMOR/fusion_project/summarize_weights.R


# pkgs <- c("parallel")

# for (p in pkgs) {
#     if (!requireNamespace(p, quietly = TRUE))
#         install.packages(p)
#     suppressMessages(library(p, character.only = TRUE))
# }

# # 設定資料夾
# wgt_dir <- "/mnt/d/Peter/gene_enrichment/code_project/data/TCGA-LUSC.TUMOR/TCGA-LUSC.TUMOR"

# output_csv <- "/mnt/d/Peter/gene_enrichment/code_project/data/TCGA-LUSC.TUMOR/fusion_project/fusion_weights_summary.csv"


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
#     # 例如：TCGA-LUSC.TUMOR.ABCA8_10351.wgt.RDat -> ABCA8_10351
#     fname <- basename(f)
#     gene_id <- gsub("^TCGA-LUSC\\.TUMOR\\.|\\.wgt\\.RDat$", "", fname)

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


# # 按以下按鍵，存檔後退出
# # Ctrl+O, Enter, Ctrl+X


# # 執行
# # Rscript /mnt/d/Peter/gene_enrichment/code_project/data/TCGA-LUSC.TUMOR/fusion_project/summarize_weights.R


# 建立腳本 ----
# # nano /mnt/d/Peter/gene_enrichment/code_project/data/TCGA-LUAD.TUMOR/fusion_project/summarize_weights.R


# pkgs <- c("parallel")

# for (p in pkgs) {
#     if (!requireNamespace(p, quietly = TRUE))
#         install.packages(p)
#     suppressMessages(library(p, character.only = TRUE))
# }

# # 設定資料夾
# wgt_dir <- "/mnt/d/Peter/gene_enrichment/code_project/data/TCGA-LUAD.TUMOR/TCGA-LUAD.TUMOR"
# output_csv <- "/mnt/d/Peter/gene_enrichment/code_project/data/TCGA-LUAD.TUMOR/fusion_project/fusion_weights_summary.csv"


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
#     # 例如：TCGA-LUAD.TUMOR.ABCA8_10351.wgt.RDat -> ABCA8_10351
#     fname <- basename(f)
#     gene_id <- gsub("^TCGA-LUAD\\.TUMOR\\.|\\.wgt\\.RDat$", "", fname)

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


# # 按以下按鍵，存檔後退出
# # Ctrl+O, Enter, Ctrl+X

# # 執行
# # Rscript  /mnt/d/Peter/gene_enrichment/code_project/data/TCGA-LUAD.TUMOR/fusion_project/summarize_weights.R


# 讀權重檔案 ----
# target_env <- new.env()
# # 2. 將檔案載入到這個特定環境中


# load("D:/Peter/gene_enrichment/code_project/data/GTEx_Heart_Left/GTExv8.ALL.Heart_Left_Ventricle/ENSG00000000419.12.wgt.RDat", envir = target_env)
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
  cancer_type,
  weight_path,
  pos_path,
  output_dir
) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  cat(
    "After remove R2<0.8 snp, run eQTL genes:",
    uniqueN(all_eQTL_gene),
    "\n"
  )
  cat(
    "After remove R2<0.8 snp, run eQTL probes:",
    uniqueN(all_eQTL_probe),
    "\n"
  )

  tcga_weight <- fread(weight_path)
  tcga_weight[, Gene := gsub("_.*$", "", Gene_ID)]
  cat(
    cancer_type, "genes:",
    length(tcga_weight$Gene),
    "\n"
  )


  tcga_weight <- tcga_weight[!is.na(Gene) & Gene != ""]
  tcga_weight <- tcga_weight[!is.na(Top1_Pval) & Top1_Pval != ""]
  cat(
    cancer_type, "After remove NULL and NA GeneID/Top1_Pval, genes:",
    uniqueN(tcga_weight$Gene),
    "\n"
  )

  tcga_pos <- fread(pos_path)
  gt_N <- fread("D:/Peter/rawData_eQTL/r2_filter_0.8/outcome/raw_N_FDR_R2_0.8.txt")

  eQTL_gene <- unique(gt_N$Gene)

  tcga_sigGene <- tcga_weight[Top1_Pval < 0.01, Gene] %>%
    unique()
  cat(
    cancer_type, "sig. genes:",
    uniqueN(tcga_sigGene),
    "\n"
  )
  cat(
    "eQTL sig. genes:",
    uniqueN(eQTL_gene),
    "\n"
  )


  sample_gene_number <- tcga_sigGene[tcga_sigGene %in% all_eQTL_gene] %>% uniqueN()
  observed_intersect <- intersect(eQTL_gene, tcga_sigGene) %>% uniqueN()

  cat(
    cancer_type, "significant genes recorded in eQTL:",
    tcga_sigGene[tcga_sigGene %in% all_eQTL_gene] %>% uniqueN(), "\n"
  )

  cat(
    "intersect of eQTL,", cancer_type, "sig gene:",
    observed_intersect, "\n"
  )


  common_idx <- match(all_eQTL_gene, eQTL_gene)
  repeat_gene_number <- c()
  random_times <- 1e6

  for (i in 1:random_times) {
    set.seed(i)
    idx <- sample(common_idx, sample_gene_number)
    repeat_gene_number[i] <- idx[!is.na(idx)] %>% length()
  }

  cat(
    "sample", sample_gene_number, "genes, mean sig in", cancer_type, "gene number:",
    mean(repeat_gene_number), "\n"
  )
  cat(
    "sample", sample_gene_number, "genes, sig in", cancer_type, "gene number range:",
    range(repeat_gene_number), "\n"
  )
  cat(
    "sample", sample_gene_number, "genes, sig in", cancer_type, "gene number >=",
    observed_intersect, "times:",
    length(which(repeat_gene_number >= observed_intersect)), "\n"
  )

  png(file.path(output_dir, sprintf("random_%sgene.png", sample_gene_number)),
    width = 8, height = 6, units = "in", res = 300
  )
  hist(repeat_gene_number,
    main = sprintf("1e4 times Gene Intersection Size between %s Sig. and eQTL Sig.", cancer_type),
    xlab = "Intersection Size each Times",
    col = "skyblue",
    breaks = c(min(repeat_gene_number):max(repeat_gene_number))
  )
  abline(v = observed_intersect, col = "red", lty = 1)
  dev.off()

  ## 5 次隨機 pval freq plot ----
  setorder(gt_N, `p-value`)
  a <- all_eQTL[Gene %in% tcga_sigGene, ][, .(Gene, pval)]
  a <- a[, .SD[1], by = Gene]
  a[, type := "eQTL"]

  random_list <- list()
  random_times <- 5

  for (i in 1:random_times) {
    set.seed(i)
    b_tmp <- all_eQTL[Gene %in% sample(all_eQTL_gene, sample_gene_number), ][, .(Gene, pval)]
    b_tmp[, type := paste0("random_", i)]
    random_list[[i]] <- b_tmp
  }

  all_random_b <- rbindlist(random_list)
  combine_all <- rbind(a[, .(Gene, pval, type)], all_random_b)

  for (x_cutoff in c(1, 0.1, 0.025)) {
    png(file.path(output_dir, sprintf("%sgene_pvalDIST_%s.png", sample_gene_number, x_cutoff)),
      width = 8, height = 6, units = "in", res = 300
    )
    print(
      density_plot(
        combine_all, "pval",
        sprintf("Density Plot for %s Gene pval vs 1e4 Random", sample_gene_number),
        x_range = c(0, x_cutoff)
      )
    )
    dev.off()
  }


  # 取平均畫 dist，加上 quantile band  ----
  random_times <- 1e6
  plot_x_range <- c(0, 1)
  density_n <- 2^12
  quantile_sample_times <- min(random_times, 10000L)
  eqtl_small <- all_eQTL[!is.na(Gene) & Gene != "" & !is.na(pval), .SD[1], by = Gene][, .(Gene, pval)]
  common_idx <- match(all_eQTL_gene, eqtl_small$Gene)
  common_idx <- common_idx[!is.na(common_idx)]
  density_x <- NULL
  mean_density <- NULL
  m2_density <- NULL
  density_quantile_sample <- matrix(NA_real_, nrow = density_n, ncol = quantile_sample_times)

  for (i in seq_len(random_times)) {
    set.seed(i)
    idx <- sample(common_idx, sample_gene_number)
    pval_tmp <- eqtl_small[idx, pval]

    density_fit <- density(
      pval_tmp,
      from = plot_x_range[1],
      to = plot_x_range[2],
      n = density_n
    )
    y <- density_fit$y

    if (is.null(mean_density)) {
      density_x <- density_fit$x
      mean_density <- y
      m2_density <- numeric(length(y))
    } else {
      delta <- y - mean_density
      mean_density <- mean_density + delta / i
      m2_density <- m2_density + delta * (y - mean_density)
    }

    if (i <= quantile_sample_times) {
      density_quantile_sample[, i] <- y
    } else {
      replace_idx <- sample.int(i, 1)
      if (replace_idx <= quantile_sample_times) {
        density_quantile_sample[, replace_idx] <- y
      }
    }
  }

  sd_density <- sqrt(m2_density / (random_times - 1))
  se_density <- sd_density / sqrt(random_times)
  density_quantile <- t(apply(
    density_quantile_sample,
    1,
    quantile,
    probs = c(0.025, 0.975),
    na.rm = TRUE
  ))
  random_density_summary <- data.table(
    x = density_x,
    mean_density = mean_density,
    ci_lower = mean_density - qt(0.975, random_times - 1) * se_density,
    ci_upper = mean_density + qt(0.975, random_times - 1) * se_density,
    q025 = density_quantile[, 1],
    q975 = density_quantile[, 2]
  )

  eqtl_density_fit <- density(
    a[, pval],
    from = plot_x_range[1],
    to = plot_x_range[2],
    n = density_n
  )
  eqtl_density <- data.table(x = eqtl_density_fit$x, density = eqtl_density_fit$y)

  all_pval <- c(a[, pval], eqtl_small[common_idx, pval])
  x_auto <- quantile(all_pval, probs = 0.9, na.rm = TRUE)
  x_auto <- max(x_auto, 1e-4)
  p_start <- ceiling(log10(x_auto))
  pow10_large <- 10^seq(0, p_start, by = -1)

  for (x_cutoff in c(pow10_large, x_auto)) {
    png(file.path(output_dir, sprintf("%sgene_confi_xlim_%s.png", sample_gene_number, format(x_cutoff, digits = 2, scientific = TRUE))),
      width = 8, height = 6, units = "in", res = 300
    )
    print(
      ggplot() +
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
          title = sprintf("Mean Density Plot for %s Gene pval vs %s Random", sample_gene_number, random_times),
          x = "p-value",
          y = "Density",
          color = "Group"
        ) +
        theme_minimal()
    )
    dev.off()
  }
}


# 取出刪掉 R2<0.8 snp 後，仍有跑出 eQTL 的 probe
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
all_eQTL <- probe_info[all_eQTL, on = .(probe)]
all_eQTL_gene <- unique(all_eQTL$Gene)
all_eQTL_probe <- unique(all_eQTL$probe)






# TCGA-HNSC ----
run_tcga_eqtl_enrichment(
  cancer_type = "TCGA-HNSC",
  weight_path = "D:/Peter/gene_enrichment/code_project/data/TCGA-HNSC.TUMOR/fusion_project/fusion_weights_summary.csv",
  pos_path = "D:/Peter/gene_enrichment/code_project/data/TCGA-HNSC.TUMOR/TCGA-HNSC.TUMOR.pos",
  output_dir = plot_output_dir[2]
)


# TCGA-LUAD ----
run_tcga_eqtl_enrichment(
  cancer_type = "TCGA-LUAD",
  weight_path = "D:/Peter/gene_enrichment/code_project/data/TCGA-LUAD.TUMOR/fusion_project/fusion_weights_summary.csv",
  pos_path = "D:/Peter/gene_enrichment/code_project/data/TCGA-LUAD.TUMOR/TCGA-LUAD.TUMOR.pos",
  output_dir = plot_output_dir[3]
)


# TCGA-LUSC ----
run_tcga_eqtl_enrichment(
  cancer_type = "TCGA-LUSC",
  weight_path = "D:/Peter/gene_enrichment/code_project/data/TCGA-LUSC.TUMOR/fusion_project/fusion_weights_summary.csv",
  pos_path = "D:/Peter/gene_enrichment/code_project/data/TCGA-LUSC.TUMOR/TCGA-LUSC.TUMOR.pos",
  output_dir = plot_output_dir[4]
)




# tcga 1e4 次
# After remove R2<0.8 snp, run eQTL genes: 15992
# After remove R2<0.8 snp, run eQTL probes: 20906
# TCGA-HNSC genes: 2767
# TCGA-HNSC After remove NULL and NA GeneID/Top1_Pval, genes: 2763
# TCGA-HNSC sig. genes: 2069
# eQTL sig. genes: 177
# TCGA-HNSC significant genes recorded in eQTL: 1540
# intersect of eQTL, TCGA-HNSC sig gene: 44
# sample 1540 genes, mean sig in TCGA-HNSC gene number: 17.1075
# sample 1540 genes, sig in TCGA-HNSC gene number range: 5 33
# sample 1540 genes, sig in TCGA-HNSC gene number >= 44 times: 0


# After remove R2<0.8 snp, run eQTL genes: 15992
# After remove R2<0.8 snp, run eQTL probes: 20906
# TCGA-LUAD genes: 2978
# TCGA-LUAD After remove NULL and NA GeneID/Top1_Pval, genes: 2974
# TCGA-LUAD sig. genes: 2248
# eQTL sig. genes: 177
# TCGA-LUAD significant genes recorded in eQTL: 1644
# intersect of eQTL, TCGA-LUAD sig gene: 47
# sample 1644 genes, mean sig in TCGA-LUAD gene number: 18.2626
# sample 1644 genes, sig in TCGA-LUAD gene number range: 6 35
# sample 1644 genes, sig in TCGA-LUAD gene number >= 47 times: 0


# After remove R2<0.8 snp, run eQTL genes: 15992
# After remove R2<0.8 snp, run eQTL probes: 20906
# TCGA-LUSC genes: 2548
# TCGA-LUSC After remove NULL and NA GeneID/Top1_Pval, genes: 2545
# TCGA-LUSC sig. genes: 1862
# eQTL sig. genes: 177
# TCGA-LUSC significant genes recorded in eQTL: 1366
# intersect of eQTL, TCGA-LUSC sig gene: 47
# sample 1366 genes, mean sig in TCGA-LUSC gene number: 15.1641
# sample 1366 genes, sig in TCGA-LUSC gene number range: 4 32
# sample 1366 genes, sig in TCGA-LUSC gene number >= 47 times: 0


# 1e6 次
# lung probe:  24216
# lung gene:  18103
# lung probe with pval (not NA):  22951
# lung gene with pval (not NA):  17157
# lung sig probe:  3137
# lung sig gene:  2875
# eQTL sig probe:  179
# eQTL sig gene:  177
# intersect of eQTL, lung total probe: 20776
# intersect of eQTL, lung total gene: 13708
# eQTL significant probes recorded in lung TWAS: 178
# eQTL significant genes recorded in lung TWAS: 154
# lung significant probes recorded in eQTL: 2937
# lung significant genes recorded in eQTL: 2357
# intersect of eQTL, lung sig probe: 127
# intersect of eQTL, lung sig gene: 111
# 抽 178 個，平均顯著 in lung probe number: 25.16189
# 抽 154 個，平均顯著 in lung gene number: 26.4774
# 抽 178 個，顯著 in lung probe number 範圍: 4, 52
# 抽 154 個，顯著 in lung gene number 範圍: 6, 51
# sample 178 genes, sig in lung probe number >= 127 times: 0
# sample 154 genes, sig in lung gene number >= 111 times: 0
# |--------------------------------------------------|
# |==================================================|
# |--------------------------------------------------|
# |==================================================|
# After remove R2<0.8 snp, run eQTL genes: 15992
# After remove R2<0.8 snp, run eQTL probes: 20906
# TCGA-HNSC genes: 2767
# TCGA-HNSC After remove NULL and NA GeneID/Top1_Pval, genes: 2763
# TCGA-HNSC sig. genes: 2069
# eQTL sig. genes: 177
# TCGA-HNSC significant genes recorded in eQTL: 1540
# intersect of eQTL, TCGA-HNSC sig gene: 44
# sample 1540 genes, mean sig in TCGA-HNSC gene number: 17.04804
# sample 1540 genes, sig in TCGA-HNSC gene number range: 2, 38
# sample 1540 genes, sig in TCGA-HNSC gene number >= 44 times: 0


# After remove R2<0.8 snp, run eQTL genes: 15992
# After remove R2<0.8 snp, run eQTL probes: 20906
# TCGA-LUAD genes: 2978
# TCGA-LUAD After remove NULL and NA GeneID/Top1_Pval, genes: 2974
# TCGA-LUAD sig. genes: 2248
# eQTL sig. genes: 177
# TCGA-LUAD significant genes recorded in eQTL: 1644
# intersect of eQTL, TCGA-LUAD sig gene: 47
# sample 1644 genes, mean sig in TCGA-LUAD gene number: 18.19917
# sample 1644 genes, sig in TCGA-LUAD gene number range: 2, 40
# sample 1644 genes, sig in TCGA-LUAD gene number >= 47 times: 0


# After remove R2<0.8 snp, run eQTL genes: 15992
# After remove R2<0.8 snp, run eQTL probes: 20906
# TCGA-LUSC genes: 2548
# TCGA-LUSC After remove NULL and NA GeneID/Top1_Pval, genes: 2545
# TCGA-LUSC sig. genes: 1862
# eQTL sig. genes: 177
# TCGA-LUSC significant genes recorded in eQTL: 1366
# intersect of eQTL, TCGA-LUSC sig gene: 47
# sample 1366 genes, mean sig in TCGA-LUSC gene number: 15.12285
# sample 1366 genes, sig in TCGA-LUSC gene number range: 1, 34
# sample 1366 genes, sig in TCGA-LUSC gene number >= 47 times: 0
