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
  "D:/Peter/gene_enrichment/code_project/outcome/rawData_eQTL_T/for_gene/%s",
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
all_eQTL <- fread("D:/Peter/rawData_eQTL/trash/raw_maf_gt_T_pvalue.txt", header = T)[
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


gt_T <- fread("D:/Peter/rawData_eQTL/r2_filter_0.8/outcome/raw_T_FDR_R2_0.8.txt")


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

eQTL_probe <- unique(gt_T$Probe)
eQTL_gene <- unique(gt_T$Gene)
cat("eQTL sig probe: ", uniqueN(eQTL_probe), "\n")
cat("eQTL sig gene: ", uniqueN(eQTL_gene), "\n")


# 挑出我們資料有紀錄、lung 紀錄的 probe/gene
common_probe <- intersect(all_eQTL_probe, lung_twas$ProbeID)
common_gene <- intersect(all_eQTL_gene, lung_twas$GeneSymbol)

cat("intersect of eQTL, lung total probe:", intersect(all_eQTL_probe, lung_twas$ProbeID) %>% uniqueN(), "\n")
cat("intersect of eQTL, lung total gene:", intersect(all_eQTL_gene, lung_twas$GeneSymbol) %>% uniqueN(), "\n")

# 在一邊顯著，同時有被紀錄在另一邊 數量
cat("eQTL significant probes recorded in lung TWAS:", uniqueN(gt_T$Probe[gt_T$Probe %in% lung_twas$ProbeID]), "\n")
cat("eQTL significant genes recorded in lung TWAS:", uniqueN(gt_T$Gene[gt_T$Gene %in% lung_twas$GeneSymbol]), "\n")
cat("lung significant probes recorded in eQTL:", lung_sig$ProbeID[lung_sig$ProbeID %in% all_eQTL_probe] %>% uniqueN(), "\n")
cat("lung significant genes recorded in eQTL:", lung_sig$GeneSymbol[lung_sig$GeneSymbol %in% all_eQTL_gene] %>% uniqueN(), "\n")


# 交集數量

cat("intersect of eQTL, lung sig probe:", intersect(eQTL_probe, lung_twas[top1_p < 0.01, ProbeID]) %>% uniqueN(), "\n")
cat("intersect of eQTL, lung sig gene:", intersect(eQTL_gene, lung_twas[top1_p < 0.01, GeneSymbol]) %>% uniqueN(), "\n")

sample_probe_number <- uniqueN(gt_T$Probe[gt_T$Probe %in% lung_twas$ProbeID])
sample_gene_number <- uniqueN(gt_T$Gene[gt_T$Gene %in% lung_twas$GeneSymbol])


## 選跟 eQTL sig. 同數量的，跟 lung_twas 交集數量 hist ----
# 隨機選 105 probe, 91 gene，多少顯著 in lung?
random_times <- 10000
repeat_probe_number <- c()
repeat_gene_number <- c()
for (i in 1:random_times) {
  set.seed(i)
  repeat_probe_number[i] <- (sample(
    common_probe,
    sample_probe_number
  ) %in% lung_sig$ProbeID[lung_sig$ProbeID %in% all_eQTL_probe]) %>%
    which() %>%
    length()

  repeat_gene_number[i] <- (sample(
    common_gene,
    sample_gene_number
  ) %in% lung_sig$GeneSymbol[lung_sig$GeneSymbol %in% all_eQTL_gene]) %>%
    which() %>%
    length()
}

cat("抽", sample_probe_number, "個，平均顯著 in lung probe number:", mean(repeat_probe_number), "\n")
cat("抽", sample_gene_number, "個，平均顯著 in lung gene number:", mean(repeat_gene_number), "\n")
cat("抽", sample_probe_number, "個，顯著 in lung probe number 範圍:", range(repeat_probe_number), "\n")
cat("抽", sample_gene_number, "個，顯著 in lung gene number 範圍:", range(repeat_gene_number), "\n")

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
setorder(gt_T, `p-value`)
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
random_times <- 10000
random_list <- list()

for (i in 1:random_times) {
  set.seed(i)
  b_tmp <- lung_twas[ProbeID %in% sample(common_probe, sample_probe_number), ][, .(ProbeID, pval = top1_p)]
  b_tmp[, type := paste0("random_", i)]
  random_list[[i]] <- b_tmp
}

all_random_b <- rbindlist(random_list)
combine_all <- rbind(a[, .(ProbeID, pval, type)], all_random_b)


plot_x_range <- c(0, 1)
# 不同 type 的 row分組，每組的 pdf 從 from 到 to 切成n份，用 預設 gaussian kernel 估計 pdf
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


for (x_cutoff in c(1, 0.05, 0.01)) {
  png(file.path(plot_output_dir[1], sprintf("%sprobe_confi_xlim_%s.png", sample_probe_number, x_cutoff)),
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
        title = sprintf("Mean Density Plot for %s probe pval vs 1e4 Random", sample_probe_number),
        x = "P-value",
        y = "Density",
        color = "Group"
      ) +
      theme_minimal()
  )
  dev.off()
}


# 對 gene ----
setorder(gt_T, `p-value`)
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
random_times <- 10000
random_list <- list()

for (i in 1:random_times) {
  set.seed(i)
  b_tmp <- lung_twas[GeneSymbol %in% sample(common_gene, sample_gene_number), ][, .(GeneSymbol, pval = top1_p)]
  b_tmp[, type := paste0("random_", i)]
  random_list[[i]] <- b_tmp
}

all_random_b <- rbindlist(random_list)
combine_all <- rbind(a[, .(GeneSymbol, pval, type)], all_random_b)


plot_x_range <- c(0, 1)
# 不同 type 的 row分組，每組的 pdf 從 from 到 to 切成n份，用 預設 gaussian kernel 估計 pdf
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


for (x_cutoff in c(1, 0.05, 0.01)) {
  png(file.path(plot_output_dir[1], sprintf("%sgene_confi_xlim_%s.png", sample_gene_number, x_cutoff)),
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
        title = sprintf("Mean Density Plot for %s Gene pval vs 1e4 Random", sample_gene_number),
        x = "P-value",
        y = "Density",
        color = "Group"
      ) +
      theme_minimal()
  )
  dev.off()
}


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
  gt_T <- fread("D:/Peter/rawData_eQTL/r2_filter_0.8/outcome/raw_T_FDR_R2_0.8.txt")

  random_times <- 10000
  eQTL_gene <- unique(gt_T$Gene)

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

  # cat(cancer_type, "genes recorded in eQTL genes:",
  #     tcga_pos[ID %in% probe_info$Gene, ID] %>% uniqueN(), "\n")

  # cat(
  #   cancer_type, "and eQTL both significant genes:",
  #   observed_intersect,
  #   "\n"
  # )


  # cat("eQTL sig genes recorded in", cancer_type, ":",
  #     eQTL_gene[eQTL_gene %in% tcga_weight$Gene] %>% uniqueN(), "\n")

  # cat("intersect of eQTL,", cancer_type, "total gene:",
  #     intersect(all_eQTL_gene, tcga_weight$Gene) %>% uniqueN(), "\n")
  cat(
    "intersect of eQTL,", cancer_type, "sig gene:",
    observed_intersect, "\n"
  )

  repeat_gene_number <- c()
  for (i in 1:random_times) {
    set.seed(i)
    repeat_gene_number[i] <- (sample(all_eQTL_gene, sample_gene_number) %in% eQTL_gene) %>%
      which() %>%
      length()
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
  setorder(gt_T, `p-value`)
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


  # 取平均畫 dist，加上 quantile band
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

  for (x_cutoff in c(1, 0.05, 0.01)) {
    png(file.path(output_dir, sprintf("%sgene_confi_xlim_%s.png", sample_gene_number, x_cutoff)),
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
          title = sprintf("Mean Density Plot for %s Gene pval vs 1e4 Random", sample_gene_number),
          x = "P-value",
          y = "Density",
          color = "Group"
        ) +
        theme_minimal()
    )
    dev.off()
  }
}


# 取出刪掉 R2<0.8 snp 後，仍有跑出 eQTL 的 probe
all_eQTL <- fread("D:/Peter/rawData_eQTL/trash/raw_maf_gt_T_pvalue.txt", header = T)[
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


# GTEx ----

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


# 下載 GTEx v8 用的 annotation ----
# 到 <https://www.gencodegenes.org/human/release_26.html> 下載，
gtf_file <- "D:/Peter/gene_enrichment/code_project/data/gencode.v26.annotation.gtf.gz"
gtf <- rtracklayer::import(gtf_file)
gene_anno <- as.data.table(gtf)
gene_anno <- gene_anno[type == "gene"]
gene_anno <- unique(gene_anno[, .(
  ensg_version = gene_id,
  query = sub("\\..*$", "", gene_id),
  symbol = gene_name,
  gene_type = gene_type,
  chr = as.character(seqnames),
  start = start,
  end = end
)])

gene_anno <- gene_anno[chr %in% sprintf("chr%s", 1:22), ]
fwrite(
  data.table(
    a = gene_anno$chr,
    b = gene_anno$start,
    c = gene_anno$end
  ),
  "D:/Peter/gene_enrichment/code_project/data/gene_hg38.txt",
  row.names = F, col.names = F, sep = "\t"
)


# liftOver 轉換 pos
cmd <- paste(
  "/mnt/d/Peter/oral_cancer/liftover/liftOver",
  "/mnt/d/Peter/gene_enrichment/code_project/data/gene_hg38.txt",
  "/mnt/d/Peter/oral_cancer/liftover/hg38ToHg18.over.chain",
  "/mnt/d/Peter/gene_enrichment/code_project/data/gene_hg38TOhg18.bed",
  "/mnt/d/Peter/gene_enrichment/code_project/data//hg38TOhg18_unmapped.bed"
)

system2("wsl", cmd)

a <- fread("D:/Peter/gene_enrichment/code_project/data/gene_hg38TOhg18.bed", header = F)
# 忽略轉失敗 .bed 檔案的 #Deleted in new 該行

unmapped_path <- "D:/Peter/gene_enrichment/code_project/data//hg38TOhg18_unmapped.bed"

if (!file.exists(unmapped_path) || file.info(unmapped_path)$size == 0) {
  a_unmapped <- data.table(V1 = character(), V2 = integer(), V3 = integer())
} else {
  a_unmapped <- fread(
    unmapped_path,
    header = FALSE,
    comment.char = "#"
  )
}

col_order <- names(gene_anno)
names(a_unmapped) <- c("chr", "start", "end")
names(a) <- c("chr_hg18", "start_hg18", "end_hg18")

# 新增轉換後 pos
a_unmapped[, c("chr_hg18", "start_hg18", "end_hg18") := "no"]
final <- merge(gene_anno, a_unmapped, all.x = T)

success_idx <- is.na(final$chr_hg18)
final[
  success_idx,
  c("chr_hg18", "start_hg18", "end_hg18") := a[, .(chr_hg18, start_hg18, end_hg18)]
]

# 轉換紀錄
# fwrite(final,
#   "D:/Peter/gene_enrichment/code_project/data/hg38_hg18.txt",
#   row.names = F, col.names = T, sep = "\t"
# )

final[, c("chr", "start", "end") := NULL]
setnames(final, old = c("chr_hg18", "start_hg18", "end_hg18"), new = c("chr", "start", "end"))
setcolorder(final, col_order)
final[which(gene_anno$start == "no"), c("chr", "start", "end") := NA_character_]
final[, start := as.numeric(start)]
final[, end := as.numeric(end)]


fwrite(final,
  "D:/Peter/gene_enrichment/code_project/data/gencode_v26_gene_anno_hg18.txt",
  row.names = F, col.names = T, sep = "\t"
)


# fusion_weights ENSG -> gene ----
# install.packages("BiocManager")
# BiocManager::install("mygene")

library(mygene)


# 讀 OC eQTL 所有的 probe
all_eQTL <- fread("D:/Peter/rawData_eQTL/trash/raw_maf_gt_T_pvalue.txt", header = T)[
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

  gt_T <- fread("D:/Peter/rawData_eQTL/r2_filter_0.8/outcome/raw_T_FDR_R2_0.8.txt")
  random_times <- 10000
  all_eQTL_gene <- unique(all_eQTL$Gene)
  eQTL_gene <- unique(all_eQTL[Gene %in% gt_T$Gene, Gene])

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
    main = sprintf("1e4 times Gene Intersection Size between %s Sig. and eQTL Sig.", tissue_name),
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
    png(file.path(output_dir, sprintf("%sgene_Top1_PvalDIST_%s.png", sample_gene_number, x_cutoff)),
      width = 8, height = 6, units = "in", res = 300
    )
    print(
      density_plot(
        combine_all, "pval",
        sprintf("Density Plot for %s Gene Top1_Pval vs 1e4 Random", sample_gene_number),
        x_range = c(0, x_cutoff)
      )
    )
    dev.off()
  }


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
    png(file.path(output_dir, sprintf("%sgene_confi_xlim_%s.png", sample_gene_number, x_cutoff)),
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
          title = sprintf("Mean Density Plot for %s Gene Top1_Pval vs 1e4 Random", sample_gene_number),
          x = "p-value",
          y = "Density",
          color = "Group"
        ) +
        theme_minimal()
    )
    dev.off()
  }
}






tissue_name_list <- c("GTEx-salivary", "GTEx-esophagus", "GTEx-thyroid", "GTEx-lung")
file_list <- sprintf(
  "D:/Peter/gene_enrichment/code_project/data/%s/fusion_project/fusion_weights_summary.csv",
  c("GTEx_salivary", "GTEx_esophagus", "GTEx_thyroid", "GTEx_lung")
)


plot_output_dir <- sprintf(
  "D:/Peter/gene_enrichment/code_project/outcome/rawData_eQTL_T/for_gene/%s",
  tissue_name_list
)

sapply(plot_output_dir, function(x) {
  if (!dir.exists(x)) dir.create(x, recursive = TRUE)
})


for (i in seq_along(tissue_name_list)) {
  gene_enrichment(
    tissue_name = tissue_name_list[i],
    weight_path = file_list[i],
    output_dir = plot_output_dir[i]
  )
}


# 跑多種組織需要的 csv ----
#  Rscript /mnt/d/Peter/gene_enrichment/code_project/data/summarize_weights_LOTStissue.R


# 設定資料夾
base_dir <- "D:/Peter/gene_enrichment/code_project/data"
out_dir <- "D:/Peter/gene_enrichment/code_project/outcome/rawData_eQTL_T/for_gene"

gtex <- data.frame(
  tissue_dir = c(
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


sapply(plot_output_dir, function(x) {
  if (!dir.exists(x)) dir.create(x, recursive = TRUE)
})


all_eQTL <- fread("D:/Peter/rawData_eQTL/trash/raw_maf_gt_T_pvalue.txt", header = T)[
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


for (i in seq_along(gtex$tissue_dir)) {
  gene_enrichment(
    tissue_name = gtex$tissue_dir[i],
    weight_path = output_csv[i],
    output_dir = plot_output_dir[i]
  )
}




# After remove R2<0.8 snp, run eQTL genes: 15992
# After remove R2<0.8 snp, run eQTL probes: 20906
# TCGA-HNSC genes: 2767
# TCGA-HNSC After remove NULL and NA GeneID/Top1_Pval, genes: 2763
# TCGA-HNSC sig. genes: 2069
# eQTL sig. genes: 129
# TCGA-HNSC significant genes recorded in eQTL: 1540
# intersect of eQTL, TCGA-HNSC sig gene: 34
# sample 1540 genes, mean sig in TCGA-HNSC gene number: 12.4399
# sample 1540 genes, sig in TCGA-HNSC gene number range: 2 27
# sample 1540 genes, sig in TCGA-HNSC gene number >= 34 times: 0
# Processed 32768 groups out of 32768. 100% done. Time elapsed: 18s. ETA: 0s.
# After remove R2<0.8 snp, run eQTL genes: 15992
# After remove R2<0.8 snp, run eQTL probes: 20906
# TCGA-LUAD genes: 2978
# TCGA-LUAD After remove NULL and NA GeneID/Top1_Pval, genes: 2974
# TCGA-LUAD sig. genes: 2248
# eQTL sig. genes: 129
# TCGA-LUAD significant genes recorded in eQTL: 1644
# intersect of eQTL, TCGA-LUAD sig gene: 37
# sample 1644 genes, mean sig in TCGA-LUAD gene number: 13.2794
# sample 1644 genes, sig in TCGA-LUAD gene number range: 3 28
# sample 1644 genes, sig in TCGA-LUAD gene number >= 37 times: 0
# Processed 32768 groups out of 32768. 100% done. Time elapsed: 13s. ETA: 0s.
# After remove R2<0.8 snp, run eQTL genes: 15992
# After remove R2<0.8 snp, run eQTL probes: 20906
# TCGA-LUSC genes: 2548
# TCGA-LUSC After remove NULL and NA GeneID/Top1_Pval, genes: 2545
# TCGA-LUSC sig. genes: 1862
# eQTL sig. genes: 129
# TCGA-LUSC significant genes recorded in eQTL: 1366
# intersect of eQTL, TCGA-LUSC sig gene: 31
# sample 1366 genes, mean sig in TCGA-LUSC gene number: 11.0279
# sample 1366 genes, sig in TCGA-LUSC gene number range: 2 27
# sample 1366 genes, sig in TCGA-LUSC gene number >= 31 times: 0
# Processed 32768 groups out of 32768. 100% done. Time elapsed: 14s. ETA: 0s.


# GTEx-lung raw data genes and ENSG genes (remove NA): 24645 24687
# GTEx-lung After remove NA Top1_Pval, genes and ENSg gene: 24645 24687
# After choose most nearest gene for repeat gene, genes and ENSG genes: 24645 24645
# GTEx-lung sig. genes: 10376
# eQTL sig. genes: 129
# intersect of eQTL, GTEx-lung total genes: 12282
# GTEx-lung significant genes recorded in eQTL: 5471
# eQTL significant genes recorded in GTEx-lung : 98
# intersect of eQTL, GTEx-lung sig genes: 64
# sample 98 genes, mean sig in GTEx-lung genes number: 43.6238
# sample 98 genes, sig in GTEx-lung genes number range: 25, 62
# sample 98 genes, sig in GTEx-lung genes number >= 64 times: 0


# GTEx-salivary raw data genes and ENSG genes (remove NA): 23814 23847
# GTEx-salivary After remove NA Top1_Pval, genes and ENSg gene: 23814 23847
# After choose most nearest gene for repeat gene, genes and ENSG genes: 23814 23814
# GTEx-salivary sig. genes: 4281
# eQTL sig. genes: 129
# intersect of eQTL, GTEx-salivary total genes: 12059
# GTEx-salivary significant genes recorded in eQTL: 2057
# eQTL significant genes recorded in GTEx-salivary : 90
# intersect of eQTL, GTEx-salivary sig genes: 30
# sample 90 genes, mean sig in GTEx-salivary genes number: 15.373
# sample 90 genes, sig in GTEx-salivary genes number range: 4, 29
# sample 90 genes, sig in GTEx-salivary genes number >= 30 times: 0


# GTEx-esophagus raw data genes and ENSG genes (remove NA): 22588 22607
# GTEx-esophagus After remove NA Top1_Pval, genes and ENSg gene: 22588 22607
# After choose most nearest gene for repeat gene, genes and ENSG genes: 22588 22588
# GTEx-esophagus sig. genes: 10345
# eQTL sig. genes: 129
# intersect of eQTL, GTEx-esophagus total genes: 12014
# GTEx-esophagus significant genes recorded in eQTL: 5901
# eQTL significant genes recorded in GTEx-esophagus : 94
# intersect of eQTL, GTEx-esophagus sig genes: 70
# sample 94 genes, mean sig in GTEx-esophagus genes number: 46.092
# sample 94 genes, sig in GTEx-esophagus genes number range: 26, 65
# sample 94 genes, sig in GTEx-esophagus genes number >= 70 times: 0


# GTEx-thyroid raw data genes and ENSG genes (remove NA): 24722 24767
# GTEx-thyroid After remove NA Top1_Pval, genes and ENSg gene: 24722 24767
# After choose most nearest gene for repeat gene, genes and ENSG genes: 24722 24722
# GTEx-thyroid sig. genes: 12638
# eQTL sig. genes: 129
# intersect of eQTL, GTEx-thyroid total genes: 12272
# GTEx-thyroid significant genes recorded in eQTL: 6608
# eQTL significant genes recorded in GTEx-thyroid : 99
# intersect of eQTL, GTEx-thyroid sig genes: 70
# sample 99 genes, mean sig in GTEx-thyroid genes number: 53.2995
# sample 99 genes, sig in GTEx-thyroid genes number range: 33, 73
# sample 99 genes, sig in GTEx-thyroid genes number >= 70 times: 3
