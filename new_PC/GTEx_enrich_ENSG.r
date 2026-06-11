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
#   "D:/Peter/gene_enrichment/code_project/outcome/rawData_eQTL",
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


## ft ----
find_nearest_interval <- function(dt, unit, query_start, query_end, for_gene = TRUE) {
  if (for_gene) {
    dt_sub <- copy(dt[gene == unit])
  } else {
    dt_sub <- copy(dt[ENSEMBL == unit])
  }

  query_dt <- data.table(
    query_start = query_start,
    query_end = query_end
  )

  result <- rbindlist(lapply(seq_len(nrow(query_dt)), function(i) {
    tmp <- copy(dt_sub)

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





# OC eQTL gene -> ENSG ----
# 下載 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE37991 我們資料有紀錄 Entrez_Gene_ID 檔案，用 Entrez_Gene_ID 轉 ENSG
ENTREZID_probe <- readLines("D:/Peter/oral_cancer/GSE37991_family.soft/GSE37991_family.soft", warn = FALSE)

# 找 platform table 的開始與結束位置
start <- grep("^!platform_table_begin", ENTREZID_probe)
end <- grep("^!platform_table_end", ENTREZID_probe)
platform_txt <- ENTREZID_probe[(start + 1):(end - 1)]

# 跟我們 probe gene 資料 "D:/Peter/oral_cancer/expression/expression_data/probe24526infoB.txt 對照過
# ID, ILMN_Gene 表示 gene, probe，但有些沒對到，可能是 probe24526infoB.txt 更新過
platform_dt <- fread(
  text = paste(platform_txt, collapse = "\n"),
  sep = "\t",
  header = TRUE,
  fill = TRUE,
  quote = ""
)[, .(ID, ILMN_Gene, Entrez_Gene_ID)]
platform_dt$Entrez_Gene_ID <- as.character(platform_dt$Entrez_Gene_ID)


# 哪些 gene 跟手邊資料不合
probe_pos <- fread("D:/Peter/oral_cancer/expression/expression_data/probe_pos_mulInterval.txt")
gene_probe <- fread("D:/Peter/oral_cancer/expression/expression_data/probe24526infoB.txt")[
  ,
  .(TargetID, PROBE_ID)
]
probe_pos <- gene_probe[probe_pos, on = .(PROBE_ID = Gene), nomatch = 0]

merge_result <- merge(
  probe_pos,
  platform_dt[, .(PROBE_ID = ID, TargetID = ILMN_Gene)],
  by = c("PROBE_ID", "TargetID"),
  all = TRUE
)

merge_result[, in_both := fifelse(
  is.na(start),
  "F",
  "T"
)]
nomatch_gene <- merge_result[is.na(start), PROBE_ID] %>% unique()
merge_result <- merge_result[PROBE_ID %in% nomatch_gene] %>% setkey(., PROBE_ID)

fwrite(merge_result,
  "D:/Peter/oral_cancer/GSE37991_family.soft/GSE37991_family_not_OCgene.txt",
  row.names = F, col.names = T, sep = "\t"
)


# 讀 OC eQTL 所有的 probe
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


# 挑出有在我們資料的 Entrez_Gene_ID row
platform_dt <- platform_dt[ID %in% all_eQTL$probe]


# OC eQTL gene -> ENSG
library(org.Hs.eg.db)
library(AnnotationDbi)
library(mygene)

# 紀錄 ENTREZID 的檔案，轉乘 ENSG
tmp <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(na.omit(as.character(platform_dt$Entrez_Gene_ID))),
  keytype = "ENTREZID",
  columns = c("ENSEMBL")
) %>% as.data.table()

# 確認過沒 na, " "
cat(
  "After QC filter, eQTL ENntrezID and probe :",
  platform_dt$Entrez_Gene_ID %>% uniqueN(), platform_dt$ID %>% uniqueN(),
  "\n"
)

tmp <- tmp[!is.na(ENTREZID) & !is.na(tmp$ENSEMBL)]
cat(
  "After transform to ENSG, eQTL probe and Entrez_Gene_ID and ENSG :",
  platform_dt[Entrez_Gene_ID %in% tmp$ENTREZID, ID] %>% uniqueN(),
  tmp$ENTREZID %>% uniqueN(),
  tmp$ENSEMBL %>% uniqueN(),
  "\n"
)

# 有時 Entrez ID 對到多個 Ensembl ID
all_eQTL_ENSG <- tmp[platform_dt, on = .(ENTREZID = Entrez_Gene_ID), nomatch = 0]

# 加上 symbol 出現次數
n_entrez_times <- all_eQTL_ENSG[, .(
  n_entrez = uniqueN(ENTREZID)
), by = ENSEMBL]
all_eQTL_ENSG <- all_eQTL_ENSG[n_entrez_times, on = .(ENSEMBL), nomatch = 0]

fwrite(all_eQTL_ENSG,
  "D:/Peter/gene_enrichment/code_project/outcome/rawData_eQTL/for_ENSG/manyGENE_to_oneENSG.txt",
  row.names = F, col.names = T, sep = "\t"
)


# 每個 ENSEMBL 選最顯著的
# setkey(all_eQTL_ENSG, pval)
# all_eQTL_ENSG <- all_eQTL_ENSG[, .SD[1], by = ENSEMBL]
# cat(
#   "After choose most sig. ENSG for repeat ENSG, ENSG_genes:",
#   all_eQTL_ENSG$Gene %>% uniqueN(), all_eQTL_ENSG$ENSEMBL %>% uniqueN(),
#   "\n"
# )

# all_eQTL_ENSGgene <- unique(all_eQTL_ENSG$ENSEMBL)


# ft
EHSG_enrichment <- function(
  tissue_name,
  weight_path,
  output_dir
) {
  gtex_weight <- fread(weight_path)
  gtex_weight[, query := sub("\\..*", "", Gene_ID)]

  gene_anno <- readRDS("D:/Peter/gene_enrichment/code_project/data/gencode_v26_gene_anno.rds")
  gene_anno <- as.data.table(gene_anno)

  # gene 名稱併入fusion weight
  gtex_weight <- gene_anno[gtex_weight, on = .(query), nomatch = 0]
  setnames(gtex_weight, old = "symbol", new = "gene")

  all_eQTL_ENSG <- fread("D:/Peter/gene_enrichment/code_project/outcome/rawData_eQTL/for_ENSG/manyGENE_to_oneENSG.txt")[
    ,
    .(ENTREZID, ENSEMBL, ID, n_entrez)
  ]

  probe_pos <- fread("D:/Peter/oral_cancer/expression/expression_data/probe_pos_mulInterval.txt")
  gene_probe <- fread("D:/Peter/oral_cancer/expression/expression_data/probe24526infoB.txt")[
    ,
    .(TargetID, PROBE_ID)
  ]
  probe_pos <- gene_probe[probe_pos, on = .(PROBE_ID = Gene), nomatch = 0]
  # 　轉乘 ENSG 的 probe 補上我們資料的 pos
  all_eQTL_ENSG <- probe_pos[all_eQTL_ENSG, on = .(PROBE_ID = ID), nomatch = 0]


  repeat_ENSG <- all_eQTL_ENSG[n_entrez > 1 & ENSEMBL %in% gtex_weight$query, ENSEMBL] %>%
    unique()

  repeat_row <- list()
  for (i in seq_along(repeat_ENSG)) {
    repeat_row[[i]] <- find_nearest_interval(
      dt = all_eQTL_ENSG,
      unit = repeat_ENSG[i],
      query_start = gtex_weight[query == repeat_ENSG[i], start],
      query_end = gtex_weight[query == repeat_ENSG[i], end],
      for_gene = F
    )
  }

  repeat_row_result <- rbindlist(repeat_row) %>% as.data.table()
  col_order <- names(all_eQTL_ENSG)
  repeat_row_result <- repeat_row_result[, ..col_order]

  # repeat gene choose one
  # 處理出現在 probe_pos 的 repeat ensg 就好，其他的不重要，因為後續跟我們資料 取交集，會被刪掉
  repeat_row_result[, choose := 1]
  all_eQTL_ENSG <- repeat_row_result[all_eQTL_ENSG, on = col_order]

  all_eQTL_ENSG[is.na(choose), choose := 2]
  setkey(all_eQTL_ENSG, choose)
  all_eQTL_ENSG <- all_eQTL_ENSG[, .SD[1], by = ENSEMBL]

  cat(
    "After choose most nearest ENSG_gene for repeat ENSG_gene, probe and entrez_ID and ENSG genes:",
    uniqueN(all_eQTL_ENSG$PROBE_ID), uniqueN(all_eQTL_ENSG$ENTREZID), uniqueN(all_eQTL_ENSG$ENSEMBL),
    "\n"
  )


  # ----------------

  gt_N <- fread("D:/Peter/rawData_eQTL/r2_filter_0.8/outcome/raw_N_FDR_R2_0.8.txt")
  random_times <- 10000
  eQTL_gene <- unique(all_eQTL_ENSG[Gene %in% gt_N$Gene, ENSEMBL])

  gtex_sigGene <- gtex_weight[Top1_Pval < 0.01, query] %>%
    unique()

  cat(
    tissue_name, "sig. ENSG_genes:",
    uniqueN(gtex_sigGene),
    "\n"
  )
  cat(
    "eQTL sig. ENSG_genes:",
    uniqueN(eQTL_gene),
    "\n"
  )

  common_gene <- intersect(all_eQTL_ENSGgene, gtex_weight$query)
  sample_gene_number <- uniqueN(eQTL_gene[eQTL_gene %in% gtex_weight$query])
  observed_intersect <- intersect(eQTL_gene, gtex_sigGene) %>% uniqueN()

  cat(
    tissue_name, "significant ENSG_genes recorded in eQTL:",
    gtex_sigGene[gtex_sigGene %in% all_eQTL_ENSGgene] %>% uniqueN(), "\n"
  )
  cat(
    "eQTL significant ENSG_genes recorded in", tissue_name, ":",
    uniqueN(eQTL_gene[eQTL_gene %in% gtex_weight$query]), "\n"
  )
  cat(
    "intersect of eQTL,", tissue_name, "total ENSG_genes:",
    intersect(all_eQTL_ENSGgene, gtex_weight$query) %>% uniqueN(), "\n"
  )
  cat(
    "intersect of eQTL,", tissue_name, "sig ENSG_genes:",
    observed_intersect, "\n"
  )

  repeat_gene_number <- c()
  for (i in 1:random_times) {
    set.seed(i)
    repeat_gene_number[i] <- (sample(
      common_gene,
      sample_gene_number
    ) %in% gtex_sigGene[gtex_sigGene %in% all_eQTL_ENSGgene]) %>%
      which() %>%
      length()
  }

  cat(
    "sample", sample_gene_number, "ENSG_genes, mean sig in", tissue_name, "ENSG_genes number:",
    mean(repeat_gene_number), "\n"
  )
  cat(
    "sample", sample_gene_number, "ENSG_genes, sig in", tissue_name, "ENSG_genes number range:",
    range(repeat_gene_number), "\n"
  )
  cat(
    "sample", sample_gene_number, "ENSG_genes, sig in", tissue_name, "ENSG_genes number >=",
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

  a <- gtex_weight[query %in% eQTL_gene, ][, .(gene_name = query, pval = Top1_Pval)]
  setkey(a, pval)
  a <- a[, .SD[1], by = gene_name]
  a[, type := "eQTL"]

  random_list <- list()
  random_times <- 5
  for (i in 1:random_times) {
    set.seed(i)
    b_tmp <- gtex_weight[query %in% sample(common_gene, sample_gene_number), ][, .(gene_name = query, pval = Top1_Pval)]
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



plot_output_dir <- sprintf(
  "D:/Peter/gene_enrichment/code_project/outcome/rawData_eQTL/for_ENSG/%s",
  c("lungTWAS", "HNSC", "TCGA-LUAD", "TCGA-LUSC", "GTEx-salivary", "GTEx-esophagus", "GTEx-thyroid", "GTEx-lung")
)

sapply(plot_output_dir, function(x) {
  if (!dir.exists(x)) dir.create(x, recursive = TRUE)
})



EHSG_enrichment(
  tissue_name = "GTEx-salivary",
  weight_path = "D:/Peter/gene_enrichment/code_project/data/GTEx_salivary/fusion_project/fusion_weights_summary.csv",
  output_dir = plot_output_dir[5]
)
