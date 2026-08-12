{
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
  #   "C:/Peter/gene_enrichment/code_project/outcome/rawData_eQTL_N",
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




  # 下載 GTEx v8 用的 annotation，得到 gene 的 ENSG ----
  # 到 <https://www.gencodegenes.org/human/release_26.html> 下載，
  # gtf_file <- "C:/Peter/gene_enrichment/code_project/data/gencode.v26.annotation.gtf.gz"
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
  #   "C:/Peter/gene_enrichment/code_project/data/gene_hg38.txt",
  #   row.names = F, col.names = F, sep = "\t"
  # )


  # liftOver 轉換 pos, hg38 轉成 hg18 ----
  # cmd <- paste(
  #   "/mnt/c/Peter/oral_cancer/liftover/liftOver",
  #   "/mnt/c/Peter/gene_enrichment/code_project/data/gene_hg38.txt",
  #   "/mnt/c/Peter/oral_cancer/liftover/hg38ToHg18.over.chain",
  #   "/mnt/c/Peter/gene_enrichment/code_project/data/gene_hg38TOhg18.bed",
  #   "/mnt/c/Peter/gene_enrichment/code_project/data//hg38TOhg18_unmapped.bed"
  # )

  # system2("wsl", cmd)

  # a <- fread("C:/Peter/gene_enrichment/code_project/data/gene_hg38TOhg18.bed", header = F)
  # # 忽略轉失敗 .bed 檔案的 #Deleted in new 該行

  # unmapped_path <- "C:/Peter/gene_enrichment/code_project/data//hg38TOhg18_unmapped.bed"

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
  # #   "C:/Peter/gene_enrichment/code_project/data/hg38_hg18.txt",
  # #   row.names = F, col.names = T, sep = "\t"
  # # )

  # final[, c("chr", "start", "end") := NULL]
  # setnames(final, old = c("chr_hg18", "start_hg18", "end_hg18"), new = c("chr", "start", "end"))
  # setcolorder(final, col_order)
  # final[which(gene_anno$start == "no"), c("chr", "start", "end") := NA_character_]
  # final[, start := as.numeric(start)]
  # final[, end := as.numeric(end)]


  # fwrite(final,
  #   "C:/Peter/gene_enrichment/code_project/data/gencode_v26_gene_anno_hg18.txt",
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



  # ft
  gene_enrichment <- function(
    tissue_name,
    weight_path,
    output_dir
  ) {
    gtex_weight <- fread(weight_path)
    gtex_weight[, query := sub("\\..*", "", Gene_ID)]

    gene_anno <- fread("C:/Peter/gene_enrichment/code_project/data/gencode_v26_gene_anno_hg18.txt")

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
    probe_pos <- fread("C:/Peter/oral_cancer/expression/expression_data/probe_pos_mulInterval.txt")
    gene_probe <- fread("C:/Peter/oral_cancer/expression/expression_data/probe24526infoB.txt")[
      ,
      .(TargetID, PROBE_ID)
    ]
    probe_pos <- gene_probe[probe_pos, on = .(PROBE_ID = Gene), nomatch = 0]
    repeat_ensg <- gtex_weight[n_ensg > 1 & gene %in% probe_pos$TargetID, gene] %>%
      unique()


    repeat_row <- list()
    # 對 gtex_weight 重複出現的 gene，找出離我們資料任一區間最近的 gtex_weight gene
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

    gt_N <- fread("C:/Peter/rawData_eQTL/r2_filter_0.8/outcome/raw_N_FDR_R2_0.8.txt")

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


    # 從共同基因隨機挑，看多少 在 gtex 顯著
    library(parallel)
    gtex_sig_sub <- gtex_sigGene[gtex_sigGene %in% all_eQTL_gene]
    is_sig <- common_gene %in% gtex_sig_sub
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
        # is_sig 判斷 common_gene 是否在 gtex_sig_sub，一串的 T,F，sum 會把T當1, F當0
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


  # 跑多種組織需要的 csv ----
  #  Rscript /mnt/c/Peter/gene_enrichment/code_project/data/summarize_weights_LOTStissue.R


  base_dir <- "C:/Peter/gene_enrichment/code_project/data"
  out_dir <- "C:/Peter/gene_enrichment/code_project/outcome/rawData_eQTL_N/for_gene"


  gtex <- data.frame(
    tissue_dir = c(
      "GTEx_salivary",
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
  all_eQTL <- probe_info[all_eQTL, on = .(probe), nomatch = 0]
  all_eQTL_gene <- unique(all_eQTL$Gene)
  all_eQTL_probe <- unique(all_eQTL$probe)


  # seq_along(gtex$tissue_dir)

  for (i in c(6, 8, 9, 10)) {
    start_time <- Sys.time()
    gene_enrichment(
      tissue_name = gtex$tissue_dir[i],
      weight_path = output_csv[i],
      output_dir = plot_output_dir[i]
    )

    end_time <- Sys.time()
    print(end_time - start_time)
  }

  # 1e8 一個 tissue 12.79 min
  # 1e9 一個 tissue 1.85 hr
}





# 有m白球，n黑球，抽k球，抽後不放回，白球數量>q的機率。lower.tail=T 則是，白球數量<=q的機率
phyper(
  # 兩邊都顯著數量-1
  q = 79 - 1,
  # GTEx sig number
  m = 4328,
  # background - GTEx sig number
  n = 11896 - 4328,
  # oral sig number
  k = 136,
  lower.tail = FALSE
) %>% sprintf("%.2e", .)


# permutation 1e9 times

# GTEx_Esophagus_Gastroesophageal raw data genes and ENSG genes (remove NA): 22679 22703
# GTEx_Esophagus_Gastroesophageal After remove NA Top1_Pval, genes and ENSg gene: 22679 22703
# After choose most nearest gene for repeat gene, genes and ENSG genes: 22679 22679
# GTEx_Esophagus_Gastroesophageal sig. genes: 8062
# eQTL sig. genes: 177
# intersect of eQTL, GTEx_Esophagus_Gastroesophageal total genes: 11896
# GTEx_Esophagus_Gastroesophageal significant genes recorded in eQTL: 4328
# eQTL significant genes recorded in GTEx_Esophagus_Gastroesophageal : 136
# intersect of eQTL, GTEx_Esophagus_Gastroesophageal sig genes: 79
# sample 136 genes, mean sig in GTEx_Esophagus_Gastroesophageal genes number: 49.47902
# sample 136 genes, sig in GTEx_Esophagus_Gastroesophageal genes number range: 17 87
# sample 136 genes, sig in GTEx_Esophagus_Gastroesophageal genes number >= 79 times: 195
# empirical p-value: 1.96e-07
# Hypergeometric dist. pval: 1.86e-07


# GTEx_salivary raw data genes and ENSG genes (remove NA): 23814 23847
# GTEx_salivary After remove NA Top1_Pval, genes and ENSg gene: 23814 23847
# After choose most nearest gene for repeat gene, genes and ENSG genes: 23814 23814
# GTEx_salivary sig. genes: 4282
# eQTL sig. genes: 177
# intersect of eQTL, GTEx_salivary total genes: 12059
# GTEx_salivary significant genes recorded in eQTL: 2058
# eQTL significant genes recorded in GTEx_salivary : 136
# intersect of eQTL, GTEx_salivary sig genes: 42
# sample 136 genes, mean sig in GTEx_salivary genes number: 23.20981
# sample 136 genes, sig in GTEx_salivary genes number range: 2 52
# sample 136 genes, sig in GTEx_salivary genes number >= 42 times: 49773
# empirical p-value: 4.98e-05
# Hypergeometric dist. pval: 4.96e-05


# GTEx_esophagus raw data genes and ENSG genes (remove NA): 22588 22607
# GTEx_esophagus After remove NA Top1_Pval, genes and ENSg gene: 22588 22607
# After choose most nearest gene for repeat gene, genes and ENSG genes: 22588 22588
# GTEx_esophagus sig. genes: 10347
# eQTL sig. genes: 177
# intersect of eQTL, GTEx_esophagus total genes: 12014
# GTEx_esophagus significant genes recorded in eQTL: 5903
# eQTL significant genes recorded in GTEx_esophagus : 140
# intersect of eQTL, GTEx_esophagus sig genes: 94
# sample 140 genes, mean sig in GTEx_esophagus genes number: 68.78782
# sample 140 genes, sig in GTEx_esophagus genes number range: 34 104
# sample 140 genes, sig in GTEx_esophagus genes number >= 94 times: 11510
# empirical p-value: 1.15e-05
# Hypergeometric dist. pval: 1.14e-05


# GTEx_thyroid raw data genes and ENSG genes (remove NA): 24722 24767
# GTEx_thyroid After remove NA Top1_Pval, genes and ENSg gene: 24722 24767
# After choose most nearest gene for repeat gene, genes and ENSG genes: 24722 24722
# GTEx_thyroid sig. genes: 12640
# eQTL sig. genes: 177
# intersect of eQTL, GTEx_thyroid total genes: 12272
# GTEx_thyroid significant genes recorded in eQTL: 6610
# eQTL significant genes recorded in GTEx_thyroid : 142
# intersect of eQTL, GTEx_thyroid sig genes: 95
# sample 142 genes, mean sig in GTEx_thyroid genes number: 76.48505
# sample 142 genes, sig in GTEx_thyroid genes number range: 42 110
# sample 142 genes, sig in GTEx_thyroid genes number >= 95 times: 1019724
# empirical p-value: 1.02e-03
# Hypergeometric dist. pval: 1.02e-03


# GTEx_lung raw data genes and ENSG genes (remove NA): 24645 24687
# GTEx_lung After remove NA Top1_Pval, genes and ENSg gene: 24645 24687
# After choose most nearest gene for repeat gene, genes and ENSG genes: 24645 24645
# GTEx_lung sig. genes: 10378
# eQTL sig. genes: 177
# intersect of eQTL, GTEx_lung total genes: 12282
# GTEx_lung significant genes recorded in eQTL: 5473
# eQTL significant genes recorded in GTEx_lung : 142
# intersect of eQTL, GTEx_lung sig genes: 82
# sample 142 genes, mean sig in GTEx_lung genes number: 63.27668
# sample 142 genes, sig in GTEx_lung genes number range: 29 99
# sample 142 genes, sig in GTEx_lung genes number >= 82 times: 1021130
# empirical p-value: 1.02e-03
# Hypergeometric dist. pval: 1.02e-03


# GTEx_Adipose_Subcutaneous raw data genes and ENSG genes (remove NA): 23361 23395
# GTEx_Adipose_Subcutaneous After remove NA Top1_Pval, genes and ENSg gene: 23361 23395
# After choose most nearest gene for repeat gene, genes and ENSG genes: 23361 23361
# GTEx_Adipose_Subcutaneous sig. genes: 10993
# eQTL sig. genes: 177
# intersect of eQTL, GTEx_Adipose_Subcutaneous total genes: 12075
# GTEx_Adipose_Subcutaneous significant genes recorded in eQTL: 5965
# eQTL significant genes recorded in GTEx_Adipose_Subcutaneous : 135
# intersect of eQTL, GTEx_Adipose_Subcutaneous sig genes: 95
# sample 135 genes, mean sig in GTEx_Adipose_Subcutaneous genes number: 66.68925
# sample 135 genes, sig in GTEx_Adipose_Subcutaneous genes number range: 33 101
# sample 135 genes, sig in GTEx_Adipose_Subcutaneous genes number >= 95 times: 560
# empirical p-value: 5.61e-07
# Hypergeometric dist. pval: 5.44e-07


# GTEx_Brain_Cerebellum raw data genes and ENSG genes (remove NA): 23892 23922
# GTEx_Brain_Cerebellum After remove NA Top1_Pval, genes and ENSg gene: 23892 23922
# After choose most nearest gene for repeat gene, genes and ENSG genes: 23892 23892
# GTEx_Brain_Cerebellum sig. genes: 8303
# eQTL sig. genes: 177
# intersect of eQTL, GTEx_Brain_Cerebellum total genes: 11911
# GTEx_Brain_Cerebellum significant genes recorded in eQTL: 4355
# eQTL significant genes recorded in GTEx_Brain_Cerebellum : 133
# intersect of eQTL, GTEx_Brain_Cerebellum sig genes: 74
# sample 133 genes, mean sig in GTEx_Brain_Cerebellum genes number: 48.62822
# sample 133 genes, sig in GTEx_Brain_Cerebellum genes number range: 17 83
# sample 133 genes, sig in GTEx_Brain_Cerebellum genes number >= 74 times: 5260
# empirical p-value: 5.26e-06
# Hypergeometric dist. pval: 5.21e-06


# GTEx_Esophagus_Muscularis raw data genes and ENSG genes (remove NA): 22481 22503
# GTEx_Esophagus_Muscularis After remove NA Top1_Pval, genes and ENSg gene: 22481 22503
# After choose most nearest gene for repeat gene, genes and ENSG genes: 22481 22481
# GTEx_Esophagus_Muscularis sig. genes: 10263
# eQTL sig. genes: 177
# intersect of eQTL, GTEx_Esophagus_Muscularis total genes: 11919
# GTEx_Esophagus_Muscularis significant genes recorded in eQTL: 5649
# eQTL significant genes recorded in GTEx_Esophagus_Muscularis : 136
# intersect of eQTL, GTEx_Esophagus_Muscularis sig genes: 87
# sample 136 genes, mean sig in GTEx_Esophagus_Muscularis genes number: 64.45737
# sample 136 genes, sig in GTEx_Esophagus_Muscularis genes number range: 28 99
# sample 136 genes, sig in GTEx_Esophagus_Muscularis genes number >= 87 times: 66650
# empirical p-value: 6.67e-05
# Hypergeometric dist. pval: 6.71e-05


# GTEx_Heart_Left raw data genes and ENSG genes (remove NA): 19870 19884
# GTEx_Heart_Left After remove NA Top1_Pval, genes and ENSg gene: 19870 19884
# After choose most nearest gene for repeat gene, genes and ENSG genes: 19870 19870
# GTEx_Heart_Left sig. genes: 7374
# eQTL sig. genes: 177
# intersect of eQTL, GTEx_Heart_Left total genes: 11126
# GTEx_Heart_Left significant genes recorded in eQTL: 4269
# eQTL significant genes recorded in GTEx_Heart_Left : 131
# intersect of eQTL, GTEx_Heart_Left sig genes: 75
# sample 131 genes, mean sig in GTEx_Heart_Left genes number: 50.26387
# sample 131 genes, sig in GTEx_Heart_Left genes number range: 20 85
# sample 131 genes, sig in GTEx_Heart_Left genes number >= 75 times: 8131
# empirical p-value: 8.13e-06
# Hypergeometric dist. pval: 8.24e-06


# GTEx_Muscle_Skeletal raw data genes and ENSG genes (remove NA): 19963 19978
# GTEx_Muscle_Skeletal After remove NA Top1_Pval, genes and ENSg gene: 19963 19978
# After choose most nearest gene for repeat gene, genes and ENSG genes: 19963 19963
# GTEx_Muscle_Skeletal sig. genes: 9673
# eQTL sig. genes: 177
# intersect of eQTL, GTEx_Muscle_Skeletal total genes: 11429
# GTEx_Muscle_Skeletal significant genes recorded in eQTL: 5725
# eQTL significant genes recorded in GTEx_Muscle_Skeletal : 133
# intersect of eQTL, GTEx_Muscle_Skeletal sig genes: 89
# sample 133 genes, mean sig in GTEx_Muscle_Skeletal genes number: 66.6224
# sample 133 genes, sig in GTEx_Muscle_Skeletal genes number range: 32 100
# sample 133 genes, sig in GTEx_Muscle_Skeletal genes number >= 89 times: 58894
# empirical p-value: 5.89e-05
# Hypergeometric dist. pval: 5.90e-05