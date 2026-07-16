rm(list = ls())
gc()

{
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
eQTL_ld <- fread("C:/Peter/gene_enrichment/code_project/outcome/ld_bootstrap/OC_gene_blocks.txt")
ld_blocks <- fread("C:/Peter/gene_enrichment/code_project/outcome/ld_bootstrap/ld_block.txt")


build_eqtl_ld <- function(all_eQTL) {
  block <- fread("C:/Peter/gene_enrichment/code_project/outcome/ld_bootstrap/gene_blocks_detail.txt")
block <- block[gene_name %in% all_eQTL_gene]
block[, LD_block := NA_character_]
block[
  ld_blocks,
  on = .(
    chr,
    tss >= start,
    tss <= end
  ),
  LD_block := i.LD_block
]
block <- block[,.(gene_name,LD_block)]

  eqtl_ld <- all_eQTL[block, on = .(Gene = gene_name), nomatch = 0]
  eqtl_ld <- eqtl_ld[!is.na(Gene) & Gene != "" & !is.na(LD_block)]
  setkey(eqtl_ld, Gene, pval)
  eqtl_ld[, .SD[1], by = "Gene"]
}

make_ld_block_n <- function(dt, sample_ids, id_col) {
  sample_dt <- dt[dt[[id_col]] %in% sample_ids]
  sample_dt <- sample_dt[!is.na(LD_block)]
  sample_dt <- sample_dt[, .SD[1], by = id_col]
  block_n <- sample_dt[, .(n_sample = .N), by = LD_block]
  if (block_n[, sum(n_sample)] == 0) {
    stop("No sampled IDs have LD_block.")
  }
  block_n
}

sample_ld_ids <- function(background_dt, block_n, id_col) {
  background_dt <- background_dt[!is.na(LD_block)]
  background_dt <- background_dt[, .SD[1], by = id_col]
  block_idx <- split(seq_len(nrow(background_dt)), as.character(background_dt$LD_block))
  block_n[, LD_block := as.character(LD_block)]

  missing_blocks <- setdiff(block_n$LD_block, names(block_idx))
  if (length(missing_blocks) > 0) {
    stop("These LD blocks are missing from sampling background: ", paste(missing_blocks, collapse = ", "))
  }

  sample_idx <- unlist(
    block_n[
      ,
      {
        pool <- block_idx[[LD_block]]
        if (length(pool) < n_sample) {
          stop(
            "LD block ", LD_block,
            " has only ", length(pool),
            " background rows, but n_sample = ", n_sample
          )
        }
        .(idx = list(pool[sample.int(length(pool), n_sample)]))
      },
      by = LD_block
    ]$idx,
    use.names = FALSE
  )

  background_dt[sample_idx, get(id_col)]
}

make_ld_block_param <- function(background_dt, block_n, hit_ids, id_col) {
  background_dt <- copy(background_dt[!is.na(LD_block)])
  background_dt <- background_dt[, .SD[1], by = id_col]
  background_dt[, LD_block := as.character(LD_block)]

  block_n <- copy(block_n)
  block_n[, LD_block := as.character(LD_block)]

  hit_ids <- unique(as.character(hit_ids))
  background_dt[, hit_tmp := as.character(get(id_col)) %chin% hit_ids]

  block_param <- background_dt[
    ,
    .(
      N_bg = .N,
      K_hit = sum(hit_tmp)
    ),
    by = LD_block
  ][
    block_n,
    on = "LD_block"
  ]

  if (anyNA(block_param$N_bg)) {
    print(block_param[is.na(N_bg)])
    stop("Some LD blocks in block_n are not found in background_dt.")
  }
  if (any(block_param$n_sample > block_param$N_bg)) {
    print(block_param[n_sample > N_bg])
    stop("Some LD blocks have n_sample larger than background size.")
  }

  block_param[]
}

simulate_overlap_hyper <- function(block_param,
                                   random_times,
                                   chunk_size = 50000,
                                   seed = 1L) {
  set.seed(seed)

  N_bg <- as.integer(block_param$N_bg)
  K_hit <- as.integer(block_param$K_hit)
  n_sample <- as.integer(block_param$n_sample)

  out <- integer(random_times)
  starts <- seq.int(1L, random_times, by = chunk_size)
  for (st in starts) {
    len <- min(chunk_size, random_times - st + 1L)
    tmp <- integer(len)

    for (b in seq_along(N_bg)) {
      if (n_sample[b] == 0L || K_hit[b] == 0L) {
        next
      }

      if (K_hit[b] == N_bg[b]) {
        tmp <- tmp + n_sample[b]
      } else {
        tmp <- tmp + rhyper(
          nn = len,
          m = K_hit[b],
          n = N_bg[b] - K_hit[b],
          k = n_sample[b]
        )
      }
    }

    out[st:(st + len - 1L)] <- tmp
  }

  out
}

make_ld_random_density <- function(value_dt,
                                   background_dt,
                                   block_n,
                                   id_col,
                                   value_id_col = id_col,
                                   value_col = "pval",
                                   random_times = 1e6,
                                   plot_x_range = c(0, 1),
                                   density_n = 2^12,
                                   quantile_sample_times = min(random_times, 10000L),
                                   seed = 1L) {
  value_small <- copy(value_dt[
    !is.na(get(value_id_col)) & !is.na(get(value_col)),
    .SD[1],
    by = value_id_col
  ])[
    ,
    .(id = as.character(get(value_id_col)), pval = get(value_col))
  ]

  pval_map <- setNames(value_small$pval, value_small$id)

  background_dt <- copy(background_dt[
    as.character(get(id_col)) %in% names(pval_map) &
      !is.na(get(id_col)) &
      !is.na(LD_block)
  ])
  background_dt <- background_dt[, .SD[1], by = id_col]
  background_dt[, LD_block := as.character(LD_block)]
  background_dt[, pval := pval_map[as.character(get(id_col))]]

  if (anyNA(background_dt$pval)) {
    print(background_dt[is.na(pval)])
    stop("Some background IDs do not have p-values.")
  }

  block_n <- copy(block_n)
  block_n[, LD_block := as.character(LD_block)]
  block_size <- background_dt[, .(block_size = .N), by = LD_block]
  block_info <- merge(block_n, block_size, by = "LD_block", all.x = TRUE, sort = FALSE)

  if (anyNA(block_info$block_size)) {
    print(block_info[is.na(block_size)])
    stop("Some LD blocks are in block_n but not in background_dt.")
  }
  if (any(block_info$block_size < block_info$n_sample)) {
    print(block_info[block_size < n_sample])
    stop("Some LD blocks have fewer background rows than required samples.")
  }

  sampling_bg <- merge(
    background_dt[, .(id = as.character(get(id_col)), LD_block, pval)],
    block_info[, .(LD_block, n_sample, block_size)],
    by = "LD_block",
    all.x = FALSE,
    sort = FALSE
  )

  sampling_bg[, weight_raw := n_sample / block_size]
  w <- sampling_bg$weight_raw / sum(sampling_bg$weight_raw)

  setorder(block_info, LD_block)
  pval_pool <- split(sampling_bg$pval, sampling_bg$LD_block)
  pval_pool <- pval_pool[block_info$LD_block]
  n_vec <- as.integer(block_info$n_sample)
  total_n <- sum(n_vec)

  sample_one_pval <- function() {
    out <- numeric(total_n)
    pos <- 1L

    for (b in seq_along(pval_pool)) {
      pool_b <- pval_pool[[b]]
      n_b <- n_vec[b]
      idx_range <- pos:(pos + n_b - 1L)

      if (length(pool_b) == n_b) {
        out[idx_range] <- pool_b
      } else {
        out[idx_range] <- pool_b[sample.int(length(pool_b), n_b)]
      }

      pos <- pos + n_b
    }

    out
  }

  set.seed(seed)
  pilot_bws <- replicate(200L, bw.nrd0(sample_one_pval()))
  pilot_bws <- pilot_bws[is.finite(pilot_bws) & pilot_bws > 0]

  if (length(pilot_bws) == 0) {
    stop("Cannot estimate a valid bandwidth.")
  }

  bw_fixed <- median(pilot_bws)
  random_mean_fit <- density(
    sampling_bg$pval,
    weights = w,
    bw = bw_fixed,
    from = plot_x_range[1],
    to = plot_x_range[2],
    n = density_n
  )

  density_quantile_sample <- matrix(
    NA_real_,
    nrow = density_n,
    ncol = quantile_sample_times
  )

  set.seed(seed)
  for (i in seq_len(quantile_sample_times)) {
    density_quantile_sample[, i] <- density(
      sample_one_pval(),
      bw = bw_fixed,
      from = plot_x_range[1],
      to = plot_x_range[2],
      n = density_n
    )$y

    if (i %% 1000L == 0L) {
      message("Finished ", i, " / ", quantile_sample_times)
    }
  }

  if (requireNamespace("matrixStats", quietly = TRUE)) {
    density_quantile <- matrixStats::rowQuantiles(
      density_quantile_sample,
      probs = c(0.025, 0.975),
      na.rm = TRUE
    )
    sd_density <- matrixStats::rowSds(density_quantile_sample, na.rm = TRUE)
  } else {
    density_quantile <- t(apply(
      density_quantile_sample,
      1,
      quantile,
      probs = c(0.025, 0.975),
      na.rm = TRUE
    ))
    sd_density <- apply(density_quantile_sample, 1, sd, na.rm = TRUE)
  }

  se_density <- sd_density / sqrt(random_times)
  random_density_summary <- data.table(
    x = random_mean_fit$x,
    mean_density = random_mean_fit$y,
    ci_lower = random_mean_fit$y - qt(0.975, random_times - 1) * se_density,
    ci_upper = random_mean_fit$y + qt(0.975, random_times - 1) * se_density,
    q025 = density_quantile[, 1],
    q975 = density_quantile[, 2]
  )

  list(
    summary = random_density_summary,
    bw = bw_fixed,
    sampling_bg = sampling_bg
  )
}


plot_output_dir <- sprintf(
  "C:/Peter/gene_enrichment/code_project/outcome/rawData_eQTL_N/LD_boostrap/%s",
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


gt_N <- fread("C:/Peter/rawData_eQTL/r2_filter_0.8/outcome/raw_N_FDR_R2_0.8.txt")


lung_twas <- fread(
  "C:/Peter/gene_enrichment/code_project/data/lung_國衛院/20260322_115samples_24216probes_ModelPerformance.txt",
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

eQTL_ld <- build_eqtl_ld(all_eQTL)
eQTL_probe_ld <- all_eQTL[eQTL_ld[, .(Gene, LD_block)], on = .(Gene), nomatch = 0]
eQTL_probe_ld <- eQTL_probe_ld[!is.na(probe) & probe != "" & !is.na(LD_block)]
eQTL_probe_ld <- eQTL_probe_ld[, .SD[1], by = "probe"]

all_eQTL_gene_ld <- unique(eQTL_ld$Gene)
all_eQTL_probe_ld <- unique(eQTL_probe_ld$probe)
eQTL_gene_ld <- eQTL_gene[eQTL_gene %in% all_eQTL_gene_ld]
eQTL_probe_ld_sig <- eQTL_probe[eQTL_probe %in% all_eQTL_probe_ld]


# 挑出我們資料有紀錄、lung 紀錄的 probe/gene
common_probe <- intersect(all_eQTL_probe_ld, lung_twas$ProbeID)
common_gene <- intersect(all_eQTL_gene_ld, lung_twas$GeneSymbol)

cat("intersect of eQTL, lung total probe:", intersect(all_eQTL_probe, lung_twas$ProbeID) %>% uniqueN(), "\n")
cat("intersect of eQTL, lung total gene:", intersect(all_eQTL_gene, lung_twas$GeneSymbol) %>% uniqueN(), "\n")

# 在一邊顯著，同時有被紀錄在另一邊 數量
cat("eQTL significant probes recorded in lung TWAS:", uniqueN(gt_N$Probe[gt_N$Probe %in% lung_twas$ProbeID]), "\n")
cat("eQTL significant genes recorded in lung TWAS:", uniqueN(gt_N$Gene[gt_N$Gene %in% lung_twas$GeneSymbol]), "\n")
cat("lung significant probes recorded in eQTL:", lung_sig$ProbeID[lung_sig$ProbeID %in% all_eQTL_probe] %>% uniqueN(), "\n")
cat("lung significant genes recorded in eQTL:", lung_sig$GeneSymbol[lung_sig$GeneSymbol %in% all_eQTL_gene] %>% uniqueN(), "\n")


# 交集數量
observed_intersect_probe <- intersect(eQTL_probe_ld_sig, lung_twas[top1_p < 0.01, ProbeID]) %>% uniqueN()
observed_intersect_gene <- intersect(eQTL_gene_ld, lung_twas[top1_p < 0.01, GeneSymbol]) %>% uniqueN()
cat("intersect of eQTL, lung sig probe:", observed_intersect_probe, "\n")
cat("intersect of eQTL, lung sig gene:", observed_intersect_gene, "\n")

sample_probe_number <- uniqueN(eQTL_probe_ld_sig[eQTL_probe_ld_sig %in% lung_twas$ProbeID])
sample_gene_number <- uniqueN(eQTL_gene_ld[eQTL_gene_ld %in% lung_twas$GeneSymbol])

lung_probe_bg <- eQTL_probe_ld[probe %in% lung_twas$ProbeID]
lung_gene_bg <- eQTL_ld[Gene %in% lung_twas$GeneSymbol]
lung_probe_block_n <- make_ld_block_n(
  lung_probe_bg,
  eQTL_probe_ld_sig[eQTL_probe_ld_sig %in% lung_twas$ProbeID],
  "probe"
)
lung_gene_block_n <- make_ld_block_n(
  lung_gene_bg,
  eQTL_gene_ld[eQTL_gene_ld %in% lung_twas$GeneSymbol],
  "Gene"
)







## 選跟 eQTL sig. 同數量的，跟 lung_twas 交集數量 hist ----
# 隨機選 105 probe, 91 gene，多少顯著 in lung?
lung_sig_probe_sub <- lung_sig$ProbeID[lung_sig$ProbeID %in% all_eQTL_probe_ld]
lung_sig_gene_sub <- lung_sig$GeneSymbol[lung_sig$GeneSymbol %in% all_eQTL_gene_ld]

random_times <- 1e6
repeat_probe_number <- simulate_overlap_hyper(
  make_ld_block_param(lung_probe_bg, lung_probe_block_n, lung_sig_probe_sub, "probe"),
  random_times = random_times,
  chunk_size = 50000,
  seed = 1L
)
repeat_gene_number <- simulate_overlap_hyper(
  make_ld_block_param(lung_gene_bg, lung_gene_block_n, lung_sig_gene_sub, "Gene"),
  random_times = random_times,
  chunk_size = 50000,
  seed = 1L
)
# aa
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
  main = sprintf("%s times Probe Intersection Size between Lung Sig. and eQTL Sig.", random_times),
  xlab = "Intersection Size each Times",
  col = "skyblue"
)
dev.off()

png(file.path(plot_output_dir[1], sprintf("random_%sgene.png", sample_gene_number)),
  width = 8, height = 6, units = "in", res = 300
)
hist(repeat_gene_number,
  main = sprintf("%s times Gene Intersection Size between Lung Sig. and eQTL Sig.", random_times),
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
  sampled_probe <- sample_ld_ids(lung_probe_bg, lung_probe_block_n, "probe")
  b_tmp <- lung_twas[ProbeID %in% sampled_probe, ][, .(ProbeID, pval = top1_p)]

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



# 畫 random dist 的平均與 quantile band ----
lung_probe_small <- lung_twas[
  ProbeID %in% common_probe & !is.na(top1_p),
  .(ProbeID, pval = top1_p)
]
setorder(lung_probe_small, pval)
lung_probe_small <- lung_probe_small[, .SD[1], by = ProbeID]
random_times <- 1e6


plot_x_range <- c(0, 1)
# 不同 type 的 row分組，每組的 pdf 從 from 到 to 切成n份，用 預設 gaussian kernel 估計 pdf
density_n <- 2^12
quantile_sample_times <- min(random_times, 10000L)
probe_density_random <- make_ld_random_density(
  value_dt = lung_probe_small,
  background_dt = lung_probe_bg,
  block_n = lung_probe_block_n,
  id_col = "probe",
  value_id_col = "ProbeID",
  random_times = random_times,
  plot_x_range = plot_x_range,
  density_n = density_n,
  quantile_sample_times = quantile_sample_times,
  seed = 1L
)
random_density_summary <- probe_density_random$summary



eqtl_density_fit <- density(
  a[, pval],
  bw = probe_density_random$bw,
  from = plot_x_range[1],
  to = plot_x_range[2],
  n = density_n
)
eqtl_density <- data.table(x = eqtl_density_fit$x, density = eqtl_density_fit$y)

all_pval <- c(a[, pval], probe_density_random$sampling_bg$pval)
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
  sampled_gene <- sample_ld_ids(lung_gene_bg, lung_gene_block_n, "Gene")
  b_tmp <- lung_twas[GeneSymbol %in% sampled_gene, ][, .(GeneSymbol, pval = top1_p)]
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




# 畫 random dist 的平均與 quantile band ----
lung_gene_small <- lung_twas[
  GeneSymbol %in% common_gene & !is.na(top1_p),
  .(GeneSymbol, pval = top1_p)
]
setorder(lung_gene_small, pval)
lung_gene_small <- lung_gene_small[, .SD[1], by = GeneSymbol]
random_times <- 1e6


plot_x_range <- c(0, 1)
# 不同 type 的 row分組，每組的 pdf 從 from 到 to 切成n份，用 預設 gaussian kernel 估計 pdf
density_n <- 2^12
quantile_sample_times <- min(random_times, 10000L)
gene_density_random <- make_ld_random_density(
  value_dt = lung_gene_small,
  background_dt = lung_gene_bg,
  block_n = lung_gene_block_n,
  id_col = "Gene",
  value_id_col = "GeneSymbol",
  random_times = random_times,
  plot_x_range = plot_x_range,
  density_n = density_n,
  quantile_sample_times = quantile_sample_times,
  seed = 1L
)
random_density_summary <- gene_density_random$summary



eqtl_density_fit <- density(
  a[, pval],
  bw = gene_density_random$bw,
  from = plot_x_range[1],
  to = plot_x_range[2],
  n = density_n
)
eqtl_density <- data.table(x = eqtl_density_fit$x, density = eqtl_density_fit$y)

all_pval <- c(a[, pval], gene_density_random$sampling_bg$pval)
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

}









{

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
  gt_N <- fread("C:/Peter/rawData_eQTL/r2_filter_0.8/outcome/raw_N_FDR_R2_0.8.txt")

  block <- fread("C:/Peter/gene_enrichment/code_project/outcome/ld_bootstrap/gene_blocks_detail.txt")
  block_gene <- unique(block[!is.na(gene_name) & gene_name != "", gene_name])

  tcga_weight <- tcga_weight[Gene %in% block_gene]
  tcga_gene_bg <- eQTL_ld[Gene %in% tcga_weight$Gene & Gene %in% block_gene]
  eQTL_gene <- unique(gt_N$Gene)
  eQTL_gene <- eQTL_gene[eQTL_gene %in% tcga_gene_bg$Gene & eQTL_gene %in% block_gene]

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


  tcga_sig_sub <- tcga_sigGene[tcga_sigGene %in% tcga_gene_bg$Gene]
  sample_gene_number <- uniqueN(eQTL_gene)
  observed_intersect <- intersect(eQTL_gene, tcga_sig_sub) %>% uniqueN()
  tcga_gene_block_n <- make_ld_block_n(tcga_gene_bg, eQTL_gene, "Gene")

  cat(
    cancer_type, "significant genes recorded in eQTL:",
    tcga_sigGene[tcga_sigGene %in% all_eQTL_gene] %>% uniqueN(), "\n"
  )

  cat(
    "intersect of eQTL,", cancer_type, "sig gene:",
    observed_intersect, "\n"
  )


  random_times <- 1e6
  repeat_gene_number <- simulate_overlap_hyper(
    make_ld_block_param(tcga_gene_bg, tcga_gene_block_n, tcga_sig_sub, "Gene"),
    random_times = random_times,
    chunk_size = 50000,
    seed = 1L
  )

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
    main = sprintf("%s times Gene Intersection Size between %s Sig. and eQTL Sig.", random_times, cancer_type),
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
    sampled_gene <- sample_ld_ids(tcga_gene_bg, tcga_gene_block_n, "Gene")
    b_tmp <- all_eQTL[Gene %in% sampled_gene, ][, .(Gene, pval)]
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


  # 畫 random dist 的平均與 quantile band ----
  random_times <- 1e6
  plot_x_range <- c(0, 1)
  density_n <- 2^12
  quantile_sample_times <- min(random_times, 10000L)
  eqtl_small <- all_eQTL[
    Gene %in% tcga_gene_bg$Gene & !is.na(Gene) & Gene != "" & !is.na(pval),
    .SD[1],
    by = Gene
  ][, .(Gene, pval)]
  tcga_density_random <- make_ld_random_density(
    value_dt = eqtl_small,
    background_dt = tcga_gene_bg,
    block_n = tcga_gene_block_n,
    id_col = "Gene",
    value_id_col = "Gene",
    random_times = random_times,
    plot_x_range = plot_x_range,
    density_n = density_n,
    quantile_sample_times = quantile_sample_times,
    seed = 1L
  )
  random_density_summary <- tcga_density_random$summary

  eqtl_density_fit <- density(
    a[, pval],
    bw = tcga_density_random$bw,
    from = plot_x_range[1],
    to = plot_x_range[2],
    n = density_n
  )
  eqtl_density <- data.table(x = eqtl_density_fit$x, density = eqtl_density_fit$y)

  all_pval <- c(a[, pval], tcga_density_random$sampling_bg$pval)
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
eQTL_ld <- build_eqtl_ld(all_eQTL)






# TCGA-HNSC ----
run_tcga_eqtl_enrichment(
  cancer_type = "TCGA-HNSC",
  weight_path = "C:/Peter/gene_enrichment/code_project/data/TCGA-HNSC.TUMOR/fusion_project/fusion_weights_summary.csv",
  pos_path = "C:/Peter/gene_enrichment/code_project/data/TCGA-HNSC.TUMOR/TCGA-HNSC.TUMOR.pos",
  output_dir = plot_output_dir[2]
)


# TCGA-LUAD ----
run_tcga_eqtl_enrichment(
  cancer_type = "TCGA-LUAD",
  weight_path = "C:/Peter/gene_enrichment/code_project/data/TCGA-LUAD.TUMOR/fusion_project/fusion_weights_summary.csv",
  pos_path = "C:/Peter/gene_enrichment/code_project/data/TCGA-LUAD.TUMOR/TCGA-LUAD.TUMOR.pos",
  output_dir = plot_output_dir[3]
)


# TCGA-LUSC ----
run_tcga_eqtl_enrichment(
  cancer_type = "TCGA-LUSC",
  weight_path = "C:/Peter/gene_enrichment/code_project/data/TCGA-LUSC.TUMOR/fusion_project/fusion_weights_summary.csv",
  pos_path = "C:/Peter/gene_enrichment/code_project/data/TCGA-LUSC.TUMOR/TCGA-LUSC.TUMOR.pos",
  output_dir = plot_output_dir[4]
)



}




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
# intersect of eQTL, lung sig probe: 105 
# intersect of eQTL, lung sig gene: 105 
# 抽 147 個，平均顯著 in lung probe number: 28.08365 
# 抽 145 個，平均顯著 in lung gene number: 33.98918 
# 抽 147 個，顯著 in lung probe number 範圍: 8 54 
# 抽 145 個，顯著 in lung gene number 範圍: 10 57 
# sample 147 genes, sig in lung probe number >= 105 times: 0 
# sample 145 genes, sig in lung gene number >= 105 times: 0 


# After remove R2<0.8 snp, run eQTL genes: 15992 
# After remove R2<0.8 snp, run eQTL probes: 20906 
# TCGA-HNSC genes: 2767 
# TCGA-HNSC After remove NULL and NA GeneID/Top1_Pval, genes: 2763 
# TCGA-HNSC sig. genes: 1689 
# eQTL sig. genes: 45 
# TCGA-HNSC significant genes recorded in eQTL: 1412 
# intersect of eQTL, TCGA-HNSC sig gene: 40 
# sample 45 genes, mean sig in TCGA-HNSC gene number: 35.29268 
# sample 45 genes, sig in TCGA-HNSC gene number range: 25 42 
# sample 45 genes, sig in TCGA-HNSC gene number >= 40 times: 12665 
 
# After remove R2<0.8 snp, run eQTL genes: 15992 
# After remove R2<0.8 snp, run eQTL probes: 20906 
# TCGA-LUAD genes: 2978 
# TCGA-LUAD After remove NULL and NA GeneID/Top1_Pval, genes: 2974 
# TCGA-LUAD sig. genes: 1822 
# eQTL sig. genes: 50 
# TCGA-LUAD significant genes recorded in eQTL: 1509 
# intersect of eQTL, TCGA-LUAD sig gene: 44 
# sample 50 genes, mean sig in TCGA-LUAD gene number: 41.45204 
# sample 50 genes, sig in TCGA-LUAD gene number range: 30 49 
# sample 50 genes, sig in TCGA-LUAD gene number >= 44 times: 164656 
 
# After remove R2<0.8 snp, run eQTL genes: 15992 
# After remove R2<0.8 snp, run eQTL probes: 20906 
# TCGA-LUSC genes: 2548 
# TCGA-LUSC After remove NULL and NA GeneID/Top1_Pval, genes: 2545 
# TCGA-LUSC sig. genes: 1493 
# eQTL sig. genes: 51 
# TCGA-LUSC significant genes recorded in eQTL: 1241 
# intersect of eQTL, TCGA-LUSC sig gene: 44 
# sample 51 genes, mean sig in TCGA-LUSC gene number: 38.11943 
# sample 51 genes, sig in TCGA-LUSC gene number range: 27 46 
# sample 51 genes, sig in TCGA-LUSC gene number >= 44 times: 2872 