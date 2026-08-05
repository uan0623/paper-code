library(data.table)
library(dplyr)
library(tidyverse)
library(stringr)
library(ggrepel)
library(matrixStats) # rowMins
library(qvalue)
library(bestNormalize)
library(ggplot2)
library(openxlsx)
library(magrittr)




gt <- fread("C:/Peter/rawData_eQTL/trash/raw_maf_gt_N_pvalue.txt")
gene_pos <- fread("C:/Peter/oral_cancer/expression/expression_data/probe24526infoB.txt", header = T)


# 按基因、pvalue 排序
setkey(gt, gene, `p-value`)
gt <- gt[, .SD[`p-value` == min(`p-value`, na.rm = TRUE)], by = "gene"]
gt_lots_asso <- gt[, if (.N > 1) .SD, by = gene]
uniqueN(gt$gene)
lead_snp <- gt[, if (.N == 1) .SD, by = gene] %>%
  unique(., by = "SNP") %>%
  dplyr::select(SNP)



lead_vec <- unique(lead_snp$SNP)
gt[, in_lead := fifelse(SNP %in% lead_vec, 1L, 0L)]

# 把 in_lead=1 的，也就是不只出現一次的 probe，snp 與只出現一次的 probe 的 snp 重疊，優先選
gt_one <- gt[
  order(gene, -in_lead, `p-value`)
][
  , .SD[1L],
  by = gene
][
  , in_lead := NULL
]


uniqueN(gt_one$gene)
uniqueN(gt_one$SNP)
gt_one[, .N, by = SNP]$N %>% table()


# liftOver 轉換 pos

gene_pos <- fread("C:/Peter/oral_cancer/expression/expression_data/probe_pos_mulInterval.txt", header = T)
gene_probe <- fread("C:/Peter/oral_cancer/expression/expression_data/probe24526infoB.txt", header = T)
gene_pos <- merge(gene_pos[, .(PROBE_ID = Gene, chr = CHROMOSOME, start, end)], gene_probe[, .(PROBE_ID, TargetID)])

fwrite(
  data.table(
    a = paste0("chr", gene_pos$chr),
    b = gene_pos$start,
    c = gene_pos$end,
    d =  gene_pos$PROBE_ID
  ),
  "C:/Peter/gene_enrichment/code_project/data/ld_bootstrap/gene_hg18.txt",
  row.names = F, col.names = F, sep = "\t"
)


cmd <- paste(
  "/mnt/c/Peter/oral_cancer/liftover/liftOver",
  "/mnt/c/Peter/gene_enrichment/code_project/data/ld_bootstrap/gene_hg18.txt",
  "/mnt/c/Peter/oral_cancer/liftover/hg18ToHg19.over.chain",
  "/mnt/c/Peter/gene_enrichment/code_project/data/ld_bootstrap/gene_hg18TOhg19.bed",
  "/mnt/c/Peter/gene_enrichment/code_project/data/ld_bootstrap/hg18TOhg19_unmapped.bed"
)

system2("wsl", cmd)

a <- fread("C:/Peter/gene_enrichment/code_project/data/ld_bootstrap/gene_hg18TOhg19.bed", header = F)
# 忽略轉失敗 .bed 檔案的 #Deleted in new 該行

unmapped_path <- "C:/Peter/gene_enrichment/code_project/data/ld_bootstrap/hg18TOhg19_unmapped.bed"

if (!file.exists(unmapped_path) || file.info(unmapped_path, extra_cols = FALSE)$size == 0) {
  a_unmapped <- data.table(V1 = character(), V2 = integer(), V3 = integer(),V4 = character())
} else {
  a_unmapped <- fread(
    unmapped_path,
    header = FALSE,
    comment.char = "#"
  )
}

col_order <- names(gene_pos)
names(a_unmapped) <- c("chr", "start", "end", "PROBE_ID")
names(a) <- c("chr_hg19", "start_hg19", "end_hg19","PROBE_ID")
a_unmapped[, chr := sub("^chr", "", chr) %>% as.numeric()]
a[, chr_hg19 := sub("^chr", "", chr_hg19) %>% as.numeric()]



# 新增轉換後 pos
a_unmapped[, c("chr_hg19", "start_hg19", "end_hg19") := 0]
final <- merge(gene_pos, a_unmapped, all.x = T)

success_idx <- is.na(final$chr_hg19)
final[
  success_idx,
  c("chr_hg19", "start_hg19", "end_hg19") := a[, .(chr_hg19, start_hg19, end_hg19)]
]

# 轉換紀錄
# fwrite(final,
#   "C:/Peter/gene_enrichment/code_project/data/hg38_hg18.txt",
#   row.names = F, col.names = T, sep = "\t"
# )


# 下載 grch37 annotation ----
library(data.table)

gtf <- fread(
  "C:/Peter/gene_enrichment/code_project/data/ld_bootstrap/gencode.v37lift37.annotation.gtf",
  sep = "\t",
  header = FALSE,
  quote = "",
  data.table = T
)

colnames(gtf) <- c(
  "chr", "source", "feature", "start", "end",
  "score", "strand", "frame", "attribute"
)

gene_gtf <- gtf[feature == "gene"]
extract_attr <- function(x, key) {
  sub(
    paste0(".*", key, " \"([^\"]+)\".*"),
    "\\1",
    x
  )
}

gene_gtf[
  ,
  `:=`(
    gene_id = extract_attr(attribute, "gene_id"),
    gene_name = extract_attr(attribute, "gene_name"),
    gene_type = extract_attr(attribute, "gene_type")
  )
][
  ,
  `:=`(
    ensembl_id = sub("\\..*", "", gene_id),
    tss = fifelse(strand == "-", end, start),
    chr = gsub("^chr", "", chr) %>% as.numeric()
  )
]

gene_gtf <- gene_gtf[chr %in% c(1:22)]

fwrite(gene_gtf,
  "C:/Peter/gene_enrichment/code_project/outcome/ld_bootstrap/gene_blocks_detail.txt",
  row.names = F, col.names = T, sep = "\t"
)


# under 0.8 threshold, 有用 cis snp 的 probe ----
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


# 挑出我們資料有 cis snp 的 gene
gene_gtf <- gene_gtf[gene_name %in% all_eQTL_gene]



# 下載 LD block https://bitbucket.org/nygcresearch/ldetect-data.git，hg19 (same as 1000 Genomes phase 1)
ld_blocks <- fread(
  "C:/Peter/gene_enrichment/code_project/data/ld_bootstrap/ldetect-data/ASN/fourier_ls-all.bed"
)
colnames(ld_blocks)[1:3] <- c("chr", "start", "end")
setDT(ld_blocks)

ld_blocks[, chr := gsub("^chr", "", chr) %>% as.numeric()]

ld_blocks[, LD_block := paste0(
  "ASN_chr", chr, "_",
  start, "_",
  end
)]

gene_gtf[, LD_block := NA_character_]

gene_gtf[
  ld_blocks,
  on = .(
    chr,
    tss >= start,
    tss <= end
  ),
  LD_block := i.LD_block
]




png("C:/Peter/gene_enrichment/code_project/outcome/ld_bootstrap/LD_block_box.png",
  width = 1000, height = 1000, res = 200
)
gene_gtf[, .N, by = LD_block]$N %>% boxplot(main = "Gene number in each LD_block")
dev.off()

fwrite(ld_blocks,
  "C:/Peter/gene_enrichment/code_project/outcome/ld_bootstrap/ld_block.txt",
  row.names = F, col.names = T, sep = "\t"
)


eQTL_ld <- all_eQTL[gene_gtf[,.(gene_name,LD_block)], on=.(Gene=gene_name), nomatch=0]
setkey(eQTL_ld, Gene, pval)
eQTL_ld <- eQTL_ld[, .SD[1], by = "Gene"]
fwrite(eQTL_ld,
  "C:/Peter/gene_enrichment/code_project/outcome/ld_bootstrap/OC_gene_blocks.txt",
  row.names = F, col.names = T, sep = "\t"
)




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

    sig_block <- eQTL_ld[
  Gene %in% eQTL_gene &
    Gene %in% gtex_weight$gene &
    !is.na(LD_block)
]

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

  # 限縮有 ld block 的 gene
  block <- fread("C:/Peter/gene_enrichment/code_project/outcome/ld_bootstrap/gene_blocks_detail.txt")
  gtex_weight <-  gtex_weight[gene %in% block$gene_name]


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

  # 限縮有 ld block 的 gene
  block <- fread("C:/Peter/gene_enrichment/code_project/outcome/ld_bootstrap/gene_blocks_detail.txt")
  all_eQTL_gene <-  all_eQTL_gene[all_eQTL_gene %in% block$gene_name]

  eQTL_gene <- unique(all_eQTL[Gene %in% gt_N$Gene, Gene])
  eQTL_gene <-  eQTL_gene[eQTL_gene %in% block$gene_name]

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


  
  if(sig_block[is.na(LD_block), .N] !=0 | eQTL_ld[is.na(LD_block), .N] !=0){
    stop("LD_block shouldn't have NA")
  }

  # 從共同基因隨機挑，看多少 在 gtex 顯著
  gtex_sig_sub <- gtex_sigGene[gtex_sigGene %in% all_eQTL_gene]


setDT(eQTL_ld)
setDT(sig_block)

random_times <- 1e6
block_n <- sig_block[, .(n_sample = .N), by = LD_block]
gtex_sig_sub <- unique(as.character(gtex_sig_sub))
eQTL_ld[, hit_tmp := as.character(Gene) %chin% gtex_sig_sub]
block_param <- eQTL_ld[
  ,
  .(
    N_bg = .N,
    K_hit = sum(hit_tmp)
  ),
  by = LD_block
]

eQTL_ld[, hit_tmp := NULL]

block_param <- block_param[
  block_n,
  on = "LD_block"
]

if (anyNA(block_param$N_bg)) {
  stop("Some LD blocks in sig_block are not found in eQTL_ld.")
}

if (any(block_param$n_sample > block_param$N_bg)) {
  print(block_param[n_sample > N_bg])
  stop("Some LD blocks have n_sample larger than background size.")
}




# simulation ----
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
        # N_bg 中抽 n_sample 個，其中 K_hit 是有重複在 GTEx 顯著的個數，模擬 len 次
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

repeat_gene_number <- simulate_overlap_hyper(
  block_param = block_param,
  random_times = random_times,
  chunk_size = 50000,
  seed = 1L
)




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

  cat(rep("\n", 2))
  png(file.path(output_dir, sprintf("random_%sgene.png", sample_gene_number)),
    width = 8, height = 6, units = "in", res = 300
  )
  hist(repeat_gene_number,
    main = sprintf("%s times Gene Intersection Size between %s Sig. and eQTL Sig.", random_times, tissue_name),
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
block_n <- sig_block[, .(n_sample = .N), by = LD_block]
block_idx <- split(seq_len(nrow(eQTL_ld)), eQTL_ld$LD_block)

for (i in seq_len(random_times)) {
  set.seed(i)
  sample_idx <- unlist(
    block_n[
      ,
      .(idx = list(sample(block_idx[[LD_block]], n_sample))),
      by = LD_block
    ]$idx
  )
  common_gene <- eQTL_ld[sample_idx, Gene]
  b_tmp <- gtex_weight[
    gene %in% common_gene,
  ][
    ,
    .(gene_name = gene, pval = Top1_Pval)
  ]
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
    png(file.path(output_dir, sprintf("%sgene_Top1_PvalDIST_%s.png", sample_gene_number, format(x_cutoff, digits = 2, scientific = TRUE))),
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








# 畫 random dist 的平均與 quantile band ----
gtex_small <- gtex_weight[, .(gene_name = gene, pval = Top1_Pval)]
gtex_small <- gtex_small[!is.na(gene_name) & !is.na(pval)]

## 如果 gene_name 有重複，match() 原本會抓第一個
## 這裡保留第一個，讓行為接近你原本的 code
gtex_small <- gtex_small[!duplicated(gene_name)]

random_times <- 1e6
plot_x_range <- c(0, 1)
density_n <- 2^12
quantile_sample_times <- min(random_times, 10000L)

gtex_pval_map <- setNames(gtex_small$pval, gtex_small$gene_name)

## 只保留 GTEx 有 p-value 的 gene 當抽樣背景
eQTL_ld_bg <- copy(eQTL_ld[
  Gene %in% gtex_small$gene_name &
    !is.na(Gene) &
    !is.na(LD_block)
])

sig_block_bg <- copy(sig_block[
  Gene %in% gtex_small$gene_name &
    !is.na(Gene) &
    !is.na(LD_block)
])

eQTL_ld_bg[, LD_block := as.character(LD_block)]
sig_block_bg[, LD_block := as.character(LD_block)]

## 直接把 GTEx p-value 加到背景資料
eQTL_ld_bg[, pval := gtex_pval_map[Gene]]

if (anyNA(eQTL_ld_bg$pval)) {
  print(eQTL_ld_bg[is.na(pval)])
  stop("Some background genes do not have GTEx p-values.")
}

## 每個 LD block 要抽幾個 row
block_n <- sig_block_bg[
  ,
  .(n_sample = .N),
  by = LD_block
]

## 每個 LD block 在 background 中有幾個 row
block_size <- eQTL_ld_bg[
  ,
  .(block_size = .N),
  by = LD_block
]

block_info <- merge(
  block_n,
  block_size,
  by = "LD_block",
  all.x = TRUE,
  sort = FALSE
)

## 檢查 block 是否存在
if (anyNA(block_info$block_size)) {
  print(block_info[is.na(block_size)])
  stop("Some LD blocks are in sig_block_bg but not in eQTL_ld_bg.")
}

## 檢查每個 block 是否夠抽
if (any(block_info$block_size < block_info$n_sample)) {
  print(block_info[block_size < n_sample])
  stop("Some LD blocks have fewer background rows than required samples.")
}

## 只保留會被抽樣的 LD block
sampling_bg <- merge(
  eQTL_ld_bg[, .(Gene, LD_block, pval)],
  block_info[, .(LD_block, n_sample, block_size)],
  by = "LD_block",
  all.x = FALSE,
  sort = FALSE
)

## ============================================================
## 1. 用 weighted density 直接算 random mean density
## ============================================================

## 每個 row 被抽到的 raw probability
sampling_bg[, weight_raw := n_sample / block_size]

## density() 的 weights 需要 sum to 1
w <- sampling_bg$weight_raw / sum(sampling_bg$weight_raw)

## ============================================================
## 2. 準備抽樣 pool，只給 quantile band 用
## ============================================================

setorder(block_info, LD_block)

pval_pool <- split(
  sampling_bg$pval,
  sampling_bg$LD_block
)

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

## ============================================================
## 3. 固定 bandwidth
## ============================================================

## 用少量 pilot random samples 估一個代表性的 bandwidth
## 這比每次 density() 自動估 bw 更快，也比較適合比較 random vs eQTL
set.seed(1)

pilot_bws <- replicate(
  200L,
  bw.nrd0(sample_one_pval())
)

pilot_bws <- pilot_bws[is.finite(pilot_bws) & pilot_bws > 0]

if (length(pilot_bws) == 0) {
  stop("Cannot estimate a valid bandwidth.")
}

bw_fixed <- median(pilot_bws)

## random mean density：不用 1e6 次 simulation，直接算 expected density
random_mean_fit <- density(
  sampling_bg$pval,
  weights = w,
  bw = bw_fixed,
  from = plot_x_range[1],
  to = plot_x_range[2],
  n = density_n
)

density_x <- random_mean_fit$x
mean_density <- random_mean_fit$y

## ============================================================
## 4. 只用 10,000 次 simulation 算 q025 / q975
## ============================================================

density_quantile_sample <- matrix(
  NA_real_,
  nrow = density_n,
  ncol = quantile_sample_times
)

set.seed(1)

for (i in seq_len(quantile_sample_times)) {
  pval_tmp <- sample_one_pval()
  
  density_quantile_sample[, i] <- density(
    pval_tmp,
    bw = bw_fixed,
    from = plot_x_range[1],
    to = plot_x_range[2],
    n = density_n
  )$y
  
  if (i %% 1000L == 0L) {
    message("Finished ", i, " / ", quantile_sample_times)
  }
}

## 算 quantile 和 sd
if (requireNamespace("matrixStats", quietly = TRUE)) {
  density_quantile <- matrixStats::rowQuantiles(
    density_quantile_sample,
    probs = c(0.025, 0.975),
    na.rm = TRUE
  )
  
  sd_density <- matrixStats::rowSds(
    density_quantile_sample,
    na.rm = TRUE
  )
} else {
  density_quantile <- t(apply(
    density_quantile_sample,
    1,
    quantile,
    probs = c(0.025, 0.975),
    na.rm = TRUE
  ))
  
  sd_density <- apply(
    density_quantile_sample,
    1,
    sd,
    na.rm = TRUE
  )
}

## 這裡的 se 是估計「如果真的跑 random_times 次」mean density 的 Monte Carlo SE
## 不過 mean_density 已經由 weighted density 直接算出，所以實務上主要看 q025/q975
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
  bw = bw_fixed,
  from = plot_x_range[1],
  to = plot_x_range[2],
  n = density_n
)

eqtl_density <- data.table(
  x = eqtl_density_fit$x,
  density = eqtl_density_fit$y
)

all_pval <- c(a[, pval], gtex_small[common_idx, pval])
x_auto <- quantile(all_pval, probs = 0.9, na.rm = TRUE)
x_auto <- max(x_auto, 1e-4)
p_start <- ceiling(log10(x_auto))
pow10_large <- 10^seq(0, p_start, by = -1)

for (x_cutoff in c(pow10_large, x_auto)) {
  png(
    file.path(
      output_dir,
      sprintf(
        "%sgene_confi_xlim_%s.png",
        sample_gene_number,
        format(x_cutoff, digits = 2, scientific = TRUE)
      )
    ),
    width = 8,
    height = 6,
    units = "in",
    res = 300
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
      scale_color_manual(
        values = c(
          "eQTL" = "red",
          "Random mean" = "gray40"
        )
      ) +
      coord_cartesian(xlim = c(0, x_cutoff)) +
      labs(
        title = sprintf(
          "Mean Density Plot for %s Gene Top1_Pval vs %s Random",
          sample_gene_number,
          random_times
        ),
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
#  Rscript /mnt/c/Peter/gene_enrichment/code_project/data/summarize_weights_LOTStissue.R



base_dir <- "C:/Peter/gene_enrichment/code_project/data"
out_dir <- "C:/Peter/gene_enrichment/code_project/outcome/rawData_eQTL_N/LD_boostrap"



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


eQTL_ld <- fread("C:/Peter/gene_enrichment/code_project/outcome/ld_bootstrap/OC_gene_blocks.txt")
ld_blocks <- fread("C:/Peter/gene_enrichment/code_project/outcome/ld_bootstrap/ld_block.txt")


for (i in seq_along(gtex$tissue_dir)) {
  gene_enrichment(
    tissue_name = gtex$tissue_dir[i],
    weight_path = output_csv[i],
    output_dir = plot_output_dir[i]
  )
}



# test aa ----
exp_N <- fread("C:/Peter/rawData_eQTL/raw_Exp_mulInterval_N.txt", header = T)
# exp_N <- t(exp_N) %>% as.data.table()
# names(exp_N) <- as.character(exp_N[1, ])
# exp_N <- exp_N[-c(1), ]
# exp_N[, (names(exp_N)) := lapply(.SD, as.numeric)]

gene_pos <- fread("D:/oral_cancer/expression/expression_data/probe_pos_mulInterval.txt", header = T)
setDT(gene_pos)

gene_pos <- gene_pos[order(CHROMOSOME, start, end)]
gene_pos[
  ,
  gap_prev := start - data.table::shift(end),
  by = CHROMOSOME
]

# 前一個 gene 距離 < 1e6，維持同一 group
gene_pos[
  ,
  group_id := cumsum(is.na(gap_prev) | gap_prev > 1e6),
  by = CHROMOSOME
]

gene_pos[
  ,
  group_id := paste0("chr", CHROMOSOME, "_group", group_id)
]
gene_pos[, gap_prev := NULL]



# test ----
library(data.table)

{

# 如果 exp_N 第一欄是 Gene/probe ID，例如欄名叫 Gene
exp_mat <- as.matrix(exp_N[, -"PROBE_ID"])
rownames(exp_mat) <- exp_N$PROBE_ID

# 如果 exp_N 已經是 matrix，且 rownames 是 Gene/probe ID
exp_mat <- as.matrix(exp_N)

calc_group_cor <- function(group_dt, exp_mat, cor_cutoff = 0.7) {
group_dt <- gene_pos[group_id %in% "chr1_group1"]
exp_mat <- exp_mat
  genes <- unique(group_dt$Gene)
  genes <- genes[genes %in% rownames(exp_mat)]

  if (length(genes) < 2) {
    return(NULL)
  }

  gene_pairs <- t(combn(genes, 2))
  result_list <- vector("list", nrow(gene_pairs))

  for (i in seq_len(nrow(gene_pairs))) {
    g1 <- gene_pairs[i, 1]
    g2 <- gene_pairs[i, 2]

    x <- as.numeric(exp_mat[g1, ])
    y <- as.numeric(exp_mat[g2, ])

    ok <- complete.cases(x, y)
    n_sample <- sum(ok)

    if (n_sample < 3) {
      next
    }

    r <- cor(x[ok], y[ok], method = "pearson")

    if (is.na(r) || abs(r) <= cor_cutoff) {
      next
    }

    pos1 <- group_dt[Gene == g1][1]
    pos2 <- group_dt[Gene == g2][1]

    distance <- max(
      0L,
      max(pos1$start, pos2$start) - min(pos1$end, pos2$end)
    )

    result_list[[i]] <- data.table(
      group_id = pos1$group_id,
      CHROMOSOME = pos1$CHROMOSOME,
      Gene_1 = g1,
      Gene_2 = g2,
      start_1 = pos1$start,
      end_1 = pos1$end,
      start_2 = pos2$start,
      end_2 = pos2$end,
      distance = distance,
      cor = r,
      abs_cor = abs(r),
      n_sample = n_sample
    )
  }

  rbindlist(result_list, fill = TRUE)
}

cor_result <- gene_pos[
  ,
  calc_group_cor(.SD, exp_mat, cor_cutoff = 0),
  by = group_id
]


cor_group_summary <- cor_result[
  ,
  .(
    n_pair_cor_gt_0.7 = .N,
    max_abs_cor = max(abs_cor, na.rm = TRUE),
    mean_abs_cor = mean(abs_cor, na.rm = TRUE)
  ),
  by = .(group_id, CHROMOSOME)
]

group_size <- gene_pos[
  ,
  .(
    n_gene = uniqueN(Gene),
    n_pair_tested = choose(uniqueN(Gene), 2)
  ),
  by = .(group_id, CHROMOSOME)
]

cor_group_summary <- group_size[
  cor_group_summary,
  on = .(group_id, CHROMOSOME)
]

}












# test
{

  make_group_id <- function(gap_prev, threshold = 1e6) {
  acc <- 0
  gid <- integer(length(gap_prev))
  current_group <- 1L

  for (i in seq_along(gap_prev)) {
    gap <- gap_prev[i]
    if (is.na(gap)) {
      acc <- 0
    } else {
      acc <- acc + gap
      if (acc > threshold) {
        current_group <- current_group + 1L
        acc <- 0
      }
    }
    gid[i] <- current_group
  }
  gid
}



cor_summary <- vector("list", length(seq(0,1,by=0.05)[-1]))
  for (j in seq(0,1,by=0.05)[-c(1,21)] %>% seq_along()) {
   
    cor_cutoff <- seq(0,1,by=0.05)[-1][j]
exp_N <- fread("C:/Peter/rawData_eQTL/raw_Exp_mulInterval_N.txt", header = T)
exp_N <- t(exp_N) %>% as.data.table()
names(exp_N) <- as.character(exp_N[1, ])
exp_N <- exp_N[-c(1), ]
exp_N[, (names(exp_N)) := lapply(.SD, as.numeric)]

gene_pos <- fread("D:/oral_cancer/expression/expression_data/probe_pos_mulInterval.txt", header = T)
gene_pos <- gene_pos[Gene %in% names(exp_N)]
setDT(gene_pos)

gene_pos <- gene_pos[order(CHROMOSOME, start, end)]
gene_pos[, gap_prev := start - data.table::shift(end),
  by = CHROMOSOME
]

# 距離累計 > 1e6，更新 group_id
gene_pos[, group_id := make_group_id(gap_prev), by = CHROMOSOME]
gene_pos[, group_id := paste0("chr", CHROMOSOME, "_group", group_id)]
gene_pos[, gap_prev := NULL]


ld_type <- gene_pos$group_id %>% unique()
cor_group <- vector("list", length(ld_type))
for (i in seq_along(ld_type)) {
  df <- gene_pos[group_id %in% ld_type[i]]
  group_probe <- df$Gene
  cor_dt <- cor(exp_N[, ..group_probe]) %>% 
  as.table() %>% 
  as.data.table()

  setnames(cor_dt, c("probe_1", "probe_2", "cor"))

  # 留下對角線下的 cor 且 cor 夠大的
  cor_dt <- cor_dt[probe_1< probe_2 & abs(cor) > cor_cutoff]
  cor_group[[i]] <- cor_dt
}

cor_dt <- rbindlist(cor_group)
cor_dt <- gene_pos[, .(probe_1=Gene, chr=CHROMOSOME,group_id)][cor_dt, on=.(probe_1)] 

gt_N <- fread("C:/Peter/rawData_eQTL/r2_filter_0.8/outcome/raw_N_FDR_R2_0.8.txt")
eQTL_gene <- unique(gt_N$Gene)
eQTL_probe <- unique(gt_N$Probe)


# 計算每個 group，有顯著 eQTL 的數量
gene_pos[,sig_FDR := 0]
gene_pos[Gene %in% eQTL_probe, sig_FDR := 1]
group_count <- gene_pos[, .(sig_FDR_sum = sum(sig_FDR, na.rm = TRUE)),
  by = group_id
]

# 挑出 cor>cor_cutoff 的 probe
cor_pass_probe <- gene_pos[Gene %in% union(cor_dt$probe_1,cor_dt$probe_2)]
record_sig <- cor_pass_probe[,.(High_Cor_probe_number=.N),by=group_id][group_count, on=.(group_id)]
record_sig[is.na(High_Cor_probe_number), High_Cor_probe_number := 0]
fit <- lm(High_Cor_probe_number ~ sig_FDR_sum, data = record_sig)
cor_result <- cor(record_sig$sig_FDR_sum, record_sig$High_Cor_probe_number)

# 算 Spearman cor
spearman <- cor.test(
    record_sig$High_Cor_probe_number,
    record_sig$sig_FDR_sum,
    method = "spearman"
)


# poisson reg
fit_poisson <- glm(
  sig_FDR_sum ~ High_Cor_probe_number,
  family = poisson(link = "log"),
  data = record_sig
)

summary(fit_poisson)$coefficients["High_Cor_probe_number", "Pr(>|z|)"]
sum(residuals(fit_poisson, type = "pearson")^2) / fit_poisson$df.residual


cor_summary[[j]] <- data.table(
  probe_cor_cutoff = cor_cutoff,
  cor = cor_result,
  reg_pval = summary(fit)$coefficients[2, 4],
  spearman_pval = spearman$p.value,
  spearman_rho = spearman$estimate %>% as.numeric(),
  poisson_pval = summary(fit_poisson)$coefficients["High_Cor_probe_number", "Pr(>|z|)"],
  poisson_beta = summary(fit_poisson)$coefficients["High_Cor_probe_number", "Estimate"],
  dispersion = sum(residuals(fit_poisson, type = "pearson")^2) / fit_poisson$df.residual
)
cor_summary[[j]][, overdispersion := ifelse(dispersion>1,1,0)]






png(sprintf("C:/Peter/gene_enrichment/code_project/outcome/ld_bootstrap/cor%s_eQTLsig_LDgene_relation.png",
cor_cutoff),
  width = 8, height = 6, units = "in", res = 300
)

plot(record_sig$sig_FDR_sum, record_sig$High_Cor_probe_number,
xlab = "eQTL sig number", 
ylab = "High Cor gene number",
main = sprintf("eQTL sig number v.s. High Cor gene number under %s cor. threshold", cor_cutoff),
pch = 19,
col="red")
abline(a = 0, b = 1, lty = 2, lwd = 2)
abline(lm(sig_FDR_sum ~ High_Cor_probe_number, data = record_sig), col = "blue")
legend(
  "topright",
  legend = c("y = x", "Reg"),
  col = c("black", "blue"),
  lty = c(2, 1),
  lwd = c(2, 2),
  bty = "n"
)

dev.off()
}

cor_final <- rbindlist(cor_summary)
fwrite(cor_final,
  "C:/Peter/gene_enrichment/code_project/outcome/ld_bootstrap/cor_final.txt",
  row.names = F, col.names = T, sep = "\t"
)


png("C:/Peter/gene_enrichment/code_project/outcome/ld_bootstrap/corCUTOFF_Cor.png",
  width = 8, height = 6, units = "in", res = 300
)

plot(cor_final$probe_cor_cutoff, cor_final$cor,
xlab = "probe_cor_cutoff", 
ylab = "Cor between LD gene and sig number",
pch = 19,
col="red")

dev.off()



}





# 紀錄 0.7 的 

 
cor_cutoff <- 0.7
exp_N <- fread("C:/Peter/rawData_eQTL/raw_Exp_mulInterval_N.txt", header = T)
exp_N <- t(exp_N) %>% as.data.table()
names(exp_N) <- as.character(exp_N[1, ])
exp_N <- exp_N[-c(1), ]
exp_N[, (names(exp_N)) := lapply(.SD, as.numeric)]

gene_pos <- fread("D:/oral_cancer/expression/expression_data/probe_pos_mulInterval.txt", header = T)
gene_pos <- gene_pos[Gene %in% names(exp_N)]
setDT(gene_pos)

gene_pos <- gene_pos[order(CHROMOSOME, start, end)]
gene_pos[, gap_prev := start - data.table::shift(end),
  by = CHROMOSOME
]

# 距離累計 > 1e6，更新 group_id
gene_pos[, group_id := make_group_id(gap_prev), by = CHROMOSOME]
gene_pos[, group_id := paste0("chr", CHROMOSOME, "_group", group_id)]
gene_pos[, gap_prev := NULL]


ld_type <- gene_pos$group_id %>% unique()
cor_group <- vector("list", length(ld_type))
for (i in seq_along(ld_type)) {
  df <- gene_pos[group_id %in% ld_type[i]]
  group_probe <- df$Gene
  cor_dt <- cor(exp_N[, ..group_probe]) %>% 
  as.table() %>% 
  as.data.table()

  setnames(cor_dt, c("probe_1", "probe_2", "cor"))

  # 留下對角線下的 cor 且 cor 夠大的
  cor_dt <- cor_dt[probe_1< probe_2 & abs(cor) > cor_cutoff]
  cor_group[[i]] <- cor_dt
}

cor_dt <- rbindlist(cor_group)
cor_dt <- gene_pos[, .(probe_1=Gene, chr=CHROMOSOME,group_id)][cor_dt, on=.(probe_1)] 

gt_N <- fread("C:/Peter/rawData_eQTL/r2_filter_0.8/outcome/raw_N_FDR_R2_0.8.txt")
eQTL_gene <- unique(gt_N$Gene)
eQTL_probe <- unique(gt_N$Probe)


# 計算每個 group，有顯著 eQTL 的數量 ----
gene_pos[,sig_FDR := 0]
gene_pos[Gene %in% eQTL_probe, sig_FDR := 1]
group_count <- gene_pos[, .(sig_FDR_sum = sum(sig_FDR, na.rm = TRUE)),
  by = group_id
]

# 挑出 cor>cor_cutoff 的 probe
cor_pass_probe <- gene_pos[Gene %in% union(cor_dt$probe_1,cor_dt$probe_2)]
record_sig <- cor_pass_probe[,.(High_Cor_probe_number=.N),by=group_id][group_count, on=.(group_id)]
record_sig[is.na(High_Cor_probe_number), High_Cor_probe_number := 0]



  fwrite(gene_pos,
    "C:/Peter/gene_enrichment/code_project/outcome/ld_bootstrap/probe_group.txt",
    row.names = F, col.names = T, sep = "\t"
  )

  fwrite(record_sig,
    "C:/Peter/gene_enrichment/code_project/outcome/ld_bootstrap/0.7number_record.txt",
    row.names = F, col.names = T, sep = "\t"
  )


