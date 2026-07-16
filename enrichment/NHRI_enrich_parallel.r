## package ----
library(ggplot2)
library(data.table)
library(dplyr)
library(stringr)
library(readxl) # read_xlsx ft


start_time <- Sys.time()


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
cat("lung gene: ", uniqueN(lung_twas$GeneSymbol), "\n")

lung_twas <- lung_twas[!is.na(top1_p), ]
cat("lung gene with pval (not NA): ", uniqueN(lung_twas$GeneSymbol), "\n")

lung_sig <- lung_twas[top1_p < 0.01]
cat("lung sig gene: ", uniqueN(lung_sig$GeneSymbol), "\n")

eQTL_gene <- unique(gt_N$Gene)
cat("eQTL sig gene: ", uniqueN(eQTL_gene), "\n")


# 挑出我們資料有紀錄、lung 紀錄的 probe/gene
common_gene <- intersect(all_eQTL_gene, lung_twas$GeneSymbol)

cat("intersect of eQTL, lung total gene:", intersect(all_eQTL_gene, lung_twas$GeneSymbol) %>% uniqueN(), "\n")

# 在一邊顯著，同時有被紀錄在另一邊 數量
cat("eQTL significant genes recorded in lung TWAS:", uniqueN(gt_N$Gene[gt_N$Gene %in% lung_twas$GeneSymbol]), "\n")
cat("lung significant genes recorded in eQTL:", lung_sig$GeneSymbol[lung_sig$GeneSymbol %in% all_eQTL_gene] %>% uniqueN(), "\n")


# 交集數量
observed_intersect <- intersect(eQTL_gene, lung_twas[top1_p < 0.01, GeneSymbol]) %>% uniqueN()
cat("intersect of eQTL, lung sig gene:", observed_intersect, "\n")

sample_gene_number <- uniqueN(gt_N$Gene[gt_N$Gene %in% lung_twas$GeneSymbol])


# 從共同基因隨機挑，看多少 在 lung_twas 顯著 ----
library(parallel)
tissue_name <- "NHRI_lung"
NHRI_sig_sub <- lung_sig$GeneSymbol[lung_sig$GeneSymbol %in% all_eQTL_gene]
is_sig <- common_gene %in% NHRI_sig_sub
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
    # is_sig 判斷 common_gene 是否在 NHRI_sig_sub，一串的 T,F，sum 會把T當1, F當0
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
png(file.path("C:/Peter/gene_enrichment/code_project/outcome/rawData_eQTL_N/for_gene/lungTWAS", sprintf("random_%sgene.png", sample_gene_number)),
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



end_time <- Sys.time()
print(end_time - start_time)


p_value <- phyper(
  O - 1,
  m = K,
  n = N - K,
  k = n,
  lower.tail = FALSE
)


# cumulative probability Hypergeometric Dist. ----
# N <- length(is_sig)
# K <- sum(is_sig)
# n <- sample_gene_number
# O <- observed_intersect

# # lower.tail = F 表示右尾，P(X> O)=P(X>= O-1)
# p_value <- phyper(
#   O - 1,
#   m = K,
#   n = N - K,
#   k = n,
#   lower.tail = FALSE
# )


# lung gene:  18103
# lung gene with pval (not NA):  17157
# lung sig gene:  2875
# eQTL sig gene:  177
# intersect of eQTL, lung total gene: 13708
# eQTL significant genes recorded in lung TWAS: 154
# lung significant genes recorded in eQTL: 2357
# intersect of eQTL, lung sig gene: 111
# sample 154 genes, mean sig in NHRI_lung genes number: 26.47938
# sample 154 genes, sig in NHRI_lung genes number range: 2 58
# sample 154 genes, sig in NHRI_lung genes number >= 111 times: 0
# empirical p-value: 1e-09