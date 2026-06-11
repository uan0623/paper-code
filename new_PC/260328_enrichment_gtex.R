# 目的: 算 GTEx data 的 enrichment，把 rawData_eQTL 改成 OQN_FIXpeople_before_eQTL, OQN_FIXgene_before_eQTL 即可
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
  "D:/Peter/gene_enrichment/code_project/outcome/rawData_eQTL/%s",
  c("lungTWAS", "HNSC", "TCGA-LUAD", "TCGA-LUSC", "GTEx-salivary", "GTEx-esophagus", "GTEx-thyroid", "GTEx-lung")
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






# run_gtex_eqtl_enrichment ----
run_gtex_eqtl_enrichment <- function(
  tissue_name,
  weight_path,
  output_dir
) {
  cat(
    "After remove R2<0.8 snp, run eQTL genes:",
    uniqueN(all_eQTL_gene),
    "\n"
  )

  gtex_eqtl <- fread(
    weight_path,
    select = c("gene_name", "qval", "pval_nominal")
  )
  cat(tissue_name, "genes:", uniqueN(gtex_eqtl$gene_name), "\n")

  gtex_eqtl <- gtex_eqtl[!is.na(qval), ]
  cat(
    tissue_name, "After remove NA qval, genes:",
    uniqueN(gtex_eqtl$gene_name),
    "\n"
  )

  gt_N <- fread("D:/Peter/rawData_eQTL/r2_filter_0.8/outcome/raw_N_FDR_R2_0.8.txt")
  random_times <- 10000
  eQTL_gene <- unique(gt_N$Gene)

  gtex_sigGene <- gtex_eqtl[qval < 0.05, gene_name] %>%
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

  common_gene <- intersect(all_eQTL_gene, gtex_eqtl$gene_name)
  sample_gene_number <- uniqueN(gt_N$Gene[gt_N$Gene %in% gtex_eqtl$gene_name])
  observed_intersect <- intersect(eQTL_gene, gtex_sigGene) %>% uniqueN()

  cat(
    tissue_name, "significant genes recorded in eQTL:",
    gtex_sigGene[gtex_sigGene %in% all_eQTL_gene] %>% uniqueN(), "\n"
  )
  cat(
    "eQTL significant genes recorded in", tissue_name, ":",
    uniqueN(gt_N$Gene[gt_N$Gene %in% gtex_eqtl$gene_name]), "\n"
  )
  cat(
    "intersect of eQTL,", tissue_name, "total gene:",
    intersect(all_eQTL_gene, gtex_eqtl$gene_name) %>% uniqueN(), "\n"
  )
  cat(
    "intersect of eQTL,", tissue_name, "sig gene:",
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
    "sample", sample_gene_number, "genes, mean sig in", tissue_name, "gene number:",
    mean(repeat_gene_number), "\n"
  )
  cat(
    "sample", sample_gene_number, "genes, sig in", tissue_name, "gene number range:",
    range(repeat_gene_number), "\n"
  )
  cat(
    "sample", sample_gene_number, "genes, sig in", tissue_name, "gene number >=",
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

  setorder(gt_N, `p-value`)
  a <- gtex_eqtl[gene_name %in% eQTL_gene, ][, .(gene_name, pval = pval_nominal)]
  setkey(a, pval)
  a <- a[, .SD[1], by = gene_name]
  a[, type := "eQTL"]

  random_list <- list()
  random_times <- 5
  for (i in 1:random_times) {
    set.seed(i)
    b_tmp <- gtex_eqtl[gene_name %in% sample(common_gene, sample_gene_number), ][, .(gene_name, pval = pval_nominal)]
    b_tmp[, type := paste0("random_", i)]
    random_list[[i]] <- b_tmp
  }

  all_random_b <- rbindlist(random_list)
  combine_all <- rbind(a[, .(gene_name, pval, type)], all_random_b)


  all_pval <- combine_all[, pval]
  x_auto <- quantile(all_pval, probs = 0.8, na.rm = TRUE)
  x_auto <- max(x_auto, 1e-4)

  for (x_cutoff in c(1, 0.1, 0.01, 0.001, x_auto)) {
    png(file.path(output_dir, sprintf("%sgene_pval_nominalDIST_%s.png", sample_gene_number, x_cutoff)),
      width = 8, height = 6, units = "in", res = 300
    )
    print(
      density_plot(
        combine_all, "pval",
        sprintf("Density Plot for %s Gene pval_nominal vs 1e4 Random", sample_gene_number),
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

  for (x_cutoff in c(1, 0.01, 0.001, x_auto)) {
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
          title = sprintf("Mean Density Plot for %s Gene pval_nominal vs 1e4 Random", sample_gene_number),
          x = "Q-value",
          y = "Density",
          color = "Group"
        ) +
        theme_minimal()
    )
    dev.off()
  }
}


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



# 使用 eQTL data ----
run_gtex_eqtl_enrichment(
  tissue_name = "GTEx-salivary",
  weight_path = "D:/Peter/GTEx_calculator/data/GTEx_Analysis_v11_eQTL/Minor_Salivary_Gland.v11.eGenes.txt",
  output_dir = plot_output_dir[5]
)

run_gtex_eqtl_enrichment(
  tissue_name = "GTEx-esophagus",
  weight_path = "D:/Peter/GTEx_calculator/data/GTEx_Analysis_v11_eQTL/Esophagus_Mucosa.v11.eGenes.txt",
  output_dir = plot_output_dir[6]
)

run_gtex_eqtl_enrichment(
  tissue_name = "GTEx-thyroid",
  weight_path = "D:/Peter/GTEx_calculator/data/GTEx_Analysis_v11_eQTL/Thyroid.v11.eGenes.txt",
  output_dir = plot_output_dir[7]
)

run_gtex_eqtl_enrichment(
  tissue_name = "GTEx-lung",
  weight_path = "D:/Peter/GTEx_calculator/data/GTEx_Analysis_v11_eQTL/Lung.v11.eGenes.txt",
  output_dir = plot_output_dir[8]
)


run_gtex_eqtl_enrichment(
  tissue_name = "GTEx-salivary",
  weight_path = "D:/Peter/gene_enrichment/code_project/data/GTEx_salivary/fusion_project/fusion_weights_summary.csv",
  output_dir = plot_output_dir[5]
)

run_gtex_eqtl_enrichment(
  tissue_name = "GTEx-esophagus",
  weight_path = "D:/Peter/gene_enrichment/code_project/data/GTEx_esophagus/fusion_project/fusion_weights_summary.csv",
  output_dir = plot_output_dir[6]
)

run_gtex_eqtl_enrichment(
  tissue_name = "GTEx-thyroid",
  weight_path = "D:/Peter/gene_enrichment/code_project/data/GTEx_thyroid/fusion_project/fusion_weights_summary.csv",
  output_dir = plot_output_dir[7]
)

run_gtex_eqtl_enrichment(
  tissue_name = "GTEx-lung",
  weight_path = "D:/Peter/gene_enrichment/code_project/data/GTEx_lung/fusion_project/fusion_weights_summary.csv",
  output_dir = plot_output_dir[8]
)






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
# sample 135 genes, sig in GTEx_Adipose_Subcutaneous genes number range: 46, 89
# sample 135 genes, sig in GTEx_Adipose_Subcutaneous genes number >= 95 times: 0

# GTEx_Brain_Cerebellum raw data genes and ENSG genes (remove NA): 23892 23922
# GTEx_Brain_Cerebellum After remove NA Top1_Pval, genes and ENSg gene: 23892 23922
# After choose most nearest gene for repeat gene, genes and ENSG genes: 23892 23892
# GTEx_Brain_Cerebellum sig. genes: 8302
# eQTL sig. genes: 177
# intersect of eQTL, GTEx_Brain_Cerebellum total genes: 11911
# GTEx_Brain_Cerebellum significant genes recorded in eQTL: 4354
# eQTL significant genes recorded in GTEx_Brain_Cerebellum : 133
# intersect of eQTL, GTEx_Brain_Cerebellum sig genes: 74
# sample 133 genes, mean sig in GTEx_Brain_Cerebellum genes number: 48.6226
# sample 133 genes, sig in GTEx_Brain_Cerebellum genes number range: 28, 73
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
# sample 136 genes, sig in GTEx_Esophagus_Gastroesophageal genes number range: 28, 74
# sample 136 genes, sig in GTEx_Esophagus_Gastroesophageal genes number >= 79 times: 0

# GTEx_Esophagus_Muscularis raw data genes and ENSG genes (remove NA): 22481 22503
# GTEx_Esophagus_Muscularis After remove NA Top1_Pval, genes and ENSg gene: 22481 22503
# After choose most nearest gene for repeat gene, genes and ENSG genes: 22481 22481
# GTEx_Esophagus_Muscularis sig. genes: 10262
# eQTL sig. genes: 177
# intersect of eQTL, GTEx_Esophagus_Muscularis total genes: 11919
# GTEx_Esophagus_Muscularis significant genes recorded in eQTL: 5648
# eQTL significant genes recorded in GTEx_Esophagus_Muscularis : 136
# intersect of eQTL, GTEx_Esophagus_Muscularis sig genes: 87
# sample 136 genes, mean sig in GTEx_Esophagus_Muscularis genes number: 64.4524
# sample 136 genes, sig in GTEx_Esophagus_Muscularis genes number range: 43, 84
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
# sample 131 genes, sig in GTEx_Heart_Left genes number range: 30, 71
# sample 131 genes, sig in GTEx_Heart_Left genes number >= 75 times: 0

# GTEx_Muscle_Skeletal raw data genes and ENSG genes (remove NA): 19963 19978
# GTEx_Muscle_Skeletal After remove NA Top1_Pval, genes and ENSg gene: 19963 19978
# After choose most nearest gene for repeat gene, genes and ENSG genes: 19963 19963
# GTEx_Muscle_Skeletal sig. genes: 9672
# eQTL sig. genes: 177
# intersect of eQTL, GTEx_Muscle_Skeletal total genes: 11429
# GTEx_Muscle_Skeletal significant genes recorded in eQTL: 5724
# eQTL significant genes recorded in GTEx_Muscle_Skeletal : 133
# intersect of eQTL, GTEx_Muscle_Skeletal sig genes: 89
# sample 133 genes, mean sig in GTEx_Muscle_Skeletal genes number: 66.6429
# sample 133 genes, sig in GTEx_Muscle_Skeletal genes number range: 41, 88
# sample 133 genes, sig in GTEx_Muscle_Skeletal genes number >= 89 times: 0
