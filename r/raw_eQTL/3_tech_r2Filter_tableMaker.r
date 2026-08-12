# FIN,EUR 造表格 ----
# 讀取檔案 C:/Peter/rawData_eQTL/r2_filter_%s/outcome/raw_MixFinngenPval_7_%s.txt 算 3 種手法剩餘 snp 數量
# 讀取檔案 C:/Peter/repeatSNP_clumping_raw/r2_filter_%s/eQTL_%s/%s/outcome/1000kb_LD0.2/%s 算 3 種手法clumping 後剩餘 snp 數量

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

library(data.table)
library(openxlsx)
library(magrittr)







for (race in c("EUR","FIN")) {
for (eqtl_threshold in c("bon","FDR")) {

read_clump_n <- function(r2_threshold, file) {
  path <- sprintf(
    "C:/Peter/repeatSNP_clumping_raw/r2_filter_%s/eQTL_%s/%s/outcome/1000kb_LD0.2/%s",
    r2_threshold, eqtl_threshold,  race, file
  )

  if (!file.exists(path)) {
    return(NA_integer_)
  }

  nrow(fread(path))
}

# mostLD_snp_nearest_hg38 不是空 or NA表示有附近的 LD snp，就用它取代 gt_N_hg38
clean_snp <- function(dt) {
  res <- dt[, {
    mostLD <- as.character(mostLD_snp_nearest_hg38)
    gt_N <- as.character(gt_N_hg38)
    fifelse(is.na(mostLD) | mostLD == "", gt_N, mostLD)
  }]
  unique(res[!is.na(res) & res != ""])
}

format_p <- function(x) {
  formatC(x, format = "e", digits = 2)
}

wb <- createWorkbook()

for (r2_threshold in c("no","0.6", "0.7", "0.8", "0.9")) {

  if (eqtl_threshold == "FDR") {
     eur <- fread(sprintf(
    "C:/Peter/rawData_eQTL/r2_filter_%s/outcome/raw_MixFinngenPval_7_%s.txt",
    r2_threshold, race
  ))
  } else if (eqtl_threshold == "bon") {
     eur <- fread(sprintf(
    "C:/Peter/rawData_eQTL/r2_filter_%s/outcome/raw_MixFinngenPval_7_%s_bon.txt",
    r2_threshold, race
  ))
  }


  eur_mostpval <- eur[, .SD[1], by = MOST_snp_nearest_hg38]
  eur_mostpval <- eur_mostpval[MOST_snp_nearest_hg38 != ""]

  eur_mostLD <- eur[, .SD[1], by = mostLD_snp_nearest_hg38]
  eur_mostLD <- eur_mostLD[mostLD_snp_nearest_hg38 != ""]

  ld_snp <- eur_mostLD$mostLD_snp_nearest_hg38[
    eur_mostLD$mostLD_snp_nearest_hg38 != ""
  ]

  eqtl <- eur[ref != ""]

  n_eqtl        <- nrow(eur)
  n_in_finngen  <- eur[alt != "", .N]
  n_not_finngen <- eur[ref == "", .N]
  n_ld_replace  <- eur[!is.na(LD), .N]
  n_no_ld_pval  <- eur[is.na(mostLD_or_finngen_FDR), .N]

  n_ld_total <- uniqueN(eur_mostLD$mostLD_snp_nearest_hg38)

  n_ld_not_repeat <- uniqueN(
    ld_snp[!ld_snp %in% eqtl$gt_N_hg38[eqtl$gt_N_hg38 != ""]]
  )

  n_after_ld <- n_in_finngen + n_ld_not_repeat
  n_mostpval <- nrow(eur_mostpval)

  ld_fdr  <- clean_snp(eur[mostLD_or_finngen_FDR < 0.05])
  ld_qval <- clean_snp(eur[mostLD_or_finngen_qvalue < 0.05])
  ld_bon  <- clean_snp(eur[mostLD_or_finngen_Bonfi == 1])

  min_p_1 <- min(eur[!is.na(gt_N_hg38_finngen_pval), gt_N_hg38_finngen_pval])
  min_p_2 <- min(eur[!is.na(mostLD_or_finngen_pval), mostLD_or_finngen_pval])

  n_1_fdr_clump  <- read_clump_n(r2_threshold, "1_fdr.clumped")
  n_1_qval_clump <- read_clump_n(r2_threshold, "1_qval.clumped")
  n_2_fdr_clump  <- read_clump_n(r2_threshold, "2_fdr.clumped")
  n_2_qval_clump <- read_clump_n(r2_threshold, "2_qval.clumped")
  n_3_fdr_clump  <- read_clump_n(r2_threshold, "3_fdr.clumped")
  n_3_bon_clump  <- read_clump_n(r2_threshold, "3_bon.clumped")

  ld_not_repeat_col <- sprintf(
    "LD SNP number\n(doesn't repeat\nin %s snp)",
    n_in_finngen
  )
  people_n <- ifelse(race == "FIN", 99, 503)

  out <- data.table(
    How = c(
      sprintf(
        "Among %s SNP ,  %s SNP recorded in finngen",
        n_eqtl, n_in_finngen
      ),
      sprintf(
        "%s-%s=%s SNP not in finngen. Among %s, %s SNP replace by most LD SNP through 1000G %s  %s people        ( %s SNP haven't pvalue)",
        n_eqtl, n_in_finngen, n_not_finngen,
        n_not_finngen, n_ld_replace, race, people_n, n_no_ld_pval
      ),
      sprintf(
        "Choose most sig. SNP surrounding %s SNP. If more than one, choose the close one. Retain %s SNP",
        n_eqtl, n_mostpval
      )
    ),

    `Number of SNP` = c(
      n_in_finngen,
      n_after_ld,
      n_mostpval
    ),

    `Number of snp with FDR(BH) <0.05 // after 1000kb r^2 0.2 clumping` = c(
      sprintf(
        "%s // %s",
        uniqueN(eur[gt_N_hg38_finngen_FDR < 0.05, gt_N_hg38]),
        n_1_fdr_clump
      ),
      sprintf(
        "%s // %s",
        length(ld_fdr),
        n_2_fdr_clump
      ),
      sprintf(
        "%s //  %s",
        uniqueN(eur_mostpval[
          MOST_snp_hg38_finngen_FDR < 0.05,
          MOST_snp_nearest_hg38
        ]),
        n_3_fdr_clump
      )
    ),

    `Number of snp with qvalue <0.05 // after 1000kb r^2 0.2 clumping` = c(
      sprintf(
        "%s // %s",
        uniqueN(eur[gt_N_hg38_finngen_qvalue < 0.05, gt_N_hg38]),
        n_1_qval_clump
      ),
      sprintf(
        "%s // %s",
        length(ld_qval),
        n_2_qval_clump
      ),
      "NA"
    ),

    `Number of snp with Bonferroni correction // after 1000kb r^2 0.2 clumping` = c(
      sprintf(
        "%s    (most sig. is %s )",
        uniqueN(eur[gt_N_hg38_finngen_Bonfi == 1, gt_N_hg38]),
        format_p(min_p_1)
      ),
      sprintf(
        "%s    (most sig. is %s )",
        length(ld_bon),
        format_p(min_p_2)
      ),
      sprintf(
        "%s //  %s",
        uniqueN(eur_mostpval[
          MOST_snp_hg38_finngen_Bonfi == 1,
          MOST_snp_nearest_hg38
        ]),
        n_3_bon_clump
      )
    ),

    `record in Finngen SNP number (repeat SNP choose the smallest pvalue one)` = c(
      n_in_finngen,
      n_in_finngen,
      n_in_finngen
    ),

    `LD SNP number` = c(
      "",
      n_ld_total,
      ""
    ),

  LD_not_repeat_tmp = c(
    "",
    n_ld_not_repeat,
    ""
  ),

    m = c(
      n_in_finngen,
      n_after_ld,
      n_mostpval
    ),

    `bon(0.05/ m)` = c(
      format_p(0.05 / n_in_finngen),
      format_p(0.05 / n_after_ld),
      format_p(0.05 / n_mostpval)
    )
  )

  sheet_name <- paste0("r2_", r2_threshold)
  addWorksheet(wb, sheet_name)
  setnames(out, "LD_not_repeat_tmp", ld_not_repeat_col)
  writeData(wb, sheet_name, out)

  header_style <- createStyle(
    fgFill = "#F4B183",
    halign = "center",
    valign = "center",
    textDecoration = "bold",
    border = "TopBottomLeftRight",
    wrapText = TRUE
  )

  body_style <- createStyle(
    halign = "center",
    valign = "center",
    border = "TopBottomLeftRight",
    wrapText = TRUE
  )

  blue_style <- createStyle(
    fgFill = "#BDD7EE",
    halign = "center",
    valign = "center",
    border = "TopBottomLeftRight",
    wrapText = TRUE
  )

  green_style <- createStyle(
    fgFill = "#C6E0B4",
    halign = "center",
    valign = "center",
    border = "TopBottomLeftRight",
    wrapText = TRUE
  )

  addStyle(wb, sheet_name, header_style, rows = 1, cols = 1:ncol(out), gridExpand = TRUE)
  addStyle(wb, sheet_name, body_style, rows = 2:4, cols = 1:ncol(out), gridExpand = TRUE)

  addStyle(wb, sheet_name, blue_style, rows = 2:3, cols = 1:2, gridExpand = TRUE)
  addStyle(wb, sheet_name, green_style, rows = 4, cols = 1:2, gridExpand = TRUE)

  setColWidths(wb, sheet_name, cols = 1, widths = 38)
  setColWidths(wb, sheet_name, cols = 2, widths = 18)
  setColWidths(wb, sheet_name, cols = 3:5, widths = 24)
  setColWidths(wb, sheet_name, cols = 6, widths = 25)
  setColWidths(wb, sheet_name, cols = 7:10, widths = 16)

  setRowHeights(wb, sheet_name, rows = 1, heights = 90)
  setRowHeights(wb, sheet_name, rows = 2, heights = 90)
  setRowHeights(wb, sheet_name, rows = 3, heights = 120)
  setRowHeights(wb, sheet_name, rows = 4, heights = 110)

  freezePane(wb, sheet_name, firstRow = TRUE)
}

saveWorkbook(
    wb,
    sprintf("C:/Peter/rawData_eQTL/outcome/raw_%s_%s_r2_filter_summary.xlsx", eqtl_threshold, race),
    overwrite = TRUE
  )



}
}


# 造表格的數字由來原理
# 3 種手法剩餘 snp 數量 ----
# for (r2_threshold in c("0.8","0.9"  )) {

#   eur <- sprintf("C:/Peter/rawData_eQTL/r2_filter_%s/outcome/raw_MixFinngenPval_7_EUR.txt",r2_threshold) %>% fread()
# eur_mostpval <- eur[,.SD[1],by=MOST_snp_nearest_hg38]
# eur_mostpval <- eur_mostpval[MOST_snp_nearest_hg38!="",]

# eur_mostLD <- eur[,.SD[1],by=mostLD_snp_nearest_hg38]
# eur_mostLD <- eur_mostLD[mostLD_snp_nearest_hg38!="",]
# ld_snp <- eur_mostLD$mostLD_snp_nearest_hg38[eur_mostLD$mostLD_snp_nearest_hg38!=""]
# eqtl <- eur[ref!="",]

# # fifelse 是高效版 ifelse
# clean_snp <- function(dt) {
#   res <- dt[, {
#     mostLD <- as.character(mostLD_snp_nearest_hg38)
#     gt_N <- as.character(gt_N_hg38)
#     fifelse(is.na(mostLD) | mostLD == "", gt_N, mostLD)
#   }]
#   unique(res[!is.na(res) & res != ""])
# }

# # 2. 直接針對不同條件進行篩選並處理
# ld_fdr  <- clean_snp(eur[mostLD_or_finngen_FDR < 0.05])
# ld_qval <- clean_snp(eur[mostLD_or_finngen_qvalue < 0.05])
# ld_bon  <- clean_snp(eur[mostLD_or_finngen_Bonfi == 1])


# cat("threshold: ", r2_threshold, "\n")
# cat("eQTL snp 數量: ",
#     nrow(eur),
#     "\n")
# cat("eQTL snp 紀錄在 fin 數量: ",
#     eur$alt[eur$alt!=""] %>% length(),
#     "\n")
# cat("過 fdr snp 數量: ",
#     (eur[gt_N_hg38_finngen_FDR<0.05, gt_N_hg38]) %>% uniqueN(),
#     "\n")
# cat("過 qval snp 數量: ",
#     (eur[gt_N_hg38_finngen_qvalue<0.05, gt_N_hg38]) %>% uniqueN(),
#     "\n")
# cat("過 bon snp 數量: ",
#     (eur[gt_N_hg38_finngen_Bonfi==1, gt_N_hg38]) %>% uniqueN(),
#     "\n")
# cat("\n")

# # most LD snp
# cat("不在finngen 的 snp 數量: ",
#     length(which((eur$ref==""))),
#     "\n")
# cat("多少 snp 用 LD snp 補數量: ",eur$LD[!is.na(eur$LD)] %>% length(),
#     "\n")
# cat("不在finngen，且附近1000 kb 沒有 finngen snp (不能用 LD snp 替代)的 snp 數量: ",length(which(is.na(eur$mostLD_or_finngen_FDR))),
#     "\n")
# cat("LD snp number: ",eur_mostLD %>% uniqueN(),
#     "\n")
# cat("LD snp number(跟 eQTL snp 不重複): ",ld_snp[!ld_snp %in% eqtl$gt_N_hg38[eqtl$gt_N_hg38!=""]] %>%
#       uniqueN(),
#     "\n")


# cat("過 fdr snp 數量: ",
#     length(ld_fdr),
#     "\n")
# cat("過 qval snp 數量: ",
#     length(ld_qval),
#     "\n")
# cat("過 bon snp 數量: ",
#     length(ld_bon),
#     "\n")
# cat("\n")


# # most pval snp
# cat("不重複 most snp 數量: ",nrow(eur_mostpval),
#     "\n")
# cat("過 fdr snp 數量: ",
#     (eur_mostpval[MOST_snp_hg38_finngen_FDR<0.05, MOST_snp_nearest_hg38]) %>% uniqueN(),
#     "\n")
# cat("過 bon snp 數量: ",
#     (eur_mostpval[MOST_snp_hg38_finngen_Bonfi==1, MOST_snp_nearest_hg38]) %>% uniqueN(),
#     "\n")


# cat(rep("\n",5))
# }


# clumping 後剩餘 snp 數量 ----
# 3 種手法clumping 後剩餘 snp 數量
# files <- c("1_fdr.clumped", "1_qval.clumped","2_fdr.clumped", "2_qval.clumped", "3_fdr.clumped", "3_bon.clumped")
# clump_record <- data.table()

# for (eqtl_threshold in c("bon","FDR")) {

# for (r2_threshold in c("0.9")) {
#   cat("r2 ",r2_threshold,"\n")
#   for (file in files) {

#   af_clump <- sprintf("C:/Peter/repeatSNP_clumping_raw/r2_filter_%s/eQTL_%s/EUR/outcome/1000kb_LD0.2/%s",
#                       r2_threshold, eqtl_threshold, file) %>%
#     fread()

#    cat(file, " ", r2_threshold, " ", nrow(af_clump),"\n")

#   }
#   cat("\n")
# }
   


# }


# clumping 後剩餘 snp 名稱 ----
# files <- c("1_fdr.clumped", "1_qval.clumped","2_fdr.clumped", "2_qval.clumped")
# populations <- c("EUR", "FIN")
# clump_record <- list()

# entry <- 1


# for (pop in populations) {
#   for (r2_threshold in c("no", "0.6", "0.7", "0.9")) {
#     for (file in files) {
#       af_clump <- sprintf(
#         "C:/Peter/repeatSNP_clumping_raw/r2_filter_%s/%s/outcome/1000kb_LD0.2/%s",
#         r2_threshold, pop, file
#       ) %>% fread()


#       clump_record[[entry]] <- data.table(
#         race = pop,
#         r2_threshold = r2_threshold,
#         file = file,
#         number = length(af_clump$SNP),
#         SNP = paste(af_clump$SNP, collapse = ",")
#       )
#       entry <- entry + 1
#     }
#   }
# }

# clump_record <- rbindlist(clump_record)
# if (!requireNamespace("writexl", quietly = TRUE)) {
#   install.packages("writexl")
# }

# writexl::write_xlsx(
#   clump_record,
#   "C:/Peter/rawData_eQTL/outcome/1000kb_LD0.2_clumping_record.xlsx"
# )


# clumping 後剩餘 snp 最小 pval ----
# for (r2_threshold in c("0.6" ,"0.7","0.9", "0.8")) {

#   eur <- sprintf("C:/Peter/rawData_eQTL/r2_filter_%s/outcome/raw_MixFinngenPval_7_FIN.txt",r2_threshold) %>%
#     fread()
# # record
# cat("threshold: ", r2_threshold, "\n")

# cat("第一種手法 most sig pval: ",
#     (eur[!is.na(gt_N_hg38_finngen_pval),gt_N_hg38_finngen_pval]) %>% min() %>% format(digits = 3,scientific = T),
#     "\n")
# cat("\n")

# # most LD snp
# cat("第二種手法 most sig pval: ",
#     (eur[!is.na(mostLD_or_finngen_pval), mostLD_or_finngen_pval]) %>% min() %>% format(digits = 3,scientific = T),
#     "\n")
# cat("\n")


# cat(rep("\n",3))
# }




### Record
# 做表格用，紀錄顯著數量

# eQTL result ----
# file_name <- c("C:/Peter/rawData_eQTL/trash/raw_maf_gt_N_pvalue.txt",
#                "C:/Peter/rawData_eQTL/trash/raw_maf_gt_T_pvalue.txt"
#                )


# for (i in file_name) {
#   a <- fread(i)
#   a[, gene2 := gsub("^([^_]+_[^_]+).*", "\\1",gene)]
#   a[,gene := NULL]

#   sig_fdr_snp <- unique(a, by="SNP") %>% nrow()
#   sig_fdr_asso <- a %>% nrow()
#   sig_fdr_probe <- unique(a, by="gene2") %>% nrow()


#   cat(i,"\n")

#   cat("eQTL :","\n")
#   cat("unique Probe number: ")
#   print(sig_fdr_probe)
#   cat("unique SNPs number: ")
#   print(sig_fdr_snp)
#   cat("association number: ")
#   print(sig_fdr_asso)

#   cat(rep("\n", 3))
# }

# 輸出結果:
# C:/Peter/rawData_eQTL/trash/raw_maf_gt_N_pvalue.txt
# eQTL :
# unique Probe number: [1] 20913
# unique SNPs number: [1] 5550566
# association number: [1] 88867137


# C:/Peter/rawData_eQTL/trash/raw_maf_gt_T_pvalue.txt
# eQTL :
# unique Probe number: [1] 20913
# unique SNPs number: [1] 5550566
# association number: [1] 88867137


















