# FIN,EUR 造表格 ----
# 讀取檔案 C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_%s/outcome/OQN_MixFinngenPval_7_%s.txt 算 3 種手法剩餘 snp 數量
# 讀取檔案 C:/Peter/repeatSNP_clumping_FIXpeople/r2_filter_%s/eQTL_%s/%s/outcome/1000kb_LD0.2/%s 算 3 種手法clumping 後剩餘 snp 數量

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
    "C:/Peter/repeatSNP_clumping_FIXpeople/r2_filter_%s/eQTL_%s/%s/outcome/1000kb_LD0.2/%s",
    r2_threshold, eqtl_threshold,  race, file
  )

  if (!file.exists(path)) {
    return(NA_integer_)
  }

  nrow(fread(path))
}

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
    "C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_%s/outcome/OQN_MixFinngenPval_7_%s.txt",
    r2_threshold, race
  ))
  } else if (eqtl_threshold == "bon") {
     eur <- fread(sprintf(
    "C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_%s/outcome/OQN_MixFinngenPval_7_%s_bon.txt",
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
    sprintf("C:/Peter/OQN_FIXpeople_before_eQTL/outcome/OQN_%s_%s_r2_filter_summary.xlsx", eqtl_threshold, race),
    overwrite = TRUE
  )



}
}


# 造表格的數字由來原理
# 3 種手法剩餘 snp 數量 ----
# for (r2_threshold in c("0.8","0.9"  )) {

#   eur <- sprintf("C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_%s/outcome/OQN_MixFinngenPval_7_EUR.txt",r2_threshold) %>% fread()
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

#   af_clump <- sprintf("C:/Peter/repeatSNP_clumping_FIXpeople/r2_filter_%s/eQTL_%s/EUR/outcome/1000kb_LD0.2/%s",
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
#         "C:/Peter/repeatSNP_clumping_FIXpeople/r2_filter_%s/%s/outcome/1000kb_LD0.2/%s",
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
#   "C:/Peter/OQN_FIXpeople_before_eQTL/outcome/1000kb_LD0.2_clumping_record.xlsx"
# )


# clumping 後剩餘 snp 最小 pval ----
# for (r2_threshold in c("0.6" ,"0.7","0.9", "0.8")) {

#   eur <- sprintf("C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_%s/outcome/OQN_MixFinngenPval_7_FIN.txt",r2_threshold) %>%
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
# file_name <- c("C:/Peter/OQN_FIXpeople_before_eQTL/trash/OQN_maf_gt_N_pvalue.txt",
#                "C:/Peter/OQN_FIXpeople_before_eQTL/trash/OQN_maf_gt_T_pvalue.txt"
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
# C:/Peter/OQN_FIXpeople_before_eQTL/trash/OQN_maf_gt_N_pvalue.txt
# eQTL :
# unique Probe number: [1] 20913
# unique SNPs number: [1] 5550566
# association number: [1] 88867137


# C:/Peter/OQN_FIXpeople_before_eQTL/trash/OQN_maf_gt_T_pvalue.txt
# eQTL :
# unique Probe number: [1] 20913
# unique SNPs number: [1] 5550566
# association number: [1] 88867137


# 顯著的probe, 有eQTL的比例 ----
# N, T EXPRESSION 顯著不同，eQTL FDR<0.05 的probe 數目

for (eqtl_threshold in c("bon","FDR")) {

message("==========================================")
message(sprintf("eQTL threshold %s", eqtl_threshold))
message("==========================================")

a <- fread("C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_0.8/outcome/OQN_exp_different_r2_0.8.txt",header=T)

# expression FDR 顯著不同，且至少一個 associatioin eQTL FDR<0.05 的probe 數目
b <- a %>% filter(sig_OQN == 1)
which(b$`0.8_N_sigCis-SNP_number` != 0) %>% length() %>% print()
which(b$`0.8_T_sigCis-SNP_number` != 0) %>% length() %>% print()
cat("\n")



# expression bonfi 顯著不同，且至少一個 associatioin eQTL bon 過門檻 的probe 數目
b <- a %>% filter(sig_OQN_Bonfi==1)
which(b$`0.8_N_sigCis-SNP_number (bon)` != 0) %>%
  length() %>%
  print()
which(b$`0.8_T_sigCis-SNP_number (bon)` != 0) %>%
  length() %>%
  print()
cat("\n")




# FDR>0.05 result ----
file_name <- sprintf(
  "C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_0.8/outcome/OQN_%s_%s_R2_0.8.txt",
   c("N", "T"),
  eqtl_threshold
)

total <- fread("C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_0.8/outcome/OQN_exp_different_r2_0.8.txt",header=T)

file_name <- sprintf(
  "C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_0.8/outcome/OQN_%s_%s_R2_0.8.txt",
   c("N", "T"),
  eqtl_threshold
)

  
# 對這些檔案生成N,T 顯著snp, probe, association number
for (i in file_name) {
  a <- fread(i)
  
  # 挑出FDR>0.05 probe
  probe <- setdiff(total$PROBE_ID, a$Probe)
  
  # N,T different FDR sig (OQN)
  sig_OQN_probe <- total %>% filter(sig_OQN ==1)
  sig_OQN_probe <- sig_OQN_probe[PROBE_ID %in% probe] %>% nrow()

  # N,T different Bonfi sig (OQN)
  Bonfi_OQN_probe <- total %>% filter(sig_OQN_Bonfi ==1) 
  Bonfi_OQN_probe <- Bonfi_OQN_probe[PROBE_ID %in% probe] %>% nrow()
  
  
  
  cat(i,"\n")
  cat("eQTL not sig: ","\n")
  cat("unique Probe number: ",probe %>% length(),"\n")
  
  cat("Among eQTL not sig,  N,T different FDR<0.05 (OQN): ","\n")
  cat("unique Probe number: ",sig_OQN_probe,"\n")
  
  cat("Among eQTL not sig,  N,T different Bonfi (OQN):","\n")
  cat("unique Probe number: ",Bonfi_OQN_probe,"\n")
  
  cat(rep("\n", 3))
}



# FDR<0.05 result ----


file_name <- sprintf(
  "C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_0.8/outcome/OQN_%s_%s_R2_0.8.txt",
   c("N", "T"),
  eqtl_threshold
)

# 對這些檔案生成N,T 顯著snp, probe, association number
for (i in file_name) {
  a <- fread(i)
  
  # FDR sig
  sig_fdr_snp <- unique(a, by="SNP") %>% nrow() 
  sig_fdr_asso <- a %>% nrow()
  sig_fdr_probe <- unique(a, by="ProbeID") %>% nrow()
  sig_fdr_gene <- unique(a, by="Gene") %>% nrow()
  
  # N,T different FDR sig (OQN)
  sig_OQN_snp <- a %>% filter(sig_OQN ==1) %>% unique( by="SNP") %>% nrow()
  sig_OQN_asso <- a %>% filter(sig_OQN ==1) %>% nrow()
  sig_OQN_probe <- a %>% filter(sig_OQN ==1) %>% unique( by="ProbeID") %>% nrow()
  sig_OQN_gene<- a %>% filter(sig_OQN ==1) %>% unique( by="Gene") %>% nrow()

  # N,T different Bonfi sig (OQN)
  Bonfi_OQN_snp <- a %>% filter(sig_OQN_Bonfi ==1) %>% unique( by="SNP") %>% nrow()
  Bonfi_OQN_asso <- a %>% filter(sig_OQN_Bonfi ==1) %>% nrow()
  Bonfi_OQN_probe <- a %>% filter(sig_OQN_Bonfi ==1) %>% unique( by="ProbeID") %>% nrow()
  Bonfi_OQN_gene <- a %>% filter(sig_OQN_Bonfi ==1) %>% unique( by="Gene") %>% nrow()
  
  
  
  cat(i,"\n")
  
  cat("eQTL sig: ","\n")
  cat("unique Probe number: ")
  print(sig_fdr_probe)
  cat("unique gene number: ")
  print(sig_fdr_gene)
  cat("unique SNPs number: ")
  print(sig_fdr_snp)
  cat("association number: ")
  print(sig_fdr_asso)
  
  

  
  cat("N,T different FDR<0.05 (OQN):","\n")
  cat("unique Probe number: ")
  print(sig_OQN_probe)
  cat("unique gene number: ")
  print(sig_OQN_gene)
  cat("unique SNPs number: ")
  print(sig_OQN_snp)
  cat("association number: ")
  print(sig_OQN_asso)
  
  
  cat("N,T different Bonfi (OQN):","\n")
  cat("unique Probe number: ")
  print(Bonfi_OQN_probe)
  cat("unique gene number: ")
  print(Bonfi_OQN_gene)
  cat("unique SNPs number: ")
  print(Bonfi_OQN_snp)
  cat("association number: ")
  print(Bonfi_OQN_asso)
  
  cat(rep("\n", 3))
}



# common probe in N,T ----
# N,T 都有cis eQTL 的probe，在N,T expression 有顯著不同的比例很高?


total <- fread("C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_0.8/outcome/OQN_exp_different_r2_0.8.txt",header=T)


# gt pvalue

a <- sprintf(
  "C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_0.8/outcome/OQN_N_%s_R2_0.8.txt",
  eqtl_threshold
) %>% fread()
b <- sprintf(
  "C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_0.8/outcome/OQN_T_%s_R2_0.8.txt",
  eqtl_threshold
) %>% fread()

both_probe <- intersect(a$Probe, b$Probe)

a_bon <- a[sig_pval_Bonfi==1,]
b_bon <- b[sig_pval_Bonfi==1,]
both_probe_bon <- intersect(a_bon$Probe,b_bon$Probe)

names(a)
# 都用 FDR(BH) 
both_probe %>% length()
sig_OQN_probe <- total %>% filter(sig_OQN ==1)
sig_OQN_probe <- sig_OQN_probe[PROBE_ID %in% both_probe] %>% nrow()
sig_OQN_probe


# 都用 Bonferroni
both_probe_bon %>% length()
probe_bon <- total %>% filter(sig_OQN_Bonfi ==1)
probe_bon <- probe_bon[PROBE_ID %in% both_probe_bon] %>% nrow()
probe_bon

cat(rep("\n",2))



# probe 分別有幾個snp ----
# 每個probe 選一個snp，很明顯每個probe 分別有1個snp

file_name <- c("OQN_N",
               "OQN_T")

plot_output_dir <- "C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_0.8/outcome"

for (i in file_name) {
  a <- sprintf(
    "C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_0.8/outcome/%s_%s_R2_0.8.txt",
    i, eqtl_threshold
  ) %>% fread()
  
  k <- a[, .N, by = "Probe"]
  if (eqtl_threshold == "FDR"){

plot_title <- ifelse( i == "OQN_N",
    sprintf("Number of significant cis-SNPs for %d probes in normal tissue (FDR < 0.05)", 
    nrow(k)),
    sprintf("Number of significant cis-SNPs for %d probes in tumor tissue (FDR < 0.05)",
     nrow(k))
  )
  } else if (eqtl_threshold == "bon"){
     plot_title <- ifelse( i == "OQN_N",
    sprintf("Number of significant cis-SNPs for %d probes in normal tissue (bon pass)", nrow(k)),
    sprintf("Number of significant cis-SNPs for %d probes in tumor tissue (bon pass)", nrow(k))
  )
  }
  
  
  png(sprintf("%s/%s_SNP_number_plot.png", plot_output_dir, i),
      width = 1600, height = 1200, res = 200)
  plot(k$N,ylab = "SNP Number", main = plot_title)
  dev.off()
  
  png(sprintf("%s/%s_SNP_number_hist.png", plot_output_dir, i),
      width = 1600, height = 1200, res = 200)
  hist(k$N,ylab = "Frequency",xlab = "SNP Number", main = plot_title, breaks=50, col = "skyblue") 
  dev.off()
  
  
  cat(i,"\n")
  cat("mean of SNP Number in each probe: ",mean(k$N), "\n") 
  cat("median of SNP Number in each probe: ",median(k$N), "\n")  
  cat("\n")
}


}







# 挑出這些 clumping 後的 snp ID，做表格 ----
maf_hg18_19 <- fread("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/cis_hg18_19_snp_MAF.txt",
  header = T
)
iteration <- 1
data_list <- vector("list")

for (pop in c("EUR", "FIN")) {
  for (eq in c("bon", "FDR")) {
    cat("\n開始整理 population =", pop, "\n")
    for (r2_threshold in c("0.6" ,"0.7","0.9", "0.8","no")) {
    base_dir <- sprintf("C:/Peter/repeatSNP_clumping_FIXpeople/r2_filter_%s/eQTL_%s/%s", r2_threshold, eq, pop)
    outcome_dir <- sprintf("%s/outcome", base_dir)
    trash_dir <- sprintf("%s/trash", base_dir)

    if (!dir.exists(outcome_dir) || !dir.exists(trash_dir)) {
      message(sprintf("Skip %s: missing outcome/trash directory under %s", pop, base_dir))
      next
    }

    outcome_prefix <- sprintf("%s/", outcome_dir)
    files <- c("2_fdr.txt", "2_qval.txt", "1_qval.txt", "1_fdr.txt")
    file_name <- paste0(outcome_prefix, files)
    
    for (i in seq_along(file_name)) {
      base <- basename(file_name[i])
      prefix <- sub("\\.txt$", "", base) # 取得 1_qval, 1_fdr 等前綴

      for (para in c(0.2)) {
        # 1000kb
        final_out <- sprintf("%s/1000kb_LD%s/%s.clumped", outcome_dir, para, prefix)
        if (!file.exists(final_out)) {
          next
        }
        a <- fread(final_out)
        k <- data.table(
          hg19_snpID = a$SNP,
          r2_threshold = r2_threshold,
          file = prefix,
          eqtl_threshold = eq,
          population = pop
        )
        data_list[[iteration]] <- k
        iteration <- iteration+1
      }
    }
    if (!any(lengths(data_list) > 0)) {
      message(sprintf("Skip %s %s: no clumped files found", pop, eq))
      next
    }
  }
  }
  }

result <- rbindlist(data_list, use.names = TRUE, fill = TRUE)
result_record <- maf_hg18_19[
  ,
  c("rsID", "hg18_snpID", "hg19_snpID")
][result, on = .(hg19_snpID)]
setkey(result_record, population, r2_threshold, eqtl_threshold, file)


# 加上是否 FDR<0.05，被紀錄再 eQTL 檔案
result_record_sub <- unique(result_record, by=c("hg18_snpID","r2_threshold"))
result_record[, FDR_sig := 0L]

for (r2 in c("0.6", "0.7", "0.9", "0.8", "no")) {

  eqtl_file <- sprintf(
    "C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_%s/outcome/OQN_N_FDR_R2_%s.txt",
    r2, r2
  ) %>% fread()
  
  af_clump_snphg18 <- result_record_sub[r2_threshold == r2, unique(na.omit(hg18_snpID))]
  matched_snp <- af_clump_snphg18[af_clump_snphg18 %in% eqtl_file$SNP]
  if (length(matched_snp) > 0) {
    result_record[r2_threshold == r2 & hg18_snpID %in% matched_snp, FDR_sig := 1L]
  }
}

fwrite(result_record,
  "C:/Peter/repeatSNP_clumping_FIXpeople/af_clumping_SNPID.txt",
  row.names = F, col.names = T, sep = "\t"
)


a <- fread("C:/Peter/repeatSNP_clumping_FIXgene/af_clumping_SNPID.txt")
a[,c("hg18_snpID", "hg19_snpID"):=NULL]
a <- a[,.(rsID = paste(rsID,collapse = ", ")), by= .(r2_threshold	,population,eqtl_threshold, file)]
fwrite(a,
  "C:/Peter/repeatSNP_clumping_FIXgene/af_clumping_SNPID_1.txt",
  row.names = F, col.names = T, sep = "\t"
)



# 挑出這些 clumping 後的 snp info，結合 eQTL info ----

result_record <- fread("C:/Peter/repeatSNP_clumping_FIXpeople/af_clumping_SNPID.txt")
result_record_sub <- unique(result_record, by=c("hg18_snpID","r2_threshold"))

wb <- createWorkbook()
has_sheet <- FALSE

for (r2 in c("0.6", "0.7", "0.9", "0.8", "no")) {
  eqtl_file <- sprintf(
    "C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_%s/outcome/OQN_N_FDR_R2_%s.txt",
    r2, r2
  ) %>% fread()
  
  af_clump_snphg18 <- result_record_sub[r2_threshold == r2, unique(na.omit(hg18_snpID))]
  af_clump_info <- eqtl_file[SNP %in% af_clump_snphg18]

  sheet_name <- sprintf("eQTL_N_%s", r2)
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, af_clump_info)
  has_sheet <- TRUE
}

if (has_sheet) {
  saveWorkbook(
    wb,
    "C:/Peter/repeatSNP_clumping_FIXpeople/af_clump_info.xlsx",
    overwrite = TRUE
  )
}






