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

library(data.table)
library(openxlsx)
library(magrittr)


# 要填 excel，所需檔案路徑
# BHplot
# C:/Peter/repeatSNP_clumping_FIXpeople/r2_filter_0.8/outcome/BHplot_r2_0.8_original.png
# C:/Peter/repeatSNP_clumping_FIXpeople/r2_filter_0.8/outcome/BHplot_r2_0.8_scale.png

# FDR_qvalue
# C:/Peter/repeatSNP_clumping_FIXpeople/r2_filter_0.8/outcome/FDR_qvalue_r2_0.8_original.png
# C:/Peter/repeatSNP_clumping_FIXpeople/r2_filter_0.8/outcome/FDR_qvalue_r2_0.8_scale_0.05.png

# N,T Comparison
# C:/Peter/rawData_eQTL/outcome/paper_table_r2_0.8_FDR.xlsx
# C:/Peter/OQN_FIXpeople_before_eQTL/outcome/paper_table_r2_0.8_bon.xlsx

# probe 對到幾個 snp 圖
# C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_0.8/outcome/OQN_N_SNP_number_plot.png

# probe_info
# C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_0.8/outcome/OQN_exp_different_r2_0.8.txt

# maf_gt_N_FDR
# C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_0.8/outcome/OQN_N_FDR_R2_0.8.txt

# Finngen_1000gEUR
# C:/Peter/OQN_FIXpeople_before_eQTL/outcome/OQN_FDR_EUR_r2_filter_summary.xlsx
# C:/Peter/OQN_FIXpeople_before_eQTL/outcome/OQN_FDR_FIN_r2_filter_summary.xlsx
# C:/Peter/OQN_FIXpeople_before_eQTL/outcome/OQN_bon_EUR_r2_filter_summary.xlsx
# C:/Peter/OQN_FIXpeople_before_eQTL/outcome/OQN_bon_FIN_r2_filter_summary.xlsx

# clumping_pass_FDR
# C:/Peter/repeatSNP_clumping_FIXpeople/af_clump_info.xlsx
# C:/Peter/repeatSNP_clumping_FIXpeople/af_clumping_SNPID.txt
# C:/Peter/repeatSNP_clumping_FIXpeople/af_clumping_SNPID_1.txt

# 5783_snp_3_tech_FIN
# C:/Peter/OQN_FIXpeople_before_eQTL/r2_filter_0.8/outcome/OQN_MixFinngenPval_7_FIN.txt

# rs4975538_eQTL
# C:/Peter/rs4975538_permutation/outcome/rawData_eQTL/rs4975538_N_eQTL.txt
# C:/Peter/rs4975538_permutation/outcome/OQN_FIXpeople_before_eQTL/rs4975538_asso_probe_correlation.png


# 要上 paper table
# C:/Peter/rawData_eQTL/outcome/paper_table_r2_0.8_FDR.xlsx
## package ----

1 + 1



maf <- fread("D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/cis_hg18_19_snp_MAF.txt", header = T)
maf_1 <- maf[MAF > 0.05]
maf_1[nchar(ALT) > 1 & nchar(REF) > 1]
maf_1[nchar(ALT) == nchar(REF) & nchar(REF) > 1]
qvalue(seq(0, 0.1, length.out = 11))$qvalues



set.seed(123)
m <- 1000
p_null <- runif(800)
p_alt <- rbeta(200, shape1 = 0.3, shape2 = 1)
p <- c(p_null, p_alt)
p <- sample(p)


# test ----
qvalue:::qvalue
function (p, fdr.level = NULL, pfdr = FALSE, lfdr.out = TRUE, 
    pi0 = NULL, ...) 
{
  p_in <- qvals_out <- lfdr_out <- p
  rm_na <- !is.na(p)
  p <- p[rm_na]
  if (min(p) < 0 || max(p) > 1) {
    stop("p-values not in valid range [0, 1].")
  } else if (!is.null(fdr.level) && (fdr.level <= 0 || fdr.level >
    1)) {
    stop("'fdr.level' must be in (0, 1].")
  }
  if (is.null(pi0)) {
    pi0s <- pi0est(p, ...)
  } else {
    if (pi0 > 0 && pi0 <= 1) {
      pi0s = list()
      pi0s$pi0 = pi0
    } else {
      stop("pi0 is not (0,1]")
    }
  }
  m <- length(p)
  i <- m:1L
  o <- order(p, decreasing = TRUE)
  ro <- order(o)
  if (pfdr) {
    qvals <- pi0s$pi0 * pmin(1, cummin(p[o] * m / (i * (1 -
      (1 - p[o])^m))))[ro]
  } else {
    qvals <- pi0s$pi0 * pmin(1, cummin(p[o] * m / i))[ro]
  }
  qvals_out[rm_na] <- qvals
  if (lfdr.out) {
    lfdr <- lfdr(p = p, pi0 = pi0s$pi0, ...)
    lfdr_out[rm_na] <- lfdr
  } else {
    lfdr_out <- NULL
  }
  if (!is.null(fdr.level)) {
    retval <- list(
      call = match.call(), pi0 = pi0s$pi0, qvalues = qvals_out,
      pvalues = p_in, lfdr = lfdr_out, fdr.level = fdr.level,
      significant = (qvals <= fdr.level), pi0.lambda = pi0s$pi0.lambda,
      lambda = pi0s$lambda, pi0.smooth = pi0s$pi0.smooth
    )
  } else {
    retval <- list(
      call = match.call(), pi0 = pi0s$pi0, qvalues = qvals_out,
      pvalues = p_in, lfdr = lfdr_out, pi0.lambda = pi0s$pi0.lambda,
      lambda = pi0s$lambda, pi0.smooth = pi0s$pi0.smooth
    )
  }
  class(retval) <- "qvalue"
  return(retval)
}



qvalue:::pi0est
function(p, lambda = seq(0.05, 0.95, 0.05), pi0.method = c(
           "smoother",
           "bootstrap"
         ), smooth.df = 3, smooth.log.pi0 = FALSE, ...) {
  rm_na <- !is.na(p)
  p <- p[rm_na]
  pi0.method <- match.arg(pi0.method)
  m <- length(p)
  lambda <- sort(lambda)
  ll <- length(lambda)
  if (min(p) < 0 || max(p) > 1) {
    stop("ERROR: p-values not in valid range [0, 1].")
  } else if (ll > 1 && ll < 4) {
    stop(sprintf(paste(
      "ERROR:", paste("length(lambda)=",
        ll, ".",
        sep = ""
      ), "If length of lambda greater than 1,",
      "you need at least 4 values."
    )))
  } else if (min(lambda) < 0 || max(lambda) >= 1) {
    stop("ERROR: Lambda must be within [0, 1).")
  }
  if (ll == 1) {
    pi0 <- mean(p >= lambda) / (1 - lambda)
    pi0.lambda <- pi0
    pi0 <- min(pi0, 1)
    pi0Smooth <- NULL
  } else {
    ind <- length(lambda):1
    pi0 <- cumsum(tabulate(findInterval(p, vec = lambda))[ind]) / (length(p) *
      (1 - lambda[ind]))
    pi0 <- pi0[ind]
    pi0.lambda <- pi0
    if (pi0.method == "smoother") {
      if (smooth.log.pi0) {
        pi0 <- log(pi0)
        spi0 <- smooth.spline(lambda, pi0, df = smooth.df)
        pi0Smooth <- exp(predict(spi0, x = lambda)$y)
        pi0 <- min(pi0Smooth[ll], 1)
      } else {
        spi0 <- smooth.spline(lambda, pi0, df = smooth.df)
        pi0Smooth <- predict(spi0, x = lambda)$y
        pi0 <- min(pi0Smooth[ll], 1)
      }
    } else if (pi0.method == "bootstrap") {
      minpi0 <- quantile(pi0, prob = 0.1)
      W <- sapply(lambda, function(l) sum(p >= l))
      mse <- (W / (m^2 * (1 - lambda)^2)) * (1 - W / m) + (pi0 -
        minpi0)^2
      pi0 <- min(pi0[mse == min(mse)], 1)
      pi0Smooth <- NULL
    } else {
      stop("ERROR: pi0.method must be one of \"smoother\" or \"bootstrap\".")
    }
  }
  if (pi0 <= 0) {
    stop("ERROR: The estimated pi0 <= 0. Check that you have valid p-values or use a different range of lambda.")
  }
  return(list(
    pi0 = pi0, pi0.lambda = pi0.lambda, lambda = lambda,
    pi0.smooth = pi0Smooth
  ))
}
