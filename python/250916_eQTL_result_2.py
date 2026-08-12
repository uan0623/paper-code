# %% Packages
import gc
import re
import numpy as np
import pandas as pd

from statsmodels.stats.multitest import multipletests
# 先用 pip install qvalue 安裝 in 系統 terminal, 到這下載 qvalue https://github.com/nfusi/qvalue/tree/master
from qvalue.qvalue import qvalue_estimate

# 檔名多加 test_
# r2_threshold =0.8
# inputname = "C:/Peter/rawData_eQTL/trash/raw_maf_gt_N_pvalue.txt"
# R2filter = "D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/cis_hg18_19_snp_MAF.txt"
# outputname = "C:/Users/user/Desktop/test.txt"
# qvalue_type = 2


def compute_FDR(inputname, R2filter, outputname, qvalue_type):

# 判斷 inputname 是否為 string
    if not isinstance(inputname, str):
        raise TypeError("inputname 必須是字串")

    if not isinstance(outputname, str):
        raise TypeError("outputname 必須是字串")

    x = pd.read_csv(inputname, sep="\t")
    R2_filter = pd.read_csv(R2filter, sep="\t")

    # R2 filter
    x = x[x["SNP"].isin(R2_filter["hg18_snpID"])].copy()

    # Bonferroni, astype(int) 把資料轉成 integer 
    x["sig_pval_Bonfi"] = (
        x["p-value"] < 0.05 / len(x)
    ).astype(int)

    # BH FDR
    x["FDR"] = multipletests(
        x["p-value"],
        method="fdr_bh"
    )[1]

    # qvalue
    if qvalue_type == 1:

        pi0_hat = np.clip(
            np.sum(x["p-value"] > 0.7)
            / ((1 - 0.7) * len(x)),
            0,
            1,
        )

        x["qvalue"] = qvalue_estimate(
            x["p-value"],
            pi0=pi0_hat
        )

    elif qvalue_type == 2:

        x["qvalue"] = qvalue_estimate(
            x["p-value"]
        )

    # FDR filter
    x = x[x["FDR"] < 0.05].copy()

    x = x.sort_values("FDR")

    x["gene2"] = (
        x["gene"]
        .str.replace(
            r"^([^_]+_[^_]+).*",
            r"\1",
            regex=True,
        )
    )

    x.to_csv(
        outputname,
        sep="\t",
        index=False,
    )

if r2_threshold == "no":

    for tissue in ["N", "T"]:

        compute_FDR(
            f"C:/Peter/rawData_eQTL/trash/raw_maf_gt_{tissue}_pvalue.txt",
            "D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/cis_hg18_19_snp_MAF.txt",
            f"C:/Peter/rawData_eQTL/r2_filter_{r2_threshold}/outcome/test_raw_maf_gt_{tissue}_pvalue_FDR_R2_{r2_threshold}.txt",
            2,
        )

else:

    for tissue in ["N", "T"]:

        compute_FDR(
            f"C:/Peter/rawData_eQTL/trash/raw_maf_gt_{tissue}_pvalue.txt",
            f"D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/cis_snp_R2_{r2_threshold}.txt",
            f"C:/Peter/rawData_eQTL/r2_filter_{r2_threshold}/outcome/test_raw_maf_gt_{tissue}_pvalue_FDR_R2_{r2_threshold}.txt",
            2,
        )


# 5j4ur,3 --

gt = pd.read_csv(
    "C:/Peter/rawData_eQTL/trash/raw_maf_gt_N_pvalue.txt",
    sep="\t",
)

if r2_threshold == "no":

    a = pd.read_csv(
        "D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/cis_hg18_19_snp_MAF.txt",
        sep="\t",
    )

else:

    a = pd.read_csv(
        f"D:/oral_cancer/expression/outcome/multiple_nucleotide_variant/cis_snp_R2_{r2_threshold}.txt",
        sep="\t",
    )

gt_cisSNP = gt[
    gt["SNP"].isin(a["hg18_snpID"])
].copy()

probe_cisSNP = (
    gt_cisSNP["gene"]
    .str.replace(
        r"^([^_]+_[^_]+).*",
        r"\1",
        regex=True,
    )
    .drop_duplicates()
)

pd.DataFrame(
    {
        "Has_association_probe": probe_cisSNP
    }
).to_csv(
    f"C:/Peter/rawData_eQTL/r2_filter_{r2_threshold}/outcome/test_Has_association_probe_R2_{r2_threshold}.txt",
    sep="\t",
    index=False,
)





gc.collect()

df = pd.read_csv(
    "C:/Peter/rawData_eQTL/outcome/raw_exp_different.txt",
    sep="\t",
)

for file in ["N", "T"]:

    a = pd.read_csv(
        f"C:/Peter/rawData_eQTL/r2_filter_{r2_threshold}/outcome/raw_maf_gt_{file}_pvalue_FDR_R2_{r2_threshold}.txt",
        sep="\t",
    )

    a["gene2"] = (
        a["gene"]
        .str.replace(
            r"^([^_]+_[^_]+).*",
            r"\1",
            regex=True,
        )
    )

    #################################
    # 所有顯著 eQTL 數量
    #################################

    sig = (
        a.groupby("gene2")
        .size()
        .reset_index(name="N")
    )

    df = df.merge(
        sig,
        left_on="PROBE_ID",
        right_on="gene2",
        how="left",
    )

    df["N"] = (
        df["N"]
        .fillna(0)
        .astype(int)
    )

    df.rename(
        columns={
            "N": f"{r2_threshold}_{file}_sigCis-SNP_number"
        },
        inplace=True,
    )

    if "gene2" in df.columns:
        df.drop(columns="gene2", inplace=True)

    #################################
    # Bonferroni
    #################################

    a2 = a[
        a["sig_pval_Bonfi"] == 1
    ]

    sig = (
        a2.groupby("gene2")
        .size()
        .reset_index(name="N")
    )

    df = df.merge(
        sig,
        left_on="PROBE_ID",
        right_on="gene2",
        how="left",
    )

    df["N"] = (
        df["N"]
        .fillna(0)
        .astype(int)
    )

    df.rename(
        columns={
            "N": f"{r2_threshold}_{file}_sigCis-SNP_number (bon)"
        },
        inplace=True,
    )

    if "gene2" in df.columns:
        df.drop(columns="gene2", inplace=True)

#################################
# format
#################################

df["pval_raw"] = df["pval_raw"].map(
    lambda x: "{:.4e}".format(x)
)

df.to_csv(
    f"C:/Peter/rawData_eQTL/r2_filter_{r2_threshold}/outcome/test_raw_exp_different_r2_{r2_threshold}.txt",
    sep="\t",
    index=False,
)