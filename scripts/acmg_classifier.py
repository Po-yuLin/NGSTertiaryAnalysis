#!/usr/bin/env python3
"""
 * =========================================================
 * WGS/WES Germline Analysis Pipeline
 * =========================================================
 * Author   : Po-Yu Lin (林伯昱)
 * Institute: Department of Neurology and
 *            Department of Genomic Medicine,
 *            National Cheng Kung University Hospital
 * Contact  : p88124019@gs.ncku.edu.tw
 *
 * Copyright (c) 2026, Po-Yu Lin (林伯昱)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * DISCLAIMER: This pipeline is provided "as is" without
 * warranty of any kind. The authors and their institution
 * make no representations or warranties regarding the
 * accuracy, completeness, or suitability of the analysis
 * results for any clinical or research purpose. Users are
 * solely responsible for validating and interpreting all
 * results.
 * =========================================================

acmg_classifier.py
==================
讀取 parse_vep_csq.py 產生的 56 欄 TSV，
依照 ACMG/AMP 2015 + ClinGen SVI 更新版規則，
對每個 variant 評分並輸出 ACMG 分類。

v3.1 更新（對照 AutoACMG 文件修正）：
  1. PM2 降為 PM2_Supporting（+1 分，符合 ClinGen SVI 2020）
  2. PM2 加入 MOI 考量（AD 用 AF 閾值；AR/XL 用 allele count 邏輯近似）
  3. PVS1 觸發時停用 PP3 和 PM4（避免 double counting）
  4. PM4 加入 stop_lost variant
  5. BP7 加入 PhyloP 保守性判斷（< 3.58 才觸發，參考 Pejaver 2022）
  6. BP3 加入 PhyloP 保守性判斷（PhyloP < 3.58 = 低保守性 → 支持 BP3）

Phase 1 實作的 criteria：
  BA1          ── Stand-alone Benign：gnomAD 全人群 AF > 5%
  BS1          ── Strong Benign：gnomAD 全人群 AF > 1%
  PM2_Supporting── Supporting Pathogenic：EAS 族群極罕見（ClinGen SVI 2020 降級）
  PP3          ── Supporting~Strong Pathogenic：計算工具預測（cascade fallback）
  BP4          ── Supporting~Strong Benign：計算工具預測（同上，互斥）
  BP7          ── Supporting Benign：Synonymous + 低 splice 影響 + 低保守性
  PM4          ── Moderate Pathogenic：Inframe indel 或 stop_lost
  PVS1         ── Very Strong Pathogenic：LOFTEE HC + ClinGen HI score=3
  BP3          ── Supporting Benign：Inframe indel 在 repeat region + 低保守性

Phase 2 criteria（需額外資料庫）：
  PS1, PM5, PP2, BP1, PM1

輸入：
  --input       parse_vep_csq.py 產生的 annotated TSV（.tsv 或 .tsv.gz）
  --clingen_hi  ClinGen Dosage Sensitivity TSV（PVS1 HI 過濾）
  --gene_moi    build_gene_moi.py 產生的 gene_moi.tsv.gz（PM2 MOI 判斷）
  --output      輸出 TSV 路徑

輸出新增欄位（在原 56 欄後面）：
  ACMG_CRITERIA   ── 觸發的所有 criteria，逗號分隔
  ACMG_SCORE      ── 數值化分數（用於排序）
  ACMG_CLASS      ── 最終分類
  ACMG_NOTES      ── 簡要說明

參考文獻：
  Richards et al. Genet Med 2015（ACMG/AMP 2015 原始規則）
  Tavtigian et al. Hum Mutat 2020（點數系統）
  ClinGen SVI Working Group 2020（PM2 降為 Supporting）
  Bergquist et al. Genet Med 2025（AlphaMissense/ESM1b/VARITY 校準閾值）
  Pejaver et al. Am J Hum Genet 2022（BayesDel 校準閾值，PhyloP 閾值 3.58）
  ClinGen Dosage Sensitivity Map（PVS1 HI gene list）
  AutoACMG（criteria 詮釋參考）https://github.com/bihealth/auto-acmg

使用方式：
  python3 acmg_classifier.py \\
      --input      NA12878_WES.snv_indel.full.annotated.tsv \\
      --clingen_hi /path/to/clingen/ClinGen_gene_curation_list_GRCh38.tsv \\
      --gene_moi   /path/to/clingen/gene_moi.tsv.gz \\
      --output     NA12878_WES.snv_indel.acmg.tsv

  # 所有資料庫都是 optional，不提供時退回簡化版
"""

import argparse
import gzip
import os
import sys


# ══════════════════════════════════════════════════════════════════
# 第零部分：外部資料庫載入
# ══════════════════════════════════════════════════════════════════

def load_clingen_hi(tsv_path: str | None) -> dict[str, int]:
    """
    載入 ClinGen Dosage Sensitivity TSV。
    回傳 gene_symbol → haploinsufficiency_score（int）。

    HI score 意義：
        3  = Sufficient evidence（支持 PVS1）
        2  = Emerging evidence（不觸發 PVS1，但標記供審閱）
        1  = Little evidence
        0  = No evidence
        30 = AR gene only（不適用 PVS1）
        40 = Dosage sensitivity unlikely

    不存在或路徑為 None → 回傳空 dict（PVS1 退回簡化版）
    """
    if not tsv_path or not os.path.exists(tsv_path):
        if tsv_path:
            print(f"[WARNING] ClinGen HI TSV 不存在：{tsv_path}", file=sys.stderr)
        return {}

    hi_table: dict[str, int] = {}
    n_hi3 = 0

    with open(tsv_path, "r", encoding="utf-8") as f:
        header_found = False
        gene_col = hi_col = 0

        for line in f:
            line = line.rstrip("\n")
            if not line.strip():
                continue

            if line.startswith("#Gene Symbol") or line.startswith("Gene Symbol"):
                cols = [c.strip() for c in line.lstrip("#").split("\t")]
                gene_col = cols.index("Gene Symbol") if "Gene Symbol" in cols else 0
                hi_col   = cols.index("Haploinsufficiency Score") if "Haploinsufficiency Score" in cols else 3
                header_found = True
                continue

            if line.startswith("#") or not header_found:
                continue

            parts = line.split("\t")
            if len(parts) <= max(gene_col, hi_col):
                continue

            gene   = parts[gene_col].strip()
            hi_raw = parts[hi_col].strip()
            if not gene or not hi_raw:
                continue

            try:
                score = int(hi_raw.split()[0])
            except (ValueError, IndexError):
                continue

            hi_table[gene] = score
            if score == 3:
                n_hi3 += 1

    print(f"[acmg_classifier] ClinGen HI：{len(hi_table):,} 個基因，HI=3：{n_hi3:,} 個", file=sys.stderr)
    return hi_table


def load_gene_moi(tsv_path: str | None) -> dict[str, str]:
    """
    載入 build_gene_moi.py 產生的 gene_moi.tsv.gz。
    回傳 gene_symbol → MOI 字串（"AD" / "AR" / "XL" / "AD/AR" / "MT"）。

    不存在或路徑為 None → 回傳空 dict（PM2 退回純 AF 閾值判斷）
    """
    if not tsv_path or not os.path.exists(tsv_path):
        if tsv_path:
            print(f"[WARNING] gene_moi.tsv.gz 不存在：{tsv_path}", file=sys.stderr)
        return {}

    moi_table: dict[str, str] = {}
    opener = gzip.open if tsv_path.endswith(".gz") else open

    with opener(tsv_path, "rt", encoding="utf-8") as f:
        for i, line in enumerate(f):
            line = line.rstrip("\n")
            if i == 0 or not line:
                continue  # 跳過 header 和空行
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            gene = parts[0].strip()
            moi  = parts[1].strip()
            if gene and moi:
                moi_table[gene] = moi

    print(f"[acmg_classifier] Gene MOI table：{len(moi_table):,} 個基因", file=sys.stderr)
    return moi_table


# 全域查表（process_tsv() 啟動時載入）
_HI_TABLE:  dict[str, int] = {}
_MOI_TABLE: dict[str, str] = {}


# ══════════════════════════════════════════════════════════════════
# 第一部分：ACMG 點數系統（Tavtigian 2020）
# ══════════════════════════════════════════════════════════════════

CRITERIA_POINTS = {
    # ── 致病性（正分）
    "PVS1":             8,   # Very Strong
    "PS1":              4,   # Strong
    "PS2":              4,
    "PS3":              4,
    "PS4":              4,
    "PM1":              2,   # Moderate
    "PM2":              2,   # 保留舊名稱相容性（不再使用，見 PM2_Supporting）
    "PM2_Supporting":   1,   # ★ ClinGen SVI 2020：PM2 降為 Supporting
    "PM3":              2,
    "PM4":              2,
    "PM5":              2,
    "PM6":              2,
    "PP1":              1,   # Supporting
    "PP2":              1,
    "PP3_Supporting":   1,
    "PP3_Moderate":     2,
    "PP3_P3":           3,
    "PP3_Strong":       4,
    "PP4":              1,
    "PP5":              1,
    # ── 良性（負分）
    "BA1":             -8,   # Stand-alone
    "BS1":             -4,   # Strong
    "BS2":             -4,
    "BS3":             -4,
    "BS4":             -4,
    "BP1":             -1,   # Supporting
    "BP2":             -1,
    "BP3":             -1,
    "BP4_Supporting":  -1,
    "BP4_Moderate":    -2,
    "BP4_M3":          -3,
    "BP4_Strong":      -4,
    "BP5":             -1,
    "BP6":             -1,
    "BP7":             -1,
}

# 分類閾值（Tavtigian 2020 Table 3）
PATHOGENIC_THRESHOLD  = 10
LIKELY_PATHOGENIC_MIN = 6
LIKELY_BENIGN_MAX     = -1
BENIGN_THRESHOLD      = -7

# PhyloP100way 保守性閾值（Pejaver et al. 2022）
# >= 3.58 → 高保守性（不觸發 BP3 / BP7）
# <  3.58 → 低保守性（支持 BP3 / BP7）
PHYLOP_CONSERVATION_THRESHOLD = 3.58


# ══════════════════════════════════════════════════════════════════
# 第二部分：TSV 欄位定義（56 欄）
# ══════════════════════════════════════════════════════════════════

EXPECTED_COLUMNS = [
    "CHROM", "POS", "REF", "ALT", "RS_ID",                          # 0-4
    "GENE", "TRANSCRIPT", "TRANSCRIPT_TYPE",                         # 5-7
    "HGVS_C", "HGVS_P", "CONSEQUENCE", "IMPACT",                    # 8-11
    "EXON", "INTRON", "MANE_ALL",                                    # 12-14
    "CALLERS", "DP_DV", "AD_DV", "VAF_DV", "DP_HC", "AD_HC",        # 15-20
    "ZYGOSITY", "GT_DV", "GT_HC",                                    # 21-23
    "GNOMAD_G_AF", "GNOMAD_G_EAS_AF",                               # 24-25
    "GNOMAD_E_AF", "GNOMAD_E_EAS_AF",                               # 26-27
    "GNOMAD_E_AF_DBNSFP", "GNOMAD_E_EAS_AF_DBNSFP", "TG_EAS_AF",   # 28-30
    "CLINVAR_SIG", "CLINVAR_STARS", "CLINVAR_DN",                    # 31-33
    "CLINVAR_SIGCONF", "CLINVAR_VARIATION_ID", "OMIM_IDS",           # 34-36
    "LOFTEE", "LOFTEE_FILTER", "LOFTEE_FLAGS", "LOFTOOL",           # 37-40
    "BAYESDEL_NOAF", "BAYESDEL_NOAF_PRED",                          # 41-42
    "ALPHAMISSENSE", "ALPHAMISSENSE_PRED",                           # 43-44
    "ESM1B", "ESM1B_PRED",                                           # 45-46
    "VARITY_R", "SIFT", "SIFT_PRED",                                # 47-49
    "DANN", "PHACTBOOST", "PHYLOP100", "GERP",                      # 50-53
    "PKNN_LLR", "PKNN_EVIDENCE",                                     # 54-55
    "PANGOLIN_SCORE", "PANGOLIN_DETAIL",                             # 56-57
    "DOMAINS", "SWISSPROT",                                          # 58-59
    "HGNC_ID",                                                       # 60  ★ v3.1
]

NEW_COLUMNS = ["ACMG_CRITERIA", "ACMG_SCORE", "ACMG_CLASS", "ACMG_NOTES"]


# ══════════════════════════════════════════════════════════════════
# 第三部分：輔助函式
# ══════════════════════════════════════════════════════════════════

def to_float(value: str) -> float | None:
    """字串 → float，缺失值（"." 或空字串）→ None"""
    if not value or value == ".":
        return None
    try:
        return float(value)
    except ValueError:
        return None


def max_af(row: dict) -> float | None:
    """取 gnomAD genome + exome 全人群 AF 的最大值（BA1/BS1 用）"""
    vals = [to_float(row["GNOMAD_G_AF"]), to_float(row["GNOMAD_E_AF"])]
    valid = [v for v in vals if v is not None]
    return max(valid) if valid else None


def min_eas_af(row: dict) -> float | None:
    """取所有 EAS AF 欄位的最小值（PM2 用）"""
    vals = [
        to_float(row["GNOMAD_G_EAS_AF"]),
        to_float(row["GNOMAD_E_EAS_AF"]),
        to_float(row["GNOMAD_E_EAS_AF_DBNSFP"]),
        to_float(row["TG_EAS_AF"]),
    ]
    valid = [v for v in vals if v is not None]
    return min(valid) if valid else None


def max_global_af(row: dict) -> float | None:
    """取所有 gnomAD 全人群 AF（genome + exome）的最大值（PM2 AD 模式用）"""
    vals = [to_float(row["GNOMAD_G_AF"]), to_float(row["GNOMAD_E_AF"])]
    valid = [v for v in vals if v is not None]
    return max(valid) if valid else None


def has_consequence(consequence_str: str, target: str) -> bool:
    """
    VEP CONSEQUENCE 欄位可能有多個值（& 分隔），
    檢查是否包含特定的 consequence。
    """
    if not consequence_str or consequence_str == ".":
        return False
    return target in consequence_str.split("&")


# ══════════════════════════════════════════════════════════════════
# 第四部分：Criteria 判斷函式
# ══════════════════════════════════════════════════════════════════

def check_BA1(row: dict) -> bool:
    """
    BA1 ── Stand-alone Benign
    gnomAD 全人群 AF > 5%。觸發時直接判 Benign，其他 criteria 不再考慮。
    """
    af = max_af(row)
    return af is not None and af > 0.05


def check_BS1(row: dict) -> bool:
    """
    BS1 ── Strong Benign
    gnomAD 全人群 AF > 1%。
    """
    af = max_af(row)
    return af is not None and af > 0.01


def check_PVS1(row: dict) -> tuple[bool, str]:
    """
    PVS1 ── Very Strong Pathogenic（LOFTEE HC + ClinGen HI=3）

    回傳 (triggered: bool, note: str)

    HI score 意義：
        3  → 觸發 PVS1
        2  → 不觸發，但標記 manual_review（emerging evidence）
        30 → AR gene only，不觸發
        其他 → 不觸發
        不在 ClinGen → 不觸發，標記 manual_review
        _HI_TABLE 為空 → 簡化版，直接觸發（向後相容）
    """
    loftee        = row.get("LOFTEE", ".")
    loftee_filter = row.get("LOFTEE_FILTER", ".")
    gene          = row.get("GENE", ".")

    # Step 1：LOFTEE HC 才進入後續判斷
    if not (loftee == "HC" and loftee_filter == "."):
        return False, ""

    # Step 2：若 HI table 未載入 → 簡化版
    if not _HI_TABLE:
        return True, f"PVS1:LOFTEE=HC(simplified)"

    # Step 3：查 HI score
    hi_score = _HI_TABLE.get(gene)

    if hi_score is None:
        return False, f"PVS1_NOT_TRIGGERED:gene={gene},HI=not_in_ClinGen(manual_review)"
    if hi_score == 3:
        return True, f"PVS1:LOFTEE=HC,gene={gene},HI=3"
    if hi_score == 2:
        return False, f"PVS1_NOT_TRIGGERED:gene={gene},HI=2(emerging,manual_review)"
    if hi_score == 30:
        return False, f"PVS1_NOT_TRIGGERED:gene={gene},HI=30(AR_only)"
    return False, f"PVS1_NOT_TRIGGERED:gene={gene},HI={hi_score}(insufficient)"


def check_PM2(row: dict) -> tuple[bool, str]:
    """
    PM2_Supporting ── Supporting Pathogenic（ClinGen SVI 2020 降級）

    判斷邏輯（參考 AutoACMG + ClinGen SVI 2020）：

    ★ 取得基因的 MOI（from _MOI_TABLE）：

    AD 基因（或 MOI 未知）：
        gnomAD 全人群 AF < 0.0001（0.01%）→ 觸發
        （AD 基因一個 haplotype 就足以致病，用全人群 AF）

    AR / XL 基因：
        EAS 族群 AF < 0.00001（0.001%）→ 觸發
        （AR 需要兩個 allele，理論上可以容忍稍高的 carrier frequency，
         但我們用 EAS 族群更保守的閾值來篩選台灣病人）

    AD/AR 基因（同時有兩種遺傳模式）：
        取 AD 的條件（較嚴格的全人群 AF < 0.0001）

    MOI 未知：
        退回 EAS AF < 0.00001（保守做法）

    回傳 (triggered: bool, note: str)
    """
    gene = row.get("GENE", ".")
    moi  = _MOI_TABLE.get(gene, "Unknown") if _MOI_TABLE else "Unknown"

    # AD 或 AD/AR 或 Unknown → 用全人群 AF 閾值
    if "AD" in moi or moi == "Unknown":
        global_af = max_global_af(row)
        if global_af is None:
            # 不在 gnomAD → 視為極罕見，觸發
            return True, f"PM2_Supporting:gene={gene},MOI={moi},gnomAD_AF=absent"
        if global_af < 0.0001:
            return True, f"PM2_Supporting:gene={gene},MOI={moi},gnomAD_AF={global_af:.6f}<0.01%"
        return False, ""

    # AR / XL → 用 EAS AF 閾值
    eas_min = min_eas_af(row)
    if eas_min is None:
        return True, f"PM2_Supporting:gene={gene},MOI={moi},EAS_AF=absent"
    if eas_min < 0.00001:
        return True, f"PM2_Supporting:gene={gene},MOI={moi},min_EAS_AF={eas_min:.6f}<0.001%"
    return False, ""


def _score_to_pp3_bp4(tool_name: str, score: float) -> str | None:
    """
    將單一工具的分數對應到 PP3/BP4 強度字串。
    閾值來源：Bergquist et al. Genet Med 2025 Table 1。

    ★ ESM1b 方向相反：分數越小（越負）越致病。
    ★ BayesDel_noAF：BP4_Strong（-4）不存在，最強良性是 BP4_M3（-3）。
    ★ AlphaMissense 和 ESM1b：BP4_Strong（-4）不存在，最強良性是 BP4_M3（-3）。
    ★ VARITY_R：唯一有 BP4_Strong（-4）的工具。

    回傳：
        "PP3_Supporting/Moderate/P3/Strong" 或
        "BP4_Supporting/Moderate/M3/Strong"
        None → indeterminate（不觸發任何 criteria）
    """
    if tool_name == "ALPHAMISSENSE":
        if score <= 0.070: return "BP4_M3"
        if score <= 0.099: return "BP4_Moderate"
        if score <= 0.169: return "BP4_Supporting"
        if score <= 0.791: return None          # indeterminate
        if score <= 0.905: return "PP3_Supporting"
        if score <= 0.971: return "PP3_Moderate"
        if score <= 0.989: return "PP3_P3"
        return "PP3_Strong"

    elif tool_name == "ESM1B":
        # ★ 方向相反：分數越小 = 越致病
        if score >= 8.8:   return "BP4_M3"
        if score >= -3.1:  return "BP4_Moderate"
        if score >= -6.3:  return "BP4_Supporting"
        if score >= -10.6: return None          # indeterminate
        if score >= -12.1: return "PP3_Supporting"
        if score >= -13.9: return "PP3_Moderate"
        if score >= -24.0: return "PP3_P3"
        return "PP3_Strong"

    elif tool_name == "VARITY_R":
        if score <= 0.036: return "BP4_Strong"
        if score <= 0.063: return "BP4_M3"
        if score <= 0.116: return "BP4_Moderate"
        if score <= 0.251: return "BP4_Supporting"
        if score <= 0.674: return None          # indeterminate
        if score <= 0.841: return "PP3_Supporting"
        if score <= 0.914: return "PP3_Moderate"
        if score <= 0.964: return "PP3_P3"
        return "PP3_Strong"

    elif tool_name == "BAYESDEL_NOAF":
        # BP4_Strong（-4）不存在，最強良性是 BP4_M3（-3）
        if score <= -0.520: return "BP4_M3"
        if score <= -0.360: return "BP4_Moderate"
        if score <= -0.180: return "BP4_Supporting"
        if score <= 0.129:  return None          # indeterminate
        if score <= 0.269:  return "PP3_Supporting"
        if score <= 0.409:  return "PP3_Moderate"
        if score <= 0.499:  return "PP3_P3"
        return "PP3_Strong"

    return None


def check_PP3(row: dict) -> tuple[str | None, str]:
    """
    PP3 / BP4 ── 計算工具預測（cascade fallback）

    優先順序：
        1. P-KNN（PKNN_EVIDENCE，已是 evidence 字串）
        2. AlphaMissense
        3. ESM1b
        4. VARITY_R
        5. BayesDel_noAF

    只要第一個有值且不在 indeterminate 範圍的工具給出結論就停止。
    PP3 和 BP4 兩個方向由同一個 cascade 決定（互斥）。

    回傳 (evidence_string | None, tool_name_used)
    """
    # Step 1：P-KNN
    pknn = row.get("PKNN_EVIDENCE", ".")
    if pknn and pknn != ".":
        if pknn in ("PP3_Supporting", "PP3_Moderate", "PP3_Strong",
                    "BP4_Supporting", "BP4_Moderate", "BP4_Strong"):
            return pknn, "PKNN"

    # Step 2-5：cascade
    for tool, col in [
        ("ALPHAMISSENSE", "ALPHAMISSENSE"),
        ("ESM1B",         "ESM1B"),
        ("VARITY_R",      "VARITY_R"),
        ("BAYESDEL_NOAF", "BAYESDEL_NOAF"),
    ]:
        score = to_float(row.get(col, "."))
        if score is not None:
            result = _score_to_pp3_bp4(tool, score)
            if result is not None:
                return result, f"{tool}({score:.3f})"

    return None, "none"


def check_PM4(row: dict) -> tuple[bool, str]:
    """
    PM4 ── Moderate Pathogenic
    條件（符合其一即觸發）：
        1. inframe_insertion 或 inframe_deletion（不在 repeat region）
        2. stop_lost（終止碼突變，蛋白質延長）

    ★ v3.1 新增：stop_lost（參考 AutoACMG + Richards 2015 原文）
    ★ 注意：PVS1 觸發時 PM4 會被停用（classify_variant 處理）

    回傳 (triggered: bool, reason: str)
    """
    csq = row.get("CONSEQUENCE", "")

    if has_consequence(csq, "stop_lost"):
        return True, f"PM4:stop_lost"

    if has_consequence(csq, "inframe_insertion") or has_consequence(csq, "inframe_deletion"):
        return True, f"PM4:{csq}"

    return False, ""


def check_BP3(row: dict) -> tuple[bool, str]:
    """
    BP3 ── Supporting Benign
    條件（全部滿足才觸發）：
        1. inframe indel（不是 stop_lost，那個是 PM4 的範疇）
        2. DOMAINS 欄位有 repeat 相關標記（在 repeat region）
        3. PhyloP100 < 3.58（低保守性）

    ★ v3.1 新增：PhyloP 保守性條件（參考 AutoACMG）
    ★ PhyloP 缺失時：保守地不觸發 BP3

    回傳 (triggered: bool, reason: str)
    """
    csq = row.get("CONSEQUENCE", "")
    if not (has_consequence(csq, "inframe_insertion") or
            has_consequence(csq, "inframe_deletion")):
        return False, ""

    # 確認在 repeat region（DOMAINS 欄位）
    domains = row.get("DOMAINS", ".")
    if not domains or domains == ".":
        return False, ""

    repeat_keywords = [
        "repeat", "Repeat", "REPEAT",
        "low_complexity", "Low_complexity",
        "coiled_coil", "COIL",
        "ANK", "ARM", "WD40", "LRR", "TPR", "HEAT", "Kelch", "HOOK",
    ]
    in_repeat = any(kw in domains for kw in repeat_keywords)
    if not in_repeat:
        return False, ""

    # PhyloP 保守性判斷（v3.1 新增）
    phylop = to_float(row.get("PHYLOP100", "."))
    if phylop is None:
        # PhyloP 缺失 → 保守地不觸發（不確定保守性）
        return False, ""
    if phylop >= PHYLOP_CONSERVATION_THRESHOLD:
        # 高保守性位置的 inframe indel → 不觸發 BP3（可能有功能）
        return False, ""

    return True, f"BP3:inframe_in_repeat,phylop={phylop:.2f}<{PHYLOP_CONSERVATION_THRESHOLD}"


def check_BP7(row: dict) -> tuple[bool, str]:
    """
    BP7 ── Supporting Benign
    條件（全部滿足才觸發）：
        1. synonymous_variant
        2. Pangolin score < 0.1（低 splice 影響）
        3. PhyloP100 < 3.58（低核苷酸保守性）

    ★ v3.1 新增：PhyloP 保守性條件（參考 AutoACMG + Richards 2015 原文）
    ★ Pangolin 缺失時：保守地不觸發 BP7
    ★ PhyloP 缺失時：保守地不觸發 BP7

    回傳 (triggered: bool, reason: str)
    """
    if not has_consequence(row.get("CONSEQUENCE", ""), "synonymous_variant"):
        return False, ""

    # Pangolin splice 判斷
    pangolin = to_float(row.get("PANGOLIN_SCORE", "."))
    if pangolin is None:
        return False, ""
    if pangolin >= 0.1:
        return False, ""

    # PhyloP 保守性判斷（v3.1 新增）
    phylop = to_float(row.get("PHYLOP100", "."))
    if phylop is None:
        return False, ""
    if phylop >= PHYLOP_CONSERVATION_THRESHOLD:
        # 高度保守的同義突變 → 不觸發 BP7（可能有 RNA 層面功能）
        return False, ""

    return True, f"BP7:synonymous,pangolin={pangolin:.3f}<0.1,phylop={phylop:.2f}<{PHYLOP_CONSERVATION_THRESHOLD}"


# ══════════════════════════════════════════════════════════════════
# 第五部分：分類邏輯
# ══════════════════════════════════════════════════════════════════

def classify_variant(row: dict) -> dict:
    """
    對單一 variant 套用所有 Phase 1 criteria，
    回傳 ACMG_CRITERIA / ACMG_SCORE / ACMG_CLASS / ACMG_NOTES。

    邏輯順序：
    1. BA1（觸發則直接 Benign，跳過所有其他）
    2. BS1
    3. PVS1（觸發則停用 PP3 和 PM4）
    4. PM2_Supporting
    5. PM4（PVS1 未觸發時才判斷）
    6. PP3 / BP4（cascade，PVS1 未觸發時才判斷）
    7. BP7
    8. BP3
    9. 計算總分 → 判斷分類
    """
    triggered_criteria: list[str] = []
    notes_parts: list[str] = []

    # ── BA1（觸發則直接返回）────────────────────────────────────────
    if check_BA1(row):
        af = max_af(row)
        return {
            "ACMG_CRITERIA": "BA1",
            "ACMG_SCORE":    str(CRITERIA_POINTS["BA1"]),
            "ACMG_CLASS":    "Benign",
            "ACMG_NOTES":    f"BA1:gnomAD_AF={af:.4f}>5%",
        }

    # ── BS1 ──────────────────────────────────────────────────────────
    if check_BS1(row):
        triggered_criteria.append("BS1")
        af = max_af(row)
        notes_parts.append(f"BS1:gnomAD_AF={af:.4f}>1%")

    # ── PVS1 ─────────────────────────────────────────────────────────
    pvs1_triggered, pvs1_note = check_PVS1(row)
    if pvs1_triggered:
        triggered_criteria.append("PVS1")
        notes_parts.append(pvs1_note)
    elif pvs1_note:
        # 未觸發但有值得關注的資訊（HI=2 或不在 ClinGen）→ 記錄
        notes_parts.append(pvs1_note)

    # ── PM2_Supporting ────────────────────────────────────────────────
    pm2_triggered, pm2_note = check_PM2(row)
    if pm2_triggered:
        triggered_criteria.append("PM2_Supporting")
        notes_parts.append(pm2_note)

    # ── PM4（PVS1 觸發時停用，避免 double counting）─────────────────
    # AutoACMG 說明：PVS1 已捕捉 LoF 效應，PM4 不應再加分
    if not pvs1_triggered:
        pm4_triggered, pm4_note = check_PM4(row)
        if pm4_triggered:
            triggered_criteria.append("PM4")
            notes_parts.append(pm4_note)

    # ── PP3 / BP4（PVS1 觸發時停用）────────────────────────────────
    # AutoACMG 說明：PVS1 代表 LoF，計算工具通常也會預測為致病，
    # 但這是同一個機制的不同證據，不應重複計分
    if not pvs1_triggered:
        pp3_or_bp4, tool_used = check_PP3(row)
        if pp3_or_bp4 is not None:
            triggered_criteria.append(pp3_or_bp4)
            notes_parts.append(f"{pp3_or_bp4}:tool={tool_used}")

    # ── BP7 ──────────────────────────────────────────────────────────
    bp7_triggered, bp7_note = check_BP7(row)
    if bp7_triggered:
        triggered_criteria.append("BP7")
        notes_parts.append(bp7_note)

    # ── BP3 ──────────────────────────────────────────────────────────
    bp3_triggered, bp3_note = check_BP3(row)
    if bp3_triggered:
        triggered_criteria.append("BP3")
        notes_parts.append(bp3_note)

    # ── 計算總分 ──────────────────────────────────────────────────────
    total_score = sum(CRITERIA_POINTS.get(c, 0) for c in triggered_criteria)

    # ── 判斷分類 ──────────────────────────────────────────────────────
    if total_score >= PATHOGENIC_THRESHOLD:
        acmg_class = "Pathogenic"
    elif total_score >= LIKELY_PATHOGENIC_MIN:
        acmg_class = "Likely_Pathogenic"
    elif total_score <= BENIGN_THRESHOLD:
        acmg_class = "Benign"
    elif total_score <= LIKELY_BENIGN_MAX:
        acmg_class = "Likely_Benign"
    else:
        acmg_class = "VUS"

    return {
        "ACMG_CRITERIA": ",".join(triggered_criteria) if triggered_criteria else ".",
        "ACMG_SCORE":    str(total_score),
        "ACMG_CLASS":    acmg_class,
        "ACMG_NOTES":    "|".join(notes_parts) if notes_parts else ".",
    }


# ══════════════════════════════════════════════════════════════════
# 第六部分：TSV 讀寫主流程
# ══════════════════════════════════════════════════════════════════

def validate_header(header_cols: list[str]) -> dict[str, int]:
    """驗證 TSV header，回傳欄位名 → index 的 mapping。"""
    col_index = {name: idx for idx, name in enumerate(header_cols)}
    required = [
        "GENE", "CONSEQUENCE", "IMPACT",
        "GNOMAD_G_AF", "GNOMAD_G_EAS_AF",
        "GNOMAD_E_AF", "GNOMAD_E_EAS_AF",
        "GNOMAD_E_EAS_AF_DBNSFP", "TG_EAS_AF",
        "LOFTEE", "LOFTEE_FILTER",
        "PKNN_EVIDENCE", "PKNN_LLR",
        "ALPHAMISSENSE", "ESM1B", "VARITY_R", "BAYESDEL_NOAF",
        "PANGOLIN_SCORE", "PHYLOP100",
        "DOMAINS",
    ]
    missing = [c for c in required if c not in col_index]
    if missing:
        print(f"[ERROR] TSV 缺少欄位：{missing}", file=sys.stderr)
        sys.exit(1)
    return col_index


def process_tsv(
    input_path: str,
    output_path: str,
    clingen_hi_path: str | None = None,
    gene_moi_path:   str | None = None,
):
    """主要處理流程：載入資料庫 → 逐行分類 → 寫出結果。"""
    global _HI_TABLE, _MOI_TABLE

    _HI_TABLE  = load_clingen_hi(clingen_hi_path)
    _MOI_TABLE = load_gene_moi(gene_moi_path)

    opener_in = gzip.open if input_path.endswith(".gz") else open

    print(f"[acmg_classifier] 開始處理：{input_path}", file=sys.stderr)

    n_total = n_p = n_lp = n_vus = n_lb = n_b = 0
    criteria_counts: dict[str, int] = {}

    with opener_in(input_path, "rt") as fin, open(output_path, "w") as fout:

        header_cols = fin.readline().rstrip("\n").split("\t")
        col_index   = validate_header(header_cols)
        print(f"  Header 驗證通過，共 {len(header_cols)} 欄", file=sys.stderr)

        fout.write("\t".join(header_cols + NEW_COLUMNS) + "\n")

        for line in fin:
            line = line.rstrip("\n")
            if not line:
                continue

            n_total += 1
            if n_total % 5000 == 0:
                print(f"  已處理 {n_total:,} 行...", file=sys.stderr)

            fields = line.split("\t")
            row = {col: (fields[idx] if idx < len(fields) else ".") for col, idx in col_index.items()}

            result = classify_variant(row)

            cls = result["ACMG_CLASS"]
            if   cls == "Pathogenic":        n_p  += 1
            elif cls == "Likely_Pathogenic":  n_lp += 1
            elif cls == "VUS":                n_vus+= 1
            elif cls == "Likely_Benign":      n_lb += 1
            elif cls == "Benign":             n_b  += 1

            if result["ACMG_CRITERIA"] != ".":
                for crit in result["ACMG_CRITERIA"].split(","):
                    criteria_counts[crit] = criteria_counts.get(crit, 0) + 1

            fout.write("\t".join(fields + [
                result["ACMG_CRITERIA"],
                result["ACMG_SCORE"],
                result["ACMG_CLASS"],
                result["ACMG_NOTES"],
            ]) + "\n")

    print(f"\n[acmg_classifier] 完成，共 {n_total:,} 個 variant", file=sys.stderr)
    print(f"  Pathogenic        : {n_p:>6,} ({n_p/n_total*100:.1f}%)", file=sys.stderr)
    print(f"  Likely_Pathogenic : {n_lp:>6,} ({n_lp/n_total*100:.1f}%)", file=sys.stderr)
    print(f"  VUS               : {n_vus:>6,} ({n_vus/n_total*100:.1f}%)", file=sys.stderr)
    print(f"  Likely_Benign     : {n_lb:>6,} ({n_lb/n_total*100:.1f}%)", file=sys.stderr)
    print(f"  Benign            : {n_b:>6,} ({n_b/n_total*100:.1f}%)", file=sys.stderr)
    print(f"\n  ── Criteria 觸發次數 ──", file=sys.stderr)
    for crit, cnt in sorted(criteria_counts.items(), key=lambda x: -x[1]):
        print(f"  {crit:<25} : {cnt:>6,}", file=sys.stderr)


# ══════════════════════════════════════════════════════════════════
# 第七部分：命令列介面
# ══════════════════════════════════════════════════════════════════

def parse_args():
    parser = argparse.ArgumentParser(
        description="ACMG/AMP variant classification v3.1",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用範例：
  python3 acmg_classifier.py \\
      --input      NA12878_WES.snv_indel.full.annotated.tsv \\
      --clingen_hi /path/to/clingen/ClinGen_gene_curation_list_GRCh38.tsv \\
      --gene_moi   /path/to/clingen/gene_moi.tsv.gz \\
      --output     NA12878_WES.snv_indel.acmg.tsv

  # 所有資料庫都是 optional（不提供時退回簡化版）
  # 快速驗證：
  cut -f57-60 NA12878_WES.snv_indel.acmg.tsv | head -20
        """
    )
    parser.add_argument("--input",  "-i", required=True,
                        help="輸入 TSV（支援 .gz）")
    parser.add_argument("--output", "-o", required=True,
                        help="輸出 TSV（新增 4 個 ACMG 欄位）")
    parser.add_argument("--clingen_hi", default=None,
                        help="ClinGen_gene_curation_list_GRCh38.tsv（PVS1 HI 過濾）")
    parser.add_argument("--gene_moi",   default=None,
                        help="gene_moi.tsv.gz（PM2 MOI 判斷，由 build_gene_moi.py 產生）")
    return parser.parse_args()


def main():
    args = parse_args()
    process_tsv(args.input, args.output,
                clingen_hi_path=args.clingen_hi,
                gene_moi_path=args.gene_moi)


if __name__ == "__main__":
    main()