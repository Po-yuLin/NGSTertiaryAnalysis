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
 * Licensed under the GNU General Public License v3.0
 *
 * DISCLAIMER: This pipeline is provided "as is" without
 * warranty of any kind. The authors and their institution
 * make no representations or warranties regarding the
 * accuracy, completeness, or suitability of the analysis
 * results for any clinical or research purpose. Users are
 * solely responsible for validating and interpreting all
 * results.
 * =========================================================

parse_str_vcf.py
================
目的：
    讀取 STR VCF（GangSTR 或 DRAGEN ExpansionHunter），
    對照 STRchive 查表，標記每個 locus 的 repeat count 是否超過 threshold，
    輸出臨床用 STR TSV。

    只輸出有對應到 STRchive 記錄的 locus（已知致病 STR）。
    其餘非疾病相關的 STR 略過不輸出。

支援兩種 pipeline：
    nckuh（GangSTR）：
        - 用位置查表（str_lookup_pos.tsv.gz）對應 STRchive
        - FORMAT/REPCN：Number=2，逗號分隔，例如 "19,15"
        - FORMAT/DP：read depth
        - FORMAT/Q：call quality（posterior probability）
        - FORMAT/REPCI：confidence interval，例如 "19-19,15-15"

    dragen（ExpansionHunter）：
        - 用 VARID 查表（str_lookup_varid.tsv.gz）對應 STRchive
        - INFO/VARID：locus ID，例如 "HTT"
        - FORMAT/REPCN：逗號分隔，例如 "19,15"
        - FORMAT/LC：locus coverage（相當於 DP）
        - FORMAT/REPCI：confidence interval

輸入：
    --vcf       STR VCF（.vcf 或 .vcf.gz）
    --sample    sample ID
    --lookup    查表檔路徑（nckuh 用 pos，dragen 用 varid）
    --pipeline  nckuh 或 dragen
    --output    輸出 TSV 路徑

輸出欄位（共 20 欄）：
    # 位置與 locus 資訊
    CHROM, POS, END, STR_ID, GENE, MOTIF, LOCUS_STRUCTURE, TYPE

    # 樣本資訊
    REPCN_A1, REPCN_A2（兩個 allele 的 repeat count）
    DP（read depth 或 locus coverage）
    REPCI（confidence interval）

    # STRchive threshold
    BENIGN_MAX, PATHOGENIC_MIN, INTERMEDIATE_MIN, INTERMEDIATE_MAX

    # 分類結果（per allele 判斷後取最嚴重）
    CLASSIFICATION（normal / intermediate / pathogenic / no_threshold）

    # 疾病資訊
    DISEASE, INHERITANCE, PIPELINE
"""

import argparse
import csv
import gzip
import sys

try:
    from cyvcf2 import VCF
except ImportError:
    print("[ERROR] 請先安裝 cyvcf2：pip install cyvcf2", file=sys.stderr)
    sys.exit(1)


# ──────────────────────────────────────────────────────────────
# 載入 STRchive 查表
# ──────────────────────────────────────────────────────────────

def load_lookup(lookup_path: str) -> dict:
    """
    讀取 str_lookup_varid.tsv.gz 或 str_lookup_pos.tsv.gz，
    建立以 KEY 為 index 的 dict。

    同一個 KEY 理論上只有一筆（STRchive locus 唯一），
    若有重複（不同 disease 同一位置），保留第一筆。

    回傳：{ key_str: row_dict, ... }
    """
    lookup = {}
    with gzip.open(lookup_path, "rt", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            key = row["KEY"]
            if key not in lookup:
                lookup[key] = dict(row)

    print(f"[parse_str_vcf] 查表載入：{len(lookup)} 個 key（{lookup_path}）",
          file=sys.stderr)
    return lookup


# ──────────────────────────────────────────────────────────────
# 解析 REPCN
# ──────────────────────────────────────────────────────────────

def parse_repcn(repcn_str: str) -> tuple:
    """
    解析 FORMAT/REPCN 欄位，回傳 (allele1, allele2) 的整數 tuple。

    兩種分隔符格式：
        GangSTR（NCKUH）：逗號分隔，例如 "19,15"
        DRAGEN（ExpansionHunter）：斜線分隔，例如 "17/22"

    特殊情況：
        "40,."   → (40, None)（hemizygous 或 missing）
        "."      → (None, None)（完全 missing）

    回傳 None 表示該 allele 無法解析。
    """
    if not repcn_str or repcn_str in (".", ""):
        return None, None

    # 自動偵測分隔符：有 "/" 用斜線，否則用逗號
    if "/" in repcn_str:
        parts = repcn_str.split("/")
    else:
        parts = repcn_str.split(",")

    def to_int(s):
        s = s.strip()
        if s in (".", ""):
            return None
        try:
            return int(float(s))  # float() 先處理 "19.0" 這種格式
        except ValueError:
            return None

    a1 = to_int(parts[0]) if len(parts) > 0 else None
    a2 = to_int(parts[1]) if len(parts) > 1 else None
    return a1, a2


# ──────────────────────────────────────────────────────────────
# 分類邏輯
# ──────────────────────────────────────────────────────────────

def classify_allele(repcn: int, benign_max: float, pathogenic_min: float,
                    intermediate_min, intermediate_max) -> str:
    """
    依照 STRchive threshold 對單一 allele 的 repeat count 分類。

    分類規則：
        pathogenic    : repcn >= pathogenic_min
        intermediate  : intermediate_min <= repcn <= intermediate_max（若有設定）
        normal        : repcn <= benign_max
        no_threshold  : threshold 不完整，無法分類
        borderline    : 介於 benign_max 和 pathogenic_min 之間，無 intermediate 定義

    回傳分類字串。
    """
    if repcn is None:
        return "."

    # threshold 不完整
    if benign_max is None or pathogenic_min is None:
        return "no_threshold"

    # pathogenic
    if repcn >= pathogenic_min:
        return "pathogenic"

    # intermediate（有定義才判斷）
    if intermediate_min is not None and intermediate_max is not None:
        if intermediate_min <= repcn <= intermediate_max:
            return "intermediate"

    # normal
    if repcn <= benign_max:
        return "normal"

    # 介於 benign_max 和 pathogenic_min 之間，沒有 intermediate 定義
    return "borderline"


def classify_locus(a1: int, a2: int, rec: dict) -> str:
    """
    對一個 locus 的兩個 allele 分別分類，取最嚴重的結果。

    嚴重程度順序：pathogenic > intermediate > borderline > normal > no_threshold > .

    特殊處理：缺失型 locus（如 VWA1）
        STRchive 的 benign_min 和 pathogenic_max 欄位用來識別這類 locus。
        VWA1 的規則：benign_min == benign_max（只有一個正常值），
        repeat count 小於 pathogenic_min 或大於 pathogenic_max 都是 pathogenic。
        直接用 benign_min/benign_max 區間判斷：在區間內才是 normal。
    """
    severity = {
        "pathogenic":   5,
        "intermediate": 4,
        "borderline":   3,
        "normal":       2,
        "no_threshold": 1,
        ".":            0,
    }

    # 從查表取 threshold（字串轉浮點數）
    def to_float(s):
        if s in (".", "", None):
            return None
        try:
            return float(s)
        except ValueError:
            return None

    benign_min      = to_float(rec.get("BENIGN_MIN"))
    benign_max      = to_float(rec.get("BENIGN_MAX"))
    pathogenic_min  = to_float(rec.get("PATHOGENIC_MIN"))
    pathogenic_max  = to_float(rec.get("PATHOGENIC_MAX"))
    intermediate_min = to_float(rec.get("INTERMEDIATE_MIN"))
    intermediate_max = to_float(rec.get("INTERMEDIATE_MAX"))

    # 偵測缺失型 locus：
    #   benign_min 有值，且 pathogenic_min < benign_min（repeat 太少也是致病）
    #   例如 VWA1：benign_min=2, benign_max=2, pathogenic_min=1, pathogenic_max=3
    is_deletion_type = (
        benign_min is not None and
        pathogenic_min is not None and
        pathogenic_min < benign_min
    )

    def classify_one(repcn):
        if repcn is None:
            return "."
        if is_deletion_type:
            # 缺失型：只有在 [benign_min, benign_max] 區間才是 normal
            if benign_min is None or benign_max is None:
                return "no_threshold"
            if benign_min <= repcn <= benign_max:
                return "normal"
            # 在 benign 區間外，用 pathogenic_max 和 pathogenic_min 進一步判斷
            if pathogenic_max is not None and repcn > pathogenic_max:
                return "pathogenic"
            if pathogenic_min is not None and repcn < pathogenic_min:
                return "pathogenic"
            # 介於 pathogenic 和 benign 之間
            return "pathogenic"
        else:
            # 擴張型（一般情況）
            return classify_allele(
                repcn, benign_max, pathogenic_min,
                intermediate_min, intermediate_max
            )

    cls_a1 = classify_one(a1)
    cls_a2 = classify_one(a2)

    # 取嚴重程度較高的
    if severity.get(cls_a1, 0) >= severity.get(cls_a2, 0):
        return cls_a1
    return cls_a2


# ──────────────────────────────────────────────────────────────
# 讀取 FORMAT 欄位
# ──────────────────────────────────────────────────────────────

def get_repcn_str_nckuh(variant, sample_idx: int) -> str:
    """
    GangSTR（NCKUH）專用：讀取 REPCN。
    GangSTR 的 REPCN 是 Number=2 Integer 型。
    cyvcf2 回傳 numpy int32 array，例如 array([16, 16])。
    回傳逗號分隔字串，例如 "16,16"。
    """
    try:
        repcn = variant.format("REPCN")
        if repcn is None:
            return "."
        val = repcn[sample_idx]  # numpy int32 array，shape=(2,)
        INT_MISSING = -2147483648
        parts = []
        for v in val:
            iv = int(v)
            parts.append("." if iv == INT_MISSING else str(iv))
        return ",".join(parts)
    except Exception:
        return "."


def get_repcn_str_dragen(variant, sample_idx: int) -> str:
    """
    DRAGEN（ExpansionHunter）專用：讀取 REPCN。
    DRAGEN 的 REPCN 是 Number=1 String 型。
    cyvcf2 回傳 numpy Unicode string array，dtype='U5' 之類。
    val 是 numpy.str_，直接轉 str 即可，例如 "17/22"。
    回傳原始字串，parse_repcn() 負責用 "/" 分隔。
    """
    try:
        repcn = variant.format("REPCN")
        if repcn is None:
            return "."
        val = repcn[sample_idx]  # numpy.str_，例如 "17/22"
        s = str(val).strip()
        return s if s else "."
    except Exception:
        return "."


def get_repcn_str(variant, sample_idx: int, pipeline: str) -> str:
    """
    依 pipeline 選擇對應的 REPCN 讀取函式。
    """
    if pipeline == "nckuh":
        return get_repcn_str_nckuh(variant, sample_idx)
    else:
        return get_repcn_str_dragen(variant, sample_idx)


def get_dp(variant, sample_idx: int, pipeline: str) -> str:
    """
    取得 read depth。
    GangSTR：FORMAT/DP（Integer）
    DRAGEN：FORMAT/LC（Float，locus coverage）

    cyvcf2 對 Float 型 FORMAT 回傳 numpy float32 array，
    需要用 [sample_idx][0] 取值，並處理 NaN。
    """
    field = "DP" if pipeline == "nckuh" else "LC"
    try:
        import math
        arr = variant.format(field)
        if arr is None:
            return "."
        val = arr[sample_idx]
        # cyvcf2 回傳 numpy array（shape: [1]），取第一個元素
        if hasattr(val, '__len__'):
            val = val[0]
        # 處理 NaN 和負數（cyvcf2 missing 值）
        if val is None:
            return "."
        try:
            f = float(val)
        except (TypeError, ValueError):
            return "."
        if math.isnan(f) or f < 0:
            return "."
        return str(int(f))
    except Exception:
        return "."


def get_repci(variant, sample_idx: int) -> str:
    """從 FORMAT/REPCI 取得 confidence interval 字串"""
    try:
        repci = variant.format("REPCI")
        if repci is None:
            return "."
        val = repci[sample_idx]
        if isinstance(val, bytes):
            val = val.decode("utf-8")
        if isinstance(val, (list, tuple)):
            return ",".join(str(v) for v in val)
        return str(val).strip()
    except Exception:
        return "."


# ──────────────────────────────────────────────────────────────
# 主流程
# ──────────────────────────────────────────────────────────────

def parse_str_vcf(vcf_path: str, sample_id: str, lookup_path: str,
                  pipeline: str, output_path: str):
    """
    主要處理流程：
    1. 載入 STRchive 查表
    2. 逐 variant 查表，只輸出有對應 STRchive 記錄的 locus
    3. 解析 REPCN，判斷分類，輸出 TSV
    """

    # 載入查表
    lookup = load_lookup(lookup_path)

    # 開啟 VCF
    vcf_in = VCF(vcf_path)
    samples = vcf_in.samples
    print(f"[parse_str_vcf] VCF sample columns：{samples}", file=sys.stderr)

    if sample_id not in samples:
        print(f"[ERROR] 找不到 sample：{sample_id}", file=sys.stderr)
        sys.exit(1)
    sample_idx = samples.index(sample_id)

    # 輸出欄位
    output_cols = [
        "CHROM", "POS", "END",
        "STR_ID", "GENE", "MOTIF", "LOCUS_STRUCTURE", "TYPE",
        "REPCN_A1", "REPCN_A2", "DP", "REPCI",
        "BENIGN_MIN", "BENIGN_MAX",
        "PATHOGENIC_MIN", "PATHOGENIC_MAX",
        "INTERMEDIATE_MIN", "INTERMEDIATE_MAX",
        "CLASSIFICATION",
        "DISEASE", "INHERITANCE", "PIPELINE",
    ]

    n_total   = 0   # VCF 裡的總 variant 數
    n_matched = 0   # 有對應到 STRchive 的 locus 數
    n_pathogenic    = 0
    n_intermediate  = 0

    with open(output_path, "w", encoding="utf-8") as fout:
        fout.write("\t".join(output_cols) + "\n")

        for variant in vcf_in:
            n_total += 1

            chrom = variant.CHROM
            pos   = variant.POS

            # ── 查 STRchive 查表 ──────────────────────────────
            if pipeline == "dragen":
                # DRAGEN：用 INFO/VARID 查表
                varid = variant.INFO.get("VARID", None)
                if not varid:
                    continue
                rec = lookup.get(str(varid))
            else:
                # NCKUH GangSTR：用 chrom:pos 查表
                key = f"{chrom}:{pos}"
                rec = lookup.get(key)

            # 沒有對應的 STRchive 記錄，略過
            if rec is None:
                continue

            n_matched += 1

            # ── 取 END 座標 ───────────────────────────────────
            end = variant.INFO.get("END", pos)

            # ── 解析 REPCN ────────────────────────────────────
            repcn_str = get_repcn_str(variant, sample_idx, pipeline)
            a1, a2 = parse_repcn(repcn_str)

            # ── DP 和 REPCI ───────────────────────────────────
            dp_str    = get_dp(variant, sample_idx, pipeline)
            repci_str = get_repci(variant, sample_idx)

            # ── 分類 ─────────────────────────────────────────
            classification = classify_locus(a1, a2, rec)

            if classification == "pathogenic":
                n_pathogenic += 1
            elif classification == "intermediate":
                n_intermediate += 1

            # ── 寫出一行 ──────────────────────────────────────
            row = [
                chrom,
                str(pos),
                str(end),
                rec.get("STR_ID",          "."),
                rec.get("GENE",            "."),
                rec.get("MOTIF",           "."),
                rec.get("LOCUS_STRUCTURE", "."),
                rec.get("TYPE",            "."),
                str(a1) if a1 is not None else ".",
                str(a2) if a2 is not None else ".",
                dp_str,
                repci_str,
                rec.get("BENIGN_MIN",        "."),
                rec.get("BENIGN_MAX",        "."),
                rec.get("PATHOGENIC_MIN",    "."),
                rec.get("PATHOGENIC_MAX",    "."),
                rec.get("INTERMEDIATE_MIN",  "."),
                rec.get("INTERMEDIATE_MAX",  "."),
                classification,
                rec.get("DISEASE",    "."),
                rec.get("INHERITANCE","."),
                pipeline,
            ]
            fout.write("\t".join(row) + "\n")

    vcf_in.close()

    # 統計摘要
    print(f"[parse_str_vcf] 完成", file=sys.stderr)
    print(f"  VCF 總 locus 數          : {n_total:>6,}", file=sys.stderr)
    print(f"  STRchive 命中（輸出）    : {n_matched:>6,}", file=sys.stderr)
    print(f"  Pathogenic               : {n_pathogenic:>6,}", file=sys.stderr)
    print(f"  Intermediate             : {n_intermediate:>6,}", file=sys.stderr)


# ──────────────────────────────────────────────────────────────
# 命令列介面
# ──────────────────────────────────────────────────────────────

def parse_args():
    parser = argparse.ArgumentParser(
        description="解析 STR VCF，對照 STRchive threshold，輸出分類 TSV",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用範例：
  # NCKUH GangSTR
  python3 parse_str_vcf.py \\
      --vcf      NA12878_WES.str.vcf.gz \\
      --sample   NA12878_WES \\
      --lookup   str_lookup_pos.tsv.gz \\
      --pipeline nckuh \\
      --output   NA12878_WES.str.tsv

  # DRAGEN ExpansionHunter
  python3 parse_str_vcf.py \\
      --vcf      VAL-10.repeats.vcf.gz \\
      --sample   VAL-10 \\
      --lookup   str_lookup_varid.tsv.gz \\
      --pipeline dragen \\
      --output   VAL-10.str.tsv
        """
    )
    parser.add_argument("--vcf",      required=True,
                        help="STR VCF 路徑（.vcf 或 .vcf.gz）")
    parser.add_argument("--sample",   required=True,
                        help="Sample ID")
    parser.add_argument("--lookup",   required=True,
                        help="STRchive 查表（pos 或 varid）")
    parser.add_argument("--pipeline", required=True,
                        choices=["nckuh", "dragen"],
                        help="pipeline 類型")
    parser.add_argument("--output",   required=True,
                        help="輸出 TSV 路徑")
    return parser.parse_args()


def main():
    args = parse_args()
    print(f"[parse_str_vcf] VCF      : {args.vcf}",      file=sys.stderr)
    print(f"[parse_str_vcf] Sample   : {args.sample}",   file=sys.stderr)
    print(f"[parse_str_vcf] Lookup   : {args.lookup}",   file=sys.stderr)
    print(f"[parse_str_vcf] Pipeline : {args.pipeline}", file=sys.stderr)
    print(f"[parse_str_vcf] Output   : {args.output}",   file=sys.stderr)
    parse_str_vcf(args.vcf, args.sample, args.lookup, args.pipeline, args.output)


if __name__ == "__main__":
    main()
