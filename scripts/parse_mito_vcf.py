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

parse_mito_vcf.py
=================
目的：
    讀取 VEP 註解後的 mito VCF，解析 CSQ 欄位，
    查詢 MITOMAP lookup 表，輸出 mito TSV。

    同時支援 NCKUH（GATK Mutect2）和 DRAGEN 的 mito VCF，
    兩者都是單一 sample column，差別只在 INFO tag 名稱：
        NCKUH：不含 CALLERS（直接讀 FORMAT/AF 和 FORMAT/DP）
        DRAGEN：含 CALLERS=DRAGEN、DP_DRAGEN、AD_DRAGEN、VAF_DRAGEN

輸入：
    --vcf        VEP 註解後的 mito VCF（.vcf.gz）
    --sample     sample ID（用來確認 sample column 名稱）
    --mitomap    mitomap_lookup.tsv.gz（build_mitomap_lookup.py 產生）
    --pipeline   nckuh 或 dragen（決定如何讀取 DP/AF）
    --output     輸出 TSV 路徑（.tsv）

輸出欄位（共 22 欄）：
    # 位置
    CHROM, POS, REF, ALT

    # Transcript（來自 VEP CSQ）
    GENE, HGVS_C, HGVS_P, CONSEQUENCE, IMPACT, BIOTYPE

    # 樣本資訊
    GENOTYPE, DP, AF_SAMPLE（heteroplasmy level）

    # gnomAD mito（來自 VEP CSQ gnomADmt）
    GNOMAD_MITO_AF, GNOMAD_MITO_AF_HOM, GNOMAD_MITO_AF_HET

    # ClinVar（來自 VEP CSQ --custom）
    CLINVAR_SIG, CLINVAR_DN

    # MITOMAP（來自 mitomap_lookup.tsv.gz）
    MITOMAP_LOCUS, MITOMAP_DISEASE, MITOMAP_STATUS,
    MITOMAP_HOMO, MITOMAP_HETERO, MITOMAP_MITOTIP
"""

import argparse
import csv
import gzip
import re
import sys

try:
    from cyvcf2 import VCF
except ImportError:
    print("[ERROR] 請先安裝 cyvcf2：pip install cyvcf2", file=sys.stderr)
    sys.exit(1)


# ──────────────────────────────────────────────────────────────
# 載入 MITOMAP lookup 表
# ──────────────────────────────────────────────────────────────

def load_mitomap(mitomap_path: str) -> dict:
    """
    讀取 mitomap_lookup.tsv.gz，建立以 POS 為 key 的 dict。
    同一個 POS 可能有多筆（不同 ALT allele），所以 value 是 list。

    回傳格式：
        { "3243": [ {row_dict}, {row_dict}, ... ], ... }
    """
    lookup = {}  # key=POS字串, value=list of dict

    with gzip.open(mitomap_path, "rt", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            pos = row["POS"].strip()
            if pos not in lookup:
                lookup[pos] = []
            lookup[pos].append(dict(row))

    print(f"[parse_mito_vcf] MITOMAP lookup 載入：{len(lookup)} 個唯一 POS",
          file=sys.stderr)
    return lookup


def query_mitomap(lookup: dict, pos: int, ref: str, alt: str) -> dict:
    """
    用 POS + REF/ALT 查詢 MITOMAP lookup。

    精確匹配邏輯：
        1. 先依 POS 找候選筆
        2. 在候選筆中，比對 NUC_CHANGE（例如 "T-C"）：
           - REF 對應 NUC_CHANGE 的第一個字母
           - ALT 對應 NUC_CHANGE 的第二個字母
        3. 若有精確匹配，回傳第一筆
        4. 若無精確匹配但 POS 有記錄，回傳 POS 層級的第一筆
           （至少知道這個位置有 MITOMAP 記錄，disease 欄位仍有參考價值）
        5. 完全沒有 POS，回傳空 dict

    回傳 dict（key 是 TSV 欄位名稱），無 match 回傳 None。
    """
    pos_str = str(pos)
    candidates = lookup.get(pos_str, [])

    if not candidates:
        return None  # 這個位置 MITOMAP 沒有記錄

    # 嘗試精確匹配（REF-ALT 對應 NUC_CHANGE）
    for rec in candidates:
        nuc = rec.get("NUC_CHANGE", ".")
        # 格式：X-Y，X=REF, Y=ALT
        if "-" in nuc:
            parts = nuc.split("-")
            if len(parts) == 2:
                mitomap_ref = parts[0].strip().upper()
                mitomap_alt = parts[1].strip().upper()
                if mitomap_ref == ref.upper() and mitomap_alt == alt.upper():
                    return rec  # 精確匹配

        # RNA 的 NUC_CHANGE 是 "."，用 ALLELE 欄位的 REF>ALT 符號比對
        allele = rec.get("ALLELE", "")
        # 格式：m.3243A>G，提取 A 和 G
        m = re.search(r'(\d+)([ACGT])>([ACGT])', allele, re.IGNORECASE)
        if m:
            a_ref = m.group(2).upper()
            a_alt = m.group(3).upper()
            if a_ref == ref.upper() and a_alt == alt.upper():
                return rec  # 從 ALLELE 欄位精確匹配

    # 精確匹配失敗：這個 variant 在 MITOMAP 沒有對應記錄
    # 不做 fallback，避免回傳不相關的同位置其他 variant 造成誤導
    return None


# ──────────────────────────────────────────────────────────────
# 解析 VEP CSQ 欄位
# ──────────────────────────────────────────────────────────────

def parse_csq_header(vcf_obj) -> list:
    """
    從 VCF header 解析 CSQ 欄位的 Format 定義，
    回傳欄位名稱的 list（對應 pipe-separated 的順序）。

    VEP 在 header 的格式：
        ##INFO=<ID=CSQ,...,Description="...Format: A|B|C|...">
    """
    for header_line in vcf_obj.raw_header.split("\n"):
        if "ID=CSQ" in header_line and "Format:" in header_line:
            # 提取 Format: 後面的欄位串
            m = re.search(r'Format: ([^"]+)', header_line)
            if m:
                fields = m.group(1).strip().split("|")
                return fields
    print("[parse_mito_vcf] [WARN] 找不到 CSQ Format 定義", file=sys.stderr)
    return []


def pick_transcript(csq_list: list, csq_fields: list) -> dict:
    """
    從多個 transcript annotation 中選出代表 transcript。

    選取優先順序（與 SNV pipeline 的 parse_vep_csq.py 相同）：
        1. PICK=1（VEP flag_pick 標記的代表 transcript）
        2. 若無 PICK=1，取第一個
    """
    if not csq_list or not csq_fields:
        return {}

    # 建立每個 transcript 的 dict
    transcripts = []
    for csq_str in csq_list:
        vals = csq_str.split("|")
        # 補齊長度（有些欄位可能是空的）
        while len(vals) < len(csq_fields):
            vals.append("")
        tx = dict(zip(csq_fields, vals))
        transcripts.append(tx)

    # 找 PICK=1 的 transcript
    for tx in transcripts:
        if tx.get("PICK", "") == "1":
            return tx

    # 沒有 PICK=1，回傳第一個
    return transcripts[0]


def get_gnomad_mito(tx: dict) -> tuple:
    """
    從 VEP CSQ 中提取 gnomAD mito 頻率。

    VEP 115 --af_gnomadg 輸出的欄位名稱是 gnomADg_AF（整個 genome 資料庫，
    包含 chrM variant）。VEP 115 不提供獨立的 homoplasmy/heteroplasmy 欄位，
    那些欄位只存在於直接解析 gnomAD mito VCF 時才有。

    回傳 (gnomADg_AF, ".", ".")：
        第一個是整體 AF，後兩個保留為 "." 以維持 TSV 欄位結構一致。
    """
    # VEP 115 的 gnomAD genome AF 欄位名稱
    af = tx.get("gnomADg_AF", "") or "."
    return af, ".", "."


def get_clinvar(tx: dict) -> tuple:
    """
    從 VEP CSQ 中提取 ClinVar custom annotation。

    VEP --custom 在 CSQ 欄位的名稱格式：
        ClinVar_CLNSIG    → 致病性
        ClinVar_CLNDN     → 疾病名稱
    """
    clnsig = tx.get("ClinVar_CLNSIG", ".") or "."
    clndn  = tx.get("ClinVar_CLNDN",  ".") or "."
    return clnsig, clndn


# ──────────────────────────────────────────────────────────────
# 讀取樣本 DP / AF（heteroplasmy level）
# ──────────────────────────────────────────────────────────────

def get_sample_metrics(variant, sample_idx: int, pipeline: str) -> tuple:
    """
    從 FORMAT 欄位讀取 GT、DP、AF（heteroplasmy level）。

    NCKUH（GATK Mutect2 mito）和 DRAGEN 都有 FORMAT/DP 和 FORMAT/AF，
    所以 pipeline 參數目前只影響 CALLERS tag 的處理，
    FORMAT 讀取邏輯相同。

    回傳 (genotype_str, dp_str, af_str)
    """
    # GT
    try:
        gt_tuple = variant.genotypes[sample_idx]
        # cyvcf2 格式：[allele1, allele2, phased]
        a1, a2 = gt_tuple[0], gt_tuple[1]
        if a1 < 0 or a2 < 0:
            gt_str = "./."
        else:
            sep = "|" if gt_tuple[2] else "/"
            gt_str = f"{a1}{sep}{a2}"
    except Exception:
        gt_str = "./."

    # DP
    try:
        dp_arr = variant.format("DP")
        dp_val = dp_arr[sample_idx][0] if dp_arr is not None else -1
        dp_str = str(dp_val) if dp_val >= 0 else "."
    except Exception:
        dp_str = "."

    # AF（heteroplasmy level）
    # GATK Mutect2 mito 用 FORMAT/AF（Number=A，per ALT）
    try:
        af_arr = variant.format("AF")
        af_val = af_arr[sample_idx][0] if af_arr is not None else -1.0
        if af_val >= 0:
            af_str = f"{af_val:.4f}"
        else:
            af_str = "."
    except Exception:
        af_str = "."

    return gt_str, dp_str, af_str


# ──────────────────────────────────────────────────────────────
# 主流程
# ──────────────────────────────────────────────────────────────

def parse_mito_vcf(vcf_path: str, sample_id: str, mitomap_path: str,
                   pipeline: str, output_path: str):
    """
    主要處理流程：
    1. 載入 MITOMAP lookup
    2. 開啟 VEP 註解後的 mito VCF
    3. 逐 variant 解析 CSQ + FORMAT，查 MITOMAP
    4. 輸出 TSV
    """

    # ── 載入 MITOMAP ──────────────────────────────────────────
    mitomap_lookup = load_mitomap(mitomap_path)

    # ── 開啟 VCF ─────────────────────────────────────────────
    vcf_in = VCF(vcf_path)

    # 確認 sample column
    samples = vcf_in.samples
    print(f"[parse_mito_vcf] VCF sample columns：{samples}", file=sys.stderr)

    if sample_id not in samples:
        print(f"[ERROR] 找不到 sample：{sample_id}", file=sys.stderr)
        sys.exit(1)
    sample_idx = samples.index(sample_id)

    # 解析 CSQ header（動態取得欄位順序）
    csq_fields = parse_csq_header(vcf_in)
    if not csq_fields:
        print("[ERROR] 無法解析 CSQ 欄位定義，請確認 VEP 有正確輸出", file=sys.stderr)
        sys.exit(1)
    print(f"[parse_mito_vcf] CSQ 欄位數：{len(csq_fields)}", file=sys.stderr)

    # ── 輸出 TSV ─────────────────────────────────────────────
    output_cols = [
        # 位置
        "CHROM", "POS", "REF", "ALT",
        # Transcript
        "GENE", "HGVS_C", "HGVS_P", "CONSEQUENCE", "IMPACT", "BIOTYPE",
        # 樣本
        "GENOTYPE", "DP", "AF_SAMPLE",
        # gnomAD mito
        "GNOMAD_MITO_AF", "GNOMAD_MITO_AF_HOM", "GNOMAD_MITO_AF_HET",
        # ClinVar
        "CLINVAR_SIG", "CLINVAR_DN",
        # MITOMAP
        "MITOMAP_LOCUS", "MITOMAP_DISEASE", "MITOMAP_STATUS",
        "MITOMAP_HOMO", "MITOMAP_HETERO", "MITOMAP_MITOTIP",
    ]

    n_total   = 0
    n_mitomap = 0  # 有 MITOMAP 記錄的 variant 數
    n_cfrm    = 0  # MITOMAP Status = Cfrm 的 variant 數

    with open(output_path, "w", encoding="utf-8") as fout:
        fout.write("\t".join(output_cols) + "\n")

        for variant in vcf_in:
            n_total += 1

            chrom = variant.CHROM  # 應該都是 chrM
            pos   = variant.POS
            ref   = variant.REF
            alt   = variant.ALT[0] if variant.ALT else "."  # 只取第一個 ALT

            # ── 解析 VEP CSQ ──────────────────────────────────
            csq_raw = variant.INFO.get("CSQ", None)
            if csq_raw:
                # CSQ 是逗號分隔的多個 transcript annotation
                csq_list = csq_raw.split(",")
                tx = pick_transcript(csq_list, csq_fields)
            else:
                tx = {}

            gene        = tx.get("SYMBOL", ".") or "."
            hgvs_c      = tx.get("HGVSc",  ".") or "."
            hgvs_p      = tx.get("HGVSp",  ".") or "."
            consequence = tx.get("Consequence", ".") or "."
            impact      = tx.get("IMPACT",      ".") or "."
            biotype     = tx.get("BIOTYPE",     ".") or "."

            # gnomAD mito AF
            gnomad_af, gnomad_hom, gnomad_het = get_gnomad_mito(tx)

            # ClinVar
            clinvar_sig, clinvar_dn = get_clinvar(tx)

            # ── 樣本 FORMAT ───────────────────────────────────
            gt_str, dp_str, af_str = get_sample_metrics(
                variant, sample_idx, pipeline
            )

            # ── 查 MITOMAP ────────────────────────────────────
            mm_rec = query_mitomap(mitomap_lookup, pos, ref, alt)

            if mm_rec:
                n_mitomap += 1
                mm_locus   = mm_rec.get("LOCUS",        ".")
                mm_disease = mm_rec.get("DISEASE",      ".")
                mm_status  = mm_rec.get("STATUS",       ".")
                mm_homo    = mm_rec.get("HOMOPLASMY",   ".")
                mm_hetero  = mm_rec.get("HETEROPLASMY", ".")
                mm_mitotip = mm_rec.get("MITOTIP",      ".")
                if mm_status.startswith("Cfrm"):
                    n_cfrm += 1
            else:
                # 這個位置 MITOMAP 沒有記錄
                mm_locus   = "."
                mm_disease = "."
                mm_status  = "."
                mm_homo    = "."
                mm_hetero  = "."
                mm_mitotip = "."

            # ── 寫出一行 ──────────────────────────────────────
            row = [
                chrom, str(pos), ref, alt,
                gene, hgvs_c, hgvs_p, consequence, impact, biotype,
                gt_str, dp_str, af_str,
                gnomad_af, gnomad_hom, gnomad_het,
                clinvar_sig, clinvar_dn,
                mm_locus, mm_disease, mm_status,
                mm_homo, mm_hetero, mm_mitotip,
            ]
            fout.write("\t".join(row) + "\n")

    vcf_in.close()

    # ── 統計摘要 ─────────────────────────────────────────────
    print(f"[parse_mito_vcf] 完成", file=sys.stderr)
    print(f"  總 variant 數        : {n_total:>6,}", file=sys.stderr)
    print(f"  有 MITOMAP 記錄      : {n_mitomap:>6,}", file=sys.stderr)
    print(f"  MITOMAP Cfrm        : {n_cfrm:>6,}", file=sys.stderr)


# ──────────────────────────────────────────────────────────────
# 命令列介面
# ──────────────────────────────────────────────────────────────

def parse_args():
    parser = argparse.ArgumentParser(
        description="解析 VEP 註解後的 mito VCF，輸出 MITO TSV",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用範例：
  python3 parse_mito_vcf.py \\
      --vcf      NA12878_WES.mito.vep.vcf.gz \\
      --sample   NA12878_WES \\
      --mitomap  mitomap_lookup.tsv.gz \\
      --pipeline nckuh \\
      --output   NA12878_WES.mito.tsv
        """
    )
    parser.add_argument("--vcf",      required=True,
                        help="VEP 註解後的 mito VCF（.vcf.gz）")
    parser.add_argument("--sample",   required=True,
                        help="Sample ID（對應 VCF 的 sample column 名稱）")
    parser.add_argument("--mitomap",  required=True,
                        help="mitomap_lookup.tsv.gz")
    parser.add_argument("--pipeline", required=True,
                        choices=["nckuh", "dragen"],
                        help="pipeline 類型（影響 CALLERS tag 處理）")
    parser.add_argument("--output",   required=True,
                        help="輸出 TSV 路徑")
    return parser.parse_args()


def main():
    args = parse_args()
    print(f"[parse_mito_vcf] 輸入 VCF  : {args.vcf}",      file=sys.stderr)
    print(f"[parse_mito_vcf] Sample ID : {args.sample}",   file=sys.stderr)
    print(f"[parse_mito_vcf] MITOMAP   : {args.mitomap}",  file=sys.stderr)
    print(f"[parse_mito_vcf] Pipeline  : {args.pipeline}", file=sys.stderr)
    print(f"[parse_mito_vcf] 輸出      : {args.output}",   file=sys.stderr)
    parse_mito_vcf(
        args.vcf, args.sample, args.mitomap, args.pipeline, args.output
    )


if __name__ == "__main__":
    main()
