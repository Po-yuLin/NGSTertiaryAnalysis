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

add_dragen_tag.py
=================
目的：
    將 DRAGEN hard-filtered VCF（單一 sample column）處理為三級分析可用的格式，
    並在 INFO 欄位新增統一的 CALLERS / DP / AD / VAF tag，
    與 ensemble pipeline（add_callers_tag.py）的輸出格式保持一致。

    同時將 VCF 依染色體分流：
        - chrM variants → 輸出到 mito VCF（供 Mito module 使用，保留所有 FILTER）
        - 非 chrM PASS  → 輸出到 snv VCF（供 SNV annotation pipeline 使用）

DRAGEN FORMAT 欄位對應：
    FORMAT/GT  → 直接使用
    FORMAT/DP  → DP_DRAGEN（read depth）
    FORMAT/AD  → AD_DRAGEN（allelic depth，REF,ALT 格式）
    FORMAT/AF  → VAF_DRAGEN（variant allele fraction，DRAGEN 用 AF 而非 VAF）
    FORMAT/GQ  → GQ_DRAGEN（genotype quality，DRAGEN 有，DV 無）

輸出的 INFO tag：
    CALLERS     = "DRAGEN"（固定，區別於 ensemble pipeline 的 DV+HC / DV / HC）
    DP_DRAGEN   = read depth
    AD_DRAGEN   = allelic depth（REF,ALT）
    VAF_DRAGEN  = variant allele fraction（四捨五入至小數點後 4 位）
    GQ_DRAGEN   = genotype quality

使用方式（由 prepare_vcf_dragen.nf 呼叫）：
    python3 add_dragen_tag.py \\
        --input   VAL-10.hard-filtered.vcf.gz \\
        --sample  VAL-10 \\
        --output_snv  VAL-10.snv_for_annotation.vcf \\
        --output_mito VAL-10.mito_for_annotation.vcf

    # 輸出再由 Nextflow module 用 bgzip + tabix 壓縮

依賴：
    pip install cyvcf2
"""

import argparse
import math
import sys

try:
    from cyvcf2 import VCF, Writer
except ImportError:
    print("[ERROR] 請先安裝 cyvcf2：pip install cyvcf2", file=sys.stderr)
    sys.exit(1)


# ──────────────────────────────────────────────────────────────
# FORMAT 欄位擷取輔助函式
# ──────────────────────────────────────────────────────────────

def get_dp(variant, sample_idx: int) -> str:
    """從 FORMAT/DP 取得 read depth，缺失回傳 '.'"""
    try:
        dp = variant.format("DP")
        if dp is None:
            return "."
        val = dp[sample_idx][0]
        return "." if val < 0 else str(val)
    except Exception:
        return "."


def get_ad(variant, sample_idx: int) -> str:
    """
    從 FORMAT/AD 取得 allelic depth。
    DRAGEN AD 格式為 Number=R（REF + 所有 ALT），
    回傳 'REF,ALT'（或 'REF,ALT1,ALT2' for multi-allelic）。
    """
    try:
        ad = variant.format("AD")
        if ad is None:
            return "."
        vals = ad[sample_idx]
        if all(v < 0 for v in vals):
            return "."
        cleaned = [str(max(int(v), 0)) for v in vals]
        return ",".join(cleaned)
    except Exception:
        return "."


def get_vaf(variant, sample_idx: int) -> str:
    """
    從 FORMAT/AF 取得 variant allele fraction。
    DRAGEN 用 AF（Number=A，per ALT allele），不是 VAF。
    multi-allelic 時取第一個 ALT 的值（主要 ALT）。
    """
    try:
        af = variant.format("AF")
        if af is None:
            return "."
        val = af[sample_idx][0]
        if math.isnan(val) or val < 0:
            return "."
        return f"{val:.4f}"
    except Exception:
        return "."


def get_gq(variant, sample_idx: int) -> str:
    """從 FORMAT/GQ 取得 genotype quality，缺失回傳 '.'"""
    try:
        gq = variant.format("GQ")
        if gq is None:
            return "."
        val = gq[sample_idx][0]
        return "." if val < 0 else str(val)
    except Exception:
        return "."


# ──────────────────────────────────────────────────────────────
# 主要處理函式
# ──────────────────────────────────────────────────────────────

def add_dragen_tag(
    input_path: str,
    sample_id: str,
    output_snv_path: str,
    output_mito_path: str,
):
    """
    主要處理流程：
    1. 開啟 DRAGEN hard-filtered VCF（單一 sample column）
    2. 新增 INFO tag：CALLERS, DP_DRAGEN, AD_DRAGEN, VAF_DRAGEN, GQ_DRAGEN
    3. 分流：
       - chrM → mito VCF（保留所有 FILTER，包含 non-PASS）
       - 非 chrM + FILTER=PASS → snv VCF
    4. 輸出兩個未壓縮 VCF（由 Nextflow 後續 bgzip + tabix）
    """

    vcf_in = VCF(input_path)

    # ── 確認 sample column ────────────────────────────────────
    samples = vcf_in.samples
    print(f"[INFO] VCF 中的 sample column：{samples}", file=sys.stderr)

    if sample_id not in samples:
        print(f"[ERROR] 找不到 sample：{sample_id}", file=sys.stderr)
        print(f"  VCF 中的 sample：{samples}", file=sys.stderr)
        sys.exit(1)

    sample_idx = samples.index(sample_id)
    print(f"[INFO] sample index：{sample_idx}", file=sys.stderr)

    # ── 新增 INFO tag 定義到 header ───────────────────────────
    new_info_fields = [
        {
            'ID': 'CALLERS',
            'Number': '1',
            'Type': 'String',
            'Description': (
                'Variant caller: DRAGEN (single caller, no ensemble)'
            )
        },
        {
            'ID': 'DP_DRAGEN',
            'Number': '1',
            'Type': 'String',
            'Description': 'Read depth from DRAGEN (FORMAT/DP). Dot if missing.'
        },
        {
            'ID': 'AD_DRAGEN',
            'Number': '1',
            'Type': 'String',
            'Description': (
                'Allelic depths from DRAGEN (FORMAT/AD), REF,ALT comma-separated. '
                'Dot if missing.'
            )
        },
        {
            'ID': 'VAF_DRAGEN',
            'Number': '1',
            'Type': 'String',
            'Description': (
                'Variant allele fraction from DRAGEN (FORMAT/AF). '
                'Rounded to 4 decimal places. Dot if missing.'
            )
        },
        {
            'ID': 'GQ_DRAGEN',
            'Number': '1',
            'Type': 'String',
            'Description': (
                'Genotype quality from DRAGEN (FORMAT/GQ). '
                'Dot if missing.'
            )
        },
    ]

    for field in new_info_fields:
        vcf_in.add_info_to_header(field)

    # ── 開啟兩個輸出 VCF ─────────────────────────────────────
    vcf_snv  = Writer(output_snv_path,  vcf_in, mode="w")
    vcf_mito = Writer(output_mito_path, vcf_in, mode="w")

    # ── 逐 variant 處理 ──────────────────────────────────────
    n_total    = 0
    n_snv_pass = 0
    n_mito     = 0
    n_filtered = 0   # 非 chrM 且非 PASS（直接丟棄）

    for variant in vcf_in:
        n_total += 1

        chrom = variant.CHROM

        # 新增 INFO tag（所有 variant 都加，不論要輸出到哪）
        variant.INFO["CALLERS"]    = "DRAGEN"
        variant.INFO["DP_DRAGEN"]  = get_dp(variant, sample_idx)
        variant.INFO["AD_DRAGEN"]  = get_ad(variant, sample_idx)
        variant.INFO["VAF_DRAGEN"] = get_vaf(variant, sample_idx)
        variant.INFO["GQ_DRAGEN"]  = get_gq(variant, sample_idx)

        # ── 分流邏輯 ─────────────────────────────────────────
        if chrom == "chrM" or chrom == "MT":
            # Mito：全部保留（包含 non-PASS），供 Mito module 使用
            vcf_mito.write_record(variant)
            n_mito += 1

        else:
            # 非 Mito：只保留 PASS
            # DRAGEN FILTER 欄位可能有多個值（如 LowDepth;PloidyConflict）
            # cyvcf2 的 variant.FILTER 回傳 None 或 list
            filters = variant.FILTER
            is_pass = (filters is None or filters == "PASS" or filters == [] or
                      (isinstance(filters, str) and filters == "PASS"))

            if is_pass:
                vcf_snv.write_record(variant)
                n_snv_pass += 1
            else:
                n_filtered += 1

        if n_total % 500_000 == 0:
            print(
                f"[INFO] 已處理 {n_total:,}（SNV PASS: {n_snv_pass:,}，"
                f"Mito: {n_mito:,}，已過濾: {n_filtered:,}）",
                file=sys.stderr
            )

    vcf_snv.close()
    vcf_mito.close()
    vcf_in.close()

    # ── 統計摘要 ─────────────────────────────────────────────
    print(f"\n[INFO] 處理完成", file=sys.stderr)
    print(f"[INFO]   總 variant 數 ：{n_total:>10,}", file=sys.stderr)
    print(f"[INFO]   SNV/indel PASS：{n_snv_pass:>10,}  → {output_snv_path}", file=sys.stderr)
    print(f"[INFO]   chrM（全部）  ：{n_mito:>10,}  → {output_mito_path}", file=sys.stderr)
    print(f"[INFO]   非 PASS 丟棄  ：{n_filtered:>10,}", file=sys.stderr)


# ──────────────────────────────────────────────────────────────
# 命令列介面
# ──────────────────────────────────────────────────────────────

def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "DRAGEN VCF 前處理：新增 CALLERS/DP/AD/VAF tag，"
            "並分流為 SNV（PASS）和 Mito 兩個輸出"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用範例：
  python3 add_dragen_tag.py \\
      --input   VAL-10.hard-filtered.vcf.gz \\
      --sample  VAL-10 \\
      --output_snv  VAL-10.snv_for_annotation.vcf \\
      --output_mito VAL-10.mito_for_annotation.vcf
        """
    )
    parser.add_argument("--input",       "-i", required=True,
                        help="DRAGEN hard-filtered VCF（.vcf.gz 或 .vcf）")
    parser.add_argument("--sample",      "-s", required=True,
                        help="Sample ID（VCF 中的 sample column 名稱，例如 VAL-10）")
    parser.add_argument("--output_snv",  required=True,
                        help="輸出 SNV/indel VCF（PASS，非 chrM）")
    parser.add_argument("--output_mito", required=True,
                        help="輸出 Mito VCF（chrM 全部，含 non-PASS）")
    return parser.parse_args()


def main():
    args = parse_args()
    print(f"[INFO] 輸入：{args.input}", file=sys.stderr)
    print(f"[INFO] Sample ID：{args.sample}", file=sys.stderr)
    print(f"[INFO] 輸出 SNV：{args.output_snv}", file=sys.stderr)
    print(f"[INFO] 輸出 Mito：{args.output_mito}", file=sys.stderr)
    add_dragen_tag(args.input, args.sample, args.output_snv, args.output_mito)


if __name__ == "__main__":
    main()
