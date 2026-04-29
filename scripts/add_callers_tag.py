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
 * Copyright (c) 2026, Po-Yu Lin
 * Licensed under the MIT License
 *
 * This pipeline was developed for clinical germline variant
 * analysis. Please cite appropriately if used in research.
 *
 * DISCLAIMER: This pipeline is provided "as is" without
 * warranty of any kind. The authors and their institution
 * make no representations or warranties regarding the
 * accuracy, completeness, or suitability of the analysis
 * results for any clinical or research purpose. Users are
 * solely responsible for validating and interpreting all
 * results. This software shall not be held liable for any
 * direct, indirect, or consequential damages arising from
 * its use.
 * =========================================================

add_callers_tag.py
==================
目的：
    將二級分析產生的 ensemble.fixed.vcf.gz（含兩個 sample column：
    {SAMPLE_ID}_DV 和 {SAMPLE_ID}_HC）處理為三級分析可用的單一樣本 VCF，
    並在 INFO 欄位新增 CALLERS tag，記錄每個 variant 是由哪些 caller 偵測到的。

輸出的 CALLERS 值：
    DV+HC  → DeepVariant 和 HaplotypeCaller 都有 call（最高信心度）
    DV     → 只有 DeepVariant call
    HC     → 只有 HaplotypeCaller call

使用方式（由 Nextflow module prepare_vcf.nf 呼叫）：
    bcftools view ensemble.fixed.vcf.gz | \\
    python3 add_callers_tag.py --sample NA12878_WES | \\
    bgzip -c > snv_for_annotation.vcf.gz
    tabix -p vcf snv_for_annotation.vcf.gz

    或直接指定輸入輸出檔：
    python3 add_callers_tag.py \\
        --input  NA12878_WES.ensemble.fixed.vcf.gz \\
        --sample NA12878_WES \\
        --output snv_for_annotation.vcf.gz

依賴：
    pip install cyvcf2   （比 pysam 更快，處理大型 VCF 效率更好）
    bgzip / tabix        （htslib，通常已在環境中）
"""

import argparse
import sys
import os

# cyvcf2 是處理 VCF 的高效 Python library，速度比 pysam 快約 3-5 倍
try:
    from cyvcf2 import VCF, Writer
except ImportError:
    print("[ERROR] 請先安裝 cyvcf2：pip install cyvcf2", file=sys.stderr)
    sys.exit(1)


# ──────────────────────────────────────────────
# GT 判斷輔助函式
# ──────────────────────────────────────────────

def is_called(gt_tuple) -> bool:
    """
    判斷一個 sample 的 GT 是否為有效的 variant call。

    cyvcf2 的 GT 格式是整數 tuple，規則如下：
        None 或 (-1, -1)  → missing（./. 或 .）
        (0, 0)            → homozygous reference（0/0）
        其他              → 有 variant call（0/1, 1/1, 1/2 等）

    Args:
        gt_tuple: cyvcf2 variant.genotypes[i] 的前兩個元素（allele1, allele2）
                  注意：cyvcf2 genotypes 格式為 [allele1, allele2, phased]
                  allele 值：-1 = missing, 0 = REF, 1 = first ALT, 2 = second ALT...

    Returns:
        True  → 此 caller 確實 call 了這個 variant
        False → missing 或 homozygous reference（視為未 call）
    """
    if gt_tuple is None:
        return False

    a1, a2 = gt_tuple[0], gt_tuple[1]

    # -1 代表 missing allele（./.）
    if a1 == -1 or a2 == -1:
        return False

    # (0, 0) 代表 homozygous reference（0/0），視為未 call
    if a1 == 0 and a2 == 0:
        return False

    # 其餘情況（0/1, 1/0, 1/1, 1/2 等）都算是有 call
    return True


def determine_callers(variant, dv_idx: int, hc_idx: int) -> str:
    """
    依據 DV 和 HC 兩個 sample column 的 GT 決定 CALLERS 值。

    Args:
        variant: cyvcf2 Variant 物件
        dv_idx:  DV sample 在 variant.genotypes 中的索引
        hc_idx:  HC sample 在 variant.genotypes 中的索引

    Returns:
        "DV+HC" / "DV" / "HC"
    """
    # cyvcf2 的 genotypes 屬性：list of [allele1, allele2, phased]
    dv_gt = variant.genotypes[dv_idx]
    hc_gt = variant.genotypes[hc_idx]

    dv_called = is_called(dv_gt)
    hc_called = is_called(hc_gt)

    if dv_called and hc_called:
        return "DV+HC"
    elif dv_called:
        return "DV"
    else:
        # 只剩 HC（不可能兩者都沒 call，因為 ensemble VCF 至少一個 caller call 了才會出現）
        return "HC"


# ──────────────────────────────────────────────
# 主要處理函式
# ──────────────────────────────────────────────

def add_callers_tag(input_path: str, sample_id: str, output_path: str):
    """
    主要處理流程：
    1. 開啟輸入 VCF
    2. 確認 DV / HC sample column 存在
    3. 在 header 新增 CALLERS INFO 定義
    4. 逐 variant 判定 CALLERS，寫出輸出 VCF

    Args:
        input_path:  輸入 VCF 路徑（支援 .vcf / .vcf.gz / "-" 代表 stdin）
        sample_id:   樣本 ID（如 "NA12878_WES"），用於找 DV/HC column
        output_path: 輸出 VCF 路徑（支援 .vcf / "-" 代表 stdout，
                     .vcf.gz 需搭配外部 bgzip）
    """

    # ── 開啟輸入 VCF ──────────────────────────────
    vcf_in = VCF(input_path)

    # ── 確認 sample column 存在 ──────────────────
    # 預期的 sample 名稱：{sample_id}_DV 和 {sample_id}_HC
    expected_dv = f"{sample_id}_DV"
    expected_hc = f"{sample_id}_HC"

    samples = vcf_in.samples  # 取得所有 sample 名稱的 list
    print(f"[INFO] VCF 中的 sample columns：{samples}", file=sys.stderr)

    if expected_dv not in samples:
        print(f"[ERROR] 找不到 DV sample column：{expected_dv}", file=sys.stderr)
        print(f"[ERROR] 可用的 columns：{samples}", file=sys.stderr)
        sys.exit(1)

    if expected_hc not in samples:
        print(f"[ERROR] 找不到 HC sample column：{expected_hc}", file=sys.stderr)
        print(f"[ERROR] 可用的 columns：{samples}", file=sys.stderr)
        sys.exit(1)

    # 取得兩個 sample 在 genotypes list 中的索引
    dv_idx = samples.index(expected_dv)
    hc_idx = samples.index(expected_hc)
    print(f"[INFO] DV sample index：{dv_idx}，HC sample index：{hc_idx}", file=sys.stderr)

    # ── 在 header 新增 CALLERS INFO 定義 ─────────
    # VCF INFO 欄位需要先在 header 宣告，否則下游工具（tabix、VEP 等）可能報錯
    vcf_in.add_info_to_header({
        'ID': 'CALLERS',
        'Number': '1',
        'Type': 'String',
        'Description': (
            'Variant callers that detected this variant: '
            'DV+HC (both DeepVariant and HaplotypeCaller), '
            'DV (DeepVariant only), '
            'HC (HaplotypeCaller only)'
        )
    })

    # ── 開啟輸出 VCF ──────────────────────────────
    # cyvcf2 的 Writer 需要傳入已修改 header 的 VCF 物件作為模板
    # mode="w" 輸出未壓縮 VCF（配合外部 bgzip 使用）
    vcf_out = Writer(output_path, vcf_in, mode="w")

    # ── 逐 variant 處理 ───────────────────────────
    n_total = 0       # 處理的 variant 總數
    n_dv_hc = 0       # DV+HC 的數量
    n_dv_only = 0     # DV only 的數量
    n_hc_only = 0     # HC only 的數量

    for variant in vcf_in:
        n_total += 1

        # 判斷此 variant 由哪些 caller 偵測
        callers = determine_callers(variant, dv_idx, hc_idx)

        # 將 CALLERS 寫入 INFO 欄位
        # cyvcf2 直接對 variant.INFO 賦值即可
        variant.INFO["CALLERS"] = callers

        # 統計計數（用於最後的 summary）
        if callers == "DV+HC":
            n_dv_hc += 1
        elif callers == "DV":
            n_dv_only += 1
        else:
            n_hc_only += 1

        # 寫出此 variant
        vcf_out.write_record(variant)

    # ── 關閉檔案 ──────────────────────────────────
    vcf_out.close()
    vcf_in.close()

    # ── 輸出統計摘要 ──────────────────────────────
    print(f"[INFO] 處理完成", file=sys.stderr)
    print(f"[INFO]   總 variant 數：{n_total:,}", file=sys.stderr)
    print(f"[INFO]   DV+HC（雙 caller）：{n_dv_hc:,} ({n_dv_hc/n_total*100:.1f}%)", file=sys.stderr)
    print(f"[INFO]   DV only：{n_dv_only:,} ({n_dv_only/n_total*100:.1f}%)", file=sys.stderr)
    print(f"[INFO]   HC only：{n_hc_only:,} ({n_hc_only/n_total*100:.1f}%)", file=sys.stderr)


# ──────────────────────────────────────────────
# 命令列介面
# ──────────────────────────────────────────────

def parse_args():
    parser = argparse.ArgumentParser(
        description="在 ensemble VCF 的 INFO 欄位新增 CALLERS tag（DV+HC / DV / HC）",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用範例：
  # 搭配 bcftools + bgzip（推薦，輸出 bgzip 壓縮 VCF）
  bcftools view NA12878_WES.ensemble.fixed.vcf.gz | \\
  python3 add_callers_tag.py --sample NA12878_WES | \\
  bgzip -c > snv_for_annotation.vcf.gz && \\
  tabix -p vcf snv_for_annotation.vcf.gz

  # 直接指定輸入輸出（輸出未壓縮 VCF）
  python3 add_callers_tag.py \\
      --input  NA12878_WES.ensemble.fixed.vcf.gz \\
      --sample NA12878_WES \\
      --output snv_for_annotation.vcf
        """
    )

    parser.add_argument(
        "--input", "-i",
        default="-",
        help="輸入 VCF 路徑（預設：stdin，用 - 表示）"
    )
    parser.add_argument(
        "--sample", "-s",
        required=True,
        help="樣本 ID，例如 NA12878_WES（腳本會自動尋找 {sample}_DV 和 {sample}_HC column）"
    )
    parser.add_argument(
        "--output", "-o",
        default="-",
        help="輸出 VCF 路徑（預設：stdout，用 - 表示；建議搭配外部 bgzip 壓縮）"
    )

    return parser.parse_args()


def main():
    args = parse_args()

    print(f"[INFO] 輸入：{args.input}", file=sys.stderr)
    print(f"[INFO] 樣本 ID：{args.sample}", file=sys.stderr)
    print(f"[INFO] 輸出：{args.output}", file=sys.stderr)

    add_callers_tag(
        input_path=args.input,
        sample_id=args.sample,
        output_path=args.output
    )


if __name__ == "__main__":
    main()
