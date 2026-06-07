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

build_str_lookup.py
===================
目的：
    把 STRchive-loci.json 轉成兩個查表檔，
    供 parse_str_vcf.py 使用：

    1. str_lookup_varid.tsv.gz
       Key = STRchive locus ID（例如 HTT_HTT）
       供 DRAGEN（ExpansionHunter）使用：
         DRAGEN VCF 的 INFO/VARID 欄位直接對應 STRchive id

    2. str_lookup_pos.tsv.gz
       Key = chrom:pos（1-based VCF 座標，例如 chr4:3074877）
       供 NCKUH（GangSTR）使用：
         GangSTR VCF 沒有 VARID，用位置對齊
         轉換規則：STRchive start_hg38（0-based）+ 1 = VCF POS（1-based）

輸出欄位（兩個查表檔格式相同，只有 KEY 欄位不同）：
    KEY             查表 key（VARID 或 chrom:pos）
    STR_ID          STRchive locus ID（例如 HTT_HTT）
    GENE            基因名稱（例如 HTT）
    CHROM           染色體
    START_HG38      0-based start（STRchive 原始座標）
    END_HG38        0-based end（STRchive 原始座標）
    MOTIF           repeat 單元（reference orientation，例如 CAG）
    DISEASE         疾病名稱
    DISEASE_ID      STRchive disease ID（例如 HD）
    INHERITANCE     遺傳模式（AD/AR/XL/...，多個以逗號分隔）
    BENIGN_MAX      正常上限（repeat count）
    PATHOGENIC_MIN  致病下限（repeat count）
    INTERMEDIATE_MIN  中間範圍下限（無則 .）
    INTERMEDIATE_MAX  中間範圍上限（無則 .）
    LOCUS_STRUCTURE repeat 結構（例如 (CAG)*）
    TYPE            locus 在基因的位置（例如 exon, 3' UTR）

使用方式：
    python3 build_str_lookup.py \\
        --json    STRchive-loci.json \\
        --output_varid  str_lookup_varid.tsv.gz \\
        --output_pos    str_lookup_pos.tsv.gz
"""

import argparse
import gzip
import json
import sys


# ──────────────────────────────────────────────────────────────
# 輔助函式
# ──────────────────────────────────────────────────────────────

def safe_str(val) -> str:
    """把任何值轉成字串，None 轉成 '.'"""
    if val is None:
        return "."
    if isinstance(val, list):
        if len(val) == 0:
            return "."
        # list 轉逗號分隔字串
        return ",".join(str(v) for v in val if v is not None)
    return str(val)


def extract_motif(locus: dict) -> str:
    """
    從 locus 取得 reference orientation 的 motif。
    優先用 reference_motif_reference_orientation（list），
    再用 locus_structure 提取（例如 '(CAG)*' → 'CAG'），
    最後 fallback 到 '.'。
    """
    # 優先：reference_motif_reference_orientation
    motif_list = locus.get("reference_motif_reference_orientation", [])
    if motif_list and len(motif_list) > 0:
        return motif_list[0]

    # 次要：從 locus_structure 提取括號內的序列
    structure = locus.get("locus_structure", "")
    if structure:
        import re
        m = re.search(r'\(([ACGTN]+)\)', structure, re.IGNORECASE)
        if m:
            return m.group(1)

    return "."


def locus_to_row(locus: dict, key: str) -> list:
    """
    把一個 STRchive locus dict 轉成查表 TSV 的一行。
    key 由呼叫端提供（VARID 或 chrom:pos）。
    """
    motif = extract_motif(locus)

    row = [
        key,
        safe_str(locus.get("id")),
        safe_str(locus.get("gene")),
        safe_str(locus.get("chrom")),
        safe_str(locus.get("start_hg38")),
        safe_str(locus.get("stop_hg38")),
        motif,
        safe_str(locus.get("disease")),
        safe_str(locus.get("disease_id")),
        safe_str(locus.get("inheritance")),
        safe_str(locus.get("benign_min")),       # 正常下限（缺失型 locus 需要）
        safe_str(locus.get("benign_max")),
        safe_str(locus.get("pathogenic_min")),
        safe_str(locus.get("pathogenic_max")),   # 致病上限（缺失型 locus 需要）
        safe_str(locus.get("intermediate_min")),
        safe_str(locus.get("intermediate_max")),
        safe_str(locus.get("locus_structure")),
        safe_str(locus.get("type")),
    ]
    return row


# ──────────────────────────────────────────────────────────────
# 主流程
# ──────────────────────────────────────────────────────────────

OUTPUT_COLS = [
    "KEY", "STR_ID", "GENE", "CHROM",
    "START_HG38", "END_HG38", "MOTIF",
    "DISEASE", "DISEASE_ID", "INHERITANCE",
    "BENIGN_MIN", "BENIGN_MAX",
    "PATHOGENIC_MIN", "PATHOGENIC_MAX",
    "INTERMEDIATE_MIN", "INTERMEDIATE_MAX",
    "LOCUS_STRUCTURE", "TYPE",
]


def build_lookup(json_path: str, output_varid: str, output_pos: str):
    """
    讀取 STRchive JSON，產生兩個查表檔。
    """

    # 讀取 JSON
    print(f"[build_str_lookup] 讀取：{json_path}", file=sys.stderr)
    with open(json_path, "r", encoding="utf-8") as f:
        data = json.load(f)

    # STRchive JSON 可能是 list 或 dict（以 id 為 key）
    if isinstance(data, dict):
        loci = list(data.values())
    elif isinstance(data, list):
        loci = data
    else:
        print("[ERROR] 無法識別 JSON 格式（非 list 也非 dict）", file=sys.stderr)
        sys.exit(1)

    print(f"[build_str_lookup] 讀取到 {len(loci)} 個 locus", file=sys.stderr)

    # 統計
    n_varid = 0   # 有 id 的 locus（DRAGEN 查表用）
    n_pos   = 0   # 有 start_hg38 的 locus（NCKUH 查表用）
    n_no_pathogenic = 0  # 沒有 pathogenic_min 的 locus（不能分類）

    header = "\t".join(OUTPUT_COLS) + "\n"

    with gzip.open(output_varid, "wt", encoding="utf-8") as f_varid, \
         gzip.open(output_pos,   "wt", encoding="utf-8") as f_pos:

        f_varid.write(header)
        f_pos.write(header)

        for locus in loci:
            locus_id = locus.get("id", "")
            chrom    = locus.get("chrom", "")
            start    = locus.get("start_hg38")

            # 沒有 pathogenic_min 的 locus 仍然輸出，只是 CLASSIFICATION 會是 "no_threshold"
            if locus.get("pathogenic_min") is None:
                n_no_pathogenic += 1

            # ── VARID 查表（DRAGEN 用）────────────────────────
            # key = STRchive id（對應 DRAGEN INFO/VARID 欄位）
            # DRAGEN 的 VARID 是 gene name（例如 VWA1），不是完整 id（VWA1_VWA1）
            # 所以同時建兩個 key：完整 id 和 gene name
            if locus_id:
                row_id = locus_to_row(locus, locus_id)
                f_varid.write("\t".join(row_id) + "\n")
                n_varid += 1

                # 同時也用 GENE 名稱建 key（DRAGEN VARID 格式）
                gene = locus.get("gene", "")
                if gene and gene != locus_id:
                    row_gene = locus_to_row(locus, gene)
                    f_varid.write("\t".join(row_gene) + "\n")

            # ── 位置查表（NCKUH GangSTR 用）──────────────────
            # key = chrom:pos（VCF 1-based）
            # STRchive start_hg38 已經是 1-based，直接使用，不需要加 1
            # 驗證：HTT start_hg38=3074877，GangSTR POS=3074877，完全一致
            if chrom and start is not None:
                key_pos = f"{chrom}:{int(start)}"
                row_pos = locus_to_row(locus, key_pos)
                f_pos.write("\t".join(row_pos) + "\n")
                n_pos += 1

    # 統計摘要
    print(f"[build_str_lookup] 完成", file=sys.stderr)
    print(f"  總 locus 數            : {len(loci):,}", file=sys.stderr)
    print(f"  VARID 查表（DRAGEN）   : {n_varid:,} 筆 → {output_varid}", file=sys.stderr)
    print(f"  位置查表（NCKUH）      : {n_pos:,} 筆 → {output_pos}", file=sys.stderr)
    print(f"  無 pathogenic_min      : {n_no_pathogenic:,} 筆（分類時標記 no_threshold）",
          file=sys.stderr)


# ──────────────────────────────────────────────────────────────
# 命令列介面
# ──────────────────────────────────────────────────────────────

def parse_args():
    parser = argparse.ArgumentParser(
        description="把 STRchive-loci.json 轉成 DRAGEN/NCKUH 查表檔",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用範例：
  python3 build_str_lookup.py \\
      --json         STRchive-loci.json \\
      --output_varid str_lookup_varid.tsv.gz \\
      --output_pos   str_lookup_pos.tsv.gz
        """
    )
    parser.add_argument("--json",         required=True,
                        help="STRchive-loci.json 路徑")
    parser.add_argument("--output_varid", required=True,
                        help="DRAGEN 用查表（key=VARID）輸出路徑")
    parser.add_argument("--output_pos",   required=True,
                        help="NCKUH 用查表（key=chrom:pos）輸出路徑")
    return parser.parse_args()


def main():
    args = parse_args()
    print(f"[build_str_lookup] JSON         : {args.json}",         file=sys.stderr)
    print(f"[build_str_lookup] output_varid : {args.output_varid}", file=sys.stderr)
    print(f"[build_str_lookup] output_pos   : {args.output_pos}",   file=sys.stderr)
    build_lookup(args.json, args.output_varid, args.output_pos)


if __name__ == "__main__":
    main()
