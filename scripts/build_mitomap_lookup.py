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

build_mitomap_lookup.py
=======================
目的：
    把 MITOMAP 下載的兩個 TSV 檔整合成一個統一的查表檔，
    供 parse_mito_vcf.py 在 annotation 時快速查詢。

輸入（兩個來自 https://mitomap.org/downloads/ 的 TSV）：
    --coding    mitomap_mutations_coding_control.tsv
                欄位：Locus, Allele, Position, Nucleotide Change,
                      Amino Acid Change, Plasmy Reports (Homo/Hetero),
                      Disease, Status, GB Freq FL(CR), GB Seqs FL(CR), References
    --rna       mitomap_mutations_rna.tsv
                欄位：Position, Locus, Disease, Allele, RNA,
                      Homoplasmy, Heteroplasmy, Status, MitoTIP,
                      GB Freq FL(CR), GB Seqs FL(CR), References

輸出：
    --output    mitomap_lookup.tsv.gz
    欄位（tab 分隔，第一行是 header）：
        POS         chrM 上的位置（整數，對應 VCF POS）
        LOCUS       MT gene 名稱（例如 MT-ND1, MT-RNR2）
        ALLELE      原始 MITOMAP allele 符號（例如 m.3243A>G）
        NUC_CHANGE  Nucleotide Change（例如 A-G），coding 檔才有
        AA_CHANGE   Amino Acid Change（例如 K→R），coding 檔才有，RNA 填 "."
        RNA_TYPE    RNA 種類（例如 tRNA, rRNA），RNA 檔才有，coding 填 "."
        DISEASE     疾病名稱（逗號分隔，若有多筆）
        STATUS      Cfrm（confirmed）/ Reported / Conflicting
        HOMOPLASMY  + / - / nr（nr = not reported）
        HETEROPLASMY + / - / nr
        MITOTIP     MitoTIP score（只有 RNA 檔有），無則 "."
        SOURCE      coding / rna（資料來源）

查表設計：
    - key = POS（整數字串）
    - 同一個 POS 可能有多筆（例如同一位置有不同的 ALT allele）
    - 此腳本輸出時保留所有筆（每行一筆），parse_mito_vcf.py 依 POS 查到後
      再比對 Nucleotide Change 做精確匹配

使用方式：
    python3 build_mitomap_lookup.py \\
        --coding  mitomap_mutations_coding_control.tsv \\
        --rna     mitomap_mutations_rna.tsv \\
        --output  mitomap_lookup.tsv.gz
"""

import argparse
import csv
import gzip
import re
import sys


# ──────────────────────────────────────────────────────────────
# 輔助函式
# ──────────────────────────────────────────────────────────────

def clean_plasmy(val: str) -> str:
    """
    統一化 Plasmy 欄位的值。
    MITOMAP 的格式不統一：有時是 "+/-"、"nr/nr"、"+" 等。
    回傳乾淨的字串，缺失或 "nr" 統一為 "nr"。
    """
    val = val.strip()
    # 把空白和 "na" 都當 nr
    if not val or val.lower() in ("na", "n/a", ""):
        return "nr"
    return val


def extract_pos_from_allele(allele: str) -> str:
    """
    從 MITOMAP allele 符號提取位置數字。
    例如：
        "m.3243A>G"  → "3243"
        "m.72T>C"    → "72"
        "3243A>G"    → "3243"（沒有 m. 前綴也能處理）

    回傳位置字串，解析失敗回傳空字串。
    """
    # 找第一個連續數字群
    match = re.search(r'(\d+)', allele)
    if match:
        return match.group(1)
    return ""


def parse_nuc_change(nuc_change: str) -> tuple:
    """
    從 Nucleotide Change 欄位提取 REF 和 ALT。
    例如：
        "A-G"   → ("A", "G")
        "T-C"   → ("T", "C")
        "noncoding" → (".", ".")
        "del"   → (".", ".")

    回傳 (ref, alt) tuple，解析失敗回傳 (".", ".")。
    """
    nuc_change = nuc_change.strip()
    # 標準格式：單字母-單字母
    match = re.match(r'^([ACGT])-([ACGT])$', nuc_change, re.IGNORECASE)
    if match:
        return match.group(1).upper(), match.group(2).upper()
    return ".", "."


# ──────────────────────────────────────────────────────────────
# 讀取 coding / control region TSV
# ──────────────────────────────────────────────────────────────

def read_coding_tsv(path: str) -> list:
    """
    讀取 mitomap_mutations_coding_control.tsv，
    回傳 list of dict，每個 dict 是一筆 lookup 紀錄。

    欄位對應（欄位名稱來自實際 TSV header）：
        Locus               → LOCUS
        Allele              → ALLELE（原始符號）
        Position            → POS（整數字串）
        Nucleotide Change   → NUC_CHANGE
        Amino Acid Change   → AA_CHANGE
        Plasmy Reports (Homo/Hetero) → HOMOPLASMY, HETEROPLASMY
        Disease             → DISEASE
        Status              → STATUS
    """
    records = []
    skipped = 0

    # errors='replace'：MITOMAP TSV 內有非 UTF-8 字元（0xa0 不斷行空格），
    # 用 replace 把無法解碼的位元組換成替代符號，不影響主要欄位內容
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            # Position 欄位直接用（比從 Allele 解析更可靠）
            pos_str = row.get("Position", "").strip()

            # 有些行 Position 可能是空白或非數字，跳過
            if not pos_str or not pos_str.isdigit():
                skipped += 1
                continue

            # 解析 Plasmy 欄位（格式是 "homo/hetero"）
            plasmy_raw = row.get("Plasmy Reports (Homo/Hetero)", "nr/nr")
            parts = plasmy_raw.split("/")
            homo  = clean_plasmy(parts[0] if len(parts) > 0 else "nr")
            hetero = clean_plasmy(parts[1] if len(parts) > 1 else "nr")

            records.append({
                "POS":          pos_str,
                "LOCUS":        row.get("Locus", ".").strip(),
                "ALLELE":       row.get("Allele", ".").strip(),
                "NUC_CHANGE":   row.get("Nucleotide Change", ".").strip(),
                "AA_CHANGE":    row.get("Amino Acid Change", ".").strip() or ".",
                "RNA_TYPE":     ".",        # coding 檔沒有 RNA 欄位
                "DISEASE":      row.get("Disease", ".").strip() or ".",
                "STATUS":       row.get("Status", ".").strip() or ".",
                "HOMOPLASMY":   homo,
                "HETEROPLASMY": hetero,
                "MITOTIP":      ".",        # coding 檔沒有 MitoTIP
                "SOURCE":       "coding",
            })

    print(f"[build_mitomap_lookup] coding TSV：讀取 {len(records)} 筆，跳過 {skipped} 筆",
          file=sys.stderr)
    return records


# ──────────────────────────────────────────────────────────────
# 讀取 RNA TSV
# ──────────────────────────────────────────────────────────────

def read_rna_tsv(path: str) -> list:
    """
    讀取 mitomap_mutations_rna.tsv，
    回傳 list of dict。

    欄位對應：
        Position    → POS
        Locus       → LOCUS
        Disease     → DISEASE
        Allele      → ALLELE
        RNA         → RNA_TYPE（例如 tRNA, rRNA）
        Homoplasmy  → HOMOPLASMY
        Heteroplasmy → HETEROPLASMY
        Status      → STATUS
        MitoTIP     → MITOTIP（tRNA in silico score）
    """
    records = []
    skipped = 0

    # 同 coding TSV，用 errors='replace' 處理非 UTF-8 字元
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            pos_str = row.get("Position", "").strip()

            if not pos_str or not pos_str.isdigit():
                skipped += 1
                continue

            records.append({
                "POS":          pos_str,
                "LOCUS":        row.get("Locus", ".").strip(),
                "ALLELE":       row.get("Allele", ".").strip(),
                "NUC_CHANGE":   ".",        # RNA 檔沒有這欄，用 Allele 補
                "AA_CHANGE":    ".",        # RNA 不適用
                "RNA_TYPE":     row.get("RNA", ".").strip() or ".",
                "DISEASE":      row.get("Disease", ".").strip() or ".",
                "STATUS":       row.get("Status", ".").strip() or ".",
                "HOMOPLASMY":   clean_plasmy(row.get("Homoplasmy", "nr")),
                "HETEROPLASMY": clean_plasmy(row.get("Heteroplasmy", "nr")),
                "MITOTIP":      row.get("MitoTIP", ".").strip() or ".",
                "SOURCE":       "rna",
            })

    print(f"[build_mitomap_lookup] RNA TSV：讀取 {len(records)} 筆，跳過 {skipped} 筆",
          file=sys.stderr)
    return records


# ──────────────────────────────────────────────────────────────
# 主流程
# ──────────────────────────────────────────────────────────────

def build_lookup(coding_path: str, rna_path: str, output_path: str):
    """
    合併兩個 TSV，按 POS 排序後輸出 mitomap_lookup.tsv.gz。
    """

    # 讀取兩個檔案
    coding_records = read_coding_tsv(coding_path)
    rna_records    = read_rna_tsv(rna_path)

    # 合併（coding 優先放前面，同 POS 的 rna 排後面）
    all_records = coding_records + rna_records

    # 按 POS 數字排序
    all_records.sort(key=lambda r: int(r["POS"]))

    # 輸出 header 和內容
    output_cols = [
        "POS", "LOCUS", "ALLELE", "NUC_CHANGE", "AA_CHANGE",
        "RNA_TYPE", "DISEASE", "STATUS",
        "HOMOPLASMY", "HETEROPLASMY", "MITOTIP", "SOURCE"
    ]

    print(f"[build_mitomap_lookup] 寫出：{output_path}", file=sys.stderr)
    with gzip.open(output_path, "wt", encoding="utf-8") as fout:
        fout.write("\t".join(output_cols) + "\n")
        for rec in all_records:
            row = [rec[col] for col in output_cols]
            fout.write("\t".join(row) + "\n")

    # 統計摘要
    n_coding = sum(1 for r in all_records if r["SOURCE"] == "coding")
    n_rna    = sum(1 for r in all_records if r["SOURCE"] == "rna")
    n_cfrm   = sum(1 for r in all_records if r["STATUS"].lower() == "cfrm")
    pos_set  = set(r["POS"] for r in all_records)

    print(f"[build_mitomap_lookup] 完成", file=sys.stderr)
    print(f"  總筆數         : {len(all_records):,}", file=sys.stderr)
    print(f"  coding/control : {n_coding:,}", file=sys.stderr)
    print(f"  RNA            : {n_rna:,}", file=sys.stderr)
    print(f"  Cfrm（確認致病）: {n_cfrm:,}", file=sys.stderr)
    print(f"  唯一 POS 數    : {len(pos_set):,}", file=sys.stderr)


# ──────────────────────────────────────────────────────────────
# 命令列介面
# ──────────────────────────────────────────────────────────────

def parse_args():
    parser = argparse.ArgumentParser(
        description="整合 MITOMAP coding/control 和 RNA TSV，建立查表檔",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用範例：
  python3 build_mitomap_lookup.py \\
      --coding  mitomap_mutations_coding_control.tsv \\
      --rna     mitomap_mutations_rna.tsv \\
      --output  mitomap_lookup.tsv.gz
        """
    )
    parser.add_argument("--coding",  required=True,
                        help="mitomap_mutations_coding_control.tsv")
    parser.add_argument("--rna",     required=True,
                        help="mitomap_mutations_rna.tsv")
    parser.add_argument("--output",  required=True,
                        help="輸出 mitomap_lookup.tsv.gz")
    return parser.parse_args()


def main():
    args = parse_args()
    print(f"[build_mitomap_lookup] coding TSV : {args.coding}", file=sys.stderr)
    print(f"[build_mitomap_lookup] RNA TSV    : {args.rna}",    file=sys.stderr)
    print(f"[build_mitomap_lookup] 輸出        : {args.output}", file=sys.stderr)
    build_lookup(args.coding, args.rna, args.output)


if __name__ == "__main__":
    main()
