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
 * =========================================================

build_gene_moi.py
=================
合併 OMIM genemap2.txt 與 ClinGen Gene-Disease Validity TSV，
產生 gene_symbol → mode_of_inheritance 的查表檔案。

輸出的 MOI 值（統一字串格式）：
  AD       → Autosomal Dominant
  AR       → Autosomal Recessive
  XL       → X-linked（dominant 或 recessive，合併為 XL）
  AD/AR    → 同一基因同時有 AD 和 AR 疾病（例如 MYH7）
  MT       → Mitochondrial
  Unknown  → 不在任何資料庫中，或無法確定

輸入：
  --omim_genemap   OMIM genemap2.txt（需向 OMIM 申請免費帳號下載）
                   下載：https://data.omim.org/downloads/{your_api_key}/genemap2.txt
  --clingen_gene   ClinGen Gene-Disease Validity TSV
                   下載：wget https://search.clinicalgenome.org/kb/gene-validity/download
                   （注意：這個跟 ClinGen_gene_curation_list_GRCh38.tsv 是不同的檔案）
  --output         輸出 gene_moi.tsv.gz

優先順序：
  1. ClinGen Gene-Disease Validity（更精確，直接標注 MOI）
  2. OMIM genemap2.txt（覆蓋範圍更廣）
  若同一基因兩個資料來源衝突，以 ClinGen 為準；
  若 ClinGen 未收錄，用 OMIM。

使用方式：
  python3 build_gene_moi.py \\
      --omim_genemap  genemap2.txt \\
      --clingen_gene  gene_disease_validity.csv \\
      --output        gene_moi.tsv.gz

  # 只有 OMIM（沒有申請 ClinGen 下載的情況）：
  python3 build_gene_moi.py \\
      --omim_genemap  genemap2.txt \\
      --output        gene_moi.tsv.gz

輸出格式（gene_moi.tsv.gz）：
  GENE    MOI     SOURCE
  MECP2   XL      ClinGen
  CFTR    AR      ClinGen
  MYH7    AD/AR   OMIM
  ...
"""

import argparse
import csv
import gzip
import re
import sys
from collections import defaultdict


# ══════════════════════════════════════════════════════════════════
# MOI 字串標準化
# ══════════════════════════════════════════════════════════════════

def normalize_moi(raw: str) -> set[str]:
    """
    將各資料庫的 MOI 原始字串統一為 {AD, AR, XL, MT} 的組合。

    OMIM genemap2 的 phenotype 欄位包含遺傳模式縮寫：
        Autosomal dominant → AD
        Autosomal recessive → AR
        X-linked dominant / X-linked recessive / X-linked → XL
        Mitochondrial → MT
        Somatic mutation / Multifactorial / ... → 忽略（不影響 ACMG）

    ClinGen 的 MOI 欄位（Gene-Disease Validity CSV）：
        "Autosomal dominant" → AD
        "Autosomal recessive" → AR
        "X-linked" / "X-linked recessive" / "X-linked dominant" → XL
        "Mitochondrial" → MT
        "Semidominant" → AD（近似）
    """
    raw_lower = raw.lower()
    result = set()

    if "autosomal dominant" in raw_lower or "semidominant" in raw_lower:
        result.add("AD")
    if "autosomal recessive" in raw_lower:
        result.add("AR")
    if "x-linked" in raw_lower or "x linked" in raw_lower:
        result.add("XL")
    if "mitochondrial" in raw_lower:
        result.add("MT")

    return result


def moi_set_to_str(moi_set: set[str]) -> str:
    """
    將 MOI set 轉為排序過的字串。
    優先順序：AD > AR > XL > MT（方便閱讀）

    範例：
        {"AD", "AR"} → "AD/AR"
        {"XL"}       → "XL"
        {}           → "Unknown"
    """
    if not moi_set:
        return "Unknown"
    order = ["AD", "AR", "XL", "MT"]
    sorted_moi = [m for m in order if m in moi_set]
    return "/".join(sorted_moi) if sorted_moi else "Unknown"


# ══════════════════════════════════════════════════════════════════
# ClinGen Gene-Disease Validity 解析
# ══════════════════════════════════════════════════════════════════

def load_clingen_gene_disease(path: str) -> dict[str, str]:
    """
    解析 ClinGen Gene-Disease Validity CSV/TSV。

    下載來源：
        https://search.clinicalgenome.org/kb/gene-validity/download
        （點「Download All」，下載的是 CSV）

    格式（CSV，有 header）：
        GENE SYMBOL, DISEASE LABEL, DISEASE MIM, MOI, SOP, CLASSIFICATION, ...

    注意：
        同一個基因可能對應多個疾病（不同 MOI），
        例如 MYH7 有 AD 的擴張性心肌病和 AR 的肌病，
        這種情況合併為 "AD/AR"。

    回傳：
        gene_symbol → MOI 字串（例如 "AD" / "AR" / "AD/AR"）
    """
    if not path:
        return {}

    import os
    if not os.path.exists(path):
        print(f"[WARNING] ClinGen Gene-Disease TSV 不存在：{path}", file=sys.stderr)
        return {}

    gene_moi: dict[str, set[str]] = defaultdict(set)
    n_loaded = 0

    with open(path, "r", encoding="utf-8-sig") as f:  # utf-8-sig 處理 BOM
        # 嘗試偵測分隔符號（CSV 或 TSV）
        sample = f.read(1024)
        f.seek(0)
        delimiter = "," if sample.count(",") > sample.count("\t") else "\t"

        reader = csv.DictReader(f, delimiter=delimiter)

        # 找 MOI 欄位（不同版本欄位名稱可能不同）
        moi_col = None
        gene_col = None

        for row in reader:
            # 第一行時確定欄位名稱
            if moi_col is None:
                fieldnames = reader.fieldnames or []
                # 找 gene symbol 欄位
                for candidate in ["GENE SYMBOL", "Gene Symbol", "gene_symbol", "HGNC gene symbol"]:
                    if candidate in fieldnames:
                        gene_col = candidate
                        break
                if gene_col is None:
                    # 用第一欄
                    gene_col = fieldnames[0] if fieldnames else None

                # 找 MOI 欄位
                for candidate in ["MOI", "Mode of Inheritance", "mode_of_inheritance",
                                   "Inheritance", "INHERITANCE"]:
                    if candidate in fieldnames:
                        moi_col = candidate
                        break

                if moi_col is None:
                    print(f"[WARNING] ClinGen Gene-Disease CSV：找不到 MOI 欄位", file=sys.stderr)
                    print(f"  現有欄位：{list(fieldnames)[:10]}", file=sys.stderr)
                    return {}

            gene = row.get(gene_col, "").strip()
            moi_raw = row.get(moi_col, "").strip()

            if not gene or not moi_raw:
                continue

            moi_set = normalize_moi(moi_raw)
            if moi_set:
                gene_moi[gene] |= moi_set
                n_loaded += 1

    result = {gene: moi_set_to_str(moi_set) for gene, moi_set in gene_moi.items()}
    print(f"[build_gene_moi] ClinGen Gene-Disease：載入 {len(result):,} 個基因", file=sys.stderr)
    return result


# ══════════════════════════════════════════════════════════════════
# OMIM genemap2.txt 解析
# ══════════════════════════════════════════════════════════════════

def load_omim_genemap2(path: str) -> dict[str, str]:
    """
    解析 OMIM genemap2.txt。

    下載來源：
        https://data.omim.org/downloads/{your_api_key}/genemap2.txt
        （需免費申請 OMIM API key）

    格式（Tab-separated，有 # 開頭的 comment header）：
        # Chromosome  Genomic Position Start  Genomic Position End
        # Cyto Location  Computed Cyto Location  MIM Number
        # Gene Symbols  Gene Name  Approved Symbol  Entrez Gene ID
        # Ensembl Gene ID  Comments  Phenotypes  Mouse Gene Symbol/ID

    關鍵欄位：
        - "Approved Symbol"（index 8）：HGNC 官方基因名稱
        - "Phenotypes"（index 12）：含遺傳模式，例如
          "Noonan syndrome 1, 163950 (3), Autosomal dominant"

    MOI 提取邏輯：
        從 Phenotypes 欄位用 regex 找 "Autosomal dominant/recessive", "X-linked" 等。
        一個基因可能有多個 phenotype（; 分隔），每個都要解析。
    """
    if not path:
        return {}

    import os
    if not os.path.exists(path):
        print(f"[WARNING] OMIM genemap2.txt 不存在：{path}", file=sys.stderr)
        return {}

    gene_moi: dict[str, set[str]] = defaultdict(set)
    n_genes = 0
    n_with_moi = 0

    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")

            # 跳過 comment 行（# 開頭）
            if line.startswith("#"):
                continue
            if not line.strip():
                continue

            parts = line.split("\t")
            if len(parts) < 13:
                continue

            # Approved Symbol（HGNC 官方名稱）
            gene = parts[8].strip()
            if not gene:
                # 備用：Gene Symbols 欄位（可能有多個，取第一個）
                gene_symbols = parts[7].strip()
                gene = gene_symbols.split(",")[0].strip() if gene_symbols else ""
            if not gene:
                continue

            n_genes += 1

            # Phenotypes 欄位（index 12）
            phenotypes = parts[12].strip()
            if not phenotypes:
                continue

            # 每個 phenotype 用 ";" 分隔
            moi_set = normalize_moi(phenotypes)
            if moi_set:
                gene_moi[gene] |= moi_set
                n_with_moi += 1

    result = {gene: moi_set_to_str(moi_set) for gene, moi_set in gene_moi.items()}
    print(f"[build_gene_moi] OMIM genemap2：掃描 {n_genes:,} 個基因，"
          f"有 MOI 記錄 {len(result):,} 個", file=sys.stderr)
    return result


# ══════════════════════════════════════════════════════════════════
# 主流程：合併兩個來源
# ══════════════════════════════════════════════════════════════════

def build_moi(omim_path: str | None, clingen_path: str | None, output_path: str):
    """
    合併 ClinGen（優先）+ OMIM（補充），輸出 gene_moi.tsv.gz。

    合併邏輯：
        - ClinGen 有記錄 → 用 ClinGen（標注 source = ClinGen）
        - ClinGen 沒有但 OMIM 有 → 用 OMIM（標注 source = OMIM）
        - 兩個都沒有 → 不寫入（classifier 遇到這種基因會用 Unknown）
    """
    print(f"[build_gene_moi] 開始建立 MOI lookup table", file=sys.stderr)

    # 載入兩個來源
    clingen_moi = load_clingen_gene_disease(clingen_path) if clingen_path else {}
    omim_moi    = load_omim_genemap2(omim_path)           if omim_path    else {}

    if not clingen_moi and not omim_moi:
        print("[ERROR] 兩個輸入都是空的，請確認至少提供一個資料庫", file=sys.stderr)
        sys.exit(1)

    # 合併：ClinGen 優先
    merged: dict[str, tuple[str, str]] = {}  # gene → (moi, source)

    # 先填 OMIM
    for gene, moi in omim_moi.items():
        if moi != "Unknown":
            merged[gene] = (moi, "OMIM")

    # 再用 ClinGen 覆蓋（優先）
    for gene, moi in clingen_moi.items():
        if moi != "Unknown":
            merged[gene] = (moi, "ClinGen")

    print(f"[build_gene_moi] 合併後共 {len(merged):,} 個基因有 MOI 資訊", file=sys.stderr)

    # 統計 MOI 分布
    moi_counts: dict[str, int] = defaultdict(int)
    source_counts: dict[str, int] = defaultdict(int)
    for gene, (moi, source) in merged.items():
        moi_counts[moi] += 1
        source_counts[source] += 1

    print(f"  MOI 分布：", file=sys.stderr)
    for moi, cnt in sorted(moi_counts.items(), key=lambda x: -x[1]):
        print(f"    {moi:<12} : {cnt:>5,}", file=sys.stderr)
    print(f"  資料來源：", file=sys.stderr)
    for src, cnt in sorted(source_counts.items(), key=lambda x: -x[1]):
        print(f"    {src:<12} : {cnt:>5,}", file=sys.stderr)

    # 寫出 TSV（gzip 壓縮，讀取快）
    print(f"[build_gene_moi] 寫出：{output_path}", file=sys.stderr)
    with gzip.open(output_path, "wt", encoding="utf-8") as fout:
        fout.write("GENE\tMOI\tSOURCE\n")
        for gene, (moi, source) in sorted(merged.items()):
            fout.write(f"{gene}\t{moi}\t{source}\n")

    print(f"[build_gene_moi] 完成", file=sys.stderr)


# ══════════════════════════════════════════════════════════════════
# 命令列介面
# ══════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="合併 OMIM + ClinGen 建立 gene MOI lookup table",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
資料下載方式：

  OMIM genemap2.txt（需免費申請 API key）：
    https://www.omim.org/downloads → 申請 → 下載 genemap2.txt

  ClinGen Gene-Disease Validity CSV：
    wget "https://search.clinicalgenome.org/kb/gene-validity/download" \\
         -O clingen_gene_disease_validity.csv

使用範例：
  # 兩個都有
  python3 build_gene_moi.py \\
      --omim_genemap  genemap2.txt \\
      --clingen_gene  clingen_gene_disease_validity.csv \\
      --output        gene_moi.tsv.gz

  # 只有 OMIM
  python3 build_gene_moi.py \\
      --omim_genemap  genemap2.txt \\
      --output        gene_moi.tsv.gz

  # 只有 ClinGen
  python3 build_gene_moi.py \\
      --clingen_gene  clingen_gene_disease_validity.csv \\
      --output        gene_moi.tsv.gz
        """
    )
    parser.add_argument(
        "--omim_genemap",
        default=None,
        help="OMIM genemap2.txt 路徑（需申請免費帳號）"
    )
    parser.add_argument(
        "--clingen_gene",
        default=None,
        help="ClinGen Gene-Disease Validity CSV/TSV 路徑"
    )
    parser.add_argument(
        "--output",
        required=True,
        help="輸出 gene_moi.tsv.gz 路徑"
    )
    args = parser.parse_args()

    if not args.omim_genemap and not args.clingen_gene:
        parser.error("至少需要提供 --omim_genemap 或 --clingen_gene 其中一個")

    build_moi(args.omim_genemap, args.clingen_gene, args.output)


if __name__ == "__main__":
    main()
