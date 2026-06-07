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
 
build_gnomad_mito_lookup.py
===========================
從 gnomAD v3.1 mito VCF 建立本機查表。

輸入：
    --vcf     gnomad.genomes.v3.1.sites.chrM.vcf.bgz
    --output  gnomad_mito_lookup.tsv.gz

輸出欄位（TSV）：
    KEY         : POS:REF:ALT（e.g. "1234:A:G"）
    AF_HOM      : 同質性 AF（heteroplasmy >= 0.95）
    AF_HET      : 異質性 AF（heteroplasmy 0.10-0.95）
    AN          : total allele number（樣本數）

License：gnomAD v3.1 是 CC0，完全開放，包括商業用途。

使用範例：
    python3 build_gnomad_mito_lookup.py \\
        --vcf    /path/to/gnomad.genomes.v3.1.sites.chrM.vcf.bgz \\
        --output /path/to/gnomad_mito_lookup.tsv.gz
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


def build_lookup(vcf_path: str, output_path: str):
    vcf_in = VCF(vcf_path)

    n_total   = 0
    n_written = 0

    with gzip.open(output_path, "wt", encoding="utf-8") as fout:
        writer = csv.writer(fout, delimiter="\t")
        writer.writerow(["KEY", "AF_HOM", "AF_HET", "AN"])

        for variant in vcf_in:
            n_total += 1

            pos = variant.POS
            ref = variant.REF
            alts = variant.ALT

            if not alts:
                continue

            # gnomAD mito VCF 通常是 bi-allelic，只取第一個 ALT
            alt = alts[0]

            # 跳過不需要的 variant
            # FILTER 不是 PASS 的跳過（artifact_prone_sites 等）
            filters = variant.FILTER
            if filters and filters != "PASS":
                continue

            af_hom = variant.INFO.get("AF_hom", None)
            af_het = variant.INFO.get("AF_het", None)
            an     = variant.INFO.get("AN",     None)

            # 兩個 AF 都是 0 且 AN 很小的話沒有意義，但仍保留（讓下游判斷）
            af_hom_str = f"{af_hom:.6f}" if af_hom is not None else "."
            af_het_str = f"{af_het:.6f}" if af_het is not None else "."
            an_str     = str(an) if an is not None else "."

            key = f"{pos}:{ref}:{alt}"
            writer.writerow([key, af_hom_str, af_het_str, an_str])
            n_written += 1

    vcf_in.close()

    print(f"[build_gnomad_mito_lookup] 完成", file=sys.stderr)
    print(f"  VCF 總 variant 數  : {n_total:>6,}", file=sys.stderr)
    print(f"  寫入 lookup 筆數   : {n_written:>6,}", file=sys.stderr)
    print(f"  輸出               : {output_path}", file=sys.stderr)


def parse_args():
    parser = argparse.ArgumentParser(
        description="從 gnomAD v3.1 mito VCF 建立本機查表（CC0）"
    )
    parser.add_argument("--vcf",    required=True,
                        help="gnomad.genomes.v3.1.sites.chrM.vcf.bgz")
    parser.add_argument("--output", required=True,
                        help="gnomad_mito_lookup.tsv.gz")
    return parser.parse_args()


def main():
    args = parse_args()
    print(f"[build_gnomad_mito_lookup] VCF    : {args.vcf}",    file=sys.stderr)
    print(f"[build_gnomad_mito_lookup] Output : {args.output}", file=sys.stderr)
    build_lookup(args.vcf, args.output)


if __name__ == "__main__":
    main()
