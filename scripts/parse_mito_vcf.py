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

parse_mito_vcf.py
=================
目的：
    讀取 VEP 註解後的 mito VCF，解析 CSQ 欄位，
    查詢 gnomAD mito v3.1 lookup 表，輸出 mito TSV。

    License 說明：
        舊版使用 MITOMAP（CC BY-NC，收費臨床檢查屬商業使用需授權）。
        此版改用 gnomAD mito v3.1（CC0 完全開放）+ ClinVar（完全開放）。

輸入：
    --vcf          VEP 註解後的 mito VCF（.vcf.gz）
    --sample       sample ID
    --gnomad_mito  gnomad_mito_lookup.tsv.gz（build_gnomad_mito_lookup.py 產生）
    --pipeline     nckuh 或 dragen
    --output       輸出 TSV 路徑（.tsv）

輸出欄位（共 21 欄）：
    CHROM, POS, REF, ALT,
    GENE, HGVS_C, HGVS_P, CONSEQUENCE, IMPACT, BIOTYPE,
    GENOTYPE, DP, AF_SAMPLE,
    GNOMAD_MITO_AF_HOM, GNOMAD_MITO_AF_HET, GNOMAD_MITO_AN,
    CLINVAR_SIG, CLINVAR_DN, CLINVAR_VARIATION_ID,
    OMIM_IDS,
    PIPELINE
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


def load_gnomad_mito(lookup_path: str) -> dict:
    """
    讀取 gnomad_mito_lookup.tsv.gz。
    key = "POS:REF:ALT"，value = {"AF_HOM": ..., "AF_HET": ..., "AN": ...}
    """
    lookup = {}
    with gzip.open(lookup_path, "rt", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            lookup[row["KEY"]] = row
    print(f"[parse_mito_vcf] gnomAD mito lookup 載入：{len(lookup)} 個 variant",
          file=sys.stderr)
    return lookup


def query_gnomad_mito(lookup: dict, pos: int, ref: str, alt: str):
    key = f"{pos}:{ref}:{alt}"
    return lookup.get(key, None)


def parse_csq_header(vcf_obj) -> list:
    for header_line in vcf_obj.raw_header.split("\n"):
        if "ID=CSQ" in header_line and "Format:" in header_line:
            m = re.search(r'Format: ([^"]+)', header_line)
            if m:
                return m.group(1).strip().split("|")
    print("[parse_mito_vcf] [WARN] 找不到 CSQ Format 定義", file=sys.stderr)
    return []


def pick_transcript(csq_list: list, csq_fields: list) -> dict:
    if not csq_list or not csq_fields:
        return {}
    transcripts = []
    for csq_str in csq_list:
        vals = csq_str.split("|")
        while len(vals) < len(csq_fields):
            vals.append("")
        transcripts.append(dict(zip(csq_fields, vals)))
    for tx in transcripts:
        if tx.get("PICK", "") == "1":
            return tx
    return transcripts[0]


def get_clinvar(tx: dict) -> tuple:
    clnsig = tx.get("ClinVar_CLNSIG", ".") or "."
    clndn  = tx.get("ClinVar_CLNDN",  ".") or "."
    # variation ID 欄位名稱依 VEP custom annotation 設定而定
    var_id = tx.get("ClinVar_AF_VARIATION_ID", ".") or "."
    if var_id == ".":
        var_id = tx.get("ClinVar_VARIATION_ID", ".") or "."
    return clnsig, clndn, var_id


def get_sample_metrics(variant, sample_idx: int, pipeline: str) -> tuple:
    # GT
    try:
        gt_tuple = variant.genotypes[sample_idx]
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
    try:
        af_arr = variant.format("AF")
        af_val = af_arr[sample_idx][0] if af_arr is not None else -1.0
        af_str = f"{af_val:.4f}" if af_val >= 0 else "."
    except Exception:
        af_str = "."

    return gt_str, dp_str, af_str


def parse_mito_vcf(vcf_path: str, sample_id: str, gnomad_mito_path: str,
                   pipeline: str, output_path: str):

    gnomad_lookup = load_gnomad_mito(gnomad_mito_path)

    vcf_in = VCF(vcf_path)
    samples = vcf_in.samples
    print(f"[parse_mito_vcf] VCF sample columns：{samples}", file=sys.stderr)

    if sample_id not in samples:
        print(f"[ERROR] 找不到 sample：{sample_id}", file=sys.stderr)
        sys.exit(1)
    sample_idx = samples.index(sample_id)

    csq_fields = parse_csq_header(vcf_in)
    if not csq_fields:
        print("[ERROR] 無法解析 CSQ 欄位定義", file=sys.stderr)
        sys.exit(1)
    print(f"[parse_mito_vcf] CSQ 欄位數：{len(csq_fields)}", file=sys.stderr)

    output_cols = [
        "CHROM", "POS", "REF", "ALT",
        "GENE", "HGVS_C", "HGVS_P", "CONSEQUENCE", "IMPACT", "BIOTYPE",
        "GENOTYPE", "DP", "AF_SAMPLE",
        "GNOMAD_MITO_AF_HOM", "GNOMAD_MITO_AF_HET", "GNOMAD_MITO_AN",
        "CLINVAR_SIG", "CLINVAR_DN", "CLINVAR_VARIATION_ID",
        "OMIM_IDS",
        "PIPELINE",
    ]

    n_total  = 0
    n_gnomad = 0

    with open(output_path, "w", encoding="utf-8") as fout:
        fout.write("\t".join(output_cols) + "\n")

        for variant in vcf_in:
            n_total += 1
            chrom = variant.CHROM
            pos   = variant.POS
            ref   = variant.REF
            alt   = variant.ALT[0] if variant.ALT else "."

            csq_raw = variant.INFO.get("CSQ", None)
            if csq_raw:
                tx = pick_transcript(csq_raw.split(","), csq_fields)
            else:
                tx = {}

            gene        = tx.get("SYMBOL",     ".") or "."
            hgvs_c      = tx.get("HGVSc",      ".") or "."
            hgvs_p      = tx.get("HGVSp",      ".") or "."
            consequence = tx.get("Consequence", ".") or "."
            impact      = tx.get("IMPACT",      ".") or "."
            biotype     = tx.get("BIOTYPE",     ".") or "."
            omim_ids    = tx.get("OMIM_IDS",    ".") or "."

            clinvar_sig, clinvar_dn, clinvar_varid = get_clinvar(tx)
            gt_str, dp_str, af_str = get_sample_metrics(variant, sample_idx, pipeline)

            gm_rec = query_gnomad_mito(gnomad_lookup, pos, ref, alt)
            if gm_rec:
                n_gnomad += 1
                gm_af_hom = gm_rec.get("AF_HOM", ".") or "."
                gm_af_het = gm_rec.get("AF_HET", ".") or "."
                gm_an     = gm_rec.get("AN",     ".") or "."
            else:
                gm_af_hom = gm_af_het = gm_an = "."

            row = [
                chrom, str(pos), ref, alt,
                gene, hgvs_c, hgvs_p, consequence, impact, biotype,
                gt_str, dp_str, af_str,
                gm_af_hom, gm_af_het, gm_an,
                clinvar_sig, clinvar_dn, clinvar_varid,
                omim_ids,
                pipeline,
            ]
            fout.write("\t".join(row) + "\n")

    vcf_in.close()

    print(f"[parse_mito_vcf] 完成", file=sys.stderr)
    print(f"  總 variant 數        : {n_total:>6,}", file=sys.stderr)
    print(f"  gnomAD mito 命中     : {n_gnomad:>6,}", file=sys.stderr)


def parse_args():
    parser = argparse.ArgumentParser(description="解析 VEP 註解後的 mito VCF，輸出 MITO TSV（21 欄）")
    parser.add_argument("--vcf",         required=True)
    parser.add_argument("--sample",      required=True)
    parser.add_argument("--gnomad_mito", required=True)
    parser.add_argument("--pipeline",    required=True, choices=["nckuh", "dragen"])
    parser.add_argument("--output",      required=True)
    return parser.parse_args()


def main():
    args = parse_args()
    print(f"[parse_mito_vcf] VCF        : {args.vcf}",         file=sys.stderr)
    print(f"[parse_mito_vcf] Sample     : {args.sample}",      file=sys.stderr)
    print(f"[parse_mito_vcf] gnomAD mito: {args.gnomad_mito}", file=sys.stderr)
    print(f"[parse_mito_vcf] Pipeline   : {args.pipeline}",    file=sys.stderr)
    print(f"[parse_mito_vcf] Output     : {args.output}",      file=sys.stderr)
    parse_mito_vcf(args.vcf, args.sample, args.gnomad_mito, args.pipeline, args.output)


if __name__ == "__main__":
    main()
