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
parse_pgx_report.py
===================
解析 PharmCAT 3.2.0 report.json，整合 mito_tsv（MT-RNR1），
輸出臨床用 PGx TSV。

PharmCAT 3.2.0 report.json 結構：
  {
    "genes": {
      "CYP2C19": {
        "sourceDiplotypes": [{"label": "*1/*2", "phenotypes": [...], "activityScore": ...}],
        "callSource": "MATCHER" | "OUTSIDE",
        ...
      }
    },
    "drugs": {
      "CPIC Guideline Annotation": {
        "clopidogrel": {
          "guidelines": [
            {
              "annotations": [
                {
                  "implications": [...],
                  "drugRecommendation": "...",
                  "classification": "Strong",
                  "phenotypes": [...],
                  "genotypes": [{"diplotypes": [...]}]
                }
              ]
            }
          ]
        }
      },
      "DPWG Guideline Annotation": { ... }
    }
  }

輸出欄位（16 欄）：
  SAMPLE_ID, PIPELINE, GENE, DIPLOTYPE, ACTIVITY_SCORE,
  PHENOTYPE, DRUG, GUIDELINE_SOURCE, RECOMMENDATION,
  IMPLICATION, CPIC_LEVEL, DPWG_LEVEL,
  OUTSIDE_CALLER, MTRN1_RISK, NOTES, EVIDENCE_STRENGTH

作者：Po-Yu Lin（林伯昱）
授權：GNU GPL v3.0
"""

import argparse
import json
import os
import re
import sys

# CPIC Level A 重點基因
CPIC_LEVEL_A_GENES = {
    "CYP2D6", "CYP2C19", "CYP2C9", "DPYD", "TPMT", "NUDT15",
    "HLA-A", "HLA-B", "SLCO1B1", "G6PD", "MT-RNR1",
    "IFNL3", "CACNA1S", "RYR1", "UGT1A1",
}

# MT-RNR1 已知致病位點
MTRN1_PATHOGENIC = {
    1555: ("m.1555A>G", "aminoglycoside-induced deafness"),
    827:  ("m.827A>G",  "aminoglycoside-induced deafness"),
    1494: ("m.1494C>T", "aminoglycoside-induced deafness"),
}

# guideline source 對應表
SOURCE_MAP = {
    "CPIC Guideline Annotation": "CPIC",
    "DPWG Guideline Annotation": "DPWG",
    "FDA Label Annotation":      "FDA",
    "FDA PGx Association":       "FDA",
}


def strip_html(text):
    if not text or text == ".":
        return text
    clean = re.sub(r"<[^>]+>", " ", str(text))
    clean = re.sub(r"\s+", " ", clean).strip()
    # HTML entity decode
    clean = clean.replace("&quot;", '"').replace("&amp;", "&")                  .replace("&lt;", "<").replace("&gt;", ">")                  .replace("&#39;", "'").replace("&apos;", "'")
    return clean or "."


def parse_gene_info(genes_dict: dict) -> dict:
    """
    從 genes dict 提取每個 gene 的 diplotype、activity score、phenotype、caller。
    回傳 {gene_symbol: {...}} dict。
    """
    gene_info = {}
    for gene, g in genes_dict.items():
        src_dips = g.get("sourceDiplotypes", [])
        if src_dips:
            sd = src_dips[0]
            diplotype     = sd.get("label", ".")
            activity      = sd.get("activityScore", ".")
            if activity is None:
                activity = "."
            phenotypes    = sd.get("phenotypes", [])
            phenotype     = "; ".join(phenotypes) if phenotypes else "."
        else:
            diplotype = "."
            activity  = "."
            phenotype = "."

        call_source = g.get("callSource", "")
        outside_caller = "PharmCAT-outside" if call_source == "OUTSIDE" else ""

        gene_info[gene] = {
            "diplotype":      diplotype,
            "activity_score": str(activity),
            "phenotype":      phenotype,
            "outside_caller": outside_caller,
        }
    return gene_info


def parse_report_json(json_path: str, gene_info: dict) -> list[dict]:
    """
    解析 PharmCAT 3.2.0 report.json 的 drugs 區段。

    drugs 結構：
      drugs[source_name][drug_name]["guidelines"][0]["annotations"][annot]
        annot.implications        → list of str（每個 gene 一條）
        annot.drugRecommendation  → str
        annot.classification      → str（"Strong" / "Moderate" 等）
        annot.phenotypes          → list of str
        annot.genotypes           → list（含 diplotypes）
    """
    rows = []

    try:
        with open(json_path, encoding="utf-8") as f:
            data = json.load(f)
    except (json.JSONDecodeError, OSError) as e:
        print(f"[PGX_PARSE] 警告：無法解析 {json_path}：{e}", file=sys.stderr)
        return rows

    drugs = data.get("drugs", {})

    for source_key, source_drugs in drugs.items():
        source_name = SOURCE_MAP.get(source_key, source_key)
        if not isinstance(source_drugs, dict):
            continue

        for drug_name, drug_data in source_drugs.items():
            if not isinstance(drug_data, dict):
                continue

            guidelines = drug_data.get("guidelines", [])
            if not guidelines:
                continue

            for guideline in guidelines:
                annotations = guideline.get("annotations", [])
                for annot in annotations:
                    # 找出這個 annotation 涉及哪些 gene
                    genotypes = annot.get("genotypes", [])
                    genes_in_annot = set()
                    for gt in genotypes:
                        for dip in gt.get("diplotypes", []):
                            g = dip.get("allele1", {}).get("gene") or dip.get("gene", "")
                            if g:
                                genes_in_annot.add(g)

                    if not genes_in_annot:
                        # 從 implications list 推導 gene
                        for impl_str in annot.get("implications", []):
                            m = re.match(r"^([A-Z0-9\-]+):", impl_str)
                            if m:
                                genes_in_annot.add(m.group(1))

                    if not genes_in_annot:
                        # 最後嘗試：找 gene_info 裡有哪些 gene 在 drug 的 relatedDrugs
                        genes_in_annot = {"Unknown"}

                    classification = annot.get("classification", ".")
                    rec_text = strip_html(annot.get("drugRecommendation", "."))
                    phenotypes = annot.get("phenotypes", [])
                    pheno_str = "; ".join(phenotypes) if phenotypes else "."
                    implications = annot.get("implications", [])
                    impl_str = strip_html("; ".join(implications)) if implications else "."

                    cpic_level = classification if source_name == "CPIC" else "."
                    dpwg_level = classification if source_name == "DPWG" else "."

                    for gene in genes_in_annot:
                        gi = gene_info.get(gene, {})
                        rows.append({
                            "gene":             gene,
                            "diplotype":        gi.get("diplotype", "."),
                            "activity_score":   gi.get("activity_score", "."),
                            "phenotype":        gi.get("phenotype", pheno_str),
                            "drug":             drug_name,
                            "guideline_source": source_name,
                            "recommendation":   rec_text,
                            "implication":      impl_str,
                            "cpic_level":       cpic_level,
                            "dpwg_level":       dpwg_level,
                            "outside_caller":   gi.get("outside_caller", ""),
                            "evidence_strength": classification,
                            "mtrn1_risk":       ".",
                            "notes":            ".",
                        })

    return rows



def parse_mtrn1_vcf(vcf_path: str) -> dict:
    """
    讀取 PGX_MTRN1 的 bcftools mpileup VCF，
    回傳每個位點的 GT 和 DP。
    {pos: {"gt": "0/0", "dp": 361, "af": 0.0}}
    """
    result = {}
    if not vcf_path or vcf_path.startswith("NO_") or not os.path.exists(vcf_path):
        return result

    try:
        import subprocess
        proc = subprocess.run(
            ["bcftools", "view", "-H", vcf_path],
            capture_output=True, text=True
        )
        for line in proc.stdout.splitlines():
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 10:
                continue
            chrom = fields[0]
            if chrom not in ("chrM", "MT", "M"):
                continue
            try:
                pos = int(fields[1])
            except ValueError:
                continue

            fmt  = fields[8].split(":")
            samp = fields[9].split(":")
            fmt_dict = dict(zip(fmt, samp))

            gt = fmt_dict.get("GT", "./.")
            dp = int(fmt_dict.get("DP", "0"))
            ad = fmt_dict.get("AD", "0")
            try:
                ad_vals = [int(x) for x in ad.split(",")]
                ref_ad = ad_vals[0]
                alt_ad = sum(ad_vals[1:]) if len(ad_vals) > 1 else 0
                af = round(alt_ad / dp, 4) if dp > 0 else 0.0
            except (ValueError, ZeroDivisionError):
                af = 0.0

            result[pos] = {"gt": gt, "dp": dp, "af": af}
    except Exception as e:
        print(f"[PGX_PARSE] 警告：讀取 mtrn1 VCF 失敗：{e}", file=sys.stderr)

    return result

def parse_mito_tsv(mito_path: str, mtrn1_vcf_data: dict | None = None) -> list[dict]:
    """讀取 mito.tsv，找 MT-RNR1 aminoglycoside 已知致病位點。

    邏輯：
    - mito.tsv 只有 variant calls，reference 位點不會出現
    - 找到 pathogenic 位點 → HIGH risk
    - 找不到任何 pathogenic 位點（但有 mito.tsv）→ Reference，LOW risk
      （代表這些位點在 mito pipeline 跑過了，只是都是 reference allele）
    - 沒有 mito.tsv → Unknown（無法判斷，不補 reference）
    """
    rows = []
    no_mito = not mito_path or mito_path.startswith("NO_") or not os.path.exists(mito_path)

    if no_mito:
        # 沒有 mito.tsv，無法判斷，回傳空（會顯示 Unknown）
        return rows

    found_pathogenic = []

    try:
        with open(mito_path, encoding="utf-8") as f:
            header = f.readline().rstrip("\n").split("\t")
            for line in f:
                line = line.strip()
                if not line:
                    continue
                vals = line.split("\t")
                row  = dict(zip(header, vals))

                chrom = row.get("CHROM", "")
                if chrom not in ("chrM", "MT", "M"):
                    continue
                try:
                    pos = int(row.get("POS", "0"))
                except ValueError:
                    continue
                if pos not in MTRN1_PATHOGENIC:
                    continue

                hgvs, drug_ctx = MTRN1_PATHOGENIC[pos]
                af     = row.get("AF_SAMPLE", row.get("HETEROPLASMY_AF", "."))
                clnsig = row.get("CLINVAR_SIG", ".")
                note   = f"Heteroplasmy AF={af}; ClinVar={clnsig}"

                found_pathogenic.append(hgvs)
                rows.append({
                    "gene":             "MT-RNR1",
                    "diplotype":        hgvs,
                    "activity_score":   ".",
                    "phenotype":        "Aminoglycoside-induced deafness risk",
                    "drug":             "aminoglycosides (gentamicin, tobramycin, streptomycin)",
                    "guideline_source": "ClinVar/CPIC",
                    "recommendation":   (
                        "Avoid aminoglycoside antibiotics. "
                        "If essential, use with extreme caution and audiological monitoring. "
                        "Counsel all maternal relatives (maternal inheritance)."
                    ),
                    "implication":      drug_ctx,
                    "cpic_level":       "A",
                    "dpwg_level":       ".",
                    "outside_caller":   "mito_pipeline",
                    "evidence_strength": "Strong",
                    "mtrn1_risk":       "HIGH",
                    "notes":            note,
                })
    except OSError as e:
        print(f"[PGX_PARSE] 警告：讀取 mito TSV 失敗：{e}", file=sys.stderr)
        return rows

    # 找不到任何 pathogenic 位點
    if not found_pathogenic:
        # 用 mtrn1_vcf_data 確認 coverage
        # 三個位點中只要一個有 coverage（DP >= 10）就補 Reference
        covered_positions = []
        if mtrn1_vcf_data:
            for pos, data in mtrn1_vcf_data.items():
                if data.get("dp", 0) >= 10:
                    covered_positions.append(pos)

        if covered_positions or no_mito is False:
            # 有 coverage 確認（或有 mito.tsv）→ Reference（LOW risk）
            coverage_note = (
                f"Coverage confirmed at chrM positions: "
                f"{','.join(str(p) for p in sorted(covered_positions))}"
                if covered_positions else
                "Coverage confirmed via mito pipeline"
            )
            rows.append({
                "gene":             "MT-RNR1",
                "diplotype":        "Reference",
                "activity_score":   ".",
                "phenotype":        "Normal - no aminoglycoside deafness risk",
                "drug":             "aminoglycosides (gentamicin, tobramycin, streptomycin)",
                "guideline_source": "ClinVar/CPIC",
                "recommendation":   (
                    "No MT-RNR1 pathogenic variants detected (m.1555A>G, m.1494C>T, m.827A>G). "
                    "Standard aminoglycoside dosing applies."
                ),
                "implication":      "No known MT-RNR1 risk variant detected",
                "cpic_level":       "A",
                "dpwg_level":       ".",
                "outside_caller":   "mito_pipeline",
                "evidence_strength": "Strong",
                "mtrn1_risk":       "LOW",
                "notes":            coverage_note,
            })
        # coverage 不足 → 不補 row，維持 Unknown

    return rows


def write_pgx_tsv(rows: list[dict], sample_id: str, pipeline: str, output: str):
    COLUMNS = [
        "SAMPLE_ID", "PIPELINE", "GENE", "DIPLOTYPE", "ACTIVITY_SCORE",
        "PHENOTYPE", "DRUG", "GUIDELINE_SOURCE", "RECOMMENDATION",
        "IMPLICATION", "CPIC_LEVEL", "DPWG_LEVEL",
        "OUTSIDE_CALLER", "MTRN1_RISK", "NOTES", "EVIDENCE_STRENGTH",
    ]

    # 只保留 CPIC Level A 基因，去重（同一 gene+drug+source 只留一筆）
    seen = set()
    filtered = []
    for r in rows:
        gene = r.get("gene", "")
        if gene not in CPIC_LEVEL_A_GENES:
            continue
        key = (gene, r.get("drug", ""), r.get("guideline_source", ""))
        if key in seen:
            continue
        seen.add(key)
        filtered.append(r)

    filtered.sort(key=lambda r: (r.get("gene", ""), r.get("drug", "")))

    with open(output, "w", encoding="utf-8") as f:
        f.write("\t".join(COLUMNS) + "\n")
        for r in filtered:
            def _v(x): return str(x) if x is not None else "."
            f.write("\t".join([
                _v(sample_id),
                _v(pipeline),
                _v(r.get("gene")),
                _v(r.get("diplotype")),
                _v(r.get("activity_score")),
                _v(r.get("phenotype")),
                _v(r.get("drug")),
                _v(r.get("guideline_source")),
                _v(r.get("recommendation")),
                _v(r.get("implication")),
                _v(r.get("cpic_level")),
                _v(r.get("dpwg_level")),
                _v(r.get("outside_caller")),
                _v(r.get("mtrn1_risk")),
                _v(r.get("notes")),
                _v(r.get("evidence_strength")),
            ]) + "\n")

    print(f"[PGX_PARSE] 輸出 {output}，共 {len(filtered)} 筆", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pharmcat_json",  required=True)
    parser.add_argument("--outside_calls",  required=True)
    parser.add_argument("--mito_tsv",       required=True)
    parser.add_argument("--mtrn1_vcf",      default=None,
                        help="bcftools mpileup VCF for MT-RNR1 coverage check")
    parser.add_argument("--sample",         required=True)
    parser.add_argument("--pipeline",       required=True)
    parser.add_argument("--output",         required=True)
    args = parser.parse_args()

    print(f"[PGX_PARSE] 開始解析 {args.sample}", file=sys.stderr)

    # 載入 JSON 取 gene_info（diplotype / phenotype / activity score）
    try:
        with open(args.pharmcat_json, encoding="utf-8") as f:
            data = json.load(f)
        gene_info = parse_gene_info(data.get("genes", {}))
    except (json.JSONDecodeError, OSError) as e:
        print(f"[PGX_PARSE] 警告：{e}", file=sys.stderr)
        gene_info = {}

    print(f"[PGX_PARSE] gene_info 基因數：{len(gene_info)}", file=sys.stderr)

    pharmcat_rows = parse_report_json(args.pharmcat_json, gene_info)
    print(f"[PGX_PARSE] PharmCAT drug 記錄數：{len(pharmcat_rows)}", file=sys.stderr)

    mtrn1_vcf_data = parse_mtrn1_vcf(args.mtrn1_vcf) if args.mtrn1_vcf else {}
    print(f"[PGX_PARSE] MT-RNR1 VCF 位點數：{len(mtrn1_vcf_data)}", file=sys.stderr)

    mito_rows = parse_mito_tsv(args.mito_tsv, mtrn1_vcf_data)
    print(f"[PGX_PARSE] MT-RNR1 記錄數：{len(mito_rows)}", file=sys.stderr)

    write_pgx_tsv(pharmcat_rows + mito_rows, args.sample, args.pipeline, args.output)


if __name__ == "__main__":
    main()
