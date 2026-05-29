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
 *  * This program is free software: you can redistribute it and/or modify
 *  * it under the terms of the GNU General Public License as published by
 *  * the Free Software Foundation, either version 3 of the License, or
 *  * (at your option) any later version.
 *  *
 *  * This program is distributed in the hope that it will be useful,
 *  * but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  * GNU General Public License for more details.
 *  *
 *  * You should have received a copy of the GNU General Public License
 *  * along with this program. If not, see <https://www.gnu.org/licenses/>.
 *  *
 *  * THIRD-PARTY TOOLS NOTICE:
 *  * This pipeline orchestrates third-party tools subject to their own licenses.
 *  * Users of main_research.nf must comply with:
 *  *   - Manta (Illumina): PolyForm Strict License 1.0.0 (non-commercial only)
 *  *   - ExpansionHunter (Illumina): PolyForm Strict License 1.0.0 (non-commercial only)
 *  * See README.md and LICENSE for details.
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
 
build_clinvar_lookup.py
=======================
從 ClinVar variant_summary.txt.gz 建立精簡的查表檔案，
供 parse_vep_csq.py 填入 CLINVAR_VARIATION_ID 和 OMIM_ID 欄位。

輸入：
  --input    ClinVar variant_summary.txt.gz（月度 archive 版本）
  --output   輸出的 clinvar_lookup.tsv.gz

輸出欄位（clinvar_lookup.tsv.gz）：
  KEY                 查表 key：chr{CHROM}:{POS}:{REF}:{ALT}
  VARIATION_ID        ClinVar Variation ID（組成連結用）
  OMIM_IDS            逗號分隔的 OMIM ID（數字，無 OMIM: 前綴）
  RS_ID               rsID（rs 前綴，無則空）

設計原則：
  - 只保留 GRCh38 且 REF/ALT 不為 "-" 的 variant（真正的 SNV/Indel）
  - 多個 OMIM ID 以逗號分隔（一個 variant 可能對應多個疾病）
  - 同一個 key 若有多筆（multi-allelic 拆分），保留 VariationID 最小的（最早收錄）
  - 輸出約 20-30MB，parse_vep_csq.py 啟動時直接載入

作者：Po-Yu Lin（林伯昱）
機構：國立成功大學醫院基因醫學部
"""

import argparse
import csv
import gzip
import re
import sys
from collections import defaultdict


def extract_omim_ids(phenotype_ids: str) -> str:
    """
    從 PhenotypeIDs 欄位提取所有 OMIM ID（數字部分）。

    輸入格式範例：
      MONDO:MONDO:0013342,MedGen:C3150901,OMIM:613647,Orphanet:306511
      OMIM:113705|OMIM:604370（多個疾病用 | 分隔）

    輸出：逗號分隔的 OMIM 數字 ID，例如 "613647" 或 "113705,604370"
    無 OMIM 則回傳空字串 ""
    """
    if not phenotype_ids or phenotype_ids in ("-", ""):
        return ""

    # 用 | 和 , 分隔後逐一找 OMIM:數字
    omim_ids = re.findall(r'OMIM:(\d+)', phenotype_ids)
    # 去重並保持順序
    seen = set()
    unique_omim = []
    for oid in omim_ids:
        if oid not in seen:
            seen.add(oid)
            unique_omim.append(oid)

    return ",".join(unique_omim)


def build_lookup(input_path: str, output_path: str):
    """
    讀取 variant_summary.txt.gz，建立 GRCh38 的查表檔案。
    """

    # key → (variation_id, omim_ids, rs_id) 的 dict
    # 若同一 key 有多筆，保留 variation_id 較小的（較早收錄）
    lookup = {}

    opener = gzip.open if input_path.endswith(".gz") else open
    total = 0
    kept = 0
    skipped_assembly = 0
    skipped_allele = 0

    print(f"[build_clinvar_lookup] 讀取：{input_path}", file=sys.stderr)

    with opener(input_path, "rt", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)

        # 確認欄位位置（以免未來 ClinVar 改欄位順序）
        col = {name.strip("#"): idx for idx, name in enumerate(header)}

        required = ["Assembly", "Chromosome", "PositionVCF",
                    "ReferenceAlleleVCF", "AlternateAlleleVCF",
                    "VariationID", "PhenotypeIDS", "RS# (dbSNP)"]
        for r in required:
            if r not in col:
                print(f"[ERROR] 找不到欄位 '{r}'，請確認 variant_summary 版本", file=sys.stderr)
                print(f"  現有欄位：{list(col.keys())}", file=sys.stderr)
                sys.exit(1)

        i_assembly    = col["Assembly"]
        i_chrom       = col["Chromosome"]
        i_pos         = col["PositionVCF"]
        i_ref         = col["ReferenceAlleleVCF"]
        i_alt         = col["AlternateAlleleVCF"]
        i_varid       = col["VariationID"]
        i_phenotype   = col["PhenotypeIDS"]
        i_rs          = col["RS# (dbSNP)"]

        for row in reader:
            total += 1
            if total % 500000 == 0:
                print(f"[build_clinvar_lookup] 已處理 {total:,} 筆，保留 {kept:,} 筆...",
                      file=sys.stderr)

            if len(row) <= max(i_assembly, i_chrom, i_pos, i_ref, i_alt, i_varid):
                continue

            # 只保留 GRCh38
            if row[i_assembly] != "GRCh38":
                skipped_assembly += 1
                continue

            chrom   = row[i_chrom].strip()
            pos     = row[i_pos].strip()
            ref     = row[i_ref].strip()
            alt     = row[i_alt].strip()
            varid   = row[i_varid].strip()
            pheno   = row[i_phenotype].strip()
            rs_raw  = row[i_rs].strip()

            # 跳過沒有座標或 REF/ALT 為 "-" 的（deletion 表示法，無法直接對位）
            if not chrom or not pos or ref in ("-", "", "na") or alt in ("-", "", "na"):
                skipped_allele += 1
                continue

            # 跳過非數字的 pos（有些是範圍或空值）
            if not pos.isdigit():
                skipped_allele += 1
                continue

            # 建立 key：加上 chr 前綴（ClinVar 用 "1"，VCF 用 "chr1"）
            chrom_vcf = chrom if chrom.startswith("chr") else f"chr{chrom}"
            key = f"{chrom_vcf}:{pos}:{ref}:{alt}"

            # 提取 OMIM ID
            omim_ids = extract_omim_ids(pheno)

            # rsID：ClinVar 存的是數字（-1 表示無），加上 rs 前綴
            rs_id = ""
            if rs_raw and rs_raw != "-1" and rs_raw.isdigit():
                rs_id = f"rs{rs_raw}"

            # 同一 key 有多筆時，保留 VariationID 較小的
            try:
                varid_int = int(varid)
            except ValueError:
                continue

            if key in lookup:
                existing_varid = int(lookup[key][0])
                if varid_int >= existing_varid:
                    # 保留舊的，但補充 OMIM（可能不同疾病）
                    existing_omim = lookup[key][1]
                    if omim_ids and omim_ids not in existing_omim:
                        merged = existing_omim + ("," if existing_omim else "") + omim_ids
                        lookup[key] = (lookup[key][0], merged, lookup[key][2])
                    continue

            lookup[key] = (varid, omim_ids, rs_id)
            kept += 1

    print(f"[build_clinvar_lookup] 讀取完成：總計 {total:,} 筆", file=sys.stderr)
    print(f"  保留（GRCh38 SNV/Indel）：{kept:,} 筆", file=sys.stderr)
    print(f"  跳過（非 GRCh38）：{skipped_assembly:,} 筆", file=sys.stderr)
    print(f"  跳過（無座標或特殊 allele）：{skipped_allele:,} 筆", file=sys.stderr)

    # 寫出 lookup TSV
    print(f"[build_clinvar_lookup] 寫出：{output_path}", file=sys.stderr)
    with gzip.open(output_path, "wt", encoding="utf-8") as fout:
        fout.write("KEY\tVARIATION_ID\tOMIM_IDS\tRS_ID\n")
        for key, (varid, omim_ids, rs_id) in sorted(lookup.items()):
            fout.write(f"{key}\t{varid}\t{omim_ids}\t{rs_id}\n")

    print(f"[build_clinvar_lookup] 完成，共 {len(lookup):,} 個唯一 key", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description="從 ClinVar variant_summary.txt.gz 建立查表檔案"
    )
    parser.add_argument("--input",  required=True,
                        help="ClinVar variant_summary.txt.gz")
    parser.add_argument("--output", required=True,
                        help="輸出 clinvar_lookup.tsv.gz")
    args = parser.parse_args()
    build_lookup(args.input, args.output)


if __name__ == "__main__":
    main()