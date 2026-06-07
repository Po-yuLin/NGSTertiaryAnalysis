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
 
parse_stellarpgx.py
===================
解析 StellarPGx .alleles 輸出檔，轉成標準格式 TSV。

StellarPGx .alleles 格式（純文字，NA12878_WGS 實測）：
  --------------------------------------------
  CYP2D6 Star Allele Calling with StellarPGx
  --------------------------------------------
  Initially computed CN = 1
  Sample core variants:
  Candidate alleles:
  1.v1_1.v1
  Result:
  *5/*1
  Activity score:
  1.0
  Metaboliser status:
  Intermediate metaboliser (IM)

解析規則：
  - "Result:" 的下一個非空行 → diplotype
  - "Activity score:" 的下一個非空行 → activity score
  - "Metaboliser status:" 的下一個非空行 → phenotype

輸出欄位：
  GENE  DIPLOTYPE  ACTIVITY_SCORE  PHENOTYPE  SOURCE

作者：Po-Yu Lin（林伯昱）
授權：GNU GPL v3.0
"""

import argparse
import os
import sys


def parse_alleles_file(input_path: str, sample_id: str) -> dict:
    result = {
        "gene":           "CYP2D6",
        "diplotype":      "Unknown",
        "activity_score": ".",
        "phenotype":      "Indeterminate",
        "source":         "StellarPGx",
    }

    if not os.path.exists(input_path):
        print(f"[PARSE_STELLARPGX] 找不到輸入檔：{input_path}", file=sys.stderr)
        return result

    try:
        with open(input_path, encoding="utf-8") as f:
            lines = [l.rstrip("\n") for l in f]
    except OSError as e:
        print(f"[PARSE_STELLARPGX] 讀取失敗：{e}", file=sys.stderr)
        return result

    # 逐行解析：找到 key 行後取下一個非空行作為值
    i = 0
    while i < len(lines):
        line = lines[i].strip()

        if line == "Result:":
            # 找下一個非空行
            j = i + 1
            while j < len(lines) and not lines[j].strip():
                j += 1
            if j < len(lines):
                result["diplotype"] = lines[j].strip()
            i = j

        elif line == "Activity score:":
            j = i + 1
            while j < len(lines) and not lines[j].strip():
                j += 1
            if j < len(lines):
                result["activity_score"] = lines[j].strip()
            i = j

        elif line == "Metaboliser status:":
            j = i + 1
            while j < len(lines) and not lines[j].strip():
                j += 1
            if j < len(lines):
                result["phenotype"] = lines[j].strip()
            i = j

        else:
            i += 1

    print(
        f"[PARSE_STELLARPGX] {sample_id}："
        f"diplotype={result['diplotype']}, "
        f"activity={result['activity_score']}, "
        f"phenotype={result['phenotype']}",
        file=sys.stderr
    )
    return result


def main():
    parser = argparse.ArgumentParser(description="解析 StellarPGx .alleles 輸出")
    parser.add_argument("--input",  required=True, help="StellarPGx .alleles 檔案路徑")
    parser.add_argument("--sample", required=True, help="樣本 ID")
    parser.add_argument("--output", required=True, help="輸出 TSV 路徑")
    args = parser.parse_args()

    result = parse_alleles_file(args.input, args.sample)

    with open(args.output, "w", encoding="utf-8") as f:
        f.write("GENE\tDIPLOTYPE\tACTIVITY_SCORE\tPHENOTYPE\tSOURCE\n")
        f.write(
            f"{result['gene']}\t{result['diplotype']}\t"
            f"{result['activity_score']}\t{result['phenotype']}\t{result['source']}\n"
        )

    print(f"[PARSE_STELLARPGX] 輸出：{args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
