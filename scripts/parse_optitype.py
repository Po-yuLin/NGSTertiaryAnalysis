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

parse_optitype.py
=================
解析 OptiType 1.3.1 *_result.tsv，轉成標準格式 TSV。

OptiType *_result.tsv 格式（tab 分隔）：
  欄位：  （空白index）  A1        A2        B1        B2        C1        C2        Reads  Objective
  範例：  0             A*01:01   A*02:01   B*07:02   B*08:01   C*07:01   C*07:02   1234   5678.9

  - 第一欄是 row index（數字），不是 header 的一部分
  - A1/A2 → HLA-A 的兩個 allele
  - B1/B2 → HLA-B 的兩個 allele
  - C1/C2 → HLA-C（PharmCAT 不用，忽略）

輸出欄位：
  GENE  ALLELE_1  ALLELE_2  SOURCE

作者：Po-Yu Lin（林伯昱）
授權：GNU GPL v3.0
"""

import argparse
import os
import sys


def parse_optitype_result(input_path, sample_id):
    """
    解析 OptiType *_result.tsv，回傳 HLA-A 和 HLA-B 的 allele 資訊。

    OptiType 輸出的 allele 格式為 "A*01:01"（含基因前綴）。
    輸出統一轉成 "HLA-A*01:01" 格式（加 "HLA-" prefix），
    讓 build_outside_calls.py 的 _strip_hla_prefix() 可以正確處理。
    """
    results = []

    if not os.path.exists(input_path):
        print("[PARSE_OPTITYPE] 找不到輸入檔：{}".format(input_path), file=sys.stderr)
        return _empty_results()

    try:
        with open(input_path, encoding="utf-8") as f:
            lines = [l.rstrip("\n") for l in f if l.strip()]
    except OSError as e:
        print("[PARSE_OPTITYPE] 讀取失敗：{}".format(e), file=sys.stderr)
        return _empty_results()

    if len(lines) < 2:
        print("[PARSE_OPTITYPE] 檔案內容不足（少於 2 行）：{}".format(input_path), file=sys.stderr)
        return _empty_results()

    # 解析 header（第一行）
    # OptiType header 格式："\tA1\tA2\tB1\tB2\tC1\tC2\tReads\tObjective"
    # 第一個 tab 前是空字串（row index 欄位），所以 split 後第一個元素是 ""
    header = lines[0].split("\t")

    # 找各欄位的 index
    col_idx = {}
    for i, col in enumerate(header):
        col_idx[col.strip()] = i

    # 確認必要欄位存在
    required = ["A1", "A2", "B1", "B2"]
    missing = [c for c in required if c not in col_idx]
    if missing:
        print("[PARSE_OPTITYPE] 缺少欄位 {}，header={}".format(missing, header), file=sys.stderr)
        return _empty_results()

    # 取第一筆資料行（通常只有一行）
    data_line = lines[1].split("\t")

    def _get(col):
        idx = col_idx.get(col)
        if idx is None or idx >= len(data_line):
            return "."
        val = data_line[idx].strip()
        return val if val else "."

    a1 = _get("A1")
    a2 = _get("A2")
    b1 = _get("B1")
    b2 = _get("B2")

    # OptiType 輸出格式是 "A*01:01"，統一加 "HLA-" prefix
    # build_outside_calls.py 的 _strip_hla_prefix() 會再把它去掉
    def _add_hla_prefix(gene_letter, allele):
        if allele == ".":
            return "."
        if allele.upper().startswith("HLA-"):
            return allele
        # "A*01:01" → "HLA-A*01:01"
        if "*" in allele:
            return "HLA-{}".format(allele)
        return allele

    hla_a1 = _add_hla_prefix("A", a1)
    hla_a2 = _add_hla_prefix("A", a2)
    hla_b1 = _add_hla_prefix("B", b1)
    hla_b2 = _add_hla_prefix("B", b2)

    results.append({"gene": "HLA-A", "allele_1": hla_a1, "allele_2": hla_a2})
    results.append({"gene": "HLA-B", "allele_1": hla_b1, "allele_2": hla_b2})

    print("[PARSE_OPTITYPE] {}：HLA-A={}/{}  HLA-B={}/{}".format(sample_id, hla_a1, hla_a2, hla_b1, hla_b2), file=sys.stderr)
    return results


def _empty_results():
    """回傳空白結果（無法解析時的 fallback）。"""
    return [
        {"gene": "HLA-A", "allele_1": ".", "allele_2": "."},
        {"gene": "HLA-B", "allele_1": ".", "allele_2": "."},
    ]


def main():
    parser = argparse.ArgumentParser(description="解析 OptiType *_result.tsv 輸出")
    parser.add_argument("--input",  required=True, help="OptiType *_result.tsv 路徑")
    parser.add_argument("--sample", required=True, help="樣本 ID")
    parser.add_argument("--output", required=True, help="輸出 TSV 路徑")
    args = parser.parse_args()

    results = parse_optitype_result(args.input, args.sample)

    with open(args.output, "w", encoding="utf-8") as f:
        f.write("GENE\tALLELE_1\tALLELE_2\tSOURCE\n")
        for r in results:
            f.write("{}\t{}\t{}\tOptiType\n".format(r['gene'], r['allele_1'], r['allele_2']))

    print("[PARSE_OPTITYPE] 輸出：{}".format(args.output), file=sys.stderr)


if __name__ == "__main__":
    main()
