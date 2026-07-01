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

build_outside_calls.py
======================
整合 StellarPGx（CYP2D6）和 OptiType（HLA-A/HLA-B）輸出，
產生 PharmCAT outside_calls.tsv。

PharmCAT outside calls 格式（tab 分隔，無 header）：
  欄1: HGNC gene symbol（必填）
  欄2: Diplotype（如 *1/*2）
  欄3: Phenotype（如 Normal Metabolizer）
  欄4: Activity score（如 2.0）
  ← 欄 2~4 可擇一或任意組合使用

CYP2D6 範例（有 activity score 最佳）：
  CYP2D6\t*1/*2\tNormal Metabolizer\t2.0

HLA-B 範例（只有 phenotype，注意兩個 tab）：
  HLA-B\t\t*57:01 positive

MT-RNR1 範例（單一 allele）：
  MT-RNR1\t1555A>G

參考：https://pharmcat.clinpgx.org/using/Outside-Call-Format/

作者：Po-Yu Lin（林伯昱）
授權：GNU GPL v3.0
"""

import argparse
import os
import sys


def is_valid_pharmcat_diplotype(value: str) -> bool:
    """
    檢查 StellarPGx 回傳的 diplotype 是否為 PharmCAT 可接受格式。

    PharmCAT 要求 CYP2D6 diplotype 必須是 *allele/*allele 格式（如 *1/*2）。
    StellarPGx 在某些情況下會回傳警告文字（如 "Possible novel allele..."）
    而非合法 diplotype，直接寫入會導致 PharmCAT BadOutsideCallException。

    合法 diplotype 特徵：
      - 包含 "/"（分隔兩個 allele）
      - 不含空格（警告文字通常很長且有空格）
      - 以 "*" 開頭（CYP2D6 allele 格式）
    """
    value = (value or "").strip()
    if not value:
        return False
    if "/" not in value:
        return False
    if " " in value:
        # 警告文字如 "Possible novel allele or suballele present: ..."
        return False
    if not value.startswith("*"):
        return False
    return True


def parse_stellarpgx(path: str) -> dict | None:
    """解析 StellarPGx TSV，回傳 CYP2D6 diplotype 資訊。"""
    if not path or path.startswith("NO_") or not os.path.exists(path):
        return None

    try:
        with open(path, encoding="utf-8") as f:
            header = f.readline().rstrip("\n").split("\t")
            for line in f:
                line = line.strip()
                if not line:
                    continue
                vals = line.split("\t")
                row = dict(zip(header, vals))
                gene      = row.get("GENE", "")
                diplotype = row.get("DIPLOTYPE", "Unknown")
                activity  = row.get("ACTIVITY_SCORE", ".")
                phenotype = row.get("PHENOTYPE", ".")

                if gene != "CYP2D6":
                    continue

                if not is_valid_pharmcat_diplotype(diplotype):
                    # StellarPGx 回傳警告文字或非法格式（如 "Possible novel allele..."）
                    # 不寫入 PharmCAT outside calls，避免 BadOutsideCallException
                    print(
                        f"[BUILD_OUTSIDE_CALLS] CYP2D6 diplotype 不合法，跳過 outside call：{diplotype!r}",
                        file=sys.stderr
                    )
                    return None

                return {
                    "diplotype": diplotype,
                    "activity":  activity if activity != "." else "",
                    "phenotype": phenotype if phenotype not in (".", "Indeterminate") else "",
                }
    except OSError as e:
        print(f"[BUILD_OUTSIDE_CALLS] 警告：讀取 StellarPGx 失敗：{e}", file=sys.stderr)

    return None


def parse_optitype(path: str) -> list[dict]:
    """解析 OptiType TSV，回傳 HLA-A / HLA-B allele 資訊。"""
    results = []
    if not path or path.startswith("NO_") or not os.path.exists(path):
        return results

    try:
        with open(path, encoding="utf-8") as f:
            header = f.readline().rstrip("\n").split("\t")
            for line in f:
                line = line.strip()
                if not line:
                    continue
                vals = line.split("\t")
                row = dict(zip(header, vals))
                gene    = row.get("GENE", "")
                allele1 = row.get("ALLELE_1", ".")
                allele2 = row.get("ALLELE_2", ".")

                if gene not in ("HLA-A", "HLA-B"):
                    continue
                if allele1 == ".":
                    continue

                # PharmCAT HLA outside call 格式：phenotype 欄（欄3），
                # 寫法如 "*57:01 positive" 或 "*57:01/*07:02"
                # PharmCAT 可接受 diplotype（欄2）格式：*07:02/*08:01
                # 但 HLA 建議用 phenotype 欄，格式見官方範例：
                #   HLA-B\t\t*57:01 positive
                # 若有兩個 allele，寫成 diplotype 欄（欄2）比較清楚：
                #   HLA-B\t*07:02/*08:01
                a1 = _strip_hla_prefix(gene, allele1)
                a2 = _strip_hla_prefix(gene, allele2)

                if a2 and a2 != ".":
                    # 有兩個 allele → 寫 diplotype 欄（欄2）
                    results.append({
                        "gene":      gene,
                        "diplotype": f"{a1}/{a2}",
                        "phenotype": "",
                    })
                else:
                    # 只有一個 allele → 寫 phenotype 欄（欄3）
                    results.append({
                        "gene":      gene,
                        "diplotype": "",
                        "phenotype": f"{a1} positive",
                    })

    except OSError as e:
        print(f"[BUILD_OUTSIDE_CALLS] 警告：讀取 OptiType 失敗：{e}", file=sys.stderr)

    return results


def _strip_hla_prefix(gene: str, allele: str) -> str:
    """把 HLA-A*01:01 或 A*01:01 統一成 *01:01（PharmCAT 不需要 gene prefix）。"""
    if not allele or allele == ".":
        return "."
    # 移除 "HLA-A" 或 "A" prefix，只保留 *XX:XX 部分
    for prefix in (f"{gene}*", f"{gene[-1]}*"):
        if allele.upper().startswith(prefix.upper()):
            return "*" + allele.split("*", 1)[1]
    # 已是 *XX:XX 格式
    if allele.startswith("*"):
        return allele
    return allele


def write_outside_calls(
    cyp2d6: dict | None,
    hla_list: list[dict],
    output_path: str,
):
    """
    寫入 PharmCAT outside_calls.tsv。
    格式：Gene<TAB>Diplotype<TAB>Phenotype<TAB>ActivityScore
    無 header 行，# 開頭為 comment。
    """
    lines = []
    lines.append("# PharmCAT outside calls\n")
    lines.append("# Gene\tDiplotype\tPhenotype\tActivityScore\n")

    if cyp2d6:
        dip      = cyp2d6.get("diplotype", "")
        pheno    = cyp2d6.get("phenotype", "")
        activity = cyp2d6.get("activity", "")
        # CYP2D6：提供 diplotype + activity score（最準確）
        lines.append(f"CYP2D6\t{dip}\t{pheno}\t{activity}\n")
        print(
            f"[BUILD_OUTSIDE_CALLS] CYP2D6：{dip} / activity={activity} / {pheno}",
            file=sys.stderr
        )
    else:
        print(
            "[BUILD_OUTSIDE_CALLS] 無 CYP2D6 outside call（PharmCAT 改用 VCF call）",
            file=sys.stderr
        )

    for hla in hla_list:
        gene  = hla["gene"]
        dip   = hla.get("diplotype", "")
        pheno = hla.get("phenotype", "")
        # HLA：有 diplotype 就用 diplotype 欄，否則用 phenotype 欄（兩個 tab）
        lines.append(f"{gene}\t{dip}\t{pheno}\t\n")
        print(f"[BUILD_OUTSIDE_CALLS] {gene}：dip={dip} pheno={pheno}", file=sys.stderr)

    with open(output_path, "w", encoding="utf-8") as f:
        f.writelines(lines)

    print(
        f"[BUILD_OUTSIDE_CALLS] 輸出：{output_path}（{len(lines) - 2} 筆 outside calls）",
        file=sys.stderr
    )


def main():
    parser = argparse.ArgumentParser(
        description="整合 StellarPGx + OptiType → PharmCAT outside_calls.tsv"
    )
    parser.add_argument("--stellarpgx", required=True)
    parser.add_argument("--optitype",   required=True)
    parser.add_argument("--sample",     required=True)
    parser.add_argument("--output",     required=True)
    args = parser.parse_args()

    print(f"[BUILD_OUTSIDE_CALLS] 開始整合 {args.sample}", file=sys.stderr)

    cyp2d6   = parse_stellarpgx(args.stellarpgx)
    hla_list = parse_optitype(args.optitype)

    write_outside_calls(cyp2d6, hla_list, args.output)


if __name__ == "__main__":
    main()
