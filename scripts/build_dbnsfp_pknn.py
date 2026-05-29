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
 
build_dbnsfp_pknn.py
====================
將 P-KNN LLR 欄位 streaming left join 進 dbNSFP 4.9c。

兩個檔案都已按 Chr+pos 排序，使用 merge-join（線性掃描），
不需要全部載入記憶體，速度快。

用法：
    python3 build_dbnsfp_pknn.py \\
        --dbnsfp  /scratch/.../dbNSFP4.9c_grch38.gz \\
        --pknn_dir /data/.../P_KNN_7 \\
        --output  /scratch/.../dbNSFP4.9c_with_pknn_grch38.gz

輸入：
  --dbnsfp    dbNSFP4.9c_grch38.gz（已排序，Chr 格式為 "1"）
  --pknn_dir  P_KNN_7/ 目錄（含 chr*.csv，Chr 格式為 "1"）
  --output    輸出 gz 路徑（未 bgzip，需後續 tabix index）

輸出：
  在 dbNSFP 所有欄位後面新增一欄 PKNN_LLR，
  無對應的 row 填 "."

注意：
  輸出是 gzip（Python gzip 模組），不是 bgzip。
  跑完後需要：
    bgzip -d output.gz && bgzip output && tabix -s 1 -b 2 -e 2 output.gz

作者：Po-Yu Lin（林伯昱）
機構：國立成功大學醫院基因醫學部
"""

import argparse
import gzip
import os
import sys
from glob import glob


def load_pknn_chrom(pknn_dir: str, chrom: str) -> dict:
    """
    載入單一染色體的 P-KNN CSV，建立 (pos, ref, alt) → LLR 的 dict。
    chrom 格式為 "1"（無 chr 前綴）。
    """
    pattern = os.path.join(pknn_dir, f"P_KNN_hg38_missense_dbNSFP_chr{chrom}.csv")
    files = glob(pattern)
    if not files:
        return {}

    pknn = {}
    with open(files[0], "r") as f:
        f.readline()  # skip header
        for line in f:
            parts = line.rstrip("\n").split(",")
            if len(parts) < 21:
                continue
            pos = parts[1].strip()
            ref = parts[3].strip()
            alt = parts[4].strip()
            llr_str = parts[20].strip()

            if not llr_str or llr_str in (".", "", "nan", "NA"):
                continue
            try:
                llr = float(llr_str)
            except ValueError:
                continue

            key = (pos, ref, alt)
            # 同一 key 多筆時取絕對值最大
            if key not in pknn or abs(llr) > abs(pknn[key]):
                pknn[key] = llr

    return pknn


def main():
    parser = argparse.ArgumentParser(
        description="將 P-KNN LLR 加入 dbNSFP 4.9c"
    )
    parser.add_argument("--dbnsfp",   required=True, help="dbNSFP4.9c_grch38.gz")
    parser.add_argument("--pknn_dir", required=True, help="P_KNN_7/ 目錄")
    parser.add_argument("--output",   required=True, help="輸出 gz 路徑")
    args = parser.parse_args()

    print(f"[build_dbnsfp_pknn] 開始處理", file=sys.stderr)
    print(f"  dbNSFP : {args.dbnsfp}", file=sys.stderr)
    print(f"  P-KNN  : {args.pknn_dir}", file=sys.stderr)
    print(f"  output : {args.output}", file=sys.stderr)

    current_chrom = None
    pknn = {}

    total = 0
    matched = 0
    unmatched = 0

    with gzip.open(args.dbnsfp, "rt") as fin, \
         gzip.open(args.output, "wt") as fout:

        # 處理 header
        header = fin.readline().rstrip("\n")
        fout.write(header + "\tPKNN_LLR\n")

        for line in fin:
            total += 1
            if total % 5_000_000 == 0:
                print(f"[build_dbnsfp_pknn] 已處理 {total:,} 行，"
                      f"matched={matched:,}，unmatched={unmatched:,}",
                      file=sys.stderr)

            line = line.rstrip("\n")
            parts = line.split("\t")

            if len(parts) < 4:
                fout.write(line + "\t.\n")
                unmatched += 1
                continue

            chrom = parts[0]   # "1", "2", ...
            pos   = parts[1]   # 1-based
            ref   = parts[2]
            alt   = parts[3]

            # 換染色體時重新載入 P-KNN
            if chrom != current_chrom:
                print(f"[build_dbnsfp_pknn] 載入 chr{chrom} P-KNN...",
                      file=sys.stderr)
                pknn = load_pknn_chrom(args.pknn_dir, chrom)
                print(f"  chr{chrom} P-KNN：{len(pknn):,} 筆", file=sys.stderr)
                current_chrom = chrom

            # 查表
            key = (pos, ref, alt)
            if key in pknn:
                llr = pknn[key]
                fout.write(f"{line}\t{llr:.6f}\n")
                matched += 1
            else:
                fout.write(f"{line}\t.\n")
                unmatched += 1

    match_rate = matched / total * 100 if total > 0 else 0
    print(f"\n[build_dbnsfp_pknn] 完成", file=sys.stderr)
    print(f"  總行數    : {total:,}", file=sys.stderr)
    print(f"  有 LLR    : {matched:,} ({match_rate:.1f}%)", file=sys.stderr)
    print(f"  無 LLR    : {unmatched:,}", file=sys.stderr)
    print(f"\n下一步：", file=sys.stderr)
    print(f"  bgzip -d {args.output}", file=sys.stderr)
    outbase = args.output.replace(".gz", "")
    print(f"  bgzip {outbase}", file=sys.stderr)
    print(f"  tabix -s 1 -b 2 -e 2 {outbase}.gz", file=sys.stderr)


if __name__ == "__main__":
    main()
