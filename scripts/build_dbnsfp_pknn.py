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
import csv
import gzip
import os
import sys
from glob import glob

# P-KNN CSV 有欄位長度極大的情況（引號內含長字串），放寬 csv 模組的欄位上限
csv.field_size_limit(10 * 1024 * 1024)


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
    n_short = n_badllr = 0
    # ⚠️ 必須用 csv.reader 而非 line.split(",")：P-KNN CSV 有些欄位是「引號包住、內含逗號」
    #    的值（例如 "Pathogenic,_no_conflicts"）。naive 切割會讓那些行多出好幾欄、欄位整體
    #    位移，第 21 欄（LLR）就會拿到像 `_no_conflicts"` 這種字串 → float() 失敗 → 整筆
    #    被靜默跳過。實測 chr1 有 314 筆因此遺失 LLR。csv.reader 會正確處理引號。
    with open(files[0], "r", newline="") as f:
        reader = csv.reader(f)
        next(reader, None)          # skip header
        for parts in reader:
            if len(parts) < 21:
                n_short += 1
                continue
            pos = parts[1].strip()
            ref = parts[3].strip()
            alt = parts[4].strip()
            llr_str = parts[20].strip()

            if not llr_str or llr_str in (".", "", "nan", "NA"):
                n_badllr += 1
                continue
            try:
                llr = float(llr_str)
            except ValueError:
                n_badllr += 1
                continue

            key = (pos, ref, alt)
            # 同一 key 多筆時取絕對值最大
            if key not in pknn or abs(llr) > abs(pknn[key]):
                pknn[key] = llr

    if n_short or n_badllr:
        print(f"  chr{chrom} P-KNN 跳過：欄數不足 {n_short:,} 筆、LLR 無法解析 {n_badllr:,} 筆",
              file=sys.stderr)
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
    used = set()          # 這條染色體「真的被用到」的 P-KNN key

    total = 0
    matched = 0
    unmatched = 0
    orphan_total = 0      # P-KNN 有、卻沒有任何 dbNSFP 列吃到（★ key 對不上的警訊）
    cat_total = {}        # 對不到的 dbNSFP 列，依 amino acid 變化分類

    # dbNSFP 的 aaref / aaalt 欄位（4.x 與 5.x 都有），用來分類「對不到」的列
    aa_i = {"aaref": None, "aaalt": None}

    def classify_unmatched(parts):
        """
        把「沒有 P-KNN」的 dbNSFP 列分類。
        P-KNN 是 missense-only（檔名即 P_KNN_hg38_missense_*），而 dbNSFP 的 variant 檔
        還包含 nonsense / stoploss / splice-site SNV，這些本來就不會有 LLR，屬預期。
        只有 **missense 卻對不到** 才代表 key 真的沒接上（座標或 allele 不一致）。
        """
        ri, ai = aa_i["aaref"], aa_i["aaalt"]
        if ri is None or ai is None or len(parts) <= max(ri, ai):
            return "unknown"
        r, a = parts[ri], parts[ai]
        if a == "X":
            return "nonsense"
        if r == "X":
            return "stoploss"
        if r in (".", "") or a in (".", ""):
            return "non-coding/splice"
        if r == a:
            return "synonymous"
        return "missense"

    def flush_chrom_stats(chrom, c_tot, c_match, c_unmatch, c_cat):
        """每處理完一條染色體，印出「雙向」統計；回傳 orphan 數。"""
        if chrom is None:
            return 0
        n_load, n_used = len(pknn), len(used)
        n_orphan = n_load - n_used
        print(f"[build_dbnsfp_pknn] === chr{chrom} 統計 ===", file=sys.stderr)
        print(f"    dbNSFP 列數      : {c_tot:,}（有 LLR {c_match:,} / 無 {c_unmatch:,}）",
              file=sys.stderr)
        print(f"    P-KNN 載入       : {n_load:,}", file=sys.stderr)
        print(f"    P-KNN 被用到     : {n_used:,}", file=sys.stderr)
        print(f"    ★P-KNN 沒被用到  : {n_orphan:,}"
              + ("  ← P-KNN 有但 dbNSFP 對不到，需檢查 key" if n_orphan else ""),
              file=sys.stderr)
        if c_cat:
            detail = "、".join(f"{k}={v:,}" for k, v in
                               sorted(c_cat.items(), key=lambda x: -x[1]))
            print(f"    無 LLR 的組成    : {detail}", file=sys.stderr)
        return n_orphan

    c_tot = c_match = c_unmatch = 0
    c_cat = {}

    with gzip.open(args.dbnsfp, "rt") as fin, \
         gzip.open(args.output, "wt") as fout:

        # 處理 header
        header = fin.readline().rstrip("\n")
        header_cols = header.split("\t")
        for name in ("aaref", "aaalt"):
            if name in header_cols:
                aa_i[name] = header_cols.index(name)
        if aa_i["aaref"] is None or aa_i["aaalt"] is None:
            print("[build_dbnsfp_pknn] 警告：header 找不到 aaref/aaalt，"
                  "無法對「對不到的列」分類（合併本身不受影響）", file=sys.stderr)
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

            # 換染色體時：先結算上一條的雙向統計，再載入新的 P-KNN
            if chrom != current_chrom:
                orphan_total += flush_chrom_stats(
                    current_chrom, c_tot, c_match, c_unmatch, c_cat)
                for k, v in c_cat.items():
                    cat_total[k] = cat_total.get(k, 0) + v
                c_tot = c_match = c_unmatch = 0
                c_cat = {}

                print(f"[build_dbnsfp_pknn] 載入 chr{chrom} P-KNN...",
                      file=sys.stderr)
                pknn = load_pknn_chrom(args.pknn_dir, chrom)
                used = set()
                print(f"  chr{chrom} P-KNN：{len(pknn):,} 筆", file=sys.stderr)
                current_chrom = chrom

            # 查表
            key = (pos, ref, alt)
            c_tot += 1
            if key in pknn:
                llr = pknn[key]
                fout.write(f"{line}\t{llr:.6f}\n")
                matched += 1
                c_match += 1
                used.add(key)
            else:
                fout.write(f"{line}\t.\n")
                unmatched += 1
                c_unmatch += 1
                cat = classify_unmatched(parts)
                c_cat[cat] = c_cat.get(cat, 0) + 1

        # 最後一條染色體也要結算
        orphan_total += flush_chrom_stats(
            current_chrom, c_tot, c_match, c_unmatch, c_cat)
        for k, v in c_cat.items():
            cat_total[k] = cat_total.get(k, 0) + v

    match_rate = matched / total * 100 if total > 0 else 0
    print(f"\n[build_dbnsfp_pknn] 完成", file=sys.stderr)
    print(f"  總行數    : {total:,}", file=sys.stderr)
    print(f"  有 LLR    : {matched:,} ({match_rate:.1f}%)", file=sys.stderr)
    print(f"  無 LLR    : {unmatched:,}", file=sys.stderr)
    if cat_total:
        print(f"\n  ── 無 LLR 的組成 ──", file=sys.stderr)
        print(f"  （P-KNN 為 missense-only；nonsense / stoploss / splice 本來就不會有 LLR）",
              file=sys.stderr)
        for k, v in sorted(cat_total.items(), key=lambda x: -x[1]):
            star = "   ★ missense 卻對不到 → 需檢查座標/allele" if k == "missense" and v else ""
            print(f"  {k:<20}: {v:>12,}{star}", file=sys.stderr)
    print(f"\n  ── 反向檢查（P-KNN 是否全數用上）──", file=sys.stderr)
    print(f"  P-KNN 沒被任何 dbNSFP 列用到：{orphan_total:,}"
          + ("（0 = P-KNN 全數對上 ✓）" if orphan_total == 0
             else "  ← 需檢查座標/allele 一致性"),
          file=sys.stderr)
    print(f"\n下一步：", file=sys.stderr)
    print(f"  bgzip -d {args.output}", file=sys.stderr)
    outbase = args.output.replace(".gz", "")
    print(f"  bgzip {outbase}", file=sys.stderr)
    print(f"  tabix -s 1 -b 2 -e 2 {outbase}.gz", file=sys.stderr)


if __name__ == "__main__":
    main()
