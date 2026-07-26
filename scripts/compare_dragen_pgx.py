#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
=========================================================
WGS/WES Germline Tertiary Analysis - compare_dragen_pgx.py
=========================================================
Author : Po-Yu Lin (林伯昱)  <p88124019@gs.ncku.edu.tw>
Licensed under the GNU General Public License v3.0
=========================================================
把 DRAGEN 原生 PGx 判讀（{sample}.targeted.json）交叉註記進我們的 pgx.tsv `NOTES` 欄。
只做「我們 pgx.tsv 有、且 DRAGEN 也有」的基因；其餘不動。**欄位不變**（只改 NOTES 值）。

每個可對照的基因，在 NOTES 後面加上 DRAGEN 的**原始** genotype ＋一致性判定：
  - 一致  ：正規化後我們的 diplotype 落在 DRAGEN 的候選解內  → "DRAGEN 一致: <raw>"
  - 不一致：同命名系統但不相等                                → "DRAGEN 不一致: <raw>"
  - 未比對：命名系統不同（如 star vs HGVS/rs），不妄下判定    → "DRAGEN 未比對: <raw>"

正規化規則（見使用者確認）：
  - reference-like（Reference / wildtype / ref / *1 / B(reference) / B(wildtype)）一律視為 "REF"
    → 解決「Reference/Reference」與「*1/*1」的跨寫法等價。
  - star allele 取排序後的 allele 集合比對。
  - DRAGEN 的模糊多重解以 ';' 分隔 → 只要我們的 diplotype 落在任一候選解就算「一致」。
  - HGVS（c.XXX，如 DPYD）兩邊同系統 → 可比對；star vs HGVS/rs 視為不同系統 → 未比對。

DRAGEN gene → genotype 來源：
  - locusAnnotations[].{gene, genotype}
  - 頂層 cyp2d6 / cyp2b6 的 .genotype（不在 locusAnnotations 內）
  - 結構型（lpa/rh/smn/cyp21a2/gba/hba）無 star/HGVS genotype、我們也不呼叫 → 略過。

相依：只用 Python 標準庫。
"""

import argparse
import json
import re
import sys

REF_TOKENS = ("reference", "wildtype", "wild-type", "ref", "normal")


def is_ref_allele(a):
    """reference/wild-type 型 allele（含 star 的 *1）→ True。"""
    a = a.strip().lower().replace(" ", "")
    if a in ("*1", "1"):
        return True
    return any(t in a for t in REF_TOKENS)


def norm_allele(a):
    a = a.strip().lower().replace(" ", "")
    return "REF" if is_ref_allele(a) else a


def norm_diplo(gt):
    """單一 diplotype 字串 → 排序後的 allele tuple（比對用）；空 → None。"""
    gt = (gt or "").strip()
    parts = [p for p in gt.split("/") if p.strip() != ""]
    if not parts:
        return None
    return tuple(sorted(norm_allele(p) for p in parts))


def dragen_candidates(dragen_gt):
    """DRAGEN 模糊多重解（';' 分隔）→ 多個正規化 diplotype tuple。"""
    return [d for d in (norm_diplo(c) for c in (dragen_gt or "").split(";")) if d is not None]


def nomenclature(gt):
    """粗判命名系統：hgvs（c./>）、rs（rs+數字）、star（*）、ref（其餘 reference-like）。"""
    g = (gt or "").lower()
    if "c." in g or ">" in g:
        return "hgvs"
    if re.search(r"rs\d", g):
        return "rs"
    if "*" in g:
        return "star"
    return "ref"


_INCOMPATIBLE = {
    ("star", "hgvs"), ("hgvs", "star"),
    ("star", "rs"), ("rs", "star"),
    ("hgvs", "rs"), ("rs", "hgvs"),
}


def verdict(our_gt, dragen_gt):
    """回傳 '一致' / '不一致' / '未比對'。"""
    our_n = norm_diplo(our_gt)
    if our_n is not None and any(our_n == c for c in dragen_candidates(dragen_gt)):
        return "一致"
    our_sys = nomenclature(our_gt)
    dra_sys = nomenclature((dragen_gt or "").split(";")[0])
    if (our_sys, dra_sys) in _INCOMPATIBLE:
        return "未比對"
    return "不一致"


def build_dragen_map(js):
    """targeted.json → {GENE(大寫): raw genotype 字串}。"""
    dmap = {}
    for loc in js.get("locusAnnotations", []) or []:
        gene = (loc.get("gene") or "").strip()
        gt = (loc.get("genotype") or "").strip()
        if gene and gt:
            dmap[gene.upper()] = gt
    for key, gene in (("cyp2d6", "CYP2D6"), ("cyp2b6", "CYP2B6")):
        node = js.get(key)
        if isinstance(node, dict):
            gt = (node.get("genotype") or "").strip()
            if gt:
                dmap[gene] = gt
    return dmap


def annotate(pgx_path, json_path, out_path):
    with open(json_path) as f:
        dmap = build_dragen_map(json.load(f))

    with open(pgx_path) as f:
        header = f.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        for req in ("GENE", "NOTES", "DIPLOTYPE"):
            if req not in idx:
                sys.stderr.write("[compare_dragen_pgx] pgx.tsv 缺欄位 %s，原樣輸出\n" % req)
                # 缺欄位就原樣複製，不動
                with open(out_path, "w") as w:
                    w.write("\t".join(header) + "\n")
                    for line in f:
                        w.write(line)
                return
        gi, ni, di = idx["GENE"], idx["NOTES"], idx["DIPLOTYPE"]

        counts = {"一致": 0, "不一致": 0, "未比對": 0, "no_dragen": 0}
        seen_genes = set()
        with open(out_path, "w") as w:
            w.write("\t".join(header) + "\n")
            for line in f:
                row = line.rstrip("\n").split("\t")
                if len(row) <= max(gi, ni, di):
                    w.write(line)
                    continue
                gene = row[gi].strip().upper()
                dragen_gt = dmap.get(gene)
                if dragen_gt:
                    v = verdict(row[di], dragen_gt)
                    tag = "DRAGEN %s: %s" % (v, dragen_gt)
                    old = row[ni].strip()
                    row[ni] = tag if old in (".", "") else "%s; %s" % (old, tag)
                    if gene not in seen_genes:
                        counts[v] += 1
                        seen_genes.add(gene)
                else:
                    if gene and gene not in seen_genes:
                        counts["no_dragen"] += 1
                        seen_genes.add(gene)
                w.write("\t".join(row) + "\n")

    sys.stderr.write(
        "[compare_dragen_pgx] 基因交叉註記：一致=%d 不一致=%d 未比對=%d（DRAGEN 無此基因=%d）\n"
        % (counts["一致"], counts["不一致"], counts["未比對"], counts["no_dragen"])
    )


def main():
    ap = argparse.ArgumentParser(description="Cross-annotate DRAGEN targeted.json PGx calls into pgx.tsv NOTES.")
    ap.add_argument("--pgx", required=True, help="our pgx.tsv")
    ap.add_argument("--dragen-json", required=True, help="DRAGEN {sample}.targeted.json")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--output", required=True, help="augmented pgx.tsv（欄位不變）")
    a = ap.parse_args()
    annotate(a.pgx, a.dragen_json, a.output)


if __name__ == "__main__":
    main()
