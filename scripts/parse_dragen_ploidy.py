#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
=========================================================
WGS/WES Germline Analysis Pipeline - parse_dragen_ploidy.py
=========================================================
Author   : Po-Yu Lin (林伯昱)
Institute: Department of Neurology and
           Department of Genomic Medicine,
           National Cheng Kung University Hospital
Contact  : p88124019@gs.ncku.edu.tw

Copyright (c) 2026, Po-Yu Lin
Licensed under the GNU General Public License v3.0

DISCLAIMER: Provided "as is" without warranty. Users are solely responsible
for validating and interpreting all results.
=========================================================
scripts/parse_dragen_ploidy.py
==============================

讀 DRAGEN 原生 `*.ploidy.vcf.gz`，輸出「與二級 ploidy_check.py 同一套」的 QC 摘要
（<sample>.ploidy_qc.txt），讓 NCKUH 與 DRAGEN 兩條路的 ploidy 呈現一致。warn-only。

⚠️ DRAGEN 的 NDC 語意（實測 chrX 男性 = 1.0 而非 0.5）：
  NDC = 觀測深度 ÷「估計核型下的期望深度」，已對**性別**正規化。所以：
    - 正常樣本每條 contig（含 chrX/chrY）NDC ≈ 1.0。
    - 偏離 1.0 = 非預期（如體染色體三體 ≈1.5；或估計核型與實際不符）。
  故 aneuploidy 判定 = 「NDC 偏離 1.0」（只信 FILTER=PASS 的 contig；chrM 多拷貝跳過）。
  性別直接取 ##estimatedSexKaryotype；與 ##referenceSexKaryotype 不符 → WARN。

（對照：二級 ploidy_check.py 的 NDC = 相對體染色體中位數的**原始**覆蓋比，男性 chrX ≈ 0.5，
  用來「推」性別；兩者 NDC 數字尺度不同，但 qc.txt 的核型/sex_check/WARNINGS 語意一致。）

相依：只用 Python 標準庫。
"""

import argparse
import gzip
import sys

AUTOSOMES = ["chr%d" % i for i in range(1, 23)]
MITO = {"chrM", "chrMT", "MT", "M"}
NDC_LO, NDC_HI = 0.75, 1.25   # NDC 偏離此區間（PASS contig）→ 疑似非整倍體


def _open(path):
    with open(path, "rb") as fh:
        magic = fh.read(2)
    return gzip.open(path, "rt") if magic == b"\x1f\x8b" else open(path, "rt")


def _f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def parse(path):
    """回傳 (estimated, reference, contigs)；contigs = [(chrom, filt, dc, ndc), ...]。"""
    est = ref = None
    contigs = []
    with _open(path) as f:
        for line in f:
            if line.startswith("##"):
                if line.startswith("##estimatedSexKaryotype="):
                    est = line.strip().split("=", 1)[1]
                elif line.startswith("##referenceSexKaryotype="):
                    ref = line.strip().split("=", 1)[1]
                continue
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 10:
                continue
            chrom, filt, fmt, smp = p[0], p[6], p[8], p[9]
            d = dict(zip(fmt.split(":"), smp.split(":")))
            contigs.append((chrom, filt, _f(d.get("DC")), _f(d.get("NDC"))))
    return est, ref, contigs


def analyze(est, ref, contigs):
    warnings = []
    if est and ref and est != ref:
        warnings.append("SEX MISMATCH：DRAGEN estimated %s vs reference %s"
                        "（可能 sample swap 或性染色體 aneuploidy，請人工確認）" % (est, ref))
    flags = []
    for chrom, filt, dc, ndc in contigs:
        if ndc is None or chrom in MITO:
            continue
        if filt not in ("PASS", ".", ""):     # 只信 PASS（如 chr19 常見 LowQual → 不判）
            continue
        if ndc < NDC_LO or ndc > NDC_HI:
            flags.append((chrom, ndc))
            warnings.append("%s NDC=%.2f → 疑似非整倍體（DRAGEN NDC 已對估計核型正規化，"
                            "偏離 1.0 = 非預期）；本染色體 SNV 基因型可能不準，考慮人工確認"
                            % (chrom, ndc))
    return {
        "estimated": est or "unknown",
        "reference": ref or "unknown",
        "contigs": contigs,
        "aneuploidy": flags,
        "warnings": warnings,
    }


def write_qc(path, sample, res):
    order = AUTOSOMES + ["chrX", "chrY", "chrM"]
    ndc_by = {c: ndc for (c, filt, dc, ndc) in res["contigs"]}
    sex_ok = (res["estimated"] in ("unknown",) or res["reference"] in ("unknown",)
              or res["estimated"] == res["reference"])
    with open(path, "w") as w:
        w.write("# Ploidy QC — %s (DRAGEN)\n" % sample)
        w.write("declared_sex_karyotype : %s\n" % res["reference"])
        w.write("estimated_sex_karyotype: %s\n" % res["estimated"])
        w.write("sex_check              : %s\n" % ("OK" if sex_ok else "MISMATCH"))
        w.write("source                 : DRAGEN native ploidy.vcf "
                "(NDC normalized to expected karyotype; ~1.0 = as-expected)\n")
        w.write("\n--- per-chromosome normalized coverage (NDC) ---\n")
        listed = [c for c in order if c in ndc_by] + \
                 [c for c in ndc_by if c not in order]
        for c in listed:
            v = ndc_by[c]
            w.write("%-6s %s\n" % (c, ("%.3f" % v) if v is not None else "NA"))
        w.write("\n--- WARNINGS ---\n")
        if res["warnings"]:
            for msg in res["warnings"]:
                w.write("WARN: %s\n" % msg)
        else:
            w.write("(none)\n")


def main():
    ap = argparse.ArgumentParser(description="Summarize DRAGEN ploidy.vcf into unified QC txt.")
    ap.add_argument("--in", dest="inp", required=True, help="DRAGEN *.ploidy.vcf(.gz)")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--out-qc", required=True)
    a = ap.parse_args()
    est, ref, contigs = parse(a.inp)
    res = analyze(est, ref, contigs)
    write_qc(a.out_qc, a.sample, res)
    for msg in res["warnings"]:
        sys.stderr.write("[parse_dragen_ploidy] WARN: %s\n" % msg)
    sys.stderr.write("[parse_dragen_ploidy] %s estimated=%s reference=%s aneuploidy=%d\n"
                     % (a.sample, res["estimated"], res["reference"], len(res["aneuploidy"])))


if __name__ == "__main__":
    main()
