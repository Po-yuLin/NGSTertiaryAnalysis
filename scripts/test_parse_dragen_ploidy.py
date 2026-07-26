#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unit tests for parse_dragen_ploidy.py — dependency-free (stdlib only).
Run:  python3 scripts/test_parse_dragen_ploidy.py

Uses real DRAGEN ploidy.vcf shapes: DRAGEN's NDC is normalized to the estimated
karyotype (a normal male's chrX/chrY read ~1.0, not 0.5), so aneuploidy = NDC
deviating from 1.0 on PASS contigs; sex comes from ##estimatedSexKaryotype.
"""
import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import parse_dragen_ploidy as D  # noqa: E402

HDR = ("##fileformat=VCFv4.2\n"
       "##source=DRAGEN_PLOIDY\n"
       "##estimatedSexKaryotype=%s\n"
       "##referenceSexKaryotype=%s\n"
       "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")


def _write(est, ref, rows):
    """rows = [(chrom, filt, dc, ndc), ...]."""
    d = tempfile.mkdtemp()
    p = os.path.join(d, "s.ploidy.vcf")
    with open(p, "w") as w:
        w.write(HDR % (est, ref))
        for chrom, filt, dc, ndc in rows:
            w.write("%s\t1\t.\tN\t.\t30\t%s\tEND=1\tDC:NDC\t%s:%s\n"
                    % (chrom, filt, dc, ndc))
    return p


def _normal_autosomes(ndc=1.0):
    return [("chr%d" % i, "PASS", 43.0, ndc) for i in range(1, 23)]


def test_male_normal():
    # real 26T00001 shape: chrX/chrY NDC ~1.0 (sex-normalized) -> no aneuploidy
    rows = _normal_autosomes() + [("chrX", "PASS", 21.8, 1.018),
                                  ("chrY", "PASS", 20.9, 0.979)]
    est, ref, contigs = D.parse(_write("XY", "XY", rows))
    res = D.analyze(est, ref, contigs)
    assert res["estimated"] == "XY" and res["reference"] == "XY"
    assert res["aneuploidy"] == [], res["aneuploidy"]
    assert res["warnings"] == [], res["warnings"]
    print("PASS test_male_normal -> XY, no aneuploidy (chrX/chrY NDC~1.0)")


def test_sex_mismatch():
    rows = _normal_autosomes() + [("chrX", "PASS", 43.0, 1.0)]
    est, ref, contigs = D.parse(_write("XY", "XX", rows))
    res = D.analyze(est, ref, contigs)
    assert any("SEX MISMATCH" in w for w in res["warnings"]), res["warnings"]
    print("PASS test_sex_mismatch -> flagged")


def test_trisomy21_flagged():
    rows = _normal_autosomes()
    rows[20] = ("chr21", "PASS", 64.0, 1.50)     # chr21 NDC 1.5
    est, ref, contigs = D.parse(_write("XX", "XX", rows))
    res = D.analyze(est, ref, contigs)
    assert any(c == "chr21" for c, _ in res["aneuploidy"]), res["aneuploidy"]
    assert any("chr21" in w and "非整倍體" in w for w in res["warnings"])
    print("PASS test_trisomy21_flagged -> chr21 NDC 1.5 flagged")


def test_lowqual_not_flagged():
    # a LowQual contig with a wild NDC must NOT be flagged (untrusted)
    rows = _normal_autosomes()
    rows[18] = ("chr19", "LowQual", 60.0, 1.60)  # deviates but LowQual
    est, ref, contigs = D.parse(_write("XX", "XX", rows))
    res = D.analyze(est, ref, contigs)
    assert res["aneuploidy"] == [], res["aneuploidy"]
    print("PASS test_lowqual_not_flagged -> LowQual contig skipped")


def test_mito_skipped():
    rows = _normal_autosomes() + [("chrM", "PASS", 3000.0, 60.0)]  # huge NDC, must skip
    est, ref, contigs = D.parse(_write("XX", "XX", rows))
    res = D.analyze(est, ref, contigs)
    assert res["aneuploidy"] == [], res["aneuploidy"]
    print("PASS test_mito_skipped -> chrM not flagged")


def test_qc_written():
    rows = _normal_autosomes() + [("chrX", "PASS", 21.8, 1.018), ("chrY", "PASS", 20.9, 0.979)]
    p = _write("XY", "XY", rows)
    est, ref, contigs = D.parse(p)
    res = D.analyze(est, ref, contigs)
    out = os.path.join(tempfile.mkdtemp(), "qc.txt")
    D.write_qc(out, "26T00001", res)
    txt = open(out).read()
    assert "estimated_sex_karyotype: XY" in txt
    assert "sex_check              : OK" in txt
    assert "chrX" in txt and "chrY" in txt
    assert "--- WARNINGS ---" in txt
    print("PASS test_qc_written -> unified qc.txt schema")


if __name__ == "__main__":
    test_male_normal()
    test_sex_mismatch()
    test_trisomy21_flagged()
    test_lowqual_not_flagged()
    test_mito_skipped()
    test_qc_written()
    print("\nALL TESTS PASSED")
