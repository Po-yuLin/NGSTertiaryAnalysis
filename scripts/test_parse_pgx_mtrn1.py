#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unit tests for parse_pgx_report.parse_mito_tsv() — dependency-free (stdlib only).
Run:  python3 scripts/test_parse_pgx_mtrn1.py

Guards the MT-RNR1 coverage gate (2026-08 fix). The bug: the gate read

    if covered_positions or no_mito is False:

but `no_mito` is necessarily False there (the True case returns earlier), so the
condition was always true and depth never gated anything. Any sample with a
mito.tsv but no MT-RNR1 pathogenic call was reported as
`Reference / MTRN1_RISK=LOW` — "standard aminoglycoside dosing applies" — even at
zero depth. Clinically that is the one direction we must never fail in: an
unmeasured m.1555A>G carrier told gentamicin is safe -> irreversible hearing loss.

Invariant these tests hold: LOW is emitted only for positions bcftools mpileup
actually covered (DP >= MTRN1_MIN_DP); everything else stays Unknown (= no row).
"""
import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import parse_pgx_report as P  # noqa: E402

MITO_HDR = "CHROM\tPOS\tREF\tALT\tAF_SAMPLE\tCLINVAR_SIG\n"

# mpileup shapes: chrM is high copy number, so real samples read in the hundreds
FULL_COV    = {1555: {"dp": 412}, 1494: {"dp": 398}, 827: {"dp": 355}}
NO_COV      = {1555: {"dp": 0},   1494: {"dp": 0},   827: {"dp": 3}}
PARTIAL_COV = {1555: {"dp": 412}, 1494: {"dp": 4},   827: {"dp": 0}}
OFF_TARGET  = {3243: {"dp": 900}}          # depth, but not at a pathogenic site


def _write_mito(rows):
    """rows = [(pos, ref, alt, af, clnsig), ...] — chrM variant calls only."""
    d = tempfile.mkdtemp()
    p = os.path.join(d, "mito.tsv")
    with open(p, "w") as w:
        w.write(MITO_HDR)
        for pos, ref, alt, af, sig in rows:
            w.write("chrM\t%d\t%s\t%s\t%s\t%s\n" % (pos, ref, alt, af, sig))
    return p


def _no_mtrn1_variant():
    # m.3243A>G (MELAS) — a real chrM call that is NOT one of the three PGx sites
    return _write_mito([(3243, "A", "G", "0.31", "Pathogenic")])


def test_full_coverage_reference():
    rows = P.parse_mito_tsv(_no_mtrn1_variant(), FULL_COV)
    assert len(rows) == 1, rows
    r = rows[0]
    assert r["mtrn1_risk"] == "LOW" and r["diplotype"] == "Reference"
    # all three assessed -> all three named, no "not assessed" caveat
    for hgvs in ("m.827A>G", "m.1494C>T", "m.1555A>G"):
        assert hgvs in r["recommendation"], r["recommendation"]
    assert "NOT assessed" not in r["notes"], r["notes"]
    print("PASS test_full_coverage_reference -> LOW, all 3 positions assessed")


def test_no_coverage_stays_unknown():
    # THE REGRESSION: mito.tsv exists, no MT-RNR1 call, but nothing was measured
    rows = P.parse_mito_tsv(_no_mtrn1_variant(), NO_COV)
    assert rows == [], "no-coverage sample must not be reported as Reference/LOW"
    print("PASS test_no_coverage_stays_unknown -> 0 rows (was LOW before the fix)")


def test_partial_coverage_names_only_assessed():
    rows = P.parse_mito_tsv(_no_mtrn1_variant(), PARTIAL_COV)
    assert len(rows) == 1, rows
    r = rows[0]
    assert r["mtrn1_risk"] == "LOW"
    # only 1555 was measured -> the other two must be declared unassessed
    assert "m.1555A>G" in r["recommendation"]
    assert "m.827A>G" in r["recommendation"] and "m.1494C>T" in r["recommendation"]
    assert "could not be assessed" in r["recommendation"], r["recommendation"]
    assert "NOT assessed" in r["notes"] and "827,1494" in r["notes"], r["notes"]
    print("PASS test_partial_coverage_names_only_assessed -> LOW + unassessed flagged")


def test_missing_mtrn1_vcf_stays_unknown():
    assert P.parse_mito_tsv(_no_mtrn1_variant(), {}) == []
    assert P.parse_mito_tsv(_no_mtrn1_variant(), None) == []
    print("PASS test_missing_mtrn1_vcf_stays_unknown -> 0 rows")


def test_off_target_depth_does_not_count():
    # depth at chrM:3243 must not be mistaken for coverage of the three PGx sites
    assert P.parse_mito_tsv(_no_mtrn1_variant(), OFF_TARGET) == []
    print("PASS test_off_target_depth_does_not_count -> 0 rows")


def test_pathogenic_call_is_high_regardless_of_gate():
    mito = _write_mito([(1555, "A", "G", "0.98", "Pathogenic")])
    rows = P.parse_mito_tsv(mito, {})          # no mpileup at all
    assert len(rows) == 1 and rows[0]["mtrn1_risk"] == "HIGH", rows
    assert rows[0]["diplotype"] == "m.1555A>G"
    assert "Avoid aminoglycoside" in rows[0]["recommendation"]
    print("PASS test_pathogenic_call_is_high_regardless_of_gate -> HIGH")


def test_all_three_pathogenic_positions_detected():
    for pos, (hgvs, _) in P.MTRN1_PATHOGENIC.items():
        rows = P.parse_mito_tsv(_write_mito([(pos, "N", "N", "0.9", "Pathogenic")]), FULL_COV)
        assert len(rows) == 1 and rows[0]["mtrn1_risk"] == "HIGH", (pos, rows)
        assert rows[0]["diplotype"] == hgvs
    print("PASS test_all_three_pathogenic_positions_detected -> 1555/1494/827 all HIGH")


def test_no_mito_tsv_stays_unknown():
    assert P.parse_mito_tsv("NO_FILE", FULL_COV) == []
    assert P.parse_mito_tsv("", FULL_COV) == []
    print("PASS test_no_mito_tsv_stays_unknown -> 0 rows")


def test_threshold_boundary():
    mito = _no_mtrn1_variant()
    assert P.parse_mito_tsv(mito, {1555: {"dp": P.MTRN1_MIN_DP - 1}}) == []
    at_threshold = P.parse_mito_tsv(mito, {1555: {"dp": P.MTRN1_MIN_DP}})
    assert len(at_threshold) == 1 and at_threshold[0]["mtrn1_risk"] == "LOW"
    print("PASS test_threshold_boundary -> DP >= %d inclusive" % P.MTRN1_MIN_DP)


if __name__ == "__main__":
    test_full_coverage_reference()
    test_no_coverage_stays_unknown()
    test_partial_coverage_names_only_assessed()
    test_missing_mtrn1_vcf_stays_unknown()
    test_off_target_depth_does_not_count()
    test_pathogenic_call_is_high_regardless_of_gate()
    test_all_three_pathogenic_positions_detected()
    test_no_mito_tsv_stays_unknown()
    test_threshold_boundary()
    print("\nALL TESTS PASSED")
