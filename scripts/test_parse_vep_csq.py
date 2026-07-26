#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unit tests for parse_vep_csq.py helpers — dependency-free (stdlib only).
Run:  python3 scripts/test_parse_vep_csq.py

Currently pins strand_bias_flag() (the STRAND_BIAS review-warning column):
germline flags strand bias for manual review, never hard-filters.
Thresholds follow GATK convention (SNV FS>60/SOR>3.0; indel FS>200/SOR>10.0).
DeepVariant-only sites lack FS/SOR -> "." (manual review).
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import parse_vep_csq as P  # noqa: E402


def _eq(desc, got, exp):
    assert got == exp, "%s: got %r, expected %r" % (desc, got, exp)
    print("PASS", desc, "->", got)


def test_snv_thresholds():
    _eq("SNV biased (FS>60 & SOR>3)", P.strand_bias_flag({"FS": "72", "SOR": "4.1"}, "A", "G"),
        "WARN(FS=72.0,SOR=4.10)")
    _eq("SNV clean", P.strand_bias_flag({"FS": "10", "SOR": "1.5"}, "A", "G"), "PASS")
    _eq("SNV FS-only over cut", P.strand_bias_flag({"FS": "70", "SOR": "1.0"}, "A", "G"),
        "WARN(FS=70.0,SOR=1.00)")
    _eq("SNV SOR-only over cut", P.strand_bias_flag({"FS": "5", "SOR": "3.5"}, "A", "G"),
        "WARN(FS=5.0,SOR=3.50)")


def test_indel_thresholds_more_lenient():
    # indel needs FS>200 / SOR>10; SNV-level bias does NOT trip an indel
    _eq("indel FS100/SOR5 clean", P.strand_bias_flag({"FS": "100", "SOR": "5"}, "A", "AT"), "PASS")
    _eq("indel FS250 biased", P.strand_bias_flag({"FS": "250", "SOR": "5"}, "A", "AT"),
        "WARN(FS=250.0,SOR=5.00)")
    _eq("del (ref len>1) SOR12 biased", P.strand_bias_flag({"FS": "10", "SOR": "12"}, "AT", "A"),
        "WARN(FS=10.0,SOR=12.00)")


def test_no_or_partial_data():
    _eq("DeepVariant-only (no FS/SOR)", P.strand_bias_flag({"CALLERS": "DV"}, "A", "G"), ".")
    _eq("both dotted", P.strand_bias_flag({"FS": ".", "SOR": "."}, "A", "G"), ".")
    _eq("only FS present, biased", P.strand_bias_flag({"FS": "80"}, "A", "G"), "WARN(FS=80.0)")
    _eq("only SOR present, clean", P.strand_bias_flag({"SOR": "2.0"}, "A", "G"), "PASS")
    _eq("unparseable FS, no SOR -> no usable data", P.strand_bias_flag({"FS": "abc"}, "A", "G"), ".")


def test_column_registered():
    assert "STRAND_BIAS" in P.OUTPUT_COLUMNS, "STRAND_BIAS must be an output column"
    print("PASS test_column_registered -> STRAND_BIAS in OUTPUT_COLUMNS")


if __name__ == "__main__":
    test_snv_thresholds()
    test_indel_thresholds_more_lenient()
    test_no_or_partial_data()
    test_column_registered()
    print("\nALL TESTS PASSED")
