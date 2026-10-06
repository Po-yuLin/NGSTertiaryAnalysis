#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unit tests for parse_vep_csq.py helpers — dependency-free (stdlib only).
Run:  python3 scripts/test_parse_vep_csq.py

Pins strand_bias_flag() (the STRAND_BIAS review-warning column):
germline flags strand bias for manual review, never hard-filters.
Thresholds follow GATK convention (SNV FS>60/SOR>3.0; indel FS>200/SOR>10.0).
DeepVariant-only sites lack FS/SOR -> "." (manual review).

Also pins infer_zygosity(): ZYGOSITY comes from the GT of the caller that actually
called an ALT (same rule as CALLERS), DV first when both did; on chrX/chrY only a
haploid GT is "hemizygous" (a diploid 1/1 there is a female, PAR or unknown-sex call).

And haploid_het_callers(): the HAPLOID_HET review column (male chrX non-PAR calls that
were heterozygous before +fixploidy) keeps only callers that called this record.
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


def test_zygosity_uses_the_caller_that_called():
    # Regression: DV 0/0 used to win over HC's real call -> ZYGOSITY "ref"
    # (VAL55: 15,422 rows; DV and HC called different alleles at one POS).
    _eq("DV 0/0 + HC 1/0 -> het (was ref)", P.infer_zygosity("0/0", "1/0", "chr1"), "het")
    _eq("DV 0|0 + HC 1/1 -> hom (was ref)", P.infer_zygosity("0|0", "1/1", "chr1"), "hom")
    _eq("male chrX DV 0 + HC 1 -> hemizygous (was ref)",
        P.infer_zygosity("0", "1", "chrX"), "hemizygous")
    _eq("half-missing DV ./1 + HC 1/1 -> hom (HC is the caller that called)",
        P.infer_zygosity("./1", "1/1", "chr2"), "hom")


def test_zygosity_unchanged_cases():
    _eq("both called: DV wins", P.infer_zygosity("1/1", "0/1", "chr1"), "hom")
    _eq("DV missing -> HC", P.infer_zygosity("./.", "0|1", "chr17"), "het")
    _eq("DV only", P.infer_zygosity("0/1", "./.", "chr1"), "het")
    _eq("DRAGEN (gt_hc '.')", P.infer_zygosity("0/1", ".", "chr1"), "het")
    _eq("DRAGEN male chrX", P.infer_zygosity("1", ".", "chrX"), "hemizygous")
    _eq("neither called, DV 0/0 -> ref (fallback)", P.infer_zygosity("0/0", "./.", "chr1"), "ref")
    _eq("neither called, all missing -> unknown", P.infer_zygosity("./.", "./.", "chr1"), "unknown")


def test_zygosity_sex_chromosomes_by_ploidy():
    # Regression: any chrX/chrY record with two ALT alleles was "hemizygous", so a female's
    # chrX 1/1 (and a male PAR 1/1) were reported as hemizygous (VAL-10: no chrX hom at all).
    _eq("female chrX 1/1 -> hom (was hemizygous)", P.infer_zygosity("1/1", ".", "chrX"), "hom")
    _eq("female chrX 1|1 -> hom (was hemizygous)", P.infer_zygosity("1|1", ".", "chrX"), "hom")
    _eq("male PAR 1/1 (diploid) -> hom", P.infer_zygosity("1/1", "1/1", "chrX"), "hom")
    _eq("female chrX 0/1 -> het", P.infer_zygosity("0/1", ".", "chrX"), "het")
    _eq("chrX 1/2 -> het (was hemizygous)", P.infer_zygosity("1/2", ".", "chrX"), "het")
    _eq("chrX half-missing 1/. -> het (was hemizygous)", P.infer_zygosity("1/.", ".", "chrX"), "het")
    _eq("male non-PAR haploid 1 -> hemizygous", P.infer_zygosity("1", "1", "chrX"), "hemizygous")
    _eq("male chrY haploid 1 -> hemizygous", P.infer_zygosity("1", ".", "chrY"), "hemizygous")
    _eq("haploid 0 with no ALT anywhere -> ref (fallback)", P.infer_zygosity("0", ".", "chrX"), "ref")


def test_haploid_het_column():
    # INFO/HAPLOID_HET (secondary haploid_het.awk) narrowed to the callers that called this record
    _eq("both flagged, both called", P.haploid_het_callers("DV,HC", "DV+HC"), "DV,HC")
    _eq("both flagged, only DV called this allele (split multiallelic)",
        P.haploid_het_callers("DV,HC", "DV"), "DV")
    _eq("flagged caller did not call this allele", P.haploid_het_callers("HC", "DV"), ".")
    _eq("no flag", P.haploid_het_callers(".", "DV+HC"), ".")
    _eq("DRAGEN (no tag)", P.haploid_het_callers(".", "DRAGEN"), ".")
    # appended last: existing column positions (index-based readers, GUI) unchanged
    assert P.OUTPUT_COLUMNS[-1] == "HAPLOID_HET", P.OUTPUT_COLUMNS[-3:]
    assert P.OUTPUT_COLUMNS[-2] == "PROTEIN_POSITION", P.OUTPUT_COLUMNS[-3:]
    print("PASS test_haploid_het_column -> HAPLOID_HET is the last output column")


if __name__ == "__main__":
    test_snv_thresholds()
    test_indel_thresholds_more_lenient()
    test_no_or_partial_data()
    test_column_registered()
    test_zygosity_uses_the_caller_that_called()
    test_zygosity_unchanged_cases()
    test_zygosity_sex_chromosomes_by_ploidy()
    test_haploid_het_column()
    print("\nALL TESTS PASSED")
