#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unit tests for add_callers_tag.py — dependency-free (stdlib only; cyvcf2 is stubbed
when absent, since these tests only exercise pure functions with fake records).
Run:  python3 scripts/test_add_callers_tag.py

Guards two bugs found at SUZ12 chr17:31998950 (2026-09).

1. determine_callers() read

       if dv_called and hc_called: "DV+HC"
       elif dv_called:             "DV"
       else:                       "HC"      # also reached when NEITHER called

   so every record with no ALT genotype in either caller was tagged "HC", passed
   FILTER_FOR_ANNOTATION and reached the ACMG table as a ZYGOSITY=ref/unknown row.
   At SUZ12 the ensemble merge had fused DV's rejected candidate GAAA>GAA with HC's
   real call GAAA>GTT; after tertiary norm the rejected allele was DV ./. + HC 0|0,
   and the report gained a bogus c.2170del — plus the two rejected A>T components
   at 952/953, which HC had no record for at all. The same else also inflated the
   "HC only" stderr count that the old prepare_vcf.nf comment quoted as 23.9%.

2. get_ad() rewrote each missing AD element as 0, turning "this caller has no data
   for this allele" (10,.) into "this caller looked and saw 0 supporting reads"
   (10,0). DV never evaluated the delinsTT allele, so a true het frameshift was
   displayed with VAF 0.

Invariants held here: CALLERS is NONE unless a caller has a fully called ALT
genotype; NONE is outside the set FILTER_FOR_ANNOTATION accepts (read from the real
prepare_vcf.nf, so a vocabulary change there breaks this test); missing AD elements
stay ".", real zeros stay 0.
"""
import os
import re
import sys
import types

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

try:                                    # real cyvcf2 if installed; stub otherwise
    import cyvcf2  # noqa: F401
except ImportError:
    _stub = types.ModuleType("cyvcf2")
    _stub.VCF = _stub.Writer = object
    sys.modules["cyvcf2"] = _stub

import add_callers_tag as A  # noqa: E402

MISS = -2147483648                      # cyvcf2's integer missing value
DV, HC = 0, 1                           # ensemble sample order: DV first, HC second


class FakeVariant:
    """Minimal stand-in for a cyvcf2.Variant: genotypes + format()."""

    def __init__(self, genotypes=None, fmt=None):
        self.genotypes = genotypes or []
        self._fmt = fmt or {}

    def format(self, key):
        return self._fmt.get(key)


def gt(a1, a2, phased=False):
    """cyvcf2 genotype entry: [allele1, allele2, phased]; -1 = missing allele."""
    return [a1, a2, phased]


NOCALL = gt(-1, -1)                     # ./.


def accepted_callers():
    """CALLERS values FILTER_FOR_ANNOTATION lets through, parsed from the real .nf."""
    nf = os.path.join(HERE, "..", "modules", "prepare_vcf.nf")
    with open(nf, encoding="utf-8") as fh:
        text = fh.read()
    m = re.search(r"-i\s+'([^']*INFO/CALLERS[^']*)'", text)
    assert m, "FILTER_FOR_ANNOTATION -i expression not found in prepare_vcf.nf"
    return set(re.findall(r'INFO/CALLERS="([^"]+)"', m.group(1)))


# ── is_called ────────────────────────────────────────────────────
def test_is_called():
    called = [gt(0, 1), gt(0, 1, True), gt(1, 1), gt(1, 0), gt(0, 2, True)]
    not_called = [gt(0, 0), gt(0, 0, True), NOCALL,
                  gt(-1, 1), gt(1, -1)]            # half-missing = not called
    for g in called:
        assert A.is_called(g) is True, g
    for g in not_called:
        assert A.is_called(g) is False, g
    assert A.is_called(None) is False
    print("PASS test_is_called -> ALT genotypes called; ./., 0/0, half-missing not")


# ── determine_callers ────────────────────────────────────────────
def _callers(dv_gt, hc_gt):
    return A.determine_callers(FakeVariant([dv_gt, hc_gt]), DV, HC)


def test_suz12_rejected_allele_is_none():
    # THE REGRESSION: the c.2170del split record — DV ./. + HC 0|0 (HC's 0|2 split)
    assert _callers(NOCALL, gt(0, 0, True)) == "NONE"
    print("PASS test_suz12_rejected_allele_is_none -> DV ./. + HC 0|0 is NONE (was HC)")


def test_neither_called_is_none():
    # 952/953 A>T: DV ./. + HC has no record at all -> ./.
    assert _callers(NOCALL, NOCALL) == "NONE"
    # classic DV RefCall 0/0 with nothing from HC
    assert _callers(gt(0, 0), NOCALL) == "NONE"
    print("PASS test_neither_called_is_none -> ./.+./. and 0/0+./. are NONE (were HC)")


def test_real_calls_unchanged():
    assert _callers(NOCALL, gt(0, 1, True)) == "HC"        # SUZ12 delinsTT
    assert _callers(gt(0, 1), gt(0, 1)) == "DV+HC"
    assert _callers(gt(0, 1), NOCALL) == "DV"
    assert _callers(gt(1, 1), gt(0, 0)) == "DV"
    print("PASS test_real_calls_unchanged -> HC / DV+HC / DV still classified as before")


def test_none_excluded_by_real_filter():
    acc = accepted_callers()
    assert acc == {"DV+HC", "DV", "HC"}, acc
    assert "NONE" not in acc
    print("PASS test_none_excluded_by_real_filter -> prepare_vcf.nf accepts %s, not NONE"
          % sorted(acc))


# ── get_ad ───────────────────────────────────────────────────────
def _ad(rows, idx=0):
    return A.get_ad(FakeVariant(fmt={"AD": rows}), idx)


def test_get_ad_keeps_missing():
    assert _ad([[10, MISS]]) == "10,."            # SUZ12 delinsTT, DV (was "10,0")
    assert _ad([[3, MISS]]) == "3,."              # SUZ12 c.2170del, HC (was "3,0")
    assert _ad([[10, 14, MISS]]) == "10,14,."     # multiallelic, one allele missing
    assert _ad([[1, 2], [10, MISS]], idx=1) == "10,."
    print("PASS test_get_ad_keeps_missing -> 10,. stays 10,. (was 10,0)")


def test_get_ad_real_values_unchanged():
    assert _ad([[3, 2]]) == "3,2"                 # SUZ12 delinsTT, HC
    assert _ad([[30, 0]]) == "30,0"               # a real zero must stay zero
    assert _ad([[MISS, MISS]]) == "."             # all missing -> single dot
    assert A.get_ad(FakeVariant(), 0) == "."      # no AD field at all
    print("PASS test_get_ad_real_values_unchanged -> 3,2 / 30,0 / . as before")


# ── the SUZ12 locus as tertiary sees it after norm -m -any ───────
def test_suz12_locus_end_to_end():
    """
    The four records at chr17:31998950-31998953 before the secondary F4 fix,
    exactly as reproduced from the ensemble the lab pasted. Only the delinsTT
    may reach the ACMG table, and it must not claim DV saw 0 alt reads.
    """
    recs = {
        "950 GA>G (c.2170del)":   (NOCALL, gt(0, 0, True), [[10, 14]], [[3, MISS]]),
        "951 AAA>TT (delinsTT)":  (NOCALL, gt(0, 1, True), [[10, MISS]], [[3, 2]]),
        "952 A>T":                (NOCALL, NOCALL,         [[11, 14]], None),
        "953 A>T":                (NOCALL, NOCALL,         [[10, 15]], None),
    }
    acc = accepted_callers()
    passed = {}
    for name, (dg, hg, dv_ad, hc_ad) in recs.items():
        callers = _callers(dg, hg)
        if callers in acc:
            passed[name] = (callers,
                            _ad(dv_ad),
                            _ad(hc_ad) if hc_ad else ".")
    assert list(passed) == ["951 AAA>TT (delinsTT)"], passed
    assert passed["951 AAA>TT (delinsTT)"] == ("HC", "10,.", "3,2"), passed
    print("PASS test_suz12_locus_end_to_end -> 4 records in, only delinsTT out "
          "(CALLERS=HC, AD_DV=10,., AD_HC=3,2); was all 4")


if __name__ == "__main__":
    test_is_called()
    test_suz12_rejected_allele_is_none()
    test_neither_called_is_none()
    test_real_calls_unchanged()
    test_none_excluded_by_real_filter()
    test_get_ad_keeps_missing()
    test_get_ad_real_values_unchanged()
    test_suz12_locus_end_to_end()
    print("\nALL TESTS PASSED")
