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

parse_vep_csq.py
================
將 VEP annotation VCF + Pangolin VCF 解析為結構化 TSV，
供後續 acmg_classifier.py 和 GUI 使用。

輸入：
  --vep_vcf         VEP annotation VCF（*.vep.vcf.gz）
  --pangolin_vcf    Pangolin splice score VCF（*.pangolin.vcf.gz）
  --clinvar_lookup  clinvar_lookup.tsv.gz（build_clinvar_lookup.py 產生）
  --sample_id       樣本 ID
  --output_full     完整輸出 TSV（archive 用）
  --output_filtered 過濾輸出 TSV（GUI 用）

輸出欄位：
  見 OUTPUT_COLUMNS

過濾規則（filtered 版本移除同時符合以下所有條件的 variant）：
  - gnomAD genome 或 exome AF > 0.01
  - ClinVar 無注釋
  - VEP IMPACT = MODIFIER
  - Alt contig（_alt / random / chrUn）

作者：Po-Yu Lin（林伯昱）
機構：國立成功大學醫院基因醫學部
"""

import argparse
import gzip
import json
import re
import sys
 
# ──────────────────────────────────────────────────────────────
# Consequence 嚴重程度排序（數字越小越嚴重）
# 用於從 MANE Select transcript 中選最嚴重的 consequence
# ──────────────────────────────────────────────────────────────
 
CONSEQUENCE_RANK = {
    "transcript_ablation":                    1,
    "splice_acceptor_variant":                2,
    "splice_donor_variant":                   3,
    "stop_gained":                            4,
    "frameshift_variant":                     5,
    "stop_lost":                              6,
    "start_lost":                             7,
    "transcript_amplification":               8,
    "inframe_insertion":                      9,
    "inframe_deletion":                      10,
    "missense_variant":                      11,
    "protein_altering_variant":              12,
    "splice_region_variant":                 13,
    "splice_donor_5th_base_variant":         14,
    "splice_donor_region_variant":           15,
    "splice_polypyrimidine_tract_variant":   16,
    "incomplete_terminal_codon_variant":     17,
    "stop_retained_variant":                 18,
    "synonymous_variant":                    19,
    "coding_sequence_variant":               20,
    "mature_miRNA_variant":                  21,
    "5_prime_UTR_variant":                   22,
    "3_prime_UTR_variant":                   23,
    "non_coding_transcript_exon_variant":    24,
    "intron_variant":                        25,
    "NMD_transcript_variant":               26,
    "non_coding_transcript_variant":         27,
    "upstream_gene_variant":                 28,
    "downstream_gene_variant":               29,
    "intergenic_variant":                    38,
}
 
 
def get_worst_consequence_rank(tx: dict) -> int:
    """取 transcript 所有 consequence 中最嚴重的 rank"""
    consequences = tx.get("Consequence", "").split("&")
    ranks = [CONSEQUENCE_RANK.get(c.strip(), 99) for c in consequences if c.strip()]
    return min(ranks) if ranks else 99
 
 
# ──────────────────────────────────────────────────────────────
# ClinVar review status → stars 轉換
# ──────────────────────────────────────────────────────────────
 
CLNREVSTAT_STARS = {
    "practice_guideline":                                    4,
    "reviewed_by_expert_panel":                              3,
    "criteria_provided_multiple_submitters_no_conflicts":    2,
    "criteria_provided_conflicting_classifications":         1,
    "criteria_provided_single_submitter":                    1,
    "no_assertion_criteria_provided":                        0,
    "no_classification_provided":                            0,
    "no_classification_for_the_single_variant":              0,
}
 
 
def clnrevstat_to_stars(revstat: str) -> int:
    if not revstat or revstat == ".":
        return 0
    normalized = revstat.replace("&_", "_").replace("&", "_").lower()
    return CLNREVSTAT_STARS.get(normalized, 0)
 
 
# ──────────────────────────────────────────────────────────────
# Zygosity 推導
# ──────────────────────────────────────────────────────────────
 
def infer_zygosity(gt_dv: str, gt_hc: str, chrom: str) -> str:
    gt = gt_dv if gt_dv not in (".", "./.", ".|.") else gt_hc
    if gt in (".", "./.", ".|."):
        return "unknown"
    gt_norm = gt.replace("|", "/")
    alleles = gt_norm.split("/")
    if len(alleles) != 2:
        return "unknown"
    ref_count = alleles.count("0")
    alt_alleles = [a for a in alleles if a not in ("0", ".")]
    if ref_count == 2:
        return "ref"
    elif ref_count == 0 and len(alt_alleles) == 2:
        if chrom in ("chrX", "chrY", "X", "Y"):
            return "hemizygous"
        return "hom"
    elif ref_count == 1:
        return "het"
    return "unknown"
 
 
# ──────────────────────────────────────────────────────────────
# GT 解析
# ──────────────────────────────────────────────────────────────
 
def parse_gt_field(format_str: str, sample_str: str, field: str) -> str:
    if not format_str or not sample_str or sample_str == ".":
        return "."
    fields = format_str.split(":")
    values = sample_str.split(":")
    if field not in fields:
        return "."
    idx = fields.index(field)
    return values[idx] if idx < len(values) else "."
 
 
# ──────────────────────────────────────────────────────────────
# rsID 提取（從 Existing_variation 欄位）
# ──────────────────────────────────────────────────────────────
 
def extract_rs_id(existing_variation: str) -> str:
    """
    從 VEP Existing_variation 欄位提取 rsID。
    格式如：rs72631890&COSV58989146
    取第一個 rs 開頭的值。
    """
    if not existing_variation or existing_variation == ".":
        return "."
    for item in existing_variation.replace(",", "&").split("&"):
        if item.startswith("rs"):
            return item
    return "."
 
 
# ──────────────────────────────────────────────────────────────
# ClinVar lookup 載入
# ──────────────────────────────────────────────────────────────
 
def load_clinvar_lookup(lookup_path: str) -> dict:
    """
    載入 clinvar_lookup.tsv.gz。
    回傳 {key: (variation_id, omim_ids, rs_id)} dict。
    key 格式：chr{CHROM}:{POS}:{REF}:{ALT}
    """
    lookup = {}
    opener = gzip.open if lookup_path.endswith(".gz") else open
 
    print(f"[parse_vep_csq] 載入 ClinVar lookup：{lookup_path}", file=sys.stderr)
    with opener(lookup_path, "rt") as f:
        header = f.readline()  # 跳過 header
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            key, varid, omim_ids, rs_id = parts[0], parts[1], parts[2], parts[3]
            lookup[key] = (varid, omim_ids, rs_id)
 
    print(f"[parse_vep_csq] ClinVar lookup 載入完成：{len(lookup):,} 筆", file=sys.stderr)
    return lookup
 
 
# ──────────────────────────────────────────────────────────────
# Pangolin VCF 解析
# ──────────────────────────────────────────────────────────────
 
def load_pangolin_scores(pangolin_vcf: str) -> dict:
    scores = {}
    opener = gzip.open if pangolin_vcf.endswith(".gz") else open
 
    with opener(pangolin_vcf, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            chrom, pos, _, ref, alt = parts[0], parts[1], parts[2], parts[3], parts[4]
            info = parts[7]
            pangolin_str = ""
            for field in info.split(";"):
                if field.startswith("Pangolin="):
                    pangolin_str = field[len("Pangolin="):]
                    break
            if not pangolin_str:
                continue
            detail = pangolin_str
            max_score = 0.0
            for seg in pangolin_str.split("|")[1:]:
                if seg.startswith("Warnings"):
                    break
                if ":" in seg:
                    try:
                        score_val = float(seg.split(":")[1])
                        if abs(score_val) > abs(max_score):
                            max_score = score_val
                    except (ValueError, IndexError):
                        pass
            key = (chrom, pos, ref, alt)
            scores[key] = (max_score, detail)
 
    return scores
 
 
# ──────────────────────────────────────────────────────────────
# VCF header 解析：CSQ 欄位順序
# ──────────────────────────────────────────────────────────────
 
def parse_csq_fields(vcf_path: str) -> dict:
    opener = gzip.open if vcf_path.endswith(".gz") else open
    with opener(vcf_path, "rt") as f:
        for line in f:
            if not line.startswith("##"):
                break
            if line.startswith("##INFO=<ID=CSQ"):
                m = re.search(r'Format: ([^"]+)"', line)
                if m:
                    fields = m.group(1).rstrip(">").split("|")
                    return {name: idx for idx, name in enumerate(fields)}
    raise ValueError("找不到 CSQ FORMAT 定義")
 
 
# ──────────────────────────────────────────────────────────────
# 安全取值
# ──────────────────────────────────────────────────────────────
 
def get(tx: dict, field: str) -> str:
    val = tx.get(field, "")
    return val if val else "."
 
 
# ──────────────────────────────────────────────────────────────
# ★ v2.9 新增：代表 transcript 選取邏輯
# 優先從 MANE Select 中選最嚴重 consequence，fallback 到 PICK=1
# ──────────────────────────────────────────────────────────────
 
def pick_representative_transcript(transcripts: list) -> dict:
    """
    選取代表 transcript。
 
    選取優先順序：
    1. MANE Select transcript 中，consequence 最嚴重的
    2. 無 MANE Select → 用 VEP 原生 PICK=1
    3. 無 PICK=1 → 第一個 transcript
    """
    # Step 1：找所有 MANE Select transcript
    mane_select_txs = [tx for tx in transcripts if tx.get("MANE_SELECT")]
 
    if mane_select_txs:
        # 在 MANE Select 中選 consequence 最嚴重的
        return min(mane_select_txs, key=get_worst_consequence_rank)
 
    # Step 2：無 MANE Select，用 VEP PICK=1
    picked = [tx for tx in transcripts if tx.get("PICK") == "1"]
    if picked:
        return picked[0]
 
    # Step 3：最後 fallback
    return transcripts[0] if transcripts else {}
 
 
# ──────────────────────────────────────────────────────────────
# TRANSCRIPT_TYPE 判斷
# ──────────────────────────────────────────────────────────────
 
def get_transcript_type(tx: dict) -> str:
    if tx.get("MANE_SELECT"):
        return "MANE_SELECT"
    elif tx.get("MANE_PLUS_CLINICAL"):
        return "MANE_PLUS_CLINICAL"
    elif tx.get("CANONICAL") == "YES":
        return "CANONICAL"
    return "OTHER"
 
 
# ──────────────────────────────────────────────────────────────
# MANE_ALL JSON 建立
# ──────────────────────────────────────────────────────────────
 
def build_mane_all(transcripts: list) -> str:
    mane_entries = []
    seen_tx = set()
 
    for tx in transcripts:
        mane_select = tx.get("MANE_SELECT", "")
        mane_plus   = tx.get("MANE_PLUS_CLINICAL", "")
        feature     = tx.get("Feature", "")
 
        if feature in seen_tx:
            continue
 
        if mane_select:
            seen_tx.add(feature)
            mane_entries.append({
                "tx":          mane_select,
                "enst":        feature,
                "type":        "MANE_SELECT",
                "consequence": tx.get("Consequence", ""),
                "hgvsc":       tx.get("HGVSc", ""),
                "hgvsp":       tx.get("HGVSp", ""),
                "impact":      tx.get("IMPACT", ""),
            })
 
        if mane_plus and feature not in seen_tx:
            seen_tx.add(feature)
            mane_entries.append({
                "tx":          mane_plus,
                "enst":        feature,
                "type":        "MANE_PLUS_CLINICAL",
                "consequence": tx.get("Consequence", ""),
                "hgvsc":       tx.get("HGVSc", ""),
                "hgvsp":       tx.get("HGVSp", ""),
                "impact":      tx.get("IMPACT", ""),
            })
 
    return json.dumps(mane_entries, ensure_ascii=False) if mane_entries else "[]"
 
 
# ──────────────────────────────────────────────────────────────
# ★ v2.9 新增：過濾邏輯
# ──────────────────────────────────────────────────────────────
 
def should_filter(row: dict) -> bool:
    """
    True = 移除（不納入 filtered GUI 版本）
    四個條件同時成立才移除，任一不符合則保留。
 
    移除條件（AND）：
      1. gnomAD genome 或 exome AF > 0.01
      2. ClinVar 無注釋
      3. VEP IMPACT = MODIFIER
      4. Alt contig
    """
    # 條件 1：高頻率
    def to_float(val):
        try:
            return float(val) if val and val != "." else 0.0
        except ValueError:
            return 0.0
 
    af_too_common = (
        to_float(row.get("GNOMAD_G_AF")) > 0.01 or
        to_float(row.get("GNOMAD_E_AF")) > 0.01
    )
 
    # 條件 2：ClinVar 無注釋
    clinvar_sig = row.get("CLINVAR_SIG", ".")
    clinvar_no_annotation = clinvar_sig in (".", "", None)
 
    # 條件 3：MODIFIER
    is_modifier = row.get("IMPACT") == "MODIFIER"
 
    # 條件 4：Alt contig
    chrom = row.get("CHROM", "")
    is_alt_contig = (
        "_alt" in chrom or
        "random" in chrom or
        "Un_" in chrom or
        chrom.startswith("chrUn")
    )
 
    return af_too_common and clinvar_no_annotation and is_modifier and is_alt_contig
 
 
# ──────────────────────────────────────────────────────────────
# ★ v2.9 新增：P-KNN LLR → ACMG evidence 轉換
# ──────────────────────────────────────────────────────────────
 
def llr_to_evidence(llr_str: str) -> str:
    """
    將 P-KNN LLR 轉換為 ACMG evidence 字串。
 
    LLR 閾值（依 ClinGen SVI Bayesian framework）：
      >= 4   → PP3_Strong
      >= 2   → PP3_Moderate
      >= 1   → PP3_Supporting
      -1~1   → .（無 evidence）
      <= -1  → BP4_Supporting
      <= -2  → BP4_Moderate
      <= -4  → BP4_Strong
    """
    if not llr_str or llr_str == ".":
        return "."
    try:
        llr = float(llr_str)
    except ValueError:
        return "."
 
    if llr >= 4:
        return "PP3_Strong"
    elif llr >= 2:
        return "PP3_Moderate"
    elif llr >= 1:
        return "PP3_Supporting"
    elif llr <= -4:
        return "BP4_Strong"
    elif llr <= -2:
        return "BP4_Moderate"
    elif llr <= -1:
        return "BP4_Supporting"
    else:
        return "."
 
 
# ──────────────────────────────────────────────────────────────
# 輸出欄位定義
# ──────────────────────────────────────────────────────────────
 
OUTPUT_COLUMNS = [
    # 位置資訊
    "CHROM", "POS", "REF", "ALT",
    "RS_ID",                        # ★ v2.9：rsID
    # Transcript 資訊
    "GENE", "TRANSCRIPT", "TRANSCRIPT_TYPE",
    "HGVS_C", "HGVS_P", "CONSEQUENCE", "IMPACT",
    "EXON", "INTRON",
    "MANE_ALL",
    # Caller 資訊
    "CALLERS", "DP_DV", "AD_DV", "VAF_DV", "DP_HC", "AD_HC",
    "ZYGOSITY", "GT_DV", "GT_HC",
    # 族群頻率
    "GNOMAD_G_AF", "GNOMAD_G_EAS_AF",
    "GNOMAD_E_AF", "GNOMAD_E_EAS_AF",
    "GNOMAD_E_AF_DBNSFP", "GNOMAD_E_EAS_AF_DBNSFP",
    "TG_EAS_AF",                    # ★ v2.9：1000 Genomes EAS AF
    # ClinVar
    "CLINVAR_SIG", "CLINVAR_STARS", "CLINVAR_DN", "CLINVAR_SIGCONF",
    "CLINVAR_VARIATION_ID",         # ★ v2.9：ClinVar Variation ID（GUI 自行組 URL）
    # OMIM
    "OMIM_IDS",                     # ★ v2.9：OMIM ID（逗號分隔，GUI 自行組 URL）
    # LOFTEE
    "LOFTEE", "LOFTEE_FILTER", "LOFTEE_FLAGS",
    "LOFTOOL",
    # In silico scores
    "BAYESDEL_NOAF", "BAYESDEL_NOAF_PRED",
    "ALPHAMISSENSE", "ALPHAMISSENSE_PRED",
    "ESM1B", "ESM1B_PRED",
    "VARITY_R",
    "SIFT", "SIFT_PRED",
    "DANN",
    "PHACTBOOST",
    "PHYLOP100",
    "GERP",
    "PKNN_LLR",                     # ★ v2.9：P-KNN log likelihood ratio
    "PKNN_EVIDENCE",                 # ★ v2.9：PP3/BP4 evidence（PP3_Strong/PP3_Moderate/...）
    # Splice
    "PANGOLIN_SCORE", "PANGOLIN_DETAIL",
    # Protein 資訊
    "DOMAINS", "SWISSPROT",
    # Gene identifier
    "HGNC_ID",                      # ★ v3.1：HGNC ID（VEP cache 內建，--symbol 旗標啟用）
]
 
 
# ──────────────────────────────────────────────────────────────
# 主解析流程
# ──────────────────────────────────────────────────────────────
 
def parse_vep_vcf(vep_vcf: str, pangolin_scores: dict,
                  clinvar_lookup: dict, sample_id: str,
                  output_full: str, output_filtered: str):
 
    csq_fields = parse_csq_fields(vep_vcf)
    opener = gzip.open if vep_vcf.endswith(".gz") else open
 
    sample_dv = f"{sample_id}_DV"
    sample_hc = f"{sample_id}_HC"
    col_dv = None
    col_hc = None
 
    written_full = 0
    written_filtered = 0
    skipped = 0
 
    with opener(vep_vcf, "rt") as fin, \
         open(output_full, "w") as fout_full, \
         open(output_filtered, "w") as fout_filtered:
 
        header_line = "\t".join(OUTPUT_COLUMNS) + "\n"
        fout_full.write(header_line)
        fout_filtered.write(header_line)
 
        for line in fin:
            line = line.rstrip("\n")
 
            if line.startswith("#CHROM"):
                cols = line.split("\t")
                if sample_dv in cols:
                    col_dv = cols.index(sample_dv)
                if sample_hc in cols:
                    col_hc = cols.index(sample_hc)
                continue
 
            if line.startswith("#"):
                continue
 
            parts = line.split("\t")
            if len(parts) < 8:
                continue
 
            chrom   = parts[0]
            pos     = parts[1]
            ref     = parts[3]
            alt     = parts[4]
            info    = parts[7]
            fmt     = parts[8] if len(parts) > 8 else ""
            smp_dv  = parts[col_dv] if col_dv and col_dv < len(parts) else "."
            smp_hc  = parts[col_hc] if col_hc and col_hc < len(parts) else "."
 
            # INFO 欄位解析
            info_dict = {}
            for field in info.split(";"):
                if "=" in field:
                    k, v = field.split("=", 1)
                    info_dict[k] = v
 
            callers = info_dict.get("CALLERS", ".")
            dp_dv   = info_dict.get("DP_DV", ".")
            ad_dv   = info_dict.get("AD_DV", ".")
            vaf_dv  = info_dict.get("VAF_DV", ".")
            dp_hc   = info_dict.get("DP_HC", ".")
            ad_hc   = info_dict.get("AD_HC", ".")
 
            gt_dv    = parse_gt_field(fmt, smp_dv, "GT")
            gt_hc    = parse_gt_field(fmt, smp_hc, "GT")
            zygosity = infer_zygosity(gt_dv, gt_hc, chrom)
 
            # ClinVar（VEP custom annotation）
            clinvar_sig     = info_dict.get("ClinVar_CLNSIG", ".")
            clinvar_revstat = info_dict.get("ClinVar_CLNREVSTAT", ".")
            clinvar_dn      = info_dict.get("ClinVar_CLNDN", ".")
            clinvar_sigconf = info_dict.get("ClinVar_CLNSIGCONF", ".")
            clinvar_stars   = clnrevstat_to_stars(clinvar_revstat)
 
            # CSQ 解析
            csq_raw = info_dict.get("CSQ", "")
            if not csq_raw:
                skipped += 1
                continue
 
            transcripts_raw = csq_raw.split(",")
            transcripts = []
            for tx_raw in transcripts_raw:
                vals = tx_raw.split("|")
                while len(vals) < len(csq_fields):
                    vals.append("")
                tx_dict = {name: vals[idx] for name, idx in csq_fields.items()
                           if idx < len(vals)}
                transcripts.append(tx_dict)
 
            # ★ v2.9：從 MANE Select 中選最嚴重 consequence 的代表 transcript
            picked_tx = pick_representative_transcript(transcripts)
 
            # 從代表 transcript 提取欄位
            gene            = get(picked_tx, "SYMBOL")
            transcript      = get(picked_tx, "Feature")
            transcript_type = get_transcript_type(picked_tx)
            hgvs_c          = get(picked_tx, "HGVSc")
            hgvs_p          = get(picked_tx, "HGVSp")
            consequence     = get(picked_tx, "Consequence")
            impact          = get(picked_tx, "IMPACT")
            exon            = get(picked_tx, "EXON")
            intron          = get(picked_tx, "INTRON")
 
            # ★ v2.9：rsID（從 Existing_variation 提取）
            existing_var = get(picked_tx, "Existing_variation")
            rs_id_vep    = extract_rs_id(existing_var)
 
            # 族群頻率
            gnomad_g_af     = get(picked_tx, "gnomADg_AF")
            gnomad_g_eas_af = get(picked_tx, "gnomADg_EAS_AF")
            gnomad_e_af     = get(picked_tx, "gnomADe_AF")
            gnomad_e_eas_af = get(picked_tx, "gnomADe_EAS_AF")
            gnomad_e_af_db      = get(picked_tx, "gnomAD_exomes_AF")
            gnomad_e_eas_af_db  = get(picked_tx, "gnomAD_exomes_EAS_AF")
            tg_eas_af       = get(picked_tx, "EAS_AF")   # ★ v2.9：1000G EAS
 
            # LOFTEE
            loftee        = get(picked_tx, "LoF")
            loftee_filter = get(picked_tx, "LoF_filter")
            loftee_flags  = get(picked_tx, "LoF_flags")
            loftool       = get(picked_tx, "LoFtool")
 
            # In silico scores
            bayesdel_noaf      = get(picked_tx, "BayesDel_noAF_score")
            bayesdel_noaf_pred = get(picked_tx, "BayesDel_noAF_pred")
            alphamissense      = get(picked_tx, "AlphaMissense_score")
            alphamissense_pred = get(picked_tx, "AlphaMissense_pred")
            esm1b              = get(picked_tx, "ESM1b_score")
            esm1b_pred         = get(picked_tx, "ESM1b_pred")
            varity_r           = get(picked_tx, "VARITY_R_score")
            sift               = get(picked_tx, "SIFT_score")
            sift_pred          = get(picked_tx, "SIFT_pred")
            dann               = get(picked_tx, "DANN_score")
            phactboost         = get(picked_tx, "PHACTboost_score")
            phylop100          = get(picked_tx, "phyloP100way_vertebrate")
            gerp               = get(picked_tx, "GERP++_RS")
            pknn_llr           = get(picked_tx, "PKNN_LLR")
            pknn_evidence      = llr_to_evidence(pknn_llr)
            domains            = get(picked_tx, "DOMAINS")
            swissprot          = get(picked_tx, "SWISSPROT")
            hgnc_id            = get(picked_tx, "HGNC_ID")   # ★ v3.1：VEP 內建，--symbol 旗標
 
            # MANE_ALL JSON
            mane_all = build_mane_all(transcripts)
 
            # Pangolin 分數
            alt_first = alt.split(",")[0]
            pang_key = (chrom, pos, ref, alt_first)
            if pang_key in pangolin_scores:
                pang_score, pang_detail = pangolin_scores[pang_key]
                pangolin_score  = f"{pang_score:.4f}"
                pangolin_detail = pang_detail
            else:
                pangolin_score  = "."
                pangolin_detail = "."
 
            # ★ v2.9：ClinVar lookup 查表
            lookup_key = f"{chrom}:{pos}:{ref}:{alt_first}"
            if lookup_key in clinvar_lookup:
                cv_varid, cv_omim, cv_rs = clinvar_lookup[lookup_key]
            else:
                cv_varid, cv_omim, cv_rs = ".", ".", "."
 
            # rsID：優先用 lookup，備用 VEP 的 Existing_variation
            rs_id = cv_rs if cv_rs != "." else rs_id_vep
 
            # 組裝 row dict（供過濾邏輯使用）
            row_dict = {
                "CHROM":                chrom,
                "POS":                  pos,
                "REF":                  ref,
                "ALT":                  alt,
                "RS_ID":                rs_id,
                "GENE":                 gene,
                "TRANSCRIPT":           transcript,
                "TRANSCRIPT_TYPE":      transcript_type,
                "HGVS_C":               hgvs_c,
                "HGVS_P":               hgvs_p,
                "CONSEQUENCE":          consequence,
                "IMPACT":               impact,
                "EXON":                 exon,
                "INTRON":               intron,
                "MANE_ALL":             mane_all,
                "CALLERS":              callers,
                "DP_DV":                dp_dv,
                "AD_DV":                ad_dv,
                "VAF_DV":               vaf_dv,
                "DP_HC":                dp_hc,
                "AD_HC":                ad_hc,
                "ZYGOSITY":             zygosity,
                "GT_DV":                gt_dv,
                "GT_HC":                gt_hc,
                "GNOMAD_G_AF":          gnomad_g_af,
                "GNOMAD_G_EAS_AF":      gnomad_g_eas_af,
                "GNOMAD_E_AF":          gnomad_e_af,
                "GNOMAD_E_EAS_AF":      gnomad_e_eas_af,
                "GNOMAD_E_AF_DBNSFP":   gnomad_e_af_db,
                "GNOMAD_E_EAS_AF_DBNSFP": gnomad_e_eas_af_db,
                "TG_EAS_AF":            tg_eas_af,
                "CLINVAR_SIG":          clinvar_sig,
                "CLINVAR_STARS":        str(clinvar_stars),
                "CLINVAR_DN":           clinvar_dn,
                "CLINVAR_SIGCONF":      clinvar_sigconf,
                "CLINVAR_VARIATION_ID": cv_varid,
                "OMIM_IDS":             cv_omim,
                "LOFTEE":               loftee,
                "LOFTEE_FILTER":        loftee_filter,
                "LOFTEE_FLAGS":         loftee_flags,
                "LOFTOOL":              loftool,
                "BAYESDEL_NOAF":        bayesdel_noaf,
                "BAYESDEL_NOAF_PRED":   bayesdel_noaf_pred,
                "ALPHAMISSENSE":        alphamissense,
                "ALPHAMISSENSE_PRED":   alphamissense_pred,
                "ESM1B":                esm1b,
                "ESM1B_PRED":           esm1b_pred,
                "VARITY_R":             varity_r,
                "SIFT":                 sift,
                "SIFT_PRED":            sift_pred,
                "DANN":                 dann,
                "PHACTBOOST":           phactboost,
                "PHYLOP100":            phylop100,
                "GERP":                 gerp,
                "PKNN_LLR":             pknn_llr,
                "PKNN_EVIDENCE":        pknn_evidence,
                "PANGOLIN_SCORE":       pangolin_score,
                "PANGOLIN_DETAIL":      pangolin_detail,
                "DOMAINS":              domains,
                "SWISSPROT":            swissprot,
                "HGNC_ID":              hgnc_id,   # ★ v3.1
            }
 
            # 輸出 row
            row_str = "\t".join(row_dict[col] for col in OUTPUT_COLUMNS) + "\n"
 
            # full 版本：全部寫入
            fout_full.write(row_str)
            written_full += 1
 
            # filtered 版本：過濾後寫入
            if not should_filter(row_dict):
                fout_filtered.write(row_str)
                written_filtered += 1
 
    filtered_out = written_full - written_filtered
    print(f"[parse_vep_csq] 完成", file=sys.stderr)
    print(f"  full：{written_full:,} variants → {output_full}", file=sys.stderr)
    print(f"  filtered：{written_filtered:,} variants（移除 {filtered_out:,}）→ {output_filtered}",
          file=sys.stderr)
    print(f"  跳過（無 CSQ）：{skipped}", file=sys.stderr)
 
 
# ──────────────────────────────────────────────────────────────
# 主程式
# ──────────────────────────────────────────────────────────────
 
def main():
    parser = argparse.ArgumentParser(
        description="解析 VEP CSQ + Pangolin + ClinVar lookup，輸出結構化 TSV"
    )
    parser.add_argument("--vep_vcf",         required=True)
    parser.add_argument("--pangolin_vcf",     required=True)
    parser.add_argument("--clinvar_lookup",   required=True,
                        help="clinvar_lookup.tsv.gz（build_clinvar_lookup.py 產生）")
    parser.add_argument("--sample_id",        required=True)
    parser.add_argument("--output_full",      required=True,
                        help="完整輸出 TSV（archive 用）")
    parser.add_argument("--output_filtered",  required=True,
                        help="過濾輸出 TSV（GUI 用）")
    args = parser.parse_args()
 
    print(f"[parse_vep_csq] 載入 Pangolin 分數：{args.pangolin_vcf}", file=sys.stderr)
    pangolin_scores = load_pangolin_scores(args.pangolin_vcf)
    print(f"[parse_vep_csq] Pangolin 載入完成：{len(pangolin_scores):,} variants",
          file=sys.stderr)
 
    clinvar_lookup = load_clinvar_lookup(args.clinvar_lookup)
 
    print(f"[parse_vep_csq] 解析 VEP VCF：{args.vep_vcf}", file=sys.stderr)
    parse_vep_vcf(
        args.vep_vcf, pangolin_scores, clinvar_lookup,
        args.sample_id, args.output_full, args.output_filtered
    )
 
 
if __name__ == "__main__":
    main()