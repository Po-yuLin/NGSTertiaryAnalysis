/*
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
 * Licensed under the GNU General Public License v3.0
 * =========================================================
 * modules/str_annotation.nf
 * =========================
 * 目的：
 *   對 STR VCF 進行 threshold 分類，輸出臨床用 STR TSV。
 *   不需要 VEP（STR 不需要 HGVS 或 consequence annotation）。
 *
 * 兩個 process：
 *
 *   Process 1 - STR_PREPARE（NCKUH only）：
 *     GangSTR 輸出的 .str.vcf 是未壓縮格式，
 *     需要先 bgzip + tabix 才能讓 cyvcf2 正常讀取。
 *     DRAGEN 的 .repeats.vcf.gz 已壓縮，直接跳過此步驟。
 *
 *   Process 2 - STR_PARSE：
 *     呼叫 parse_str_vcf.py：
 *       - DRAGEN：用 INFO/VARID 查 str_lookup_varid.tsv.gz
 *       - NCKUH：用 chrom:pos 查 str_lookup_pos.tsv.gz
 *       - 對照 STRchive threshold 分類（normal/intermediate/pathogenic）
 *       - 只輸出有對應 STRchive 記錄的 locus（已知致病 STR）
 *       - 輸出 {SAMPLE_ID}.str.tsv（22 欄）
 *
 * 包含三個獨立 process，由 main_tertiary.nf 依 pipeline_type 條件呼叫：
 *   STR_PREPARE_NCKUH → NCKUH GangSTR bgzip + tabix
 *   STR_PARSE_NCKUH   → NCKUH GangSTR 解析
 *   STR_PARSE_DRAGEN  → DRAGEN ExpansionHunter 解析
 *
 * 使用容器：
 *   tertiary_python_1.0.0.sif（含 bgzip + tabix + cyvcf2）
 *
 * 踩雷記錄：
 *   - GangSTR REPCN 是 Number=2 Integer 型，cyvcf2 回傳 numpy int32 array
 *   - DRAGEN REPCN 是 Number=1 String 型，cyvcf2 回傳 numpy Unicode string（dtype='U*'）
 *   - 兩者不能共用同一套讀取邏輯，parse_str_vcf.py 內分別實作
 */

// ──────────────────────────────────────────────────────────────
// Process 1：NCKUH GangSTR 前處理（bgzip + tabix）
// ──────────────────────────────────────────────────────────────

process STR_PREPARE_NCKUH {

    label 'process_low'

    container "${params.sif_dir}/tertiary_python_1.0.0.sif"

    containerOptions "${params.apptainer_base_opts}"

    input:
    // NCKUH GangSTR 輸出：未壓縮 .str.vcf
    tuple val(sample_id), path(str_vcf)

    output:
    // 輸出壓縮後的 vcf.gz 和 tbi
    tuple val(sample_id),
          path("${sample_id}.str.vcf.gz"),
          path("${sample_id}.str.vcf.gz.tbi"),
          emit: str_prepared_ch

    script:
    """
    echo "[STR_PREPARE_NCKUH] ${sample_id}：bgzip + tabix" >&2

    # GangSTR 輸出未壓縮 VCF，需先壓縮才能讓 cyvcf2 正常讀取
    bgzip -c ${str_vcf} > ${sample_id}.str.vcf.gz
    tabix -p vcf ${sample_id}.str.vcf.gz

    echo "[STR_PREPARE_NCKUH] ${sample_id} 完成" >&2
    bcftools stats ${sample_id}.str.vcf.gz | grep "^SN" >&2
    """
}

// ──────────────────────────────────────────────────────────────
// Process 2a：NCKUH STR 解析
// ──────────────────────────────────────────────────────────────

process STR_PARSE_NCKUH {

    label 'process_low'

    container "${params.sif_dir}/tertiary_python_1.0.0.sif"

    containerOptions "${params.apptainer_base_opts}"

    publishDir "${params.out_dir}/${sample_id}/05_str", mode: 'copy'

    input:
    tuple val(sample_id), path(str_vcf_gz), path(str_tbi)

    output:
    tuple val(sample_id),
          path("${sample_id}.str.tsv"),
          emit: str_tsv_ch

    script:
    """
    echo "[STR_PARSE_NCKUH] ${sample_id}：解析 GangSTR STR VCF" >&2

    python3 ${params.scripts_dir}/parse_str_vcf.py \\
        --vcf      ${str_vcf_gz} \\
        --sample   ${sample_id} \\
        --lookup   ${params.str_lookup_pos} \\
        --pipeline nckuh \\
        --output   ${sample_id}.str.tsv

    echo "[STR_PARSE_NCKUH] ${sample_id} 完成" >&2
    echo "--- STR TSV 總行數（含 header）---" >&2
    wc -l ${sample_id}.str.tsv >&2
    echo "--- pathogenic / intermediate 筆數 ---" >&2
    awk -F'\\t' 'NR>1 && (\$19=="pathogenic" || \$19=="intermediate")' \\
        ${sample_id}.str.tsv | wc -l >&2
    """
}

// ──────────────────────────────────────────────────────────────
// Process 2b：DRAGEN STR 解析
// ──────────────────────────────────────────────────────────────

process STR_PARSE_DRAGEN {

    label 'process_low'

    container "${params.sif_dir}/tertiary_python_1.0.0.sif"

    containerOptions "${params.apptainer_base_opts}"

    publishDir "${params.out_dir}/${sample_id}/05_str", mode: 'copy'

    input:
    // DRAGEN 已壓縮，直接傳入 vcf.gz（不需要 tbi，cyvcf2 可自行處理）
    tuple val(sample_id), path(str_vcf_gz)

    output:
    tuple val(sample_id),
          path("${sample_id}.str.tsv"),
          emit: str_tsv_ch

    script:
    """
    echo "[STR_PARSE_DRAGEN] ${sample_id}：解析 DRAGEN STR VCF" >&2

    python3 ${params.scripts_dir}/parse_str_vcf.py \\
        --vcf      ${str_vcf_gz} \\
        --sample   ${sample_id} \\
        --lookup   ${params.str_lookup_varid} \\
        --pipeline dragen \\
        --output   ${sample_id}.str.tsv

    echo "[STR_PARSE_DRAGEN] ${sample_id} 完成" >&2
    echo "--- STR TSV 總行數（含 header）---" >&2
    wc -l ${sample_id}.str.tsv >&2
    echo "--- pathogenic / intermediate 筆數 ---" >&2
    awk -F'\\t' 'NR>1 && (\$19=="pathogenic" || \$19=="intermediate")' \\
        ${sample_id}.str.tsv | wc -l >&2
    """
}
