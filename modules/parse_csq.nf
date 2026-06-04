/*
 * =========================================================
 * WGS/WES Germline Analysis Pipeline - Parse CSQ Module
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
 * modules/parse_csq.nf
 * ====================
 * 目的：
 *   執行 parse_vep_csq.py，將 VEP annotation VCF 和 Pangolin VCF
 *   解析為結構化 TSV（snv_indel.annotated.tsv），供後續
 *   acmg_classifier.py 和 GUI 使用。
 *
 * 輸入：
 *   vep_ch:      tuple val(sample_id), path(vep_vcf), path(vep_tbi)
 *   pangolin_ch: tuple val(sample_id), path(pangolin_vcf), path(pangolin_tbi)
 *
 * 輸出：
 *   tsv_ch: tuple val(sample_id), path("*.snv_indel.annotated.tsv")
 *
 * 使用的容器：
 *   tertiary_python_1.0.0.sif（只需要 Python 3 標準函式庫）
 */

process PARSE_CSQ {
    label 'process_medium'

    container "${params.sif_dir}/tertiary_python_1.0.0.sif"

    // full TSV 和 filtered TSV 都不 publish：
    //   - full TSV 是 ACMG classifier 的中間輸入，最終產物是 03_acmg/
    //   - filtered TSV 目前只供除錯用，GUI 會直接讀 acmg.tsv
    //   - 需要時去 Nextflow work/ 目錄撈即可
    // publishDir "${params.out_dir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(vep_vcf), path(vep_tbi)
    tuple val(sample_id2), path(pangolin_vcf), path(pangolin_tbi)

    output:
    tuple val(sample_id),
          path("${sample_id}.snv_indel.full.annotated.tsv"),
          emit: full_tsv
    tuple val(sample_id),
          path("${sample_id}.snv_indel.annotated.tsv"),
          emit: filtered_tsv

    script:
    """
    python3 ${params.scripts_dir}/parse_vep_csq.py \\
        --vep_vcf         ${vep_vcf} \\
        --pangolin_vcf    ${pangolin_vcf} \\
        --clinvar_lookup  ${params.clinvar_lookup_tsv} \\
        --sample_id       ${sample_id} \\
        --input_type      ${params.input_type ?: (params.pipeline_type == 'dragen' ? 'dragen' : 'ensemble')} \\
        --output_full     ${sample_id}.snv_indel.full.annotated.tsv \\
        --output_filtered ${sample_id}.snv_indel.annotated.tsv

    echo "[PARSE_CSQ] ${sample_id} 完成" >&2
    FULL=\$(wc -l < ${sample_id}.snv_indel.full.annotated.tsv)
    FILT=\$(wc -l < ${sample_id}.snv_indel.annotated.tsv)
    echo "[PARSE_CSQ] full: \$(( FULL - 1 )) variants, filtered: \$(( FILT - 1 )) variants" >&2
    """
}

// ──────────────────────────────────────────────────────────────
// 組合 workflow（供 main_tertiary.nf 呼叫）
// ──────────────────────────────────────────────────────────────

workflow PARSE_VEP_CSQ {
    take:
    vep_ch
    pangolin_ch

    main:
    PARSE_CSQ(vep_ch, pangolin_ch)

    emit:
    full_tsv_ch     = PARSE_CSQ.out.full_tsv
    filtered_tsv_ch = PARSE_CSQ.out.filtered_tsv
}
