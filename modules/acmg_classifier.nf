/*
 * =========================================================
 * modules/acmg_classifier.nf  ── v3.1
 * =========================================================
 * 輸入：
 *   annotated_tsv  ── parse_vep_csq.py 產生的 56 欄 TSV
 *   clingen_hi_tsv ── ClinGen Dosage Sensitivity TSV（PVS1 HI 過濾）
 *   gene_moi_tsv   ── gene_moi.tsv.gz（PM2 MOI 判斷）
 *
 * 兩個資料庫都是 optional：
 *   不傳入 → 退回簡化版（PVS1 只看 LOFTEE HC；PM2 只看 EAS AF）
 */

process ACMG_CLASSIFY {

    label 'process_medium'

    container "${params.sif_dir}/tertiary_python_1.0.0.sif"

    publishDir "${params.out_dir}/${sample_id}/03_acmg", mode: 'copy'

    input:
        tuple val(sample_id), path(annotated_tsv)
        path clingen_hi_tsv   // optional：傳入 file("NO_FILE") 表示不使用
        path gene_moi_tsv     // optional：同上

    output:
        tuple val(sample_id), path("${sample_id}.snv_indel.acmg.tsv"), emit: acmg_tsv

    script:
    def hi_arg  = (clingen_hi_tsv.name != 'NO_FILE') ? "--clingen_hi ${clingen_hi_tsv}" : ""
    def moi_arg = (gene_moi_tsv.name   != 'NO_FILE') ? "--gene_moi ${gene_moi_tsv}"     : ""

    """
    python3 ${params.scripts_dir}/acmg_classifier.py \\
        --input  ${annotated_tsv} \\
        --output ${sample_id}.snv_indel.acmg.tsv \\
        ${hi_arg} \\
        ${moi_arg}
    """
}
