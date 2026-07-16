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
 * modules/prepare_vcf_dragen.nf
 * =============================
 * DRAGEN hard-filtered VCF 前處理。
 *
 * 輸入：tuple(sample_id, dragen_vcf)（不需要 tbi）
 * 輸出：
 *   snv_ch  → PASS SNV/indel（非 chrM），bgzip + tabix
 *   mito_ch → 所有 chrM variant（含 non-PASS），bgzip + tabix
 *
 * 使用容器：tertiary_python_1.0.0.sif
 *   需包含：cyvcf2、bcftools、bgzip、tabix
 *
 * Script 執行順序：
 *   1. tabix：若 .tbi 不存在則自動建立
 *   2. add_dragen_tag.py：加 INFO tag，分流 SNV / Mito
 *   3. bgzip + tabix：壓縮並建立 index
 *   4. bcftools stats：輸出統計
 */

process ADD_DRAGEN_TAG {

    label 'process_medium'

    container "${params.sif_dir}/tertiary_python_1.0.0.sif"

    publishDir "${params.out_dir}/${sample_id}/00_prepare", mode: 'copy'

    input:
    tuple val(sample_id), path(dragen_vcf)

    output:
    tuple val(sample_id),
          path("${sample_id}.snv_for_annotation.vcf.gz"),
          path("${sample_id}.snv_for_annotation.vcf.gz.tbi"),
          emit: snv_ch

    tuple val(sample_id),
          path("${sample_id}.mito_for_annotation.vcf.gz"),
          path("${sample_id}.mito_for_annotation.vcf.gz.tbi"),
          emit: mito_ch

    script:
    """
    # Step 1：正規化 + 拆分多等位基因（left-align + split multiallelics）
    #   讓每個 ALT 各自一列，使 AD/VAF 與後續 VEP / gnomAD / ClinVar 註解 per-allele 正確。
    #   DRAGEN VCF 多半已正規化，再跑一次為冪等、無害。
    #   -c w：REF 與參考不符時只警告不中斷。
    bcftools norm -m -any -f ${params.ref_fasta} -c w \\
        ${dragen_vcf} \\
        -Oz -o ${sample_id}.norm.vcf.gz
    tabix -p vcf ${sample_id}.norm.vcf.gz

    # Step 2：add_dragen_tag.py（吃正規化後的 biallelic VCF）
    #   - 新增 INFO tag：CALLERS, DP_DRAGEN, AD_DRAGEN, VAF_DRAGEN, GQ_DRAGEN
    #   - 分流：SNV（PASS，非 chrM）→ snv_raw.vcf；Mito（chrM 全部）→ mito_raw.vcf
    python3 ${params.scripts_dir}/add_dragen_tag.py \\
        --input       ${sample_id}.norm.vcf.gz \\
        --sample      ${sample_id} \\
        --output_snv  ${sample_id}.snv_raw.vcf \\
        --output_mito ${sample_id}.mito_raw.vcf

    # Step 3：bgzip + tabix
    bgzip -c ${sample_id}.snv_raw.vcf  > ${sample_id}.snv_for_annotation.vcf.gz
    tabix -p vcf ${sample_id}.snv_for_annotation.vcf.gz

    bgzip -c ${sample_id}.mito_raw.vcf > ${sample_id}.mito_for_annotation.vcf.gz
    tabix -p vcf ${sample_id}.mito_for_annotation.vcf.gz

    # Step 4：統計
    echo "[ADD_DRAGEN_TAG] ${sample_id} SNV stats：" >&2
    bcftools stats ${sample_id}.snv_for_annotation.vcf.gz | grep "^SN" >&2

    echo "[ADD_DRAGEN_TAG] ${sample_id} Mito stats：" >&2
    bcftools stats ${sample_id}.mito_for_annotation.vcf.gz | grep "^SN" >&2

    # 清理暫時檔
    rm -f ${sample_id}.snv_raw.vcf ${sample_id}.mito_raw.vcf
    """
}

// ──────────────────────────────────────────────────────────────
// (選用) compound 合成：把相鄰/重疊的 cis 變異合成單一 canonical MNV
//   DRAGEN 自帶 PS（phase set），直接拿它做 combine，不需 whatshap（不同於二級 NCKUH
//   要先 whatshap）。這樣即使該版 DRAGEN 把 compound 拆開（如 SUZ12 delAAAinsTT），
//   也能在進 VEP 前合回一筆，讓 HGVS p. 正確（p.Glu723_Thr724delinsAla）。
//   combine_phased.py 與二級同一支（scripts/combine_phased.py，只用 Python 標準庫）；
//   在 ADD_DRAGEN_TAG 的 norm -m -any「之前」做（此時仍是 DRAGEN 原始表示、帶 PS）。
//   由 params.combine_phased 開關（預設 true）。chrM 多半無 PS，實質不受影響。
// ──────────────────────────────────────────────────────────────
process COMBINE_DRAGEN {

    label 'process_low'

    container "${params.sif_dir}/tertiary_python_1.0.0.sif"

    input:
    tuple val(sample_id), path(dragen_vcf)

    output:
    tuple val(sample_id), path("${sample_id}.dragen.combined.vcf.gz"), emit: vcf

    script:
    """
    # combine_phased.py 只用 Python 標準庫，讀 (bgzip) VCF、自帶 faidx（讀 \${ref_fasta}.fai）。
    python3 ${params.scripts_dir}/combine_phased.py \\
        --in ${dragen_vcf} \\
        --out ${sample_id}.dragen.combined.vcf \\
        --fasta ${params.ref_fasta} \\
        --max-gap ${params.combine_max_gap}
    # 用 bcftools 排序 + 索引（-Oz 自帶 bgzip、index -t 免 tabix）；tertiary_python 已含 bcftools。
    bcftools sort ${sample_id}.dragen.combined.vcf \\
        -Oz -o ${sample_id}.dragen.combined.vcf.gz
    bcftools index -t ${sample_id}.dragen.combined.vcf.gz
    rm -f ${sample_id}.dragen.combined.vcf
    """
}

workflow PREPARE_VCF_DRAGEN {
    take:
    dragen_ch   // tuple val(sample_id), path(dragen_vcf)（不需要 tbi）

    main:
    // 先 combine（用 DRAGEN 原生 PS）再進 norm/tag；--combine_phased false 可關閉。
    if (params.combine_phased) {
        COMBINE_DRAGEN(dragen_ch)
        ADD_DRAGEN_TAG(COMBINE_DRAGEN.out.vcf)
    } else {
        ADD_DRAGEN_TAG(dragen_ch)
    }

    emit:
    snv_ch  = ADD_DRAGEN_TAG.out.snv_ch
    mito_ch = ADD_DRAGEN_TAG.out.mito_ch
}
