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
    # Step 1：確保 .tbi 存在（DRAGEN 輸出有時不含 index）
    if [ ! -f "${dragen_vcf}.tbi" ]; then
        echo "[ADD_DRAGEN_TAG] 建立 tabix index..." >&2
        tabix -p vcf ${dragen_vcf}
    fi

    # Step 2：add_dragen_tag.py
    #   - 新增 INFO tag：CALLERS, DP_DRAGEN, AD_DRAGEN, VAF_DRAGEN, GQ_DRAGEN
    #   - 分流：SNV（PASS，非 chrM）→ snv_raw.vcf；Mito（chrM 全部）→ mito_raw.vcf
    python3 ${params.scripts_dir}/add_dragen_tag.py \\
        --input       ${dragen_vcf} \\
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

workflow PREPARE_VCF_DRAGEN {
    take:
    dragen_ch   // tuple val(sample_id), path(dragen_vcf)（不需要 tbi）

    main:
    ADD_DRAGEN_TAG(dragen_ch)

    emit:
    snv_ch  = ADD_DRAGEN_TAG.out.snv_ch
    mito_ch = ADD_DRAGEN_TAG.out.mito_ch
}
