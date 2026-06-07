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
 * modules/mito_annotation.nf
 * ==========================
 * 目的：
 *   對 mito VCF 進行輕量 VEP 註解，再用 parse_mito_vcf.py
 *   整合 gnomAD mito v3.1（CC0）和 ClinVar 資訊，輸出臨床用 mito TSV。
 *
 * 兩個 process：
 *
 *   Process 1 - MITO_VEP：
 *     Step 1：依 pipeline_type 決定是否過濾 PASS
 *       - NCKUH（GATK Mutect2）：只保留 PASS，過濾 weak_evidence 等
 *       - DRAGEN：保留所有 chrM variant，臨床端依 FILTER 欄位自行篩選
 *     Step 2：VEP 115 輕量 annotation
 *       - --symbol / --hgvs / --numbers / --canonical / --biotype
 *       - --af_gnomadg（VEP cache 內建 gnomAD mito AF）
 *       - ClinVar --custom（與 SNV pipeline 共用同一個 ClinVar VCF）
 *       - 不帶 dbNSFP（chrM variants 在 nuclear genome 資料庫裡幾乎沒有資料）
 *       - 不帶 LOFTEE（不適用於 mtDNA）
 *       - 不帶 Pangolin（mito splice 機制不同）
 *
 *   Process 2 - MITO_PARSE：
 *     呼叫 parse_mito_vcf.py：
 *       - 解析 VEP CSQ 欄位
 *       - 查詢 gnomad_mito_lookup.tsv.gz（gnomAD mito v3.1，CC0）
 *       - 讀取 FORMAT/AF（heteroplasmy level）和 FORMAT/DP
 *       - 輸出 {SAMPLE_ID}.mito.tsv（22 欄）
 *
 * 輸入（來自 main_tertiary.nf 的 mito_ch）：
 *   tuple val(sample_id), val(pipeline_type),
 *         path(mito_vcf), path(mito_tbi)
 *
 * 輸出：
 *   mito_tsv_ch → tuple val(sample_id), path("{sample_id}.mito.tsv")
 *
 * 使用容器：
 *   MITO_VEP   → vep_115.sif（含 bcftools，與 SNV pipeline 共用）
 *   MITO_PARSE → tertiary_python_1.0.0.sif（含 cyvcf2）
 *
 * 注意事項：
 *   - NCKUH mito：GATK Mutect2，單一 sample column，FORMAT/AF = heteroplasmy
 *   - DRAGEN mito：已由 add_dragen_tag.py 加 CALLERS=DRAGEN tag，
 *                  FORMAT/AF 同樣是 heteroplasmy level
 *   - 兩者 FORMAT 欄位相同（GT, AD, AF, DP），parse_mito_vcf.py 共用同一邏輯
 */

// ──────────────────────────────────────────────────────────────
// Process 1：過濾 + VEP 輕量 Annotation
// ──────────────────────────────────────────────────────────────

process MITO_VEP {

    label 'process_medium'

    container "${params.sif_dir}/vep_115.sif"

    // apptainer_base_opts 提供 --bind /scratch,/data 基礎掛載
    // mito 不需要 loftee bind，所以只帶基礎選項
    containerOptions "${params.apptainer_base_opts}"

    publishDir "${params.out_dir}/${sample_id}/04_mito", mode: 'copy'

    input:
    tuple val(sample_id), val(pipeline_type), path(mito_vcf), path(mito_tbi)

    output:
    tuple val(sample_id), val(pipeline_type),
          path("${sample_id}.mito.vep.vcf.gz"),
          path("${sample_id}.mito.vep.vcf.gz.tbi"),
          emit: mito_vep_ch

    script:
    """
    echo "[MITO_VEP] ${sample_id}（${pipeline_type}）開始" >&2

    # ── Step 1：過濾 + 去除重複行 ─────────────────────────────
    #
    # NCKUH（GATK Mutect2 mito）：
    #   FILTER 有 PASS / weak_evidence / base_qual / strand_bias 等。
    #   臨床只報告 PASS，先過濾減少 VEP 工作量。
    #
    # DRAGEN：
    #   add_dragen_tag.py 輸出時已保留所有 chrM variant（含 non-PASS）。
    #   保留全部讓 TSV 完整輸出，臨床端可依 FILTER 欄位自行篩選。
    #
    # 兩者都需要 bcftools norm -d exact：
    #   GATK Mutect2 mito pipeline 有時會輸出 CHROM/POS/REF/ALT 完全相同
    #   的重複行（同一位置被多個 evidence group call 到），
    #   不去重的話 TSV 會出現重複的 variant 行。

    if [ "${pipeline_type}" = "nckuh" ]; then
        echo "[MITO_VEP] NCKUH：過濾 PASS + 去除重複行" >&2
        bcftools view -f PASS ${mito_vcf} \
            | bcftools norm -d exact -Oz -o mito_input.vcf.gz
        tabix -p vcf mito_input.vcf.gz
    else
        echo "[MITO_VEP] DRAGEN：保留所有 chrM variant + 去除重複行" >&2
        bcftools norm -d exact ${mito_vcf} -Oz -o mito_input.vcf.gz
        tabix -p vcf mito_input.vcf.gz
    fi

    # 確認輸入 variant 數
    echo "[MITO_VEP] 送 VEP 的 variant 數：" >&2
    bcftools stats mito_input.vcf.gz | grep "^SN" >&2

    # ── Step 2：VEP 輕量 annotation ───────────────────────────
    #   --af_gnomadg   → gnomAD v3 genome AF，其中包含 chrM（gnomAD mito）
    #   --custom ClinVar → 重用 SNV pipeline 的 ClinVar，取得 CLNSIG
    #   不帶 dbNSFP、LOFTEE、Pangolin（不適用 chrM）
    vep \\
        --input_file mito_input.vcf.gz \\
        --output_file ${sample_id}.mito.vep.vcf.gz \\
        --vcf \\
        --compress_output bgzip \\
        \\
        --offline \\
        --cache \\
        --dir_cache ${params.vep_cache} \\
        --dir_plugins /opt/vep/Plugins \\
        --assembly GRCh38 \\
        --fasta ${params.ref_fasta} \\
        --fork ${task.cpus} \\
        \\
        --hgvs \\
        --symbol \\
        --numbers \\
        --canonical \\
        --biotype \\
        \\
        --mane \\
        --flag_pick \\
        --pick_order mane_select,mane_plus_clinical,canonical,appris,tsl,biotype,ccds,rank,length \\
        \\
        --custom file=${params.clinvar},short_name=ClinVar,format=vcf,type=exact,coords=0,fields=CLNSIG%CLNREVSTAT%CLNDN%CLNSIGCONF \\
        \\
        --af_gnomadg \\
        \\
        --force_overwrite \\
        --no_stats \\
        --safe

    # tabix index
    tabix -p vcf ${sample_id}.mito.vep.vcf.gz

    echo "[MITO_VEP] ${sample_id} 完成" >&2
    bcftools stats ${sample_id}.mito.vep.vcf.gz | grep "^SN" >&2
    """
}

// ──────────────────────────────────────────────────────────────
// Process 2：解析 VEP CSQ + 查 gnomAD mito → 輸出 TSV
// ──────────────────────────────────────────────────────────────

process MITO_PARSE {

    label 'process_low'

    container "${params.sif_dir}/tertiary_python_1.0.0.sif"

    publishDir "${params.out_dir}/${sample_id}/04_mito", mode: 'copy'

    input:
    tuple val(sample_id), val(pipeline_type),
          path(mito_vep_vcf), path(mito_vep_tbi)

    output:
    tuple val(sample_id),
          path("${sample_id}.mito.tsv"),
          emit: mito_tsv_ch

    script:
    """
    echo "[MITO_PARSE] ${sample_id}（${pipeline_type}）開始解析" >&2

    python3 ${params.scripts_dir}/parse_mito_vcf.py \\
        --vcf      ${mito_vep_vcf} \\
        --sample   ${sample_id} \\
        --gnomad_mito ${params.gnomad_mito_lookup} \\
        --pipeline ${pipeline_type} \\
        --output   ${sample_id}.mito.tsv

    echo "[MITO_PARSE] ${sample_id} 完成" >&2

    # 輸出統計供確認
    echo "--- mito TSV 總行數（含 header）---" >&2
    wc -l ${sample_id}.mito.tsv >&2
    echo "--- 前 3 行（前 10 欄）---" >&2
    head -3 ${sample_id}.mito.tsv | cut -f1-10 >&2
    """
}

// ──────────────────────────────────────────────────────────────
// 組合 workflow（供 main_tertiary.nf 呼叫）
// ──────────────────────────────────────────────────────────────

workflow MITO_ANNOTATE {

    take:
    // tuple val(sample_id), val(pipeline_type), path(mito_vcf), path(mito_tbi)
    mito_ch

    main:
    // Step 1：過濾 PASS（NCKUH only）+ VEP 輕量 annotation
    MITO_VEP(mito_ch)

    // Step 2：解析 CSQ + 查 gnomAD mito → TSV
    MITO_PARSE(MITO_VEP.out.mito_vep_ch)

    emit:
    // 最終輸出：tuple val(sample_id), path(mito.tsv)
    mito_tsv_ch = MITO_PARSE.out.mito_tsv_ch
}
