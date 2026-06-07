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
 * modules/cnv_sv_annotation.nf
 * ============================
 * 目的：
 *   對 CNV 和 SV VCF 執行 AnnotSV annotation，
 *   輸出臨床用 CNV/SV TSV。
 *
 * 輸入來源對照：
 *
 *   NCKUH WES：
 *     CNV：{input_dir}/05_cnv_sv/{sample_id}.gcnv.vcf.gz（GATK gCNV）
 *     SV： {input_dir}/05_cnv_sv/{sample_id}.delly.vcf.gz（Delly）
 *
 *   NCKUH WGS：
 *     CNV：{input_dir}/05_cnv_sv/{sample_id}.call.cns（CNVkit，TSV 格式，需先轉 BED）
 *     SV： {input_dir}/05_cnv_sv/{sample_id}.delly.vcf.gz（Delly）
 *
 *   DRAGEN：
 *     CNV：{input_dir}/vcf.gz/{sample_id}.cnv.vcf.gz
 *     SV： {input_dir}/vcf.gz/{sample_id}.sv.vcf.gz
 *     （不用 cnv_sv.vcf.gz，避免 CNV 和 SV event 重疊）
 *
 * 包含的 Process：
 *
 *   CNVKIT_TO_BED（NCKUH WGS only）：
 *     把 CNVkit .call.cns 轉成 AnnotSV 可讀的 BED 格式
 *     輸出：{sample_id}.cnvkit.bed
 *
 *   ANNOTSV_CNV_NCKUH_WES：NCKUH WES gCNV → AnnotSV TSV
 *   ANNOTSV_CNV_NCKUH_WGS：NCKUH WGS CNVkit BED → AnnotSV TSV
 *   ANNOTSV_CNV_DRAGEN：   DRAGEN CNV VCF → AnnotSV TSV
 *   ANNOTSV_SV_NCKUH：     NCKUH Delly VCF → AnnotSV TSV（WES + WGS 共用）
 *   ANNOTSV_SV_DRAGEN：    DRAGEN SV VCF → AnnotSV TSV
 *
 * 輸出：
 *   {out_dir}/{sample_id}/06_cnv_sv/
 *     {sample_id}.cnv.annotated.tsv     ← CNV AnnotSV 輸出（both 模式）
 *     {sample_id}.cnv.unannotated.tsv   ← CNV 未標注的 variant
 *     {sample_id}.sv.annotated.tsv      ← SV AnnotSV 輸出（both 模式）
 *     {sample_id}.sv.unannotated.tsv    ← SV 未標注的 variant
 *
 * 使用容器：
 *   annotsv_3.5.10.sif（含 AnnotSV + bedtools + bcftools）
 *
 * 踩雷記錄：
 *   - annotationsDir 的路徑要指向 share/AnnotSV/（不是上層目錄）
 *   - AnnotSV 在容器內跑時 /tmp/spring.log 需要寫入權限（Exomiser 用）
 *     我們不用 Exomiser，但 AnnotSV 還是會嘗試寫，加 --bind /tmp 即可
 *   - SVminSize 預設已是 50bp，不需要額外設定
 *   - -overlap 參數只影響 user custom BED，不影響 pathogenic/benign SV overlap
 */

// ──────────────────────────────────────────────────────────────
// 輔助函式：AnnotSV 核心指令（所有 process 共用）
// ──────────────────────────────────────────────────────────────
// AnnotSV 參數說明：
//   -SVinputInfo 1  → 保留原始 VCF 的 REF/ALT/FORMAT/INFO 欄位進 TSV
//   -annotationMode both（預設）→ 輸出 full 行（整體 SV）+ split 行（per gene）
//   -genomeBuild GRCh38 → 使用 hg38 annotation databases

// ──────────────────────────────────────────────────────────────
// Process 1：CNVkit .call.cns → BED（NCKUH WGS only）
// ──────────────────────────────────────────────────────────────

process CNVKIT_TO_BED {

    label 'process_low'

    // tertiary_python sif 有 python3，用來做格式轉換
    container "${params.sif_dir}/tertiary_python_1.0.0.sif"

    containerOptions "${params.apptainer_base_opts}"

    input:
    tuple val(sample_id), path(call_cns)

    output:
    tuple val(sample_id), path("${sample_id}.cnvkit.bed"),
          emit: cnvkit_bed_ch

    script:
    """
    echo "[CNVKIT_TO_BED] ${sample_id}：轉換 CNVkit .call.cns → BED" >&2

    # CNVkit .call.cns 格式（tab 分隔）：
    #   chromosome  start  end  gene  log2  baf  cn  cn1  cn2  depth  probes  weight
    # 第一行是 header，start 是 0-based（BED 格式）
    # AnnotSV BED 需要：chrom  start  end  SVTYPE
    #   SVTYPE：cn < 2 → DEL；cn > 2 → DUP；cn == 2 → 跳過（copy-neutral）

    python3 - << 'EOF'
import sys

input_file = "${call_cns}"
output_file = "${sample_id}.cnvkit.bed"

n_del = 0
n_dup = 0
n_neutral = 0
n_no_cn = 0

with open(input_file, 'r') as fin, open(output_file, 'w') as fout:
    # 寫 BED header（AnnotSV 支援 # 開頭的 header 行）
    fout.write("#chrom\\tstart\\tend\\tSVTYPE\\tSamples_ID\\n")

    for i, line in enumerate(fin):
        if i == 0:
            continue  # 跳過 header

        parts = line.strip().split('\\t')
        if len(parts) < 7:
            continue

        chrom = parts[0]
        start = parts[1]  # 已是 0-based（BED）
        end   = parts[2]
        cn_str = parts[6]  # cn 欄位

        # cn 欄位可能是空的（CNVkit 有時無法估計）
        if not cn_str or cn_str == '':
            n_no_cn += 1
            continue

        try:
            cn = int(float(cn_str))
        except ValueError:
            n_no_cn += 1
            continue

        if cn < 2:
            svtype = 'DEL'
            n_del += 1
        elif cn > 2:
            svtype = 'DUP'
            n_dup += 1
        else:
            # cn == 2：copy-neutral，不輸出
            n_neutral += 1
            continue

        fout.write(f"{chrom}\\t{start}\\t{end}\\t{svtype}\\t${sample_id}\\n")

print(f"[CNVKIT_TO_BED] 完成：DEL={n_del}, DUP={n_dup}, copy-neutral 跳過={n_neutral}, 無 CN={n_no_cn}", file=sys.stderr)
EOF

    echo "[CNVKIT_TO_BED] 輸出行數：\$(wc -l < ${sample_id}.cnvkit.bed)" >&2
    """
}

// ──────────────────────────────────────────────────────────────
// Process 2a：NCKUH WES CNV（gCNV VCF）→ AnnotSV TSV
// ──────────────────────────────────────────────────────────────

process ANNOTSV_CNV_NCKUH_WES {

    label 'process_medium'

    container "${params.sif_dir}/annotsv_3.5.10.sif"

    containerOptions "${params.apptainer_base_opts}"

    publishDir "${params.out_dir}/${sample_id}/06_cnv_sv", mode: 'copy'

    input:
    tuple val(sample_id), path(gcnv_vcf), path(gcnv_tbi)

    output:
    tuple val(sample_id),
          path("${sample_id}.cnv.annotated.tsv"),
          path("${sample_id}.cnv.unannotated.tsv"),
          emit: cnv_tsv_ch

    script:
    """
    echo "[ANNOTSV_CNV_NCKUH_WES] ${sample_id}：gCNV VCF → AnnotSV" >&2

    AnnotSV \\
        -SVinputFile ${gcnv_vcf} \\
        -outputDir . \\
        -outputFile ${sample_id}.cnv.annotated.tsv \\
        -genomeBuild GRCh38 \\
        -annotationsDir ${params.annotsv_annotations} \\
        -SVinputInfo 1

    # AnnotSV 有時輸出名稱會有細微差異，確保檔名正確
    if [ ! -f "${sample_id}.cnv.annotated.tsv" ]; then
        mv *.annotated.tsv ${sample_id}.cnv.annotated.tsv || true
    fi
    if [ ! -f "${sample_id}.cnv.unannotated.tsv" ]; then
        touch ${sample_id}.cnv.unannotated.tsv
    fi

    echo "[ANNOTSV_CNV_NCKUH_WES] 完成" >&2
    wc -l ${sample_id}.cnv.annotated.tsv >&2
    """
}

// ──────────────────────────────────────────────────────────────
// Process 2b：NCKUH WGS CNV（CNVkit BED）→ AnnotSV TSV
// ──────────────────────────────────────────────────────────────

process ANNOTSV_CNV_NCKUH_WGS {

    label 'process_medium'

    container "${params.sif_dir}/annotsv_3.5.10.sif"

    containerOptions "${params.apptainer_base_opts}"

    publishDir "${params.out_dir}/${sample_id}/06_cnv_sv", mode: 'copy'

    input:
    tuple val(sample_id), path(cnvkit_bed)

    output:
    tuple val(sample_id),
          path("${sample_id}.cnv.annotated.tsv"),
          path("${sample_id}.cnv.unannotated.tsv"),
          emit: cnv_tsv_ch

    script:
    """
    echo "[ANNOTSV_CNV_NCKUH_WGS] ${sample_id}：CNVkit BED → AnnotSV" >&2

    # BED 輸入需要指定 -svtBEDcol（第 4 欄是 SVTYPE）和 -samplesidBEDcol（第 5 欄是 sample ID）
    AnnotSV \\
        -SVinputFile ${cnvkit_bed} \\
        -svtBEDcol 4 \\
        -samplesidBEDcol 5 \\
        -outputDir . \\
        -outputFile ${sample_id}.cnv.annotated.tsv \\
        -genomeBuild GRCh38 \\
        -annotationsDir ${params.annotsv_annotations} \\
        -SVinputInfo 1

    if [ ! -f "${sample_id}.cnv.annotated.tsv" ]; then
        mv *.annotated.tsv ${sample_id}.cnv.annotated.tsv || true
    fi
    if [ ! -f "${sample_id}.cnv.unannotated.tsv" ]; then
        touch ${sample_id}.cnv.unannotated.tsv
    fi

    echo "[ANNOTSV_CNV_NCKUH_WGS] 完成" >&2
    wc -l ${sample_id}.cnv.annotated.tsv >&2
    """
}

// ──────────────────────────────────────────────────────────────
// Process 2c-pre：DRAGEN CNV 前處理（過濾 non-PASS 和 copy-neutral）
// ──────────────────────────────────────────────────────────────

process PREPARE_CNV_DRAGEN {

    label 'process_low'

    container "${params.sif_dir}/tertiary_python_1.0.0.sif"

    containerOptions "${params.apptainer_base_opts}"

    input:
    tuple val(sample_id), path(cnv_vcf)

    output:
    tuple val(sample_id),
          path("${sample_id}.cnv_filtered.vcf.gz"),
          emit: cnv_filtered_ch

    script:
    """
    echo "[PREPARE_CNV_DRAGEN] ${sample_id}：過濾 DRAGEN CNV VCF" >&2

    # 過濾邏輯：
    #   1. 只保留 PASS（移除 cnvLength 等 filter）
    #   2. 移除 copy-neutral segment（ALT=.，DRAGEN 會輸出 REF segment 作為背景）
    #      這些 REF segment 對 AnnotSV 沒有意義，且數量龐大會拖慢速度
    bcftools view -f PASS ${cnv_vcf} \
        | bcftools view -e 'ALT=="."' \
        -Oz -o ${sample_id}.cnv_filtered.vcf.gz

    tabix -p vcf ${sample_id}.cnv_filtered.vcf.gz

    echo "[PREPARE_CNV_DRAGEN] 過濾後 variant 數：" >&2
    bcftools stats ${sample_id}.cnv_filtered.vcf.gz | grep "^SN" >&2
    """
}

// ──────────────────────────────────────────────────────────────
// Process 3a-pre：DRAGEN SV 前處理（過濾 + INS symbolic allele 轉換）
// ──────────────────────────────────────────────────────────────

process PREPARE_SV_DRAGEN {

    label 'process_low'

    container "${params.sif_dir}/tertiary_python_1.0.0.sif"

    containerOptions "${params.apptainer_base_opts}"

    input:
    tuple val(sample_id), path(sv_vcf)

    output:
    tuple val(sample_id),
          path("${sample_id}.sv_filtered.vcf.gz"),
          emit: sv_filtered_ch

    script:
    """
    echo "[PREPARE_SV_DRAGEN] ${sample_id}：過濾 DRAGEN SV VCF" >&2

    # Step 1：過濾 PASS（移除 MaxMQ0Frac 等 filter）
    bcftools view -f PASS ${sv_vcf} -Oz -o sv_pass.vcf.gz
    tabix -p vcf sv_pass.vcf.gz

    # Step 2：把 INS 真實序列 ALT 換成 <INS> symbolic allele
    # 用獨立 script 避免 Nextflow heredoc 的跳脫字元問題
    python3 ${params.scripts_dir}/prepare_sv_dragen.py \
        --input  sv_pass.vcf.gz \
        --output ${sample_id}.sv_filtered.vcf

    bgzip ${sample_id}.sv_filtered.vcf
    tabix -p vcf ${sample_id}.sv_filtered.vcf.gz
    rm -f sv_pass.vcf.gz sv_pass.vcf.gz.tbi

    echo "[PREPARE_SV_DRAGEN] 完成" >&2
    bcftools stats ${sample_id}.sv_filtered.vcf.gz | grep "^SN" >&2
    """
}

// ──────────────────────────────────────────────────────────────
// Process 2c：DRAGEN CNV → AnnotSV TSV
// ──────────────────────────────────────────────────────────────

process ANNOTSV_CNV_DRAGEN {

    label 'process_medium'

    container "${params.sif_dir}/annotsv_3.5.10.sif"

    containerOptions "${params.apptainer_base_opts}"

    publishDir "${params.out_dir}/${sample_id}/06_cnv_sv", mode: 'copy'

    input:
    // 接收 PREPARE_CNV_DRAGEN 輸出的 filtered VCF
    tuple val(sample_id), path(cnv_vcf)

    output:
    tuple val(sample_id),
          path("${sample_id}.cnv.annotated.tsv"),
          path("${sample_id}.cnv.unannotated.tsv"),
          emit: cnv_tsv_ch

    script:
    """
    echo "[ANNOTSV_CNV_DRAGEN] ${sample_id}：DRAGEN CNV VCF（filtered）→ AnnotSV" >&2

    AnnotSV \\
        -SVinputFile ${cnv_vcf} \\
        -outputDir . \\
        -outputFile ${sample_id}.cnv.annotated.tsv \\
        -genomeBuild GRCh38 \\
        -annotationsDir ${params.annotsv_annotations} \\
        -SVinputInfo 1

    if [ ! -f "${sample_id}.cnv.annotated.tsv" ]; then
        mv *.annotated.tsv ${sample_id}.cnv.annotated.tsv || true
    fi
    if [ ! -f "${sample_id}.cnv.unannotated.tsv" ]; then
        touch ${sample_id}.cnv.unannotated.tsv
    fi

    echo "[ANNOTSV_CNV_DRAGEN] 完成" >&2
    wc -l ${sample_id}.cnv.annotated.tsv >&2
    """
}

// ──────────────────────────────────────────────────────────────
// Process 3a：NCKUH SV（Delly VCF）→ AnnotSV TSV（WES + WGS 共用）
// ──────────────────────────────────────────────────────────────

process ANNOTSV_SV_NCKUH {

    label 'process_medium'

    container "${params.sif_dir}/annotsv_3.5.10.sif"

    containerOptions "${params.apptainer_base_opts}"

    publishDir "${params.out_dir}/${sample_id}/06_cnv_sv", mode: 'copy'

    input:
    tuple val(sample_id), path(delly_vcf), path(delly_tbi)

    output:
    tuple val(sample_id),
          path("${sample_id}.sv.annotated.tsv"),
          path("${sample_id}.sv.unannotated.tsv"),
          emit: sv_tsv_ch

    script:
    """
    echo "[ANNOTSV_SV_NCKUH] ${sample_id}：Delly VCF → AnnotSV" >&2

    AnnotSV \\
        -SVinputFile ${delly_vcf} \\
        -outputDir . \\
        -outputFile ${sample_id}.sv.annotated.tsv \\
        -genomeBuild GRCh38 \\
        -annotationsDir ${params.annotsv_annotations} \\
        -SVinputInfo 1

    if [ ! -f "${sample_id}.sv.annotated.tsv" ]; then
        mv *.annotated.tsv ${sample_id}.sv.annotated.tsv || true
    fi
    if [ ! -f "${sample_id}.sv.unannotated.tsv" ]; then
        touch ${sample_id}.sv.unannotated.tsv
    fi

    echo "[ANNOTSV_SV_NCKUH] 完成" >&2
    wc -l ${sample_id}.sv.annotated.tsv >&2
    """
}

// ──────────────────────────────────────────────────────────────
// Process 3b：DRAGEN SV → AnnotSV TSV
// ──────────────────────────────────────────────────────────────

process ANNOTSV_SV_DRAGEN {

    label 'process_medium'

    container "${params.sif_dir}/annotsv_3.5.10.sif"

    containerOptions "${params.apptainer_base_opts}"

    publishDir "${params.out_dir}/${sample_id}/06_cnv_sv", mode: 'copy'

    input:
    // 接收 PREPARE_SV_DRAGEN 輸出的 filtered + INS converted VCF
    tuple val(sample_id), path(sv_vcf)

    output:
    tuple val(sample_id),
          path("${sample_id}.sv.annotated.tsv"),
          path("${sample_id}.sv.unannotated.tsv"),
          emit: sv_tsv_ch

    script:
    """
    echo "[ANNOTSV_SV_DRAGEN] ${sample_id}：DRAGEN SV VCF（filtered）→ AnnotSV" >&2

    AnnotSV \\
        -SVinputFile ${sv_vcf} \\
        -outputDir . \\
        -outputFile ${sample_id}.sv.annotated.tsv \\
        -genomeBuild GRCh38 \\
        -annotationsDir ${params.annotsv_annotations} \\
        -SVinputInfo 1

    if [ ! -f "${sample_id}.sv.annotated.tsv" ]; then
        mv *.annotated.tsv ${sample_id}.sv.annotated.tsv || true
    fi
    if [ ! -f "${sample_id}.sv.unannotated.tsv" ]; then
        touch ${sample_id}.sv.unannotated.tsv
    fi

    echo "[ANNOTSV_SV_DRAGEN] 完成" >&2
    wc -l ${sample_id}.sv.annotated.tsv >&2
    """
}
