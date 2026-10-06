/*
 * =========================================================
 * WGS/WES Germline Analysis Pipeline - CNV/SV Module
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
 * modules/prepare_vcf.nf
 * ======================
 * 目的：
 *   將二級分析產生的 ensemble.fixed.vcf.gz（雙 sample column：_DV + _HC）
 *   處理為三級分析 VEP annotation 的輸入，共兩個步驟：
 *
 *   Step 1 - ADD_CALLERS_TAG：
 *     執行 add_callers_tag.py，在 INFO 欄位新增 CALLERS tag（DV+HC/DV/HC/NONE；
 *     NONE = 兩個 caller 都沒有 ALT genotype，Step 2 會擋掉）
 *
 *   Step 2 - FILTER_FOR_ANNOTATION：
 *     用 bcftools 依 CALLERS 過濾（不看 FILTER 欄）：只收 DV+HC / DV / HC，
 *     擋掉 CALLERS=NONE（兩個 sample column 都沒有 ALT genotype，./. 或 0/0）。
 *     ⚠️ 這個「擋掉兩邊都沒 call」的設計原本就寫在這裡，但舊版 add_callers_tag.py
 *        把這種紀錄標成 HC，所以實際上從未生效（見下方 FILTER_FOR_ANNOTATION 註解）。
 *     同時 bgzip 壓縮 + tabix index，產生標準的 .vcf.gz + .tbi
 *
 * 輸入（來自 main_tertiary.nf）：
 *   tuple val(sample_id), path(ensemble_vcf), path(ensemble_tbi)
 *
 * 輸出：
 *   tuple val(sample_id), path("*.snv_for_annotation.vcf.gz"), path("*.snv_for_annotation.vcf.gz.tbi")
 *
 * 使用的容器：
 *   Step 1：tertiary_python_1.0.0.sif（含 cyvcf2）
 *   Step 2：bcftools_1.23.1.sif
 */

// ──────────────────────────────────────────────────────────────
// Process 1：新增 CALLERS tag
// ──────────────────────────────────────────────────────────────

process ADD_CALLERS_TAG {
    // 標籤：對應 nextflow_tertiary.config 的資源設定
    label 'process_low'

    // 使用 tertiary_python sif（含 cyvcf2）
    container "${params.sif_dir}/tertiary_python_1.0.0.sif"

    // 輸入：ensemble VCF + 其 index
    input:
    tuple val(sample_id), path(ensemble_vcf), path(ensemble_tbi)

    // 輸出：加上 CALLERS tag 的未壓縮 VCF（暫時檔，交給下一個 process）
    output:
    tuple val(sample_id), path("${sample_id}.callers_tagged.vcf")

    script:
    """
    # 先正規化 + 拆分多等位基因（left-align + split multiallelics）。
    # 讓每個 ALT 各自一列，使 CALLERS/AD/VAF 以及後續 VEP / gnomAD / ClinVar 註解
    # 都是 per-allele 正確（避免多等位基因位點把不同 allele 的 AF / 致病性混在一起）。
    # -c w：REF 與參考不符時只警告不中斷（同一 hg38 參考下通常不會發生）。
    bcftools norm -m -any -f ${params.ref_fasta} -c w \\
        ${ensemble_vcf} \\
        -Oz -o ${sample_id}.norm.vcf.gz

    # 執行 add_callers_tag.py（吃正規化後的 biallelic VCF）
    # --sample 傳入 sample_id，腳本會自動尋找 {sample_id}_DV 和 {sample_id}_HC column
    python3 ${params.scripts_dir}/add_callers_tag.py \\
        --input  ${sample_id}.norm.vcf.gz \\
        --sample ${sample_id} \\
        --output ${sample_id}.callers_tagged.vcf
    """
}

// ──────────────────────────────────────────────────────────────
// Process 2：過濾 + bgzip + tabix
// ──────────────────────────────────────────────────────────────

process FILTER_FOR_ANNOTATION {
    label 'process_low'

    // 使用既有的 bcftools sif
    container "${params.sif_dir}/bcftools_1.23.1.sif"

    // publishDir：將最終輸出複製到三級分析輸出目錄
    // mode: 'copy' 確保輸出目錄有獨立的檔案（不是 symlink）
    publishDir "${params.out_dir}/${sample_id}/00_prepare", mode: 'copy'

    input:
    tuple val(sample_id), path(callers_tagged_vcf)

    // 輸出：bgzip 壓縮的 VCF + tabix index
    output:
    tuple val(sample_id),
          path("${sample_id}.snv_for_annotation.vcf.gz"),
          path("${sample_id}.snv_for_annotation.vcf.gz.tbi")

    script:
    """
    # 過濾策略：依據 CALLERS tag 過濾，不使用 FILTER 欄位。
    #
    # 背景：
    #   ensemble 的 FILTER 是 bcftools merge 從各 caller 的紀錄合併而來，不代表
    #   「哪個 caller 有 call」—— 用 FILTER="PASS" 過濾會把 HC-only 的 call 丟掉，
    #   違反 ensemble「只要有一個 caller call 到就保留」的設計。
    #   （二級 BCFTOOLS_ENSEMBLE 已在 merge 前丟掉 DV 沒有 ALT 的紀錄，所以
    #    ensemble 裡不會再出現 DV 的 FILTER=RefCall；見二級 postprocessing.nf。）
    #
    # CALLERS 由 add_callers_tag.py 依 GT 判斷，四種值：
    #   CALLERS=DV+HC → 兩個都有 ALT call
    #   CALLERS=DV    → 只有 DV 有 ALT call
    #   CALLERS=HC    → 只有 HC 有 ALT call
    #   CALLERS=NONE  → 兩個都沒有 ALT call（./. 或 0/0）→ 下面的 -i 條件會擋掉
    #
    # ⚠️ 舊版 add_callers_tag.py 的 determine_callers() 沒有 NONE：兩邊都沒 call
    #   會掉進 else 被標成 "HC"，於是這些非變異紀錄通過本步驟、以
    #   ZYGOSITY=ref/unknown 出現在 ACMG 表（實例：SUZ12 多出錯誤的 c.2170del）。
    #   這段註解舊版還寫著「這種 case 不會有 CALLERS tag」—— 那是錯的，tag 一律會寫。
    #   同一個 else 也讓 stderr 的「HC only」統計被灌水：舊註解引用的
    #   「NA12878_WES HC-only 23.9%（8,897/37,198）」包含了這些 no-call 紀錄。
    #   修正後實測（VAL55 WGS，2026-09）：DV+HC 89.9%、DV only 2.6%、HC only 7.4%、NONE 0.2%。

    bcftools view \\
        -i 'INFO/CALLERS="DV+HC" || INFO/CALLERS="DV" || INFO/CALLERS="HC"' \\
        ${callers_tagged_vcf} \\
        -Oz -o ${sample_id}.snv_for_annotation.vcf.gz

    # 建立 tabix index（VEP 和後續工具都需要）
    tabix -p vcf ${sample_id}.snv_for_annotation.vcf.gz

    # 輸出統計（寫進 log，方便 debug）
    echo "[FILTER_FOR_ANNOTATION] ${sample_id}" >&2
    bcftools stats ${sample_id}.snv_for_annotation.vcf.gz | \\
        grep "^SN" >&2
    """
}

// ──────────────────────────────────────────────────────────────
// 組合 workflow（供 main_tertiary.nf 呼叫）
// ──────────────────────────────────────────────────────────────

workflow PREPARE_VCF {
    // 輸入 channel：tuple(sample_id, ensemble_vcf, ensemble_tbi)
    take:
    ensemble_ch

    // 執行兩個 process，串接輸出
    main:
    ADD_CALLERS_TAG(ensemble_ch)
    FILTER_FOR_ANNOTATION(ADD_CALLERS_TAG.out)

    // 輸出 channel：tuple(sample_id, snv_vcf, snv_tbi)
    emit:
    snv_ch = FILTER_FOR_ANNOTATION.out
}
