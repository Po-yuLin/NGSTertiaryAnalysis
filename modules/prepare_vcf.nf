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
 *     執行 add_callers_tag.py，在 INFO 欄位新增 CALLERS tag（DV+HC/DV/HC）
 *
 *   Step 2 - FILTER_FOR_ANNOTATION：
 *     用 bcftools 過濾掉不適合送進 VEP 的 variant：
 *       - RefCall（FILTER=RefCall，DV 叫 0/0 但 HC 叫 ./. 的情況）
 *       - 兩個 sample column 都是 0/0 或 ./. 的 variant
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
    # 執行 add_callers_tag.py
    # --sample 傳入 sample_id，腳本會自動尋找 {sample_id}_DV 和 {sample_id}_HC column
    python3 ${params.scripts_dir}/add_callers_tag.py \\
        --input  ${ensemble_vcf} \\
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
    #   ensemble VCF 的 FILTER 欄位由二級分析設定，語義如下：
    #     FILTER=PASS    → DV 有 call（DV+HC 或 DV-only）
    #     FILTER=RefCall → DV 叫 0/0（HC 可能有 call 也可能是 ./.）
    #     FILTER=.       → 兩個 caller 都沒有 call
    #
    # 問題：
    #   用 FILTER="PASS" 過濾會把所有 HC-only variant（FILTER=RefCall）
    #   全部丟掉，違反 joint calling「只要有一個 caller call 到就保留」的規則。
    #   NA12878_WES 測試確認：HC-only 佔 23.9%（8,897/37,198 個），不應被丟棄。
    #
    # 修正：
    #   add_callers_tag.py 已根據 GT 欄位正確判斷每個 variant 的 call 狀態：
    #     CALLERS=DV+HC → 兩個都有 ALT call，來自 FILTER=PASS
    #     CALLERS=DV    → 只有 DV 有 ALT call，來自 FILTER=PASS
    #     CALLERS=HC    → 只有 HC 有 ALT call，來自 FILTER=RefCall（HC=0/1）
    #   RefCall 中 DV=0/0, HC=./. 的 variant，CALLERS 被判為 HC（但 is_called 回傳 False）
    #   → 實際上這種 case 兩個 caller 都沒有 call，不會有 CALLERS tag。
    #   所以只要 CALLERS 有值，就代表至少一個 caller 有有效的 ALT call。

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
