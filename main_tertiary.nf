#!/usr/bin/env nextflow
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
 *
 * DISCLAIMER: This pipeline is provided "as is" without
 * warranty of any kind. The authors and their institution
 * make no representations or warranties regarding the
 * accuracy, completeness, or suitability of the analysis
 * results for any clinical or research purpose. Users are
 * solely responsible for validating and interpreting all
 * results.
 * =========================================================
 * main_tertiary.nf
 * ================
 * 臨床 WGS/WES 三級分析主 workflow（支援 sample sheet 批次輸入）
 *
 * 支援兩種輸入來源（sample sheet 的 pipeline_type 欄位切換）：
 *
 *   nckuh：NCKUH 二級分析 ensemble VCF
 *     路徑規則：{input_dir}/04_snv_indel/{sample_id}.ensemble.fixed.vcf.gz
 *
 *   dragen：Illumina DRAGEN hard-filtered VCF
 *     路徑規則：{input_dir}/vcf.gz/{sample_id}.hard-filtered.vcf.gz
 *
 * --pipeline_type（選填）：過濾 sample sheet 中 pipeline_type 符合的 row。
 *   不傳入 → 要求 sample sheet 內所有 row 的 pipeline_type 一致（否則 error）
 *   傳入   → 只跑符合的 row，不符合的 warn 後跳過
 *
 * Sample sheet 格式（CSV）：
 *   sample_id,pipeline_type,input_dir,seq_type,hpo
 *   NA12878_WES,nckuh,/path/to/nckuh_output,WES,
 *   VAL-10,dragen,/path/to/dragen_output,WGS,HP:0001250|HP:0001263
 *   VAL-11,dragen,/path/to/dragen_output,WGS,
 *
 * 使用方式：
 *   # 不傳 --pipeline_type（sample sheet 裡必須全部同一種）
 *   nextflow -c nextflow_tertiary.config run main_tertiary.nf \
 *       -profile local \
 *       --samplesheet samplesheet_tertiary.csv \
 *       --out_dir /scratch/pylin1991/tertiary_test \
 *       -resume
 *
 *   # 傳 --pipeline_type 過濾（sample sheet 可以混放，只跑指定的那種）
 *   nextflow -c nextflow_tertiary.config run main_tertiary.nf \
 *       -profile dgm \
 *       --pipeline_type dragen \
 *       --samplesheet all_samples.csv \
 *       --out_dir /home/pipeline/tertiary_output \
 *       -resume
 *
 * 開發進度：
 *   ✅ PREPARE_VCF        - NCKUH ensemble VCF 前處理
 *   ✅ PREPARE_VCF_DRAGEN - DRAGEN VCF 前處理（chrM 分流）
 *   ✅ SNV_ANNOTATE       - VEP 115 annotation
 *   ✅ PARSE_VEP_CSQ      - transcript 選取 + TSV 產生
 *   ✅ ACMG_CLASSIFY      - ACMG evidence 分類
 *   ✅ MITO_ANNOTATE      - mtDNA annotation（VEP + MITOMAP）
 *   🔲 CNV_SV            - AnnotSV
 *   🔲 STR               - STRchive
 *   🔲 ROH               - consanguinity + UPD
 *   🔲 PHENOTYPE         - Exomiser + LIRICAL
 *   🔲 PGX               - Aldy
 *   🔲 SECONDARY         - ACMG SF v3.2
 */

nextflow.enable.dsl = 2

// ──────────────────────────────────────────────────────────────
// 匯入 modules
// ──────────────────────────────────────────────────────────────

include { PREPARE_VCF        } from './modules/prepare_vcf.nf'
include { PREPARE_VCF_DRAGEN } from './modules/prepare_vcf_dragen.nf'
include { SNV_ANNOTATE       } from './modules/snv_annotation.nf'
include { PARSE_VEP_CSQ      } from './modules/parse_csq.nf'
include { ACMG_CLASSIFY      } from './modules/acmg_classifier.nf'
include { MITO_ANNOTATE      } from './modules/mito_annotation.nf'

// ──────────────────────────────────────────────────────────────
// Sample sheet 解析
// ──────────────────────────────────────────────────────────────

def parse_samplesheet(csv_path, pipeline_type_filter = null) {
    /*
     * 解析三級分析 sample sheet（CSV）。
     * 回傳 list of maps，每個 map 對應一個樣本。
     *
     * 必要欄位：sample_id, pipeline_type, input_dir, seq_type
     * 選填欄位：hpo（空白表示無）
     *
     * pipeline_type_filter（來自 --pipeline_type 參數）：
     *   null   → 不過濾，但要求 sample sheet 內所有 row 的 pipeline_type 一致
     *   有值   → 只保留 pipeline_type 符合的 row，不符合的 warn 後跳過
     *
     * 驗證規則：
     *   - 必要欄位不得為空
     *   - pipeline_type 只接受 nckuh 或 dragen
     *   - seq_type 只接受 WES 或 WGS
     *   - sample_id 不得重複
     *   - 過濾後至少要有一個樣本
     */

    def csv_file = file(csv_path)
    if (!csv_file.exists()) {
        error "[ERROR] Sample sheet 不存在：${csv_path}"
    }

    def lines    = csv_file.readLines()
    def header   = lines[0].split(',').collect { it.trim() }
    def samples  = []
    def seen_ids = [] as Set

    // 確認必要欄位存在
    def required_cols = ['sample_id', 'pipeline_type', 'input_dir', 'seq_type']
    for (col in required_cols) {
        if (!header.contains(col)) {
            error "[ERROR] Sample sheet 缺少必要欄位：${col}\n       現有欄位：${header}"
        }
    }

    // 解析每一行
    lines[1..-1].eachWithIndex { line, idx ->
        line = line.trim()
        if (!line || line.startsWith('#')) return  // 跳過空行和 comment

        def row_num = idx + 2  // 給使用者看的行號（從 2 開始，1 是 header）
        def values  = line.split(',', -1).collect { it.trim() }

        // 欄位數量檢查
        if (values.size() < header.size()) {
            // 補齊空欄位
            while (values.size() < header.size()) values << ''
        }

        def row = [header, values].transpose().collectEntries { k, v -> [(k): v] }

        // 必要欄位不得為空
        for (col in required_cols) {
            if (!row[col]) {
                error "[ERROR] Sample sheet 第 ${row_num} 行，欄位 '${col}' 不得為空"
            }
        }

        // pipeline_type 驗證
        if (!['nckuh', 'dragen'].contains(row.pipeline_type)) {
            error "[ERROR] Sample sheet 第 ${row_num} 行，pipeline_type '${row.pipeline_type}' 無效\n       只接受：nckuh 或 dragen"
        }

        // seq_type 驗證
        if (!['WES', 'WGS'].contains(row.seq_type)) {
            error "[ERROR] Sample sheet 第 ${row_num} 行，seq_type '${row.seq_type}' 無效\n       只接受：WES 或 WGS"
        }

        // sample_id 重複檢查
        if (seen_ids.contains(row.sample_id)) {
            error "[ERROR] Sample sheet 第 ${row_num} 行，sample_id '${row.sample_id}' 重複"
        }
        seen_ids << row.sample_id

        // input_dir 存在檢查
        def input_dir = file(row.input_dir)
        if (!input_dir.exists()) {
            error "[ERROR] Sample sheet 第 ${row_num} 行，input_dir 不存在：${row.input_dir}"
        }

        // hpo 預設空字串
        row.hpo = row.hpo ?: ''

        // pipeline_type_filter 過濾（來自 --pipeline_type 參數）
        if (pipeline_type_filter && row.pipeline_type != pipeline_type_filter) {
            log.warn "[WARN] 樣本 ${row.sample_id} 的 pipeline_type=${row.pipeline_type}" +
                     " 不符合 --pipeline_type=${pipeline_type_filter}，跳過"
            return
        }

        samples << row
    }

    if (samples.isEmpty()) {
        error "[ERROR] Sample sheet 沒有任何有效樣本：${csv_path}"
    }

    return samples
}

// ──────────────────────────────────────────────────────────────
// 共用參數驗證（資料庫路徑）
// ──────────────────────────────────────────────────────────────

def validate_databases() {
    def checks = [
        ['VEP cache',  params.vep_cache],
        ['dbNSFP',     params.dbnsfp],
        ['LOFTEE dir', params.loftee_dir],
        ['ClinVar',    params.clinvar],
    ]
    for (chk in checks) {
        def f = file(chk[1])
        if (!f.exists()) {
            error "[ERROR] ${chk[0]} 不存在：${chk[1]}"
        }
    }

    // Optional 資料庫（只 warn）
    if (params.clingen_hi_tsv && !file(params.clingen_hi_tsv).exists()) {
        log.warn "[WARN] ClinGen HI TSV 不存在：${params.clingen_hi_tsv}\n       PVS1 將使用簡化版"
    }
    if (params.gene_moi_tsv && !file(params.gene_moi_tsv).exists()) {
        log.warn "[WARN] gene_moi.tsv.gz 不存在：${params.gene_moi_tsv}\n       PM2 將使用純 AF 閾值判斷"
    }
}

// ──────────────────────────────────────────────────────────────
// 主 workflow
// ──────────────────────────────────────────────────────────────

workflow {

    // ── 基本參數檢查 ──────────────────────────────────────────
    if (!params.samplesheet) {
        error "[ERROR] --samplesheet 為必填參數\n       例如：--samplesheet samplesheet_tertiary.csv"
    }
    if (!params.out_dir) {
        error "[ERROR] --out_dir 為必填參數"
    }

    // ── 解析 sample sheet ─────────────────────────────────────
    def samples = parse_samplesheet(params.samplesheet, params.pipeline_type ?: null)

    // ── pipeline_type 決定與一致性檢查 ──────────────────────────
    def pipeline_types = samples.collect { it.pipeline_type }.unique()

    // --pipeline_type 傳入時：已在 parse_samplesheet 過濾，pipeline_types 必定只有一種
    // --pipeline_type 未傳入時：要求 sample sheet 內所有 row 一致
    if (!params.pipeline_type && pipeline_types.size() > 1) {
        error "[ERROR] Sample sheet 含有多種 pipeline_type：${pipeline_types}\n" +
              "       請拆成兩個 sample sheet 分別執行，\n" +
              "       或用 --pipeline_type 指定要跑的那種（其他 row 會被跳過）"
    }
    def pipeline_type = pipeline_types[0]

    // ── 資料庫路徑驗證 ────────────────────────────────────────
    validate_databases()

    // ── 印出執行資訊 ──────────────────────────────────────────
    log.info """
    ╔══════════════════════════════════════════════════════╗
    ║         臨床三級分析 Pipeline  v1.0.0                ║
    ╚══════════════════════════════════════════════════════╝
    Pipeline 類型 : ${pipeline_type.toUpperCase()}
    樣本數        : ${samples.size()}
    Sample sheet  : ${params.samplesheet}
    輸出目錄      : ${params.out_dir}
    容器目錄      : ${params.sif_dir}
    VEP cache     : ${params.vep_cache}
    dbNSFP        : ${params.dbnsfp}
    ClinVar       : ${params.clinvar}
    ClinGen HI    : ${params.clingen_hi_tsv ?: '（未提供，PVS1 簡化版）'}
    Gene MOI      : ${params.gene_moi_tsv   ?: '（未提供，PM2 純 AF 模式）'}
    """.stripIndent()

    // 印出樣本清單
    samples.each { s ->
        log.info "  樣本：${s.sample_id}  seq_type=${s.seq_type}  input=${s.input_dir}"
    }

    // ── 建立 channel（依 pipeline_type 推導 VCF 路徑）────────
    if (pipeline_type == 'nckuh') {

        // NCKUH：{input_dir}/04_snv_indel/{sample_id}.ensemble.fixed.vcf.gz
        input_ch = Channel.fromList(samples).map { s ->
            def vcf = file("${s.input_dir}/04_snv_indel/${s.sample_id}.ensemble.fixed.vcf.gz")
            def tbi = file("${s.input_dir}/04_snv_indel/${s.sample_id}.ensemble.fixed.vcf.gz.tbi")
            if (!vcf.exists()) {
                error "[ERROR] 找不到 ensemble VCF：${vcf}"
            }
            if (!tbi.exists()) {
                error "[ERROR] 找不到 ensemble VCF index：${tbi}"
            }
            tuple(s.sample_id, vcf, tbi)
        }

        PREPARE_VCF(input_ch)
        snv_ch = PREPARE_VCF.out.snv_ch

        // NCKUH mito：{input_dir}/07_mitochondria/{sample_id}.mito.vcf.gz
        // 格式統一為 tuple(sample_id, pipeline_type, mito_vcf, mito_tbi)
        // 與 DRAGEN 的 mito_ch 格式相同，供 MITO_ANNOTATE 共用
        nckuh_mito_ch = Channel.fromList(samples).map { s ->
            def mito_vcf = file("${s.input_dir}/07_mitochondria/${s.sample_id}.mito.vcf.gz")
            def mito_tbi = file("${s.input_dir}/07_mitochondria/${s.sample_id}.mito.vcf.gz.tbi")
            if (!mito_vcf.exists()) {
                // mito 不是所有樣本都有，找不到只 warn 不中斷
                log.warn "[WARN] 找不到 mito VCF，跳過：${mito_vcf}"
                return null
            }
            tuple(s.sample_id, "nckuh", mito_vcf, mito_tbi)
        }
        // 過濾掉 null（找不到 mito VCF 的樣本）
        .filter { it != null }

        MITO_ANNOTATE(nckuh_mito_ch)

    } else {

        // DRAGEN：{input_dir}/vcf.gz/{sample_id}.hard-filtered.vcf.gz
        input_ch = Channel.fromList(samples).map { s ->
            def vcf = file("${s.input_dir}/vcf.gz/${s.sample_id}.hard-filtered.vcf.gz")
            if (!vcf.exists()) {
                error "[ERROR] 找不到 DRAGEN VCF：${vcf}"
            }
            // tbi 不在這裡檢查：由 ENSURE_DRAGEN_TBI process 自動建立
            tuple(s.sample_id, vcf)
        }

        PREPARE_VCF_DRAGEN(input_ch)
        snv_ch = PREPARE_VCF_DRAGEN.out.snv_ch

        // DRAGEN mito_ch 已由 PREPARE_VCF_DRAGEN 備妥
        // 格式：tuple(sample_id, mito_vcf.gz, mito_tbi)
        // 加上 pipeline_type，統一成 tuple(sample_id, pipeline_type, mito_vcf, mito_tbi)
        dragen_mito_ch = PREPARE_VCF_DRAGEN.out.mito_ch.map { sample_id, mito_vcf, mito_tbi ->
            tuple(sample_id, "dragen", mito_vcf, mito_tbi)
        }

        MITO_ANNOTATE(dragen_mito_ch)
    }

    // ── 以下完全共用（兩種 pipeline 相同）────────────────────

    SNV_ANNOTATE(snv_ch)

    PARSE_VEP_CSQ(
        SNV_ANNOTATE.out.vep_ch,
        SNV_ANNOTATE.out.pangolin_ch
    )

    def clingen_hi_file = (params.clingen_hi_tsv && file(params.clingen_hi_tsv).exists())
        ? file(params.clingen_hi_tsv)
        : file("NO_FILE")

    def gene_moi_file = (params.gene_moi_tsv && file(params.gene_moi_tsv).exists())
        ? file(params.gene_moi_tsv)
        : file("NO_FILE")

    ACMG_CLASSIFY(
        PARSE_VEP_CSQ.out.full_tsv_ch,
        clingen_hi_file,
        gene_moi_file
    )
}
