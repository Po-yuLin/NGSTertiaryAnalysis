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
 *   ✅ STR_ANNOTATE       - STRchive threshold 分類
 *   ✅ CNV_SV_ANNOTATE   - AnnotSV（NCKUH WES/WGS + DRAGEN）
 *   🔲 ROH               - consanguinity + UPD
 *   🔲 PHENOTYPE         - Exomiser + LIRICAL
 *   ✅ PGX               - PharmCAT + StellarPGx（CYP2D6）+ OptiType（HLA）
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
include { STR_PREPARE_NCKUH  } from './modules/str_annotation.nf'
include { STR_PARSE_NCKUH    } from './modules/str_annotation.nf'
include { STR_PARSE_DRAGEN   } from './modules/str_annotation.nf'
include { CNVKIT_TO_BED          } from './modules/cnv_sv_annotation.nf'
include { ANNOTSV_CNV_NCKUH_WES  } from './modules/cnv_sv_annotation.nf'
include { ANNOTSV_CNV_NCKUH_WGS  } from './modules/cnv_sv_annotation.nf'
include { ANNOTSV_SV_NCKUH       } from './modules/cnv_sv_annotation.nf'
include { PREPARE_CNV_DRAGEN     } from './modules/cnv_sv_annotation.nf'
include { PREPARE_SV_DRAGEN      } from './modules/cnv_sv_annotation.nf'
include { ANNOTSV_CNV_DRAGEN     } from './modules/cnv_sv_annotation.nf'
include { ANNOTSV_SV_DRAGEN      } from './modules/cnv_sv_annotation.nf'
include { PGX_ANNOTATE           } from './modules/pgx_annotation.nf'

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
    PGx           : ${params.run_pgx ? '啟用' : '停用'}（--run_pgx）
    PGx CYP2D6    : ${params.run_pgx_cyp2d6 ? 'StellarPGx BAM-based' : 'VCF only（準確度較低）'}
    PGx HLA       : ${params.run_pgx_hla ? 'OptiType BAM-based' : '停用'}
    """.stripIndent()

    // 印出樣本清單
    samples.each { s ->
        log.info "  樣本：${s.sample_id}  seq_type=${s.seq_type}  input=${s.input_dir}"
    }

    // ── BAM channel 預先定義（scope 跨 if/else）────────────────
    def bam_ch_pgx = Channel.empty()

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
        nckuh_mito_ch = Channel.fromList(samples).map { s ->
            def mito_vcf = file("${s.input_dir}/07_mitochondria/${s.sample_id}.mito.vcf.gz")
            def mito_tbi = file("${s.input_dir}/07_mitochondria/${s.sample_id}.mito.vcf.gz.tbi")
            if (!mito_vcf.exists()) {
                log.warn "[WARN] 找不到 mito VCF，跳過：${mito_vcf}"
                return null
            }
            tuple(s.sample_id, "nckuh", mito_vcf, mito_tbi)
        }
        .filter { it != null }

        MITO_ANNOTATE(nckuh_mito_ch)

        // NCKUH STR：{input_dir}/06_repeat/{sample_id}.str.vcf（未壓縮）
        nckuh_str_ch = Channel.fromList(samples).map { s ->
            def str_vcf = file("${s.input_dir}/06_repeat/${s.sample_id}.str.vcf")
            if (!str_vcf.exists()) {
                log.warn "[WARN] 找不到 STR VCF，跳過：${str_vcf}"
                return null
            }
            tuple(s.sample_id, str_vcf)
        }
        .filter { it != null }

        STR_PREPARE_NCKUH(nckuh_str_ch)
        STR_PARSE_NCKUH(STR_PREPARE_NCKUH.out.str_prepared_ch)

        // ── NCKUH CNV：依 seq_type 分流 WES / WGS ────────────────
        def seq_types = samples.collect { it.seq_type }.unique()

        if (seq_types.contains("WES")) {
            nckuh_cnv_wes_ch = Channel.fromList(samples)
                .filter { s -> s.seq_type == "WES" }
                .map { s ->
                    def vcf = file("${s.input_dir}/05_cnv_sv/${s.sample_id}.gcnv.vcf.gz")
                    def tbi = file("${s.input_dir}/05_cnv_sv/${s.sample_id}.gcnv.vcf.gz.tbi")
                    if (!vcf.exists()) {
                        log.warn "[WARN] 找不到 gCNV VCF，跳過：${vcf}"
                        return null
                    }
                    tuple(s.sample_id, vcf, tbi)
                }
                .filter { it != null }

            ANNOTSV_CNV_NCKUH_WES(nckuh_cnv_wes_ch)
        }

        if (seq_types.contains("WGS")) {
            nckuh_cnvkit_ch = Channel.fromList(samples)
                .filter { s -> s.seq_type == "WGS" }
                .map { s ->
                    def cns = file("${s.input_dir}/05_cnv_sv/${s.sample_id}.call.cns")
                    if (!cns.exists()) {
                        log.warn "[WARN] 找不到 CNVkit .call.cns，跳過：${cns}"
                        return null
                    }
                    tuple(s.sample_id, cns)
                }
                .filter { it != null }

            CNVKIT_TO_BED(nckuh_cnvkit_ch)
            ANNOTSV_CNV_NCKUH_WGS(CNVKIT_TO_BED.out.cnvkit_bed_ch)
        }

        // ── NCKUH SV（Delly，WES + WGS 共用）────────────────────
        nckuh_sv_ch = Channel.fromList(samples).map { s ->
            def vcf = file("${s.input_dir}/05_cnv_sv/${s.sample_id}.delly.vcf.gz")
            def tbi = file("${s.input_dir}/05_cnv_sv/${s.sample_id}.delly.vcf.gz.tbi")
            if (!vcf.exists()) {
                log.warn "[WARN] 找不到 Delly VCF，跳過：${vcf}"
                return null
            }
            tuple(s.sample_id, vcf, tbi)
        }
        .filter { it != null }

        ANNOTSV_SV_NCKUH(nckuh_sv_ch)

        // ── NCKUH BAM channel（PGx CYP2D6/HLA 用）────────────────
        // StellarPGx 只支援 WGS，WES 樣本直接跳過
        // BAM 路徑：{input_dir}/02_alignment/{sample_id}.aligned.sorted.bam
        bam_ch_pgx = Channel.fromList(samples)
            .filter { s -> s.seq_type == "WGS" }
            .map { s ->
                def bam = file("${s.input_dir}/02_alignment/${s.sample_id}.aligned.sorted.bam")
                def bai = file("${s.input_dir}/02_alignment/${s.sample_id}.aligned.sorted.bam.bai")
                if (!bam.exists()) {
                    log.warn "[WARN] 找不到 BAM，CYP2D6/HLA outside call 跳過（PGx 改用 VCF 模式）：${bam}"
                    return null
                }
                tuple(s.sample_id, bam, bai)
            }
            .filter { it != null }

    } else {

        // DRAGEN：{input_dir}/vcf.gz/{sample_id}.hard-filtered.vcf.gz
        input_ch = Channel.fromList(samples).map { s ->
            def vcf = file("${s.input_dir}/vcf.gz/${s.sample_id}.hard-filtered.vcf.gz")
            if (!vcf.exists()) {
                error "[ERROR] 找不到 DRAGEN VCF：${vcf}"
            }
            tuple(s.sample_id, vcf)
        }

        PREPARE_VCF_DRAGEN(input_ch)
        snv_ch = PREPARE_VCF_DRAGEN.out.snv_ch

        dragen_mito_ch = PREPARE_VCF_DRAGEN.out.mito_ch.map { sample_id, mito_vcf, mito_tbi ->
            tuple(sample_id, "dragen", mito_vcf, mito_tbi)
        }

        MITO_ANNOTATE(dragen_mito_ch)

        dragen_str_ch = Channel.fromList(samples).map { s ->
            def str_vcf = file("${s.input_dir}/vcf.gz/${s.sample_id}.repeats.vcf.gz")
            if (!str_vcf.exists()) {
                log.warn "[WARN] 找不到 STR VCF，跳過：${str_vcf}"
                return null
            }
            tuple(s.sample_id, str_vcf)
        }
        .filter { it != null }

        STR_PARSE_DRAGEN(dragen_str_ch)

        // ── DRAGEN CNV ────────────────────────────────────────
        dragen_cnv_ch = Channel.fromList(samples).map { s ->
            def vcf = file("${s.input_dir}/vcf.gz/${s.sample_id}.cnv.vcf.gz")
            if (!vcf.exists()) {
                log.warn "[WARN] 找不到 DRAGEN CNV VCF，跳過：${vcf}"
                return null
            }
            tuple(s.sample_id, vcf)
        }
        .filter { it != null }

        PREPARE_CNV_DRAGEN(dragen_cnv_ch)
        ANNOTSV_CNV_DRAGEN(PREPARE_CNV_DRAGEN.out.cnv_filtered_ch)

        // ── DRAGEN SV ─────────────────────────────────────────
        dragen_sv_ch = Channel.fromList(samples).map { s ->
            def vcf = file("${s.input_dir}/vcf.gz/${s.sample_id}.sv.vcf.gz")
            if (!vcf.exists()) {
                log.warn "[WARN] 找不到 DRAGEN SV VCF，跳過：${vcf}"
                return null
            }
            tuple(s.sample_id, vcf)
        }
        .filter { it != null }

        PREPARE_SV_DRAGEN(dragen_sv_ch)
        ANNOTSV_SV_DRAGEN(PREPARE_SV_DRAGEN.out.sv_filtered_ch)

        // ── DRAGEN BAM channel（PGx CYP2D6/HLA 用）──────────────
        // DRAGEN 輸出 BAM 路徑依版本可能不同，嘗試兩個慣用路徑
        bam_ch_pgx = Channel.fromList(samples).map { s ->
            def bam = file("${s.input_dir}/bam/${s.sample_id}.bam")
            def bai = file("${s.input_dir}/bam/${s.sample_id}.bam.bai")
            if (!bam.exists()) {
                // 備用路徑（DRAGEN 直接放上層目錄）
                bam = file("${s.input_dir}/${s.sample_id}.bam")
                bai = file("${s.input_dir}/${s.sample_id}.bam.bai")
            }
            if (!bam.exists()) {
                log.warn "[WARN] 找不到 DRAGEN BAM，CYP2D6/HLA outside call 跳過（PGx 改用 VCF 模式）：${bam}"
                return null
            }
            tuple(s.sample_id, bam, bai)
        }
        .filter { it != null }
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

    // ── PGx annotation（PharmCAT + StellarPGx）────────────────
    // --run_pgx false（預設）→ 跳過
    // --run_pgx true  → 啟用；WGS 走 StellarPGx（CYP2D6 BAM-based），WES 直接 VCF
    if (params.run_pgx) {

        // pipeline_type lookup channel（從 samples 建立）
        ptype_ch = Channel.fromList(samples)
            .map { s -> tuple(s.sample_id, s.pipeline_type) }

        // WGS：帶 BAM 進去，StellarPGx + OptiType 在 PGX_ANNOTATE 內部跑
        pgx_wgs_ch = snv_ch
            .join(bam_ch_pgx)
            .join(ptype_ch)
            .map { sid, vcf, tbi, bam, bai, ptype ->
                tuple(sid, ptype, vcf, tbi, bam, bai)
            }

        // WES：沒有 BAM，直接進 PharmCAT（純 VCF 模式）
        pgx_wes_ch = snv_ch
            .join(bam_ch_pgx, remainder: true)
            .filter { vals -> vals[3] == null }   // 沒有 BAM 的樣本
            .join(ptype_ch)
            .map { vals ->
                def sid   = vals[0]
                def vcf   = vals[1]
                def tbi   = vals[2]
                def ptype = vals[4]
                tuple(sid, ptype, vcf, tbi)
            }

        PGX_ANNOTATE(
            pgx_wgs_ch,
            pgx_wes_ch,
            MITO_ANNOTATE.out.mito_tsv_ch.ifEmpty(Channel.empty())
        )
    }
}
