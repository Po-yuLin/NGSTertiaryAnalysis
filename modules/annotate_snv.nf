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
 * modules/annotate_snv.nf
 * =======================
 * SNV/indel annotation 尾段 sub-workflow（NCKUH / DRAGEN 兩 pipeline 共用）：
 *   SNV_ANNOTATE（VEP 115 + Pangolin）→ PARSE_VEP_CSQ（transcript 選取 + TSV）
 *   → ACMG_CLASSIFY（ACMG evidence 分類）。
 *
 * 對外吃 snv_ch（前處理後的 SNV VCF）＋ ClinGen HI / gene MOI 兩個資料庫檔（main 依 params
 * 決定是否為 NO_FILE）。跨三個 module，故獨立成一個 module 檔（與二級 call_snv.nf 同理）。
 *
 * emit：
 *   acmg     → 最終 ACMG 分類 TSV（發布；報告彙整用）
 *   full_tsv → VEP CSQ 解析後的完整 TSV
 */

include { SNV_ANNOTATE  } from './snv_annotation.nf'
include { PARSE_VEP_CSQ } from './parse_csq.nf'
include { ACMG_CLASSIFY } from './acmg_classifier.nf'

workflow ANNOTATE_SNV {
    take:
    snv_ch              // tuple(sample_id, vcf, tbi)：前處理後的 SNV VCF
    clingen_hi_file     // ClinGen HI TSV 或 NO_FILE
    gene_moi_file       // gene MOI TSV 或 NO_FILE

    main:
    SNV_ANNOTATE(snv_ch)
    PARSE_VEP_CSQ(
        SNV_ANNOTATE.out.vep_ch,
        SNV_ANNOTATE.out.pangolin_ch
    )
    ACMG_CLASSIFY(
        PARSE_VEP_CSQ.out.full_tsv_ch,
        clingen_hi_file,
        gene_moi_file
    )

    emit:
    acmg     = ACMG_CLASSIFY.out.acmg_tsv
    full_tsv = PARSE_VEP_CSQ.out.full_tsv_ch
}
