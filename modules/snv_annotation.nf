/*
 * =========================================================
 * WGS/WES Germline Analysis Pipeline - SNV Annotation Module
 * =========================================================
 * Author   : Po-Yu Lin (林伯昱)
 * Institute: Department of Neurology and
 *            Department of Genomic Medicine,
 *            National Cheng Kung University Hospital
 * Contact  : p88124019@gs.ncku.edu.tw
 *
 * Copyright (c) 2026, Po-Yu Lin
 * Licensed under the MIT License
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
 * modules/snv_annotation.nf
 * =========================
 * 目的：
 *   對 prepare_vcf 輸出的 snv_for_annotation.vcf.gz 進行兩階段 annotation：
 *
 *   Process 1 - VEP_ANNOTATE：
 *     執行 Ensembl VEP 115，帶入以下 plugins 和 custom annotation：
 *       - dbNSFP 4.9c（BayesDel、AlphaMissense、ESM1b 等）
 *       - LOFTEE GRCh38（LoF HC/LC 判定，PVS1 基礎）
 *         * loftee 資料目錄透過 --bind 掛載至 /opt/vep/Plugins/loftee_data
 *         * gerp BigWig 透過 --bind 直接掛載至 /opt/vep/Plugins/（LOFTEE bug workaround）
 *       - LoFtool（基因 LoF 容忍度）
 *       - ClinVar 20260418（CLNSIG、CLNREVSTAT、CLNDN、CLNSIGCONF）
 *       - gnomAD（VEP cache 內建，--af_gnomadg + --af_gnomade）
 *       - 1000 Genomes（VEP cache 內建，--af_1kg，提供 EAS_AF 等次族群頻率）
 *
 *     annotation 旗標（CSQ 欄位內容）：
 *       - --hgvs    ：輸出 HGVSc / HGVSp（TSV 主表必要欄位，無此旗標欄位為空）
 *       - --symbol  ：輸出 HGNC gene symbol（GENE 欄位，無此旗標只有 Ensembl ID）
 *       - --numbers ：輸出 EXON / INTRON 編號（格式 2/15，PVS1 NMD escape 判斷需要）
 *       - --canonical：輸出 CANONICAL 旗標（transcript 備用選取依據）
 *       - --biotype ：輸出 BIOTYPE（PVS1 需確認 protein_coding）
 *       - --tsl     ：輸出 TSL（--pick_order 的 tsl 層才有資料可比較）
 *       - --appris  ：輸出 APPRIS（--pick_order 的 appris 層才有資料可比較）
 *       - --uniprot ：輸出 SWISSPROT / TREMBL（PM1 UniProt domain 判斷依據）
 *       - --domains ：輸出 DOMAINS（PM1 蛋白質 domain 資訊）
 *       - --safe    ：plugin 報錯時 VEP 正確退出（臨床 pipeline 必要）
 *
 *     transcript 策略：
 *       - --mane：標記 MANE Select / MANE Plus Clinical
 *       - --flag_pick：標記代表 transcript（PICK 欄位值為 1）
 *       - --pick_order：mane_select > mane_plus_clinical > canonical > appris > tsl > ...
 *       - 注意：不用 --pick（--pick 只保留一個 transcript，無法展開所有 MANE）
 *               只用 --flag_pick，保留所有 transcript，以 PICK=1 標記代表 transcript
 *     PICK 欄位說明：
 *       VCF 的 CSQ 欄位中，PICK 是第 22 個 pipe-separated 值（從 1 開始數）
 *       值為 "1" 表示此 transcript 是代表 transcript，值為空表示非代表
 *       加入 --hgvs/--symbol 等旗標後，PICK 的欄位位置會往後移，
 *       parse_vep_csq.py 應從 VCF header 的 CSQ 欄位定義動態解析位置，不要硬編碼
 *
 *   Process 2 - PANGOLIN_SCORE：
 *     對 VEP 輸出中帶有 splice consequence 的 variant 執行 Pangolin inference。
 *     取代 SpliceAI（SpliceAI 商用需付費授權，Pangolin 為 GPL-3.0）。
 *     GPU inference（--nv 由 config 的 process_gpu label containerOptions 設定）。
 *
 * 輸入：
 *   tuple val(sample_id), path(snv_vcf), path(snv_tbi)
 *
 * 輸出：
 *   vep_ch:      tuple val(sample_id), path("*.vep.vcf.gz"), path("*.vep.vcf.gz.tbi")
 *   pangolin_ch: tuple val(sample_id), path("*.pangolin.vcf.gz"), path("*.pangolin.vcf.gz.tbi")
 *
 * 已知踩雷：
 *   - --stats_file 和 --no_stats 不能同時使用，本 module 統一用 --no_stats（較快）
 *   - LOFTEE gerp_dist.pl bug：路徑拼接缺少 /，需要透過 --bind 直接掛載 gerp bw
 *   - tabix 必須在 VEP 輸出之後才能建 index（VEP 加 --compress_output bgzip 自動壓縮）
 *   - 加入 --hgvs/--symbol/--numbers 等旗標後，CSQ 欄位的 column 數增加，
 *     PICK 的欄位位置不再是固定第 22 欄，parse_vep_csq.py 必須從 VCF header
 *     動態解析 CSQ format（##INFO=<ID=CSQ,...,Format="..."> 那行）
 */

// ──────────────────────────────────────────────────────────────
// Process 1：VEP Annotation
// ──────────────────────────────────────────────────────────────

process VEP_ANNOTATE {
    label 'process_high'

    container "${params.sif_dir}/vep_115.sif"

    // LOFTEE bind mount 設定：
    //   loftee_data → 整個 loftee 資料目錄掛載至容器內 /opt/vep/Plugins/loftee_data
    //                 （含 human_ancestor.fa.gz、loftee.sql、gerp_...GRCh38.bw）
    //   gerp bw     → 另外直接掛載至 /opt/vep/Plugins/（保留舊 gerp_dist.pl workaround）
    //
    // GERP 現在以 gerp_bigwig: 參數明確傳入 LoF plugin（指向 loftee_data 內的 bw），
    // 不再只依賴檔名拼接的 workaround。若無此 bw，LoF 的 GERP-based flag 不會計算。
    //
    // 注意：apptainer_base_opts 提供 --bind /scratch,/data 基礎掛載
    //       這裡在基礎上追加 loftee 相關的 bind
    //   NMD.pm      → Ensembl 官方 NMD plugin（VEP_plugins repo），掛進 plugin 目錄。
    //                 用來預測「提前終止密碼子是否會逃過 NMD」，是 ClinGen SVI PVS1
    //                 決策樹的第一個分支。容器內 /opt/vep/Plugins 唯讀，故比照 gerp bw
    //                 以單檔 bind 方式掛入，不需重建容器。
    containerOptions "${params.apptainer_base_opts} \
        --bind ${params.loftee_dir}:/opt/vep/Plugins/loftee_data \
        --bind ${params.loftee_dir}/gerp_conservation_scores.homo_sapiens.GRCh38.bw:/opt/vep/Plugins/gerp_conservation_scores.homo_sapiens.GRCh38.bw \
        --bind ${params.loftee_dir}/NMD.pm:/opt/vep/Plugins/NMD.pm"

    publishDir "${params.out_dir}/${sample_id}/01_vep", mode: 'copy'

    input:
    tuple val(sample_id), path(snv_vcf), path(snv_tbi)

    output:
    tuple val(sample_id),
          path("${sample_id}.vep.vcf.gz"),
          path("${sample_id}.vep.vcf.gz.tbi"),
          emit: vep_out

    script:
    // ── dbNSFP 版本切換（--academic_dbnsfp，預設 false）──────────────────────
    //   false → dbNSFP 4.9c（預設路徑，維持全部工具可商用）
    //   true  → dbNSFP 5.3a，並額外抓 REVEL / MutPred2 / VEST4 / CADD_phred
    //           ★ 這些多為「學術免費、商業需另行授權」（CADD 尤其明確），故只在此模式取用。
    //   兩版都已把 P-KNN 合併進去（PKNN_LLR 為最後一欄），GUI 排序訊號不受影響；
    //   5.3a 的 P-KNN 覆蓋更完整（P-KNN 本來就是以 dbNSFP 5.3 產生）。
    //   族群頻率欄名在 5.3a 改了：gnomAD_exomes_* → gnomAD4.1_joint_*（gnomAD 4.1，
    //   exomes+genomes 合併），parse_vep_csq.py 以 get_any() 相容兩種欄名。
    // ⚠️ 不要用 `as boolean`：Groovy 對非空字串一律為 true，指令列傳 "--academic_dbnsfp false"
    //    會被誤判成開啟。一律轉小寫字串比對。
    def academic   = params.academic_dbnsfp.toString().toLowerCase() == 'true'
    def dbnsfp_f   = academic ? params.dbnsfp_academic : params.dbnsfp
    // ACMG 用的族群頻率：兩版都取 gnomAD 2.1.1 exomes，讓 PM2 的判讀基準不因換版而變。
    //   4.9c 有「整體」欄位；5.3a 只保留 controls / non_neuro / non_cancer 三個子集，
    //   取 non_cancer（約 118k，最接近整體的 ~125k；controls 僅約 60k 差距過大）。
    def dbnsfp_af  = academic
        ? "gnomAD2.1.1_exomes_non_cancer_AF,gnomAD2.1.1_exomes_non_cancer_EAS_AF"
        : "gnomAD_exomes_AF,gnomAD_exomes_EAS_AF"
    // 5.3a 額外抓的「參考用」欄位（不進 ACMG 計分）：新 in-silico 工具 + gnomAD 4.1。
    def dbnsfp_extra = academic
        ? ",REVEL_score,MutPred2_score,MutPred2_pred,VEST4_score,CADD_phred" +
          ",gnomAD4.1_joint_AF,gnomAD4.1_joint_EAS_AF"
        : ""
    """
    echo "[VEP_ANNOTATE] dbNSFP = ${dbnsfp_f}" >&2

    vep \\
        --input_file ${snv_vcf} \\
        --output_file ${sample_id}.vep.vcf.gz \\
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
        --total_length \\
        --canonical \\
        --biotype \\
        --tsl \\
        --appris \\
        --uniprot \\
        --domains \\
        \\
        --mane \\
        --flag_pick \\
        --pick_order mane_select,mane_plus_clinical,canonical,appris,tsl,biotype,ccds,rank,length \\
        \\
        --plugin dbNSFP,${dbnsfp_f},\\
BayesDel_noAF_score,BayesDel_noAF_pred,\\
AlphaMissense_score,AlphaMissense_pred,\\
ESM1b_score,ESM1b_pred,\\
VARITY_R_score,\\
SIFT_score,SIFT_pred,\\
DANN_score,\\
PHACTboost_score,\\
phyloP100way_vertebrate,\\
GERP++_RS,\\
${dbnsfp_af},\\
PKNN_LLR${dbnsfp_extra} \\
        \\
        --plugin LoF,\\
loftee_path:/opt/vep/Plugins/,\\
human_ancestor_fa:/opt/vep/Plugins/loftee_data/human_ancestor.fa.gz,\\
conservation_file:/opt/vep/Plugins/loftee_data/loftee.sql,\\
gerp_bigwig:/opt/vep/Plugins/loftee_data/gerp_conservation_scores.homo_sapiens.GRCh38.bw \\
        \\
        --plugin LoFtool,/opt/vep/Plugins/loftee_data/LoFtool_scores.txt \\
        \\
        --plugin NMD \\
        \\
        --custom file=${params.clinvar},short_name=ClinVar,format=vcf,type=exact,coords=0,fields=CLNSIG%CLNREVSTAT%CLNDN%CLNSIGCONF \\
        \\
        --af_gnomadg \\
        --af_gnomade \\
        --af_1kg \\
        --force_overwrite \\
        --no_stats \\
        --safe

    # 建立 tabix index（VEP 已用 --compress_output bgzip 壓縮，直接 tabix）
    tabix -p vcf ${sample_id}.vep.vcf.gz

    echo "[VEP_ANNOTATE] ${sample_id} 完成" >&2
    bcftools stats ${sample_id}.vep.vcf.gz | grep "^SN" >&2
    """
}

// ──────────────────────────────────────────────────────────────
// Process 2：Pangolin Splice Scoring
// ──────────────────────────────────────────────────────────────

process PANGOLIN_SCORE {
    label 'process_gpu'

    // ── 容器依機器的 GPU 世代切換（params.pangolin_sif，每個 profile 都要宣告）──
    //   Pangolin 用 PyTorch，而 PyTorch 的 prebuilt wheel 只包特定 compute
    //   capability 的 kernel，跨 GPU 世代不相容（cubin 只在同一個 major 版本內
    //   向前相容）。實測官方 wheel 的 arch list：
    //     cu121 → sm_50…sm_90   （含 sm_70 = V100，無 sm_120）
    //     cu128 → sm_75…sm_120  （含 sm_120 = Blackwell，**已砍 sm_70**；cu130 同理）
    //   → 沒有任何一顆 prebuilt wheel 同時涵蓋兩者，只能分成兩顆容器：
    //     pangolin_cu121_1.0.0.sif → production：dgm + dgx（V100 = sm_70）
    //     pangolin_cu130_1.0.0.sif → 開發機 local（RTX PRO 6000 Blackwell = sm_120）
    //   檔名標的是 CUDA 變體而非機器，因為差異的本質就是 wheel 的 CUDA 變體。
    //   建置方式與 arch 守門員見 DEVELOPMENT_NOTES.md「容器建立 → Pangolin」。
    container "${params.sif_dir}/${params.pangolin_sif}"

    publishDir "${params.out_dir}/${sample_id}/02_pangolin", mode: 'copy'

    input:
    tuple val(sample_id), path(vep_vcf), path(vep_tbi)

    output:
    tuple val(sample_id),
          path("${sample_id}.pangolin.vcf.gz"),
          path("${sample_id}.pangolin.vcf.gz.tbi"),
          emit: pangolin_out

    script:
    // ── GPU lock（DGX-2 共用機器）─────────────────────────────────────────
    //   與二級分析同一套機制：gpu_lock.sh 搶下 N 張空閒 V100 並輸出可 eval 的環境變數，
    //   trap EXIT 時以 gpu_unlock.sh 歸還（即使 process 失敗也會釋放）。
    //   Apptainer 預設會把 host 環境變數帶進容器，故 CUDA_VISIBLE_DEVICES 會被 Pangolin
    //   看見；這裡再明確 export 一次，避免 lock script 只設了 MY_GPUS。
    //   非 DGX（local / dgm）use_gpu_lock=false → 走 config 既有的 GPU 設定，不插入此段。
    // ── GPU / CPU 切換（--use_gpu_pangolin）──────────────────────────────
    //   Pangolin 以 PyTorch 實作，會自動偵測 CUDA；把 CUDA_VISIBLE_DEVICES 設為空字串
    //   即可強制走 CPU（慢但結果相同），讓 pipeline 能部署到沒有 GPU 的環境。
    //   容器的 --nv 由 config 的 process_gpu label 依同一參數決定是否加上。
    def use_gpu       = params.use_gpu_pangolin == null ? true
                        : params.use_gpu_pangolin.toString().toLowerCase() == 'true'
    // ⚠️ use_gpu_lock / pangolin_num_gpus 只宣告在 profile 內（不能放全域 params 區塊，
    //    那個區塊在 profiles 之後會蓋掉 profile 的值）。但 Nextflow 只要「讀到」未宣告的
    //    param 就會噴 "WARN: Access to undefined parameter" —— 即使後面接了 ?: 預設值也一樣。
    //    所以用 containsKey 先問再讀：三個 profile 都已宣告，這層只是保險，
    //    讓「沒帶 -profile」或「將來新增 profile 忘了宣告」時不會噴 warning 也不會壞。
    // 不用 GPU 就不必搶卡
    def use_lock      = use_gpu &&
                        (params.containsKey('use_gpu_lock') ? params.use_gpu_lock : false)
    def lock_script   = params.gpu_lock_script
    def unlock_script = params.gpu_unlock_script
    def num_gpus      = params.containsKey('pangolin_num_gpus') ? params.pangolin_num_gpus : 1
    def gpu_block = !use_gpu ? """
    export CUDA_VISIBLE_DEVICES=""
    echo "[PANGOLIN] ${sample_id} CPU 模式（use_gpu_pangolin=false）" >&2
    """ : use_lock ? """
    eval \$(bash ${lock_script} ${num_gpus})
    export CUDA_VISIBLE_DEVICES=\${MY_GPUS}
    echo "[PANGOLIN] ${sample_id} 取得 GPU \${MY_GPUS}" >&2
    trap "bash ${unlock_script} \${MY_GPUS}; echo '[PANGOLIN] ${sample_id} 釋放 GPU \${MY_GPUS}' >&2" EXIT
    """ : """
    echo "[PANGOLIN] ${sample_id} 使用 config 指定的 GPU（未啟用 GPU lock）" >&2
    """
    """
    ${gpu_block}

    # Step 1：從 VEP 輸出中篩選 splice candidate
    # bcftools view -h → 只取 header 行
    # bcftools view -H → 只取 variant 行，awk 篩 INFO 欄含 "splice" 字眼
    bcftools view -h ${vep_vcf} > splice_header.vcf
    bcftools view -H ${vep_vcf} | \\
        awk -F'\\t' '\$8 ~ /splice/' > splice_body.vcf

    SPLICE_COUNT=\$(wc -l < splice_body.vcf)
    echo "[PANGOLIN] ${sample_id}：splice candidates = \${SPLICE_COUNT}" >&2

    if [ "\${SPLICE_COUNT}" -eq 0 ]; then
        # 無 splice candidate，產生空輸出
        echo "[PANGOLIN] 無 splice candidate，產生空輸出" >&2
        bgzip -c splice_header.vcf > ${sample_id}.pangolin.vcf.gz

    else
        # Step 2：移除 CSQ + 過濾 alt/random contig，再送 Pangolin
        # 原因 1：WGS CSQ 可能長達 270KB，Pangolin parse 時 segfault
        # 原因 2：Pangolin gencode DB 沒有 alt/random contig 的 gene model，
        #         碰到 chr*_alt / chr*_random / chrUn_* 也會 segfault
        # bcftools view --regions 需要 bgzip+tabix 輸入，無法直接接 pipe
        # 改用 grep + awk：先移除 CSQ，再只保留標準染色體（chr1-22,X,Y,M）
        # alt/random contig 的 Pangolin gencode DB 沒有 gene model，會 segfault
        cat splice_header.vcf splice_body.vcf \
            | bcftools annotate -x INFO/CSQ \
            | grep -E '^#|^chr([0-9]+|[XYM])\t' \
            > splice_for_pangolin.vcf

        # Step 3：執行 Pangolin（GPU）
        # output_file 是 prefix，Pangolin 自動加 .vcf
        # 不要給 ${sample_id}.pangolin.vcf，否則輸出變成 .vcf.vcf
        # IUPAC ambiguity code（Y/R/W 等）已在 sif build 時 patch 處理，
        # 碰到非 ACGT 字元一律當 N（全零 encoding），不會 crash
        pangolin \\
            splice_for_pangolin.vcf \\
            ${params.ref_fasta} \\
            ${params.pangolin_db} \\
            ${sample_id}.pangolin \\
            --distance 50 \\
            -m True

        bgzip -c ${sample_id}.pangolin.vcf > ${sample_id}.pangolin.vcf.gz
    fi

    tabix -p vcf ${sample_id}.pangolin.vcf.gz
    echo "[PANGOLIN] ${sample_id} 完成" >&2
    """
}

// ──────────────────────────────────────────────────────────────
// 組合 workflow（供 main_tertiary.nf 呼叫）
// ──────────────────────────────────────────────────────────────

workflow SNV_ANNOTATE {
    take:
    snv_ch   // tuple val(sample_id), path(snv_vcf), path(snv_tbi)

    main:
    VEP_ANNOTATE(snv_ch)
    PANGOLIN_SCORE(VEP_ANNOTATE.out.vep_out)

    emit:
    vep_ch      = VEP_ANNOTATE.out.vep_out
    pangolin_ch = PANGOLIN_SCORE.out.pangolin_out
}
