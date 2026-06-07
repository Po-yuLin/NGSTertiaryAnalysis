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
 * modules/pgx_annotation.nf
 * =========================
 * 目的：
 *   Pharmacogenomics（PGx）分析：
 *     - 大部分基因（CYP2C19、CYP2C9、DPYD、TPMT、NUDT15、SLCO1B1 等）
 *       直接由 PharmCAT 從 WGS/WES VCF 呼叫，不需要外部工具。
 *     - CYP2D6：大量 structural variant → 用 StellarPGx（BAM-based）
 *     - HLA-A/HLA-B：高度多型區域 → 用 OptiType（BAM-based）
 *     - MT-RNR1：直接用現有 mito pipeline 輸出（mito_tsv）
 *
 * License：
 *   PharmCAT   MPL 2.0        ✅ 商業可用
 *   StellarPGx Open source    ✅ PharmCAT 官方推薦 CYP2D6 caller
 *   OptiType   BSD-3-Clause   ✅ 商業可用
 *   ❌ BCyrius（PolyForm Strict）、Aldy（non-commercial）禁止使用
 *
 * 包含的 Process：
 *   PGX_STELLARPGX  BAM → StellarPGx → CYP2D6 diplotype TSV（選用）
 *   PGX_OPTITYPE    BAM → OptiType   → HLA-A/HLA-B typing TSV（選用）
 *   PGX_PHARMCAT    VCF + outside_calls.tsv → pharmcat_vcf_preprocessor + pharmcat.jar -po → JSON + TSV
 *   PGX_PARSE       PharmCAT JSON + mito_tsv → parse_pgx_report.py → pgx.tsv
 *
 * PharmCAT 3.2.0 容器內已備妥：
 *   /pharmcat/pharmcat_pipeline          ← 整合 pipeline（preprocessor + matcher + phenotyper + reporter）
 *   /pharmcat/pharmcat_vcf_preprocessor  ← 單獨 preprocessor（已是執行檔，非 Python script）
 *   /pharmcat/pharmcat.jar               ← 主程式
 *   /pharmcat/reference.fna.bgz          ← 內建 GRCh38 reference（不需外掛）
 *   /pharmcat/pharmcat_positions.vcf.bgz ← 內建 PGx 位點定義
 *
 * 輸出（07_pgx/）：
 *   {sample_id}.pgx.tsv                  ← 臨床報告用（17 欄，CPIC Level A 基因）★
 *   {sample_id}.pharmcat.report.json     ← PharmCAT 完整輸出（歸檔用）
 *   {sample_id}.outside_calls.tsv        ← 整合的 outside calls（歸檔用）
 *   {sample_id}.stellarpgx.tsv           ← CYP2D6 diplotype（歸檔用，有 BAM 才有）
 *   {sample_id}.optitype.tsv             ← HLA typing（歸檔用，有 BAM 才有）
 *
 * 踩雷記錄：
 *   - pharmcat_vcf_preprocessor 在 3.x 已是編譯好的執行檔（非 Python script）
 *   - pharmcat_pipeline 沒有 -po 參數，必須拆成 preprocessor + pharmcat.jar 兩步
 *   - 容器內建 /pharmcat/reference.fna.bgz，不需要外掛 ref_fasta
 *   - outside calls TSV 格式（tab 分隔）：Gene\tDiplotype\tFunction\tSource
 *   - -reporterJson 才會輸出 JSON；預設只輸出 HTML
 *   - -rs CPIC,DPWG 限制 recommendation 來源
 *   - outside calls 的 Gene 欄位要與 PharmCAT 內部名稱一致（CYP2D6、HLA-A、HLA-B）
 */

// ──────────────────────────────────────────────────────────────
// Process 1：StellarPGx — CYP2D6 outside caller（BAM-based）
// ──────────────────────────────────────────────────────────────

process PGX_STELLARPGX {

    label 'process_medium'

    container "${params.sif_dir}/stellarpgx_graphtyper2.5.1.sif"

    containerOptions "${params.apptainer_base_opts}"

    publishDir "${params.out_dir}/${sample_id}/07_pgx", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id),
          path("${sample_id}.stellarpgx.tsv"),
          emit: stellarpgx_ch

    script:
    // StellarPGx 執行流程（移植自 StellarPGx/main.nf，hg38 模式）：
    //   call_snvs1/2 → format_snvs → get_core_var → analyse_1/2/3 → call_stars
    // 所有步驟在同一個 process 內串接，避免 nested Nextflow 問題。
    // 資料庫路徑：
    //   params.stellarpgx_db      → StellarPGx/database/cyp2d6/hg38/
    //   params.stellarpgx_res     → StellarPGx/resources/cyp2d6/res_hg38/
    //   params.stellarpgx_scripts → StellarPGx/scripts/cyp2d6/hg38/bin/
    //   params.stellarpgx_res_base → StellarPGx/resources/（annotation/ 在這一層）
    //
    // CYP2D6 hg38 座標：
    //   region_a = chr22:42126000-42137500（graphtyper 輸入範圍）
    //   region_b = chr22:42126300-42132400（最終變異過濾範圍）
    """
    echo "[PGX_STELLARPGX] ${sample_id}：CYP2D6 star allele calling（StellarPGx 1.2.8）" >&2

    REF_DIR=\$(dirname ${params.ref_fasta})
    REF_NAME=\$(basename ${params.ref_fasta})
    # (例: NA12878_WGS.aligned.sorted.bam)
    STAGED_BAM=${bam}
    STAGED_BAI=${bai}

    DB=${params.stellarpgx_db}
    RES=${params.stellarpgx_res}
    RES_BASE=${params.stellarpgx_res_base}
    CALLER=${params.stellarpgx_scripts}

    # hg38 CYP2D6 座標
    CHROM="chr22"
    REGION_A1="chr22:42126000-42137500"
    REGION_A2="042126000-042137500"
    REGION_B1="chr22:42126300-42132400"
    REGION_B2="042126300-042132400"
    TRANSCRIPT="ENST00000645361"
    DEBUG38="--minimum_extract_score_over_homref=0"

    # ── Step 1a：call_snvs1（with prior VCF）────────────────
    echo "[PGX_STELLARPGX] Step 1a：graphtyper genotype（prior VCF）" >&2
    graphtyper genotype \
        \${REF_DIR}/\${REF_NAME} \
        --sam=\${STAGED_BAM} \
        --sams_index=<(echo \${STAGED_BAI}) \
        --region=\${REGION_A1} \
        --output=var_1 \
        --prior_vcf=\${RES}/common_plus_core_var.vcf.gz \
        -a \${DEBUG38}

    bcftools concat var_1/\${CHROM}/*.vcf.gz > var_1/\${CHROM}/\${REGION_A2}.vcf
    bgzip -f var_1/\${CHROM}/\${REGION_A2}.vcf
    tabix -f var_1/\${CHROM}/\${REGION_A2}.vcf.gz

    # ── Step 1b：call_snvs2（without prior VCF）──────────────
    echo "[PGX_STELLARPGX] Step 1b：graphtyper genotype（no prior）" >&2
    graphtyper genotype \
        \${REF_DIR}/\${REF_NAME} \
        --sam=\${STAGED_BAM} \
        --sams_index=<(echo \${STAGED_BAI}) \
        --region=\${REGION_A1} \
        --output=var_2 \
        -a \${DEBUG38}

    bcftools concat var_2/\${CHROM}/*.vcf.gz > var_2/\${CHROM}/\${REGION_A2}.vcf
    bgzip -f var_2/\${CHROM}/\${REGION_A2}.vcf
    tabix -f var_2/\${CHROM}/\${REGION_A2}.vcf.gz

    # ── Step 2：call_sv_del + call_sv_dup ────────────────────
    echo "[PGX_STELLARPGX] Step 2：graphtyper genotype_sv" >&2
    graphtyper genotype_sv \
        \${REF_DIR}/\${REF_NAME} \
        --sam=\${STAGED_BAM} \
        --region=\${REGION_A1} \
        --output=sv_del \
        \${RES}/sv_test.vcf.gz

    graphtyper genotype_sv \
        \${REF_DIR}/\${REF_NAME} \
        --sam=\${STAGED_BAM} \
        --region=\${REGION_A1} \
        --output=sv_dup \
        \${RES}/sv_test3.vcf.gz

    # ── Step 3：get_depth ────────────────────────────────────
    echo "[PGX_STELLARPGX] Step 3：samtools bedcov" >&2
    samtools bedcov \
        --reference \${REF_DIR}/\${REF_NAME} \
        \${RES}/test3.bed \
        \${STAGED_BAM} \
        > ${sample_id}_cyp2d6_ctrl.depth

    # ── Step 4：format_snvs ──────────────────────────────────
    echo "[PGX_STELLARPGX] Step 4：format_snvs" >&2
    mkdir -p all_var
    bcftools isec -p all_var -Oz \
        var_1/\${CHROM}/\${REGION_A2}.vcf.gz \
        var_2/\${CHROM}/\${REGION_A2}.vcf.gz

    bcftools concat -a -D -r \${REGION_B1} \
        all_var/0000.vcf.gz all_var/0001.vcf.gz all_var/0002.vcf.gz \
        -Oz -o all_var/${sample_id}_\${REGION_B2}.vcf.gz
    tabix all_var/${sample_id}_\${REGION_B2}.vcf.gz

    bcftools norm -m - all_var/${sample_id}_\${REGION_B2}.vcf.gz \
        | bcftools view -e 'GT="1/0"' \
        | bcftools view -e 'GT="0/0"' \
        | bcftools view -e 'FILTER="PASS" & INFO/QD<10 || 0<ABHet<0.25' \
        | bgzip -c > all_var/${sample_id}_all_norm.vcf.gz
    tabix all_var/${sample_id}_all_norm.vcf.gz

    # ── Step 5：get_core_var（hg38：bcftools csq + isec）─────
    echo "[PGX_STELLARPGX] Step 5：get_core_var" >&2
    mkdir -p core_int

    bcftools csq -p m -v 0 \
        -f \${REF_DIR}/\${REF_NAME} \
        -g \${RES_BASE}/annotation/Homo_sapiens.GRCh38.110.gff3.gz \
        all_var/${sample_id}_all_norm.vcf.gz \
        -o all_var/${sample_id}_all_norm_annot.vcf
    bgzip all_var/${sample_id}_all_norm_annot.vcf
    tabix all_var/${sample_id}_all_norm_annot.vcf.gz

    bcftools isec \
        all_var/${sample_id}_all_norm_annot.vcf.gz \
        \${RES}/allele_def_var.vcf.gz \
        -p core_int -Oz

    bcftools norm -m - core_int/0002.vcf.gz \
        | bcftools view -e 'GT="1/0"' \
        | bcftools view -e 'GT="0/0"' \
        > core_int/${sample_id}_core_int1.vcf

    bgzip -d core_int/0000.vcf.gz
    python3 \${CALLER}/../../../novel/core_var.py \
        core_int/0000.vcf CYP2D6 \${TRANSCRIPT} \
        >> core_int/${sample_id}_core_int1.vcf

    bcftools sort core_int/${sample_id}_core_int1.vcf -T core_int \
        | bgzip -c > core_int/${sample_id}_core.vcf.gz
    tabix core_int/${sample_id}_core.vcf.gz

    # ── Step 6：analyse_1/2/3 ────────────────────────────────
    echo "[PGX_STELLARPGX] Step 6：analyse" >&2
    bcftools query \
        -f'%ID\t%ALT\t[%GT\t%DP]\t%INFO/ABHet\t%INFO/ABHom\n' \
        sv_del/\${CHROM}/\${REGION_A2}.vcf.gz \
        > ${sample_id}_gene_del_summary.txt

    bcftools query \
        -f'%POS~%REF>%ALT\t[%GT\t%DP]\t%INFO/ABHet\t%INFO/ABHom\n' \
        -i'GT="alt"' \
        sv_dup/\${CHROM}/\${REGION_A2}.vcf.gz \
        > ${sample_id}_gene_dup_summary.txt
    bcftools query \
        -f'%POS~%REF>%ALT\t[%GT\t%DP]\t%INFO/ABHet\t%INFO/ABHom\n' \
        -i'GT="alt"' \
        core_int/${sample_id}_core.vcf.gz \
        >> ${sample_id}_gene_dup_summary.txt

    bcftools query \
        -f'[%POS~%REF>%ALT~%GT\n]' \
        core_int/${sample_id}_core.vcf.gz \
        > ${sample_id}_core_snvs.dip
    bcftools query \
        -f '%POS~%REF>%ALT\n' \
        all_var/${sample_id}_all_norm.vcf.gz \
        > ${sample_id}_full.dip
    bcftools query \
        -f'[%POS~%REF>%ALT~%GT\n]' \
        all_var/${sample_id}_all_norm.vcf.gz \
        > ${sample_id}_gt.dip

    # ── Step 7：call_stars → .alleles ────────────────────────
    echo "[PGX_STELLARPGX] Step 7：call_stars（stellarpgx.py）" >&2
    python3 \${CALLER}/stellarpgx.py \
        \${DB}/diplo_db_debugged2.dbs \
        ${sample_id}_core_snvs.dip \
        ${sample_id}_full.dip \
        ${sample_id}_gt.dip \
        \${DB}/genotypes4.dbs \
        ${sample_id}_gene_del_summary.txt \
        ${sample_id}_gene_dup_summary.txt \
        ${sample_id}_cyp2d6_ctrl.depth \
        \${DB}/haps_var_new.dbs \
        \${DB}/a_scores.dbs \
        > ${sample_id}_cyp2d6.alleles

    echo "[PGX_STELLARPGX] .alleles 內容：" >&2
    cat ${sample_id}_cyp2d6.alleles >&2

    # ── Step 8：parse .alleles → TSV ─────────────────────────
    python3 ${params.scripts_dir}/parse_stellarpgx.py \
        --input  ${sample_id}_cyp2d6.alleles \
        --sample ${sample_id} \
        --output ${sample_id}.stellarpgx.tsv

    echo "[PGX_STELLARPGX] ${sample_id} 完成" >&2
    cat ${sample_id}.stellarpgx.tsv >&2
    """
}

// ──────────────────────────────────────────────────────────────
// Process 2：OptiType — HLA-A/HLA-B outside caller（BAM-based）
// ──────────────────────────────────────────────────────────────

process PGX_OPTITYPE {

    label 'process_medium'

    container "${params.sif_dir}/optitype_1.0.0.sif"

    containerOptions "${params.apptainer_base_opts}"

    publishDir "${params.out_dir}/${sample_id}/07_pgx", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id),
          path("${sample_id}.optitype.tsv"),
          emit: optitype_ch

    script:
    """
    echo "[PGX_OPTITYPE] ${sample_id}：HLA-A/HLA-B typing" >&2

    # ── Step 1：HLA reads 擷取（razers3）─────────────────────
    HLA_REF=\$(find /opt /usr -name "hla_reference_dna.fasta" 2>/dev/null | head -1)
    if [ -z "\$HLA_REF" ]; then
        echo "[PGX_OPTITYPE] 錯誤：找不到 hla_reference_dna.fasta" >&2
        exit 1
    fi
    echo "[PGX_OPTITYPE] HLA reference：\$HLA_REF" >&2

    samtools view -h ${bam} \
        | razers3 \
            --percent-identity 90 \
            --max-hits 1 \
            --distance-range 0 \
            --output hla_reads.bam \
            \$HLA_REF /dev/stdin

    samtools fastq hla_reads.bam > hla_reads.fastq

    READ_COUNT=\$(wc -l < hla_reads.fastq)
    echo "[PGX_OPTITYPE] HLA reads 數：\$((READ_COUNT / 4))" >&2

    if [ "\$((READ_COUNT / 4))" -eq 0 ]; then
        echo "[PGX_OPTITYPE] 警告：沒有 HLA reads，產生空白結果" >&2
        printf "GENE\tALLELE_1\tALLELE_2\tSOURCE\n" > ${sample_id}.optitype.tsv
        printf "HLA-A\t.\t.\tOptiType\n"             >> ${sample_id}.optitype.tsv
        printf "HLA-B\t.\t.\tOptiType\n"             >> ${sample_id}.optitype.tsv
    else
        # ── Step 2：OptiType typing ─────────────────────────
        OptiTypePipeline.py \
            --input hla_reads.fastq \
            --dna \
            --outdir optitype_out \
            --prefix ${sample_id}

        RESULT_FILE=\$(find optitype_out -name "*_result.tsv" | head -1)

        if [ -z "\$RESULT_FILE" ]; then
            echo "[PGX_OPTITYPE] 警告：找不到 OptiType 輸出，產生空白結果" >&2
            printf "GENE\tALLELE_1\tALLELE_2\tSOURCE\n" > ${sample_id}.optitype.tsv
            printf "HLA-A\t.\t.\tOptiType\n"             >> ${sample_id}.optitype.tsv
            printf "HLA-B\t.\t.\tOptiType\n"             >> ${sample_id}.optitype.tsv
        else
            echo "[PGX_OPTITYPE] 找到結果：\$RESULT_FILE" >&2
            python3 ${params.scripts_dir}/parse_optitype.py \
                --input  "\$RESULT_FILE" \
                --sample ${sample_id} \
                --output ${sample_id}.optitype.tsv
        fi
    fi

    echo "[PGX_OPTITYPE] ${sample_id} 完成" >&2
    cat ${sample_id}.optitype.tsv >&2
    """
}

// ──────────────────────────────────────────────────────────────
// Process 3：pharmcat_vcf_preprocessor + pharmcat.jar（matcher + phenotyper + reporter）
// ──────────────────────────────────────────────────────────────
// pharmcat_pipeline 用法（3.2.0 確認）：
//   pharmcat_pipeline <input.vcf.gz> \
//     -s <sample_id> \
//     -po outside_calls.tsv \         ← outside calls（-po flag，與 pharmcat.jar 相同）
//     -o <output_dir> \
//     -bf <basename> \
//     -rs CPIC,DPWG \                 ← 只輸出 CPIC + DPWG recommendation
//     -reporterJson \                 ← 輸出 JSON（預設只有 HTML）
//     -reporterCallsOnlyTsv           ← 額外輸出 calls-only TSV

process PGX_PHARMCAT {

    label 'process_medium'

    container "${params.sif_dir}/pharmcat_3.2.0.sif"

    containerOptions "${params.apptainer_base_opts}"

    publishDir "${params.out_dir}/${sample_id}/07_pgx", mode: 'copy'

    input:
    tuple val(sample_id),
          path(snv_vcf), path(snv_tbi),
          val(stellarpgx_path),
          val(optitype_path)

    output:
    tuple val(sample_id),
          path("${sample_id}.pharmcat.report.json"),
          path("${sample_id}.outside_calls.tsv"),
          emit: pharmcat_ch

    script:
    """
    echo "[PGX_PHARMCAT] ${sample_id}：PharmCAT 3.2.0" >&2

    # ── Step 1：整合 outside calls ────────────────────────────
    # 格式（tab 分隔，# 開頭為 comment）：
    #   欄1: Gene   欄2: Diplotype   欄3: Phenotype   欄4: ActivityScore
    # CYP2D6: CYP2D6\t*1/*2\tNormal Metabolizer\t2.0
    # HLA-B:  HLA-B\t\t*57:01 positive     ← 兩個 tab（Diplotype 欄空白）
    # 參考：https://pharmcat.clinpgx.org/using/Outside-Call-Format/
    echo "[PGX_PHARMCAT] Step 1：整合 outside calls" >&2

    python3 ${params.scripts_dir}/build_outside_calls.py \
        --stellarpgx  ${stellarpgx_path} \
        --optitype    ${optitype_path} \
        --sample      ${sample_id} \
        --output      ${sample_id}.outside_calls.tsv

    echo "[PGX_PHARMCAT] outside_calls.tsv 內容：" >&2
    cat ${sample_id}.outside_calls.tsv >&2

    # ── Step 2：VCF Preprocessor ──────────────────────────────
    # NCKUH ensemble VCF 有雙 sample column（{sample_id}_DV, {sample_id}_HC）
    # PharmCAT preprocessor 的 -s 需要用 VCF 裡實際的 sample name
    # 取第一個 sample column（_DV）給 PharmCAT 用
    # PharmCAT 只需要其中一個 column 的 genotype
    echo "[PGX_PHARMCAT] Step 2：VCF preprocessor" >&2

    # 抓 VCF 裡的第一個 sample name（處理 NCKUH 雙 column 和 DRAGEN 單 column）
    VCF_SAMPLE=\$(bcftools query -l ${snv_vcf} | head -1)
    echo "[PGX_PHARMCAT] VCF sample column：\$VCF_SAMPLE" >&2

    mkdir -p preproc_out pharmcat_out

    /pharmcat/pharmcat_vcf_preprocessor \
        -vcf ${snv_vcf} \
        -s   \$VCF_SAMPLE \
        -o   preproc_out

    PREPROC_VCF=\$(find preproc_out -name "*.preprocessed.vcf.bgz" | head -1)
    if [ -z "\$PREPROC_VCF" ]; then
        echo "[PGX_PHARMCAT] 錯誤：找不到 preprocessor 輸出" >&2
        exit 1
    fi
    echo "[PGX_PHARMCAT] Preprocessed VCF：\$PREPROC_VCF（\$(wc -c < \$PREPROC_VCF) bytes）" >&2

    # ── Step 3：pharmcat.jar（matcher + phenotyper + reporter）─
    # -vcf   → preprocessed VCF（自動跑完 matcher + phenotyper + reporter）
    # -po    → outside calls TSV（CYP2D6 diplotype + HLA typing）
    # -s     → sample ID
    # -rs    → 只輸出 CPIC + DPWG recommendation
    # -reporterJson         → 輸出 report.json
    # -reporterCallsOnlyTsv → 輸出 calls_only.tsv（可直接讀，不需 parse JSON）
    echo "[PGX_PHARMCAT] Step 3：pharmcat.jar" >&2

    java -jar /pharmcat/pharmcat.jar \
        -vcf "\$PREPROC_VCF" \
        -po  ${sample_id}.outside_calls.tsv \
        -s   \$VCF_SAMPLE \
        -rs  CPIC,DPWG \
        -reporterJson \
        -reporterCallsOnlyTsv \
        -o   pharmcat_out \
        -bf  ${sample_id}

    # ── Step 4：整理輸出 ────────────────────────────────────
    # pharmcat.jar 輸出（-bf {sample_id}）：
    #   pharmcat_out/{sample_id}.report.json
    #   pharmcat_out/{sample_id}.report.html
    #   pharmcat_out/{sample_id}.report.calls_only.tsv
    REPORT_JSON=\$(find pharmcat_out -name "*.report.json" | head -1)
    if [ -n "\$REPORT_JSON" ]; then
        cp "\$REPORT_JSON" ${sample_id}.pharmcat.report.json
        echo "[PGX_PHARMCAT] report.json 大小：\$(wc -c < \$REPORT_JSON) bytes" >&2
    else
        echo "[PGX_PHARMCAT] 警告：找不到 report.json，產生空白輸出" >&2
        echo '{}' > ${sample_id}.pharmcat.report.json
    fi

    CALLS_TSV=\$(find pharmcat_out -name "*.calls_only.tsv" | head -1)
    if [ -n "\$CALLS_TSV" ]; then
        echo "[PGX_PHARMCAT] calls_only.tsv 預覽（前 5 行）：" >&2
        head -5 "\$CALLS_TSV" >&2
    fi

    echo "[PGX_PHARMCAT] ${sample_id} 完成" >&2
    """
}

// ──────────────────────────────────────────────────────────────
// Process 4：解析 PharmCAT JSON → 臨床用 TSV
// ──────────────────────────────────────────────────────────────

process PGX_PARSE {

    label 'process_low'

    container "${params.sif_dir}/tertiary_python_1.0.0.sif"

    containerOptions "${params.apptainer_base_opts}"

    publishDir "${params.out_dir}/${sample_id}/07_pgx", mode: 'copy'

    input:
    tuple val(sample_id), val(pipeline_type),
          path(pharmcat_json),
          path(outside_calls_tsv),
          path(mito_tsv, stageAs: "mito_input.tsv")

    output:
    tuple val(sample_id),
          path("${sample_id}.pgx.tsv"),
          emit: pgx_tsv_ch

    script:
    """
    echo "[PGX_PARSE] ${sample_id}（${pipeline_type}）：解析 PharmCAT 輸出" >&2

    python3 ${params.scripts_dir}/parse_pgx_report.py \
        --pharmcat_json   ${pharmcat_json} \
        --outside_calls   ${outside_calls_tsv} \
        --mito_tsv        mito_input.tsv \
        --sample          ${sample_id} \
        --pipeline        ${pipeline_type} \
        --output          ${sample_id}.pgx.tsv

    echo "[PGX_PARSE] ${sample_id} 完成" >&2
    echo "--- PGx TSV 總行數（含 header）---" >&2
    wc -l ${sample_id}.pgx.tsv >&2
    echo "--- 前 5 行（前 8 欄）---" >&2
    head -5 ${sample_id}.pgx.tsv | cut -f1-8 >&2
    """
}

// ──────────────────────────────────────────────────────────────
// 組合 workflow（供 main_tertiary.nf 呼叫）
// ──────────────────────────────────────────────────────────────

workflow PGX_ANNOTATE {

    take:
    pgx_wgs_vcf_ch  // tuple(sample_id, pipeline_type, snv_vcf, snv_tbi, bam, bai)  ← WGS 樣本
    pgx_wes_vcf_ch  // tuple(sample_id, pipeline_type, snv_vcf, snv_tbi)             ← WES 樣本
    mito_tsv_ch     // tuple(sample_id, mito_tsv)                                    ← 全樣本

    main:

    def no_file = file("NO_FILE")

    // ── WGS：先跑 StellarPGx，再進 PharmCAT ─────────────────
    // PGX_PHARMCAT input：tuple(sample_id, vcf, tbi, stellarpgx_tsv, optitype_tsv)
    // WGS 有 stellarpgx，WES 兩個都是 no_file
    if (params.run_pgx_cyp2d6) {
        bam_only_ch = pgx_wgs_vcf_ch
            .map { sid, ptype, vcf, tbi, bam, bai -> tuple(sid, bam, bai) }
        PGX_STELLARPGX(bam_only_ch)

        wgs_pharmcat_ch = pgx_wgs_vcf_ch
            .map { sid, ptype, vcf, tbi, bam, bai -> tuple(sid, vcf, tbi) }
            .join(PGX_STELLARPGX.out.stellarpgx_ch)
            .map { sid, vcf, tbi, spgx -> tuple(sid, vcf, tbi, spgx, no_file) }
    } else {
        wgs_pharmcat_ch = pgx_wgs_vcf_ch
            .map { sid, ptype, vcf, tbi, bam, bai ->
                tuple(sid, vcf, tbi, no_file, no_file)
            }
    }

    // ── WES：直接進 PharmCAT（no outside call）───────────────
    wes_pharmcat_ch = pgx_wes_vcf_ch
        .map { sid, ptype, vcf, tbi -> tuple(sid, vcf, tbi, no_file, no_file) }

    // ── 合併 WGS + WES → PGX_PHARMCAT ───────────────────────
    PGX_PHARMCAT(wgs_pharmcat_ch.mix(wes_pharmcat_ch))

    // ── PGX_PARSE：加上 pipeline_type + mito_tsv ────────────
    all_vcf_ch = pgx_wgs_vcf_ch
        .map { sid, ptype, vcf, tbi, bam, bai -> tuple(sid, ptype) }
        .mix(pgx_wes_vcf_ch.map { sid, ptype, vcf, tbi -> tuple(sid, ptype) })

    parse_input_ch = all_vcf_ch
        .join(PGX_PHARMCAT.out.pharmcat_ch)
        .join(mito_tsv_ch.ifEmpty(Channel.empty()), remainder: true)
        .map { vals ->
            def sid   = vals[0]
            def ptype = vals[1]
            def json  = vals[2]
            def ocall = vals[3]
            def mito  = (vals.size() > 4 && vals[4] != null) ? vals[4] : no_file
            tuple(sid, ptype, json, ocall, mito)
        }

    PGX_PARSE(parse_input_ch)

    emit:
    pgx_tsv_ch = PGX_PARSE.out.pgx_tsv_ch
}
