# WGS/WES Germline Tertiary Analysis Pipeline

A clinical-grade Nextflow DSL2 pipeline for tertiary analysis of germline variants from whole-genome and whole-exome sequencing, developed for the Department of Genomic Medicine and Neurology, National Cheng Kung University Hospital.

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A523.x-brightgreen)](https://www.nextflow.io/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE)
[![Container: Apptainer](https://img.shields.io/badge/container-Apptainer-blue)](https://apptainer.org/)

---

## Overview

This pipeline takes VCF output from secondary analysis pipelines (NCKUH or DRAGEN) and performs clinical annotation, classification, and pharmacogenomics analysis. It supports both WGS and WES modes, and is designed for clinical deployment in a regulated environment.

Two input modes are supported:

| Input type | Source | SNV/Indel input | BAM required |
|------------|--------|-----------------|--------------|
| `nckuh` | NCKUH secondary pipeline | `{sample}.ensemble.fixed.vcf.gz` (DeepVariant + HaplotypeCaller ensemble) | Optional (PGx WGS) |
| `dragen` | Illumina DRAGEN | `{sample}.hard-filtered.vcf.gz` | Optional (PGx WGS) |

---

## Pipeline Flowchart

```
VCF (nckuh / dragen)                    BAM (WGS only, optional)
        │                                         │
        ▼                                         │
┌───────────────────────────────────────────────────────────────┐
│  Step 0 · Prepare VCF                                         │
│  NCKUH: Add CALLERS tag, PASS filter, split chrM              │
│  DRAGEN: Add pipeline tag, auto-build tabix index, split chrM │
└───────────────────────┬───────────────────────────────────────┘
                        │
          ┌─────────────┴──────────────────────────────┐
          ▼                                             ▼
┌─────────────────────────────────┐     ┌──────────────────────────────────┐
│  SNV / Indel Track              │     │  Parallel Tracks                 │
│                                 │     │                                  │
│  VEP 115                        │     │  MITO: VEP (light) + gnomAD mito │
│  (dbNSFP, LOFTEE, ClinVar,      │     │         + ClinVar → mito.tsv     │
│   gnomAD, 1000G)                │     │                                  │
│       ↓                         │     │  STR:  GangSTR/ExpansionHunter   │
│  Pangolin (splice, GPU)         │     │        + STRchive → str.tsv      │
│       ↓                         │     │                                  │
│  PARSE_CSQ (61 columns)         │     │  CNV/SV: AnnotSV 3.5.10          │
│       ↓                         │     │   NCKUH WES  → cnv.annotated.tsv │
│  ACMG classifier                │     │   NCKUH WGS  → sv.annotated.tsv  │
│  (ClinGen SVI 2022)             │     │   DRAGEN CNV → cnv.annotated.tsv │
│       ↓                         │     │   DRAGEN SV  → sv.annotated.tsv  │
│  snv_indel.acmg.tsv (65 cols)   │     └──────────────────────────────────┘
└─────────────────────────────────┘

          ┌──────────────────────────────────────────────┐
          │  PGx Track (WGS + WES, BAM optional)         │
          │                                              │
          │  WGS:                                        │
          │    BAM → PGX_HLA_EXTRACT (samtools)          │
          │        → PGX_OPTITYPE (razers3 + OptiType)   │
          │          → HLA-A/B allele                    │
          │    BAM → PGX_STELLARPGX (StellarPGx)         │
          │          → CYP2D6 diplotype (outside call)   │
          │    VCF → pharmcat_vcf_preprocessor           │
          │        → pharmcat.jar -po outside_calls.tsv  │
          │          → report.json                       │
          │                                              │
          │  WES:                                        │
          │    VCF → pharmcat_vcf_preprocessor           │
          │        → pharmcat.jar → report.json          │
          │                                              │
          │  report.json + mito.tsv                      │
          │        → parse_pgx_report.py                 │
          │          → pgx.tsv (16 cols, CPIC Level A)   │
          └──────────────────────────────────────────────┘
```

---

## Hardware Requirements

| Environment | CPU | GPU | RAM | Role |
|-------------|-----|-----|-----|------|
| Local (dev) | R9 9950X 16c | RTX PRO 6000 96GB | 128GB | Development & testing |
| DGM Server | Xeon w7-3565X 32c | RTX 2000 Ada 16GB | 125GB | Clinical deployment |

> GPU is only required for Pangolin splice scoring (`use_gpu_pangolin = true`). All other steps are CPU-only.

---

## Installation

### Step 1 — Install dependencies

```bash
# Install Apptainer
sudo apt update && sudo apt install -y apptainer

# Install Miniforge + Nextflow
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
bash Miniforge3-Linux-x86_64.sh -b -p /opt/miniforge
source /opt/miniforge/bin/activate
mamba create -n nextflow openjdk=17 nextflow procps-ng -y
mamba create -n genome python=3.11 bcftools tabix -y
```

### Step 2 — Clone pipeline

```bash
git clone https://github.com/your-org/NGStertiary.git
cd NGStertiary
# Pipeline code lives at NGStertiary/1_0_0/
```

### Step 3 — Build containers

> ⚠️ **Critical:** Always add `--disable-cache` and `APPTAINER_SQUASH_OPTIONS="-processors 1"` to every `apptainer build/pull` command. mksquashfs 4.7.5 has a bug that causes segfault (exit 139) without these flags.

**tertiary_python_1.0.0.sif** — Python 3.11 + cyvcf2 + bcftools + tabix

```bash
cat > /tmp/tertiary_python.def << 'EOF'
Bootstrap: docker
From: python:3.11-slim

%post
    apt-get update -qq && apt-get install -y --no-install-recommends \
        gcc g++ make \
        zlib1g-dev libbz2-dev liblzma-dev \
        libcurl4-openssl-dev libssl-dev \
        procps bcftools tabix \
        && rm -rf /var/lib/apt/lists/*
    pip install --no-cache-dir cyvcf2 pandas numpy

%test
    python3 -c "import cyvcf2; print('cyvcf2:', cyvcf2.__version__)"
    bcftools --version | head -1
    bgzip --version | head -1

%labels
    Version 1.0.0
    Description "Tertiary pipeline Python tools: cyvcf2 + pandas + bcftools + bgzip + tabix"
EOF

APPTAINER_SQUASH_OPTIONS="-processors 1" \
apptainer build --disable-cache $SIF_DIR/tertiary_python_1.0.0.sif /tmp/tertiary_python.def
```

> ⚠️ `tabix` apt package includes bgzip and tabix. Installing only `bcftools` is not enough.

---

**vep_115.sif** — VEP 115 + LOFTEE GRCh38 + samtools + bcftools

```bash
cat > /tmp/vep_115.def << 'EOF'
Bootstrap: docker
From: ensemblorg/ensembl-vep:release_115.0

%post
    apt-get update -qq && apt-get install -y --no-install-recommends \
        git procps samtools bcftools \
        && rm -rf /var/lib/apt/lists/*

    mkdir -p /opt/vep/Plugins

    cd /opt/vep/src/ensembl-vep
    perl INSTALL.pl \
        --AUTO p \
        --PLUGINS dbNSFP,LoFtool \
        --PLUGINSDIR /opt/vep/Plugins \
        --NO_HTSLIB \
        --NO_UPDATE \
        2>&1 | tail -10

    cd /tmp
    git clone --depth 1 --branch grch38 \
        https://github.com/konradjk/loftee.git loftee_grch38
    cp -r /tmp/loftee_grch38/* /opt/vep/Plugins/
    rm -rf /tmp/loftee_grch38

    cpanm --quiet --notest List::MoreUtils

%test
    vep --help 2>&1 | grep "ensembl-vep"
    samtools --version | head -1
    bcftools --version | head -1
    perl -e "use DBD::SQLite; print 'DBD::SQLite OK\n'"

%labels
    Version 1.0.5
    VEP_release 115
    Description "VEP 115 + LOFTEE GRCh38 + samtools + bcftools + dbNSFP + LoFtool"
EOF

APPTAINER_SQUASH_OPTIONS="-processors 1" \
apptainer build --disable-cache $SIF_DIR/vep_115.sif /tmp/vep_115.def
```

---

**annotsv_3.5.10.sif** — AnnotSV 3.5.10

```bash
cat > /tmp/annotsv_3.5.10.def << 'EOF'
Bootstrap: docker
From: ubuntu:22.04

%post
    apt-get update && apt-get install -y \
        tcl tcllib tclx wget curl git bedtools vcftools \
        && apt-get clean

    cd /opt
    wget https://github.com/lgmgeo/AnnotSV/archive/refs/tags/v3.5.10.tar.gz
    tar xzf v3.5.10.tar.gz
    cd AnnotSV-3.5.10 && make PREFIX=/usr/local install
    rm /opt/v3.5.10.tar.gz

%test
    AnnotSV --help 2>&1 | head -3

%labels
    Version 1.0.0
    AnnotSV 3.5.10
EOF

APPTAINER_SQUASH_OPTIONS="-processors 1" \
apptainer build --disable-cache $SIF_DIR/annotsv_3.5.10.sif /tmp/annotsv_3.5.10.def
```

---

**optitype_1.3.5.sif** — OptiType 1.3.5 + razers3 3.5.12 + samtools 1.21 (custom build)

Why custom: `fred2/optitype:latest` uses Python 2.7/3.5 (f-strings unsupported, razers3 lacks fastq output); `biocontainers/optitype:1.5.0` has a Pyomo 6.10 constraint infeasible bug.

```bash
cat > /tmp/optitype_1.3.5.def << 'EOF'
Bootstrap: docker
From: ubuntu:22.04

%labels
    maintainer p88124019@gs.ncku.edu.tw
    version 1.3.5
    description OptiType 1.3.5 HLA typing with razers3 3.5.12 + samtools 1.21

%environment
    export PATH="/opt/conda/bin:$PATH"
    export DEBIAN_FRONTEND=noninteractive

%post
    export DEBIAN_FRONTEND=noninteractive
    export PATH="/opt/conda/bin:$PATH"

    apt-get update && apt-get install -y wget bzip2 ca-certificates \
        && apt-get clean && rm -rf /var/lib/apt/lists/*

    wget -q https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh \
        -O /tmp/miniforge.sh
    bash /tmp/miniforge.sh -b -p /opt/conda
    rm /tmp/miniforge.sh

    /opt/conda/bin/conda config --add channels defaults
    /opt/conda/bin/conda config --add channels bioconda
    /opt/conda/bin/conda config --add channels conda-forge
    /opt/conda/bin/conda config --set channel_priority strict

    /opt/conda/bin/conda install -y \
        python=3.11 optitype=1.3.5 samtools=1.21 coincbc \
        && /opt/conda/bin/conda clean -afy

%test
    OptiTypePipeline.py --help 2>&1 | head -3
    razers3 --version 2>&1 | head -2
    samtools --version | head -1

%runscript
    exec OptiTypePipeline.py "$@"
EOF

APPTAINER_SQUASH_OPTIONS="-processors 1" \
apptainer build --disable-cache $SIF_DIR/optitype_1.3.5.sif /tmp/optitype_1.3.5.def
```

### Step 4 — Clone StellarPGx repo

StellarPGx requires the full git repo (database, resources, scripts) in addition to the container:

```bash
cd /path/to/NGStertiary/1_0_0/
git clone https://github.com/SBIMB/StellarPGx stellarpgx_repo
```

Expected structure:
```
stellarpgx_repo/
├── database/cyp2d6/hg38/       ← CYP2D6 allele database
├── resources/cyp2d6/res_hg38/  ← HLA + annotation resources
└── scripts/cyp2d6/hg38/bin/    ← calling scripts
```

### Step 5 — Download reference data

```bash
REF_DIR="/path/to/reference/hg38"
TERTIARY_DIR="${REF_DIR}/tertiary"
mkdir -p ${TERTIARY_DIR}

# ── VEP cache ──────────────────────────────────────────────────
mkdir -p ${TERTIARY_DIR}/vep_cache
cd ${TERTIARY_DIR}/vep_cache
curl -O https://ftp.ensembl.org/pub/release-115/variation/indexed_vep_cache/homo_sapiens_vep_115_GRCh38.tar.gz
tar -xzf homo_sapiens_vep_115_GRCh38.tar.gz

# ── dbNSFP 4.9c ────────────────────────────────────────────────
mkdir -p ${TERTIARY_DIR}/dbnsfp
# Download from https://sites.google.com/site/jpopgen/dbNSFP (registration required)
# File: dbNSFP4.9c_with_pknn_grch38.gz + .tbi

# ── LOFTEE ─────────────────────────────────────────────────────
mkdir -p ${TERTIARY_DIR}/loftee
# Follow: https://github.com/konradjk/loftee
# Required: loftee plugin + gerp_conservation_scores.homo_sapiens.GRCh38.bw

# ── ClinVar ────────────────────────────────────────────────────
mkdir -p ${TERTIARY_DIR}/clinvar
cd ${TERTIARY_DIR}/clinvar
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi
# Rename chromosomes to chr-prefixed:
bcftools annotate --rename-chrs chr_name_conv.txt clinvar.vcf.gz \
    -O z -o clinvar_20260510.vcf.gz
tabix -p vcf clinvar_20260510.vcf.gz
# Build lookup TSV (for fast ClinVar lookup in annotation step):
python3 scripts/build_clinvar_lookup.py \
    --input clinvar_20260510.vcf.gz \
    --output clinvar_lookup.tsv.gz

# ── gnomAD ─────────────────────────────────────────────────────
mkdir -p ${TERTIARY_DIR}/gnomad
# gnomAD v4 genome + exome (large files — download only needed chromosomes)
# https://gnomad.broadinstitute.org/downloads

# ── Pangolin gene annotation DB ────────────────────────────────
mkdir -p ${TERTIARY_DIR}/pangolin
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gtf.gz
# Convert to Pangolin .db format (see Pangolin docs)

# ── STRchive ───────────────────────────────────────────────────
mkdir -p ${TERTIARY_DIR}/strchive
wget https://raw.githubusercontent.com/dashnowlab/STRchive/main/data/STRchive-loci.json
# Build lookup tables:
python3 scripts/build_str_lookup.py \
    --input STRchive-loci.json \
    --out_varid str_lookup_varid.tsv.gz \
    --out_pos str_lookup_pos.tsv.gz

# ── gnomAD mito v3.1（CC0）─────────────────────────────────────
mkdir -p ${TERTIARY_DIR}/gnomad_mito
# Download from gnomAD browser: gnomAD v3.1 mitochondrial variants
python3 scripts/build_gnomad_mito_lookup.py \
    --input gnomad.genomes.v3.1.sites.chrM.vcf.bgz \
    --output gnomad_mito_lookup.tsv.gz

# ── AnnotSV annotation databases ───────────────────────────────
mkdir -p ${TERTIARY_DIR}/annotsv_annotations
# Follow AnnotSV documentation to install annotation databases
# https://lbgi.fr/AnnotSV/Documentation

# ── ClinGen ────────────────────────────────────────────────────
mkdir -p ${TERTIARY_DIR}/clingen
wget https://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh38.tsv
wget "https://search.clinicalgenome.org/kb/gene-validity/download" \
    -O clingen_gene_disease_validity.csv
python3 scripts/build_gene_moi.py \
    --clingen_gene clingen_gene_disease_validity.csv \
    --output gene_moi.tsv.gz
```

Expected directory structure:
```
${REF_DIR}/tertiary/
├── vep_cache/homo_sapiens/115_GRCh38/
├── dbnsfp/dbNSFP4.9c_with_pknn_grch38.gz(.tbi)
├── loftee/
├── clinvar/clinvar_20260510.vcf.gz(.tbi) + clinvar_lookup.tsv.gz
├── gnomad/
├── pangolin/gencode.v47.annotation.db
├── strchive/STRchive-loci.json + str_lookup_varid.tsv.gz + str_lookup_pos.tsv.gz
├── gnomad_mito/gnomad_mito_lookup.tsv.gz
├── annotsv_annotations/share/AnnotSV/
└── clingen/ClinGen_gene_curation_list_GRCh38.tsv + gene_moi.tsv.gz
```

### Step 6 — Configure pipeline

Edit `nextflow_tertiary.config` and set paths for your environment under the appropriate profile (`local` or `dgm`):

```groovy
local {
    params {
        ref_dir     = "/path/to/reference/hg38"
        sif_dir     = "/path/to/containers"
        scripts_dir = "/path/to/NGStertiary/1_0_0/scripts"
        out_dir     = "/path/to/output"
    }
}
```

All other paths are derived automatically from `ref_dir`.

---

## Quick Start

### 1. Load environment

```bash
conda activate nextflow
cd /path/to/NGStertiary/1_0_0/
```

### 2. Prepare samplesheet

```csv
sample_id,pipeline_type,input_dir,seq_type,hpo
NA12878_WES,nckuh,/path/to/secondary_output/NA12878_WES/NA12878_WES,WES,
NA12878_WGS,nckuh,/path/to/secondary_output/NA12878_WGS/NA12878_WGS,WGS,HP:0001250
VAL-10,dragen,/path/to/dragen_output/20260428_run,WGS,
```

| Field | Required | Description |
|-------|----------|-------------|
| `sample_id` | ✅ | Unique sample identifier |
| `pipeline_type` | ✅ | `nckuh` or `dragen` |
| `input_dir` | ✅ | Secondary analysis output directory |
| `seq_type` | ✅ | `WES` or `WGS` |
| `hpo` | ❌ | HPO terms, `\|`-separated (reserved for future phenotype matching) |

**Input path rules (auto-resolved from `input_dir` + `sample_id`):**
```
nckuh: {input_dir}/04_snv_indel/{sample_id}.ensemble.fixed.vcf.gz
       {input_dir}/02_alignment/{sample_id}.aligned.sorted.bam  (WGS PGx)
dragen: {input_dir}/vcf.gz/{sample_id}.hard-filtered.vcf.gz
        {input_dir}/bam/{sample_id}.bam  (WGS PGx)
```

### 3. Run

```bash
# NCKUH samples (WES + WGS mixed, PGx enabled by default for WGS)
nextflow -c nextflow_tertiary.config run main_tertiary.nf \
    -profile local \
    --pipeline_type nckuh \
    --samplesheet /path/to/samplesheet.csv \
    --out_dir /path/to/output \
    -resume

# DRAGEN samples
nextflow -c nextflow_tertiary.config run main_tertiary.nf \
    -profile dgm \
    --pipeline_type dragen \
    --samplesheet /path/to/samplesheet.csv \
    --out_dir /path/to/output \
    -resume

# Disable PGx (faster, for SNV/CNV annotation only)
nextflow -c nextflow_tertiary.config run main_tertiary.nf \
    -profile local \
    --pipeline_type nckuh \
    --samplesheet /path/to/samplesheet.csv \
    --out_dir /path/to/output \
    --run_pgx false \
    -resume
```

### Key parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--pipeline_type` | `null` | Filter samplesheet by pipeline type (`nckuh` / `dragen`) |
| `--run_pgx` | `true` | Enable PGx module (PharmCAT) |
| `--run_pgx_cyp2d6` | `true` | Enable StellarPGx CYP2D6 caller (requires BAM, WGS only) |
| `--run_pgx_hla` | `true` | Enable OptiType HLA-A/B typing (requires BAM, WGS only) |

---

## Output Structure

```
{out_dir}/{SAMPLE_ID}/
├── 00_prepare/          Preprocessed VCF (intermediate)
├── 01_vep/              VEP annotation output (intermediate)
├── 02_pangolin/         Pangolin splice scores (intermediate)
├── 03_acmg/
│   └── {SAMPLE_ID}.snv_indel.acmg.tsv     ★ SNV/Indel (65 columns)
├── 04_mito/
│   └── {SAMPLE_ID}.mito.tsv               ★ mtDNA variants (21 columns)
├── 05_str/
│   └── {SAMPLE_ID}.str.tsv                ★ STR (22 columns)
├── 06_cnv_sv/
│   ├── {SAMPLE_ID}.cnv.annotated.tsv      ★ CNV (AnnotSV)
│   └── {SAMPLE_ID}.sv.annotated.tsv       ★ SV (AnnotSV)
├── 07_pgx/
│   ├── {SAMPLE_ID}.pgx.tsv                ★ PGx report (16 columns, CPIC Level A)
│   ├── {SAMPLE_ID}.pharmcat.report.json   PharmCAT full report (archive)
│   ├── {SAMPLE_ID}.outside_calls.tsv      Outside calls (archive)
│   ├── {SAMPLE_ID}.stellarpgx.tsv         CYP2D6 diplotype (WGS only)
│   └── {SAMPLE_ID}.optitype.tsv           HLA-A/B alleles (WGS only)
└── pipeline_info/       Execution reports (HTML + trace)
```

### pgx.tsv columns (16)

`SAMPLE_ID · PIPELINE · GENE · DIPLOTYPE · ACTIVITY_SCORE · PHENOTYPE · DRUG · GUIDELINE_SOURCE · RECOMMENDATION · IMPLICATION · CPIC_LEVEL · DPWG_LEVEL · OUTSIDE_CALLER · MTRN1_RISK · NOTES · EVIDENCE_STRENGTH`

Covers CPIC Level A genes: CYP2D6, CYP2C19, CYP2C9, DPYD, TPMT, NUDT15, SLCO1B1, HLA-A, HLA-B, UGT1A1, G6PD, MT-RNR1 (via mito pipeline), IFNL3, CACNA1S, RYR1.

---

## Profiles

| Profile | Target | Notes |
|---------|--------|-------|
| `local` | Development machine (16c) | `process_high` = 16 CPUs |
| `dgm` | DGM Server (32c) | `process_high` = 32 CPUs |

---

## Validation with NA12878

```bash
# Download WES test data
wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R1_001.fastq.gz
wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R2_001.fastq.gz

# Download WGS test data
wget https://ftp.ebi.ac.uk/vol1/fastq/ERR194/ERR194147/ERR194147_1.fastq.gz
wget https://ftp.ebi.ac.uk/vol1/fastq/ERR194/ERR194147/ERR194147_2.fastq.gz
```

Expected results after running secondary → tertiary pipeline:

| Sample | SNV/Indel | P/LP | CYP2D6 | HLA-B |
|--------|-----------|------|--------|-------|
| NA12878 WES | ~37,199 variants | ~41 LP (low VAF artifact) | VCF-based only | N/A (WES) |
| NA12878 WGS | ~5,729,808 variants | — | `*1/*5`, activity=1.0 (StellarPGx) | `*08:01/*08:01`⚠️, *57:01 negative ✅ |


> ⚠️ NA12878 HLA-B ground truth is `*07:02/*40:02`（heterozygous），the current pipeline call  `*08:01/*08:01`, which is close but still incorrect.

---

## License and Third-party Tools

This pipeline is released under the [GNU General Public License v3](LICENSE) (GPL v3).

> ⚠️ **Clinical use warning:** This pipeline is designed for clinical use. All tools included are compatible with commercial/clinical use. Do **not** substitute with BCyrius (PolyForm Strict), Aldy (non-commercial), or MITOMAP (CC BY-NC).

| Tool | Version | License | Notes |
|------|---------|---------|-------|
| [Nextflow](https://github.com/nextflow-io/nextflow) | ≥ 23.x | Apache 2.0 | |
| [Apptainer](https://github.com/apptainer/apptainer) | ≥ 1.x | BSD 3-Clause | |
| [VEP](https://github.com/Ensembl/ensembl-vep) | 115 | Apache 2.0 | |
| [Pangolin](https://github.com/tkzeng/Pangolin) | custom | BSD | |
| [AnnotSV](https://github.com/lgmgeo/AnnotSV) | 3.5.10 | GNU GPL v3 | |
| [PharmCAT](https://github.com/PharmGKB/PharmCAT) | 3.2.0 | MPL 2.0 | |
| [StellarPGx](https://github.com/SBIMB/StellarPGx) | 1.2.8 | Open source | Graphtyper 2.5.1 |
| [OptiType](https://github.com/FRED-2/OptiType) | 1.3.5 | BSD 3-Clause | 自建 sif（Ubuntu 22.04 + Miniforge）|
| [SAMtools](https://github.com/samtools/samtools) | 1.23.1 | MIT | |
| [BCFtools](https://github.com/samtools/bcftools) | 1.23.1 | MIT | |
| [gnomAD v3.1 mito](https://gnomad.broadinstitute.org) | 3.1 | CC0 | Replaces MITOMAP (CC BY-NC) |
| [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) | — | Public domain | |
| [STRchive](https://strchive.org) | — | CC BY 4.0 | |
| [PharmVar / CPIC](https://www.cpicpgx.org) | — | CC0 | |
| [ClinGen](https://clinicalgenome.org) | — | CC0 | |
| Aldy | — | Non-commercial only | ❌ **Not used** |
| BCyrius | — | PolyForm Strict | ❌ **Not used** |
| MITOMAP | — | CC BY-NC | ❌ **Replaced by gnomAD mito** |

Users are responsible for compliance with each tool's license terms.

---

## Citation

If you use this pipeline in your research, please cite the relevant tools listed above.

**Key references:**
- PharmCAT: Sangkuhl et al., *Clinical Pharmacology & Therapeutics* (2020)
- StellarPGx: Twesigomwe et al., *npj Genomic Medicine* (2021)
- OptiType: Szolek et al., *Bioinformatics* (2014)
- AnnotSV: Geoffroy et al., *Bioinformatics* (2018)
- VEP: McLaren et al., *Genome Biology* (2016)
