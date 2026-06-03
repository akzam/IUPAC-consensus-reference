# IUPAC Consensus Reference — Command Line Workflow

This document describes the complete command-line workflow used in:

> **IUPAC consensus references improve variant detection in clinically challenging genomic regions**  
> Saidin A, Ricos MG, Dibbens LM. *Cell Genomics*. 2026.

---

## Table of Contents

1. [Overview](#overview)
2. [Prerequisites](#prerequisites)
3. [Software Versions](#software-versions)
4. [Step 1: Prepare Variant Databases](#step-1-prepare-variant-databases)
5. [Step 2: Build IUPAC Consensus References](#step-2-build-iupac-consensus-references)
6. [Step 3: Align Reads](#step-3-align-reads)
7. [Step 4: Variant Calling](#step-4-variant-calling)
8. [Step 5: Benchmarking with hap.py](#step-5-benchmarking-with-happy)
9. [Complete Workflow Script](#complete-workflow-script)
10. [Output File Tree](#output-file-tree)
11. [Troubleshooting](#troubleshooting)

---

## Overview

The workflow consists of five stages:

```
Variant Database → IUPAC Consensus Reference → Alignment → Variant Calling → Benchmarking
     (VCF)              (FASTA)                (BAM)        (VCF)           (hap.py CSV)
```

Three experimental arms are compared:

| Arm | Aligner | Reference | Purpose |
|-----|---------|-----------|---------|
| **Baseline** | BWA-MEM | GRCh38.p13 | Standard clinical pipeline |
| **Base** | novoAlign | GRCh38.p13 | Isolate aligner effect |
| **Consensus** | novoAlign | IUPAC consensus | Full IUPAC effect |

---

## Prerequisites

### Required software

| Tool | Purpose | License |
|------|---------|---------|
| bcftools | VCF filtering | MIT/Free |
| novoUtil / novoAlign | IUPAC consensus generation & ambiguity-aware alignment | Commercial (free academic license) |
| BWA-MEM | Baseline alignment | MIT |
| samtools | BAM manipulation | MIT/Free |
| FreeBayes | Variant calling | MIT |
| GATK | Variant calling | BSD-3 |
| hap.py | Benchmarking | MIT |

### Input data

- **Reference genome:** GRCh38.p13 (`GRCh38.p13.fa`)
- **Population variant databases:**
  - 1000 Genomes (release 20190312)
  - dbSNP build 151
  - gnomAD v3 (`gnomADr3.vcf.gz`)
- **Sequencing data:** PCR-free 30x WGS (NovaSeq 6000), paired-end 150 bp
- **Truth sets:** GIAB v4.2.1 (HG001, HG002, HG005)
- **Stratifications:** GIAB GRCh38 stratification BEDs (v3.0)

### Environment setup

```bash
# Recommended: Conda environment
conda create -n iupac python=3.8 bcftools=1.9 samtools=1.17
conda activate iupac

# Install GATK (via Conda or manual download)
conda install -c bioconda gatk4=4.0

# Install FreeBayes
conda install -c bioconda freebayes=1.3.1

# NovoAlign / novoUtil: request academic license from Novocraft Technologies
# Download and add to PATH
export PATH="/path/to/novocraft/bin:$PATH"

# Verify installations
bcftools --version    # 1.9
novoalign --version   # 4.04.02
bwa                    # 0.7.18
samtools --version     # 1.17
freebayes --version    # 1.3.1
gatk --version         # 4.0
```

---

## Software Versions

All commands were run with the following versions. Use of Conda or Docker is strongly recommended to replicate the environment.

| Tool | Version | Notes |
|------|---------|-------|
| bcftools | 1.9 | VCF filtering and SNP extraction |
| novoUtil (Novocraft) | 4.04.02 | IUPAC code substitution |
| novoAlign | 4.04.02 | Ambiguity-aware alignment |
| BWA-MEM | 0.7.18 | Baseline aligner |
| samtools | 1.17 | BAM sorting, indexing |
| FreeBayes | 1.3.1 | Haplotype-based variant caller |
| GATK | 4.0 | HaplotypeCaller |
| hap.py | 0.3.15 | vcfeval benchmarking engine |

---

## Step 1: Prepare Variant Databases

Extract single-nucleotide variants (SNVs) from population databases and filter by allele frequency. This step is repeated for each database and AF threshold.

### 1.1 gnomAD v3

```bash
# Input: gnomAD v3 VCF (downloaded from https://gnomad.broadinstitute.org/)
GNOMAD_VCF="gnomADr3.vcf.gz"
REF="GRCh38.p13.fa"

# Extract SNVs only (exclude indels)
bcftools view --exclude-types indels -Oz -o gnomADr3_snp.vcf.gz ${GNOMAD_VCF}
bcftools index gnomADr3_snp.vcf.gz

# Filter by allele frequency thresholds
echo "Filtering gnomAD by AF thresholds..."

# AF ≥ 10%
bcftools filter -i 'INFO/AF[0] >= 0.10' -Oz -o gnomADr3_SNP_AF10pc.vcf.gz gnomADr3_snp.vcf.gz
bcftools index gnomADr3_SNP_AF10pc.vcf.gz

# AF ≥ 30%
bcftools filter -i 'INFO/AF[0] >= 0.30' -Oz -o gnomADr3_SNP_AF30pc.vcf.gz gnomADr3_snp.vcf.gz
bcftools index gnomADr3_SNP_AF30pc.vcf.gz

# Population-specific filters (gnomAD subpopulations)
for POP in AFR ASJ EAS NFE; do
    bcftools filter -i "INFO/AF_${POP}[0] >= 0.10" -Oz -o gnomADr3_${POP}_SNP_AF10pc.vcf.gz gnomADr3_snp.vcf.gz
    bcftools index gnomADr3_${POP}_SNP_AF10pc.vcf.gz

    bcftools filter -i "INFO/AF_${POP}[0] >= 0.30" -Oz -o gnomADr3_${POP}_SNP_AF30pc.vcf.gz gnomADr3_snp.vcf.gz
    bcftools index gnomADr3_${POP}_SNP_AF30pc.vcf.gz
done
```

### 1.2 1000 Genomes

```bash
# Input: 1000 Genomes phase 3 VCF
KG_VCF="1000G_phase3.vcf.gz"

# Extract SNVs and filter by AF
bcftools view --exclude-types indels -Oz -o 1KG_snp.vcf.gz ${KG_VCF}
bcftools index 1KG_snp.vcf.gz

# AF ≥ 10%
bcftools filter -i 'INFO/AF >= 0.10' -Oz -o 1KG_SNP_AF10pc.vcf.gz 1KG_snp.vcf.gz
bcftools index 1KG_SNP_AF10pc.vcf.gz

# AF ≥ 30%
bcftools filter -i 'INFO/AF >= 0.30' -Oz -o 1KG_SNP_AF30pc.vcf.gz 1KG_snp.vcf.gz
bcftools index 1KG_SNP_AF30pc.vcf.gz
```

### 1.3 dbSNP

```bash
# Input: dbSNP build 151 VCF
DBSNP_VCF="dbSNP151.vcf.gz"

# Extract SNVs
bcftools view --exclude-types indels -Oz -o dbSNP151_snp.vcf.gz ${DBSNP_VCF}
bcftools index dbSNP151_snp.vcf.gz

# CAF (Common Allele Frequency) ≥ 70%
bcftools filter -i 'INFO/CAF[0] >= 0.70' -Oz -o dbSNP151_SNP_CAF70pc.vcf.gz dbSNP151_snp.vcf.gz
bcftools index dbSNP151_SNP_CAF70pc.vcf.gz

# CAF ≥ 90%
bcftools filter -i 'INFO/CAF[0] >= 0.90' -Oz -o dbSNP151_SNP_CAF90pc.vcf.gz dbSNP151_snp.vcf.gz
bcftools index dbSNP151_SNP_CAF90pc.vcf.gz
```

### 1.4 Summary of variant database files

| Database | AF Threshold | Output VCF | Description |
|----------|-------------|------------|-------------|
| gnomAD v3 | AF ≥ 10% | `gnomADr3_SNP_AF10pc.vcf.gz` | Pan-human, common variants |
| gnomAD v3 | AF ≥ 30% | `gnomADr3_SNP_AF30pc.vcf.gz` | Pan-human, high-frequency variants |
| gnomAD v3 | AFR AF ≥ 10% | `gnomADr3_AFR_SNP_AF10pc.vcf.gz` | African population-specific |
| gnomAD v3 | AFR AF ≥ 30% | `gnomADr3_AFR_SNP_AF30pc.vcf.gz` | African population-specific |
| gnomAD v3 | ASJ AF ≥ 10% | `gnomADr3_ASJ_SNP_AF10pc.vcf.gz` | Ashkenazi Jewish population-specific |
| gnomAD v3 | ASJ AF ≥ 30% | `gnomADr3_ASJ_SNP_AF30pc.vcf.gz` | Ashkenazi Jewish population-specific |
| gnomAD v3 | EAS AF ≥ 10% | `gnomADr3_EAS_SNP_AF10pc.vcf.gz` | East Asian population-specific |
| gnomAD v3 | EAS AF ≥ 30% | `gnomADr3_EAS_SNP_AF30pc.vcf.gz` | East Asian population-specific |
| gnomAD v3 | NFE AF ≥ 10% | `gnomADr3_NFE_SNP_AF10pc.vcf.gz` | Non-Finnish European population-specific |
| gnomAD v3 | NFE AF ≥ 30% | `gnomADr3_NFE_SNP_AF30pc.vcf.gz` | Non-Finnish European population-specific |
| 1000 Genomes | AF ≥ 10% | `1KG_SNP_AF10pc.vcf.gz` | 1000 Genomes common variants |
| 1000 Genomes | AF ≥ 30% | `1KG_SNP_AF30pc.vcf.gz` | 1000 Genomes high-frequency variants |
| dbSNP 151 | CAF ≥ 70% | `dbSNP151_SNP_CAF70pc.vcf.gz` | dbSNP common allele frequency |
| dbSNP 151 | CAF ≥ 90% | `dbSNP151_SNP_CAF90pc.vcf.gz` | dbSNP very common allele frequency |

---

## Step 2: Build IUPAC Consensus References

Substitute GRCh38 reference alleles with IUPAC degenerate codes at common variant sites using novoUtil. Multi-allelic SNPs are represented using IUPAC codes to preserve allelic diversity while maintaining linear coordinates.

### 2.1 Generate consensus FASTA files

```bash
REF="GRCh38.p13.fa"

# gnomAD pan-human consensuses
echo "Building gnomAD consensus references..."
novoutil iupac gnomADr3_SNP_AF10pc.vcf.gz ${REF} > GRCh38_gnomAD_AF10.fa
novoutil iupac gnomADr3_SNP_AF30pc.vcf.gz ${REF} > GRCh38_gnomAD_AF30.fa

# gnomAD population-specific consensuses
for POP in AFR ASJ EAS NFE; do
    novoutil iupac gnomADr3_${POP}_SNP_AF10pc.vcf.gz ${REF} > GRCh38_gnomAD_${POP}_AF10.fa
    novoutil iupac gnomADr3_${POP}_SNP_AF30pc.vcf.gz ${REF} > GRCh38_gnomAD_${POP}_AF30.fa
done

# 1000 Genomes consensuses
novoutil iupac 1KG_SNP_AF10pc.vcf.gz ${REF} > GRCh38_1KG_AF10.fa
novoutil iupac 1KG_SNP_AF30pc.vcf.gz ${REF} > GRCh38_1KG_AF30.fa

# dbSNP consensuses
novoutil iupac dbSNP151_SNP_CAF70pc.vcf.gz ${REF} > GRCh38_dbSNP_CAF70.fa
novoutil iupac dbSNP151_SNP_CAF90pc.vcf.gz ${REF} > GRCh38_dbSNP_CAF90.fa
```

### 2.2 Index consensus references for alignment

```bash
# Build novoAlign indices for all consensus references
echo "Indexing consensus references..."

for FA in GRCh38_gnomAD_*.fa GRCh38_1KG_*.fa GRCh38_dbSNP_*.fa; do
    echo "Indexing ${FA}..."
    novoindex ${FA%.fa}.ndx ${FA}
done

# Also index baseline GRCh38 for novoAlign base condition
novoindex GRCh38.p13.ndx GRCh38.p13.fa
```

### 2.3 Verify consensus reference integrity

```bash
# Check that consensus FASTA has same sequence lengths as original
# (IUPAC substitution preserves coordinate structure)
echo "Verifying consensus reference lengths..."
python3 << 'EOF'
from Bio import SeqIO
import sys

ref = SeqIO.to_dict(SeqIO.parse("GRCh38.p13.fa", "fasta"))
cons = SeqIO.to_dict(SeqIO.parse("GRCh38_gnomAD_AF10.fa", "fasta"))

mismatches = 0
for chrom in ref:
    if len(ref[chrom]) != len(cons[chrom]):
        print(f"LENGTH MISMATCH: {chrom}: {len(ref[chrom])} vs {len(cons[chrom])}", file=sys.stderr)
        mismatches += 1

if mismatches == 0:
    print("All chromosomes: lengths match ✓")
else:
    print(f"{mismatches} chromosomes with length mismatches ✗", file=sys.stderr)
    sys.exit(1)
EOF
```

> **Note:** Consensus references preserve GRCh38 coordinate structure (length and position), differing only in base substitution via IUPAC codes. Variant calling and benchmarking are therefore performed directly in the native GRCh38 coordinate system without liftover.

---

## Step 3: Align Reads

Align 30x WGS reads from GIAB samples (HG001, HG002, HG005) to each reference.

### 3.1 Baseline: BWA-MEM / GRCh38

```bash
SAMPLE="HG001"
R1="${SAMPLE}_R1.fastq.gz"
R2="${SAMPLE}_R2.fastq.gz"
REF="GRCh38.p13.fa"

# Align with BWA-MEM
bwa mem -t 32 ${REF} ${R1} ${R2} |     samtools view -bS - > ${SAMPLE}_bwa_GRCh38.bam

# Sort and index
samtools sort -@ 32 ${SAMPLE}_bwa_GRCh38.bam -o ${SAMPLE}_bwa_GRCh38.sort.bam
samtools index ${SAMPLE}_bwa_GRCh38.sort.bam

# Remove unsorted BAM to save space
rm ${SAMPLE}_bwa_GRCh38.bam
```

### 3.2 Base condition: novoAlign / GRCh38

```bash
# Align unmodified GRCh38 with novoAlign (to isolate aligner effect)
novoalign -a --tune NOVASEQ -r Random -o BAM -c 32     -d GRCh38.p13.ndx     -f ${R1} ${R2}     > ${SAMPLE}_novo_GRCh38.bam

novosort --pcrFree --md --kt -i -o ${SAMPLE}_novo_GRCh38.sort.bam ${SAMPLE}_novo_GRCh38.bam
samtools index ${SAMPLE}_novo_GRCh38.sort.bam
rm ${SAMPLE}_novo_GRCh38.bam
```

### 3.3 Consensus conditions: novoAlign / IUPAC

```bash
# Example: gnomAD AF10% consensus
novoalign -a --tune NOVASEQ -r Random -o BAM -c 32     -d GRCh38_gnomAD_AF10.ndx     -f ${R1} ${R2}     > ${SAMPLE}_novo_gnomAD_AF10.bam

novosort --pcrFree --md --kt -i -o ${SAMPLE}_novo_gnomAD_AF10.sort.bam ${SAMPLE}_novo_gnomAD_AF10.bam
samtools index ${SAMPLE}_novo_gnomAD_AF10.sort.bam
rm ${SAMPLE}_novo_gnomAD_AF10.bam

# Repeat for all consensus references:
# GRCh38_gnomAD_AF30, GRCh38_1KG_AF10, GRCh38_1KG_AF30,
# GRCh38_dbSNP_CAF70, GRCh38_dbSNP_CAF90,
# GRCh38_gnomAD_AFR_AF10, GRCh38_gnomAD_AFR_AF30, etc.
```

### 3.4 Alignment parameters explained

| Parameter | Value | Purpose |
|-----------|-------|---------|
| `-a` | — | Output all alignments (for ambiguous reads) |
| `--tune NOVASEQ` | — | Optimise for NovaSeq chemistry |
| `-r Random` | — | Randomly assign ambiguously mapped reads |
| `-o BAM` | — | Output BAM format |
| `-c 32` | 32 threads | Parallel alignment |
| `--pcrFree` | — | No duplicate marking (PCR-free library) |
| `--md` | — | Calculate MD tags |
| `--kt` | — | Keep temporary files (for debugging) |

> **Important:** BWA-MEM does not support IUPAC ambiguous bases in the reference, preventing use of IUPAC consensus with this aligner. This is a current limitation of the method. The IUPAC effect is therefore inseparable from the aligner change in this design.

---

## Step 4: Variant Calling

Call SNVs and INDELs against GRCh38.p13 using three independent callers. To focus evaluation on detection performance rather than post-calling filtration, all tools report variants with QUAL ≥ 1, without caller-specific filtering.

### 4.1 BCFtools

```bash
SAMPLE="HG001"
BAM="${SAMPLE}_bwa_GRCh38.sort.bam"  # or consensus BAM
REF="GRCh38.p13.fa"

bcftools mpileup -Ou -f ${REF} ${BAM} |     bcftools call -f GQ -mv -Oz -o ${SAMPLE}_bcftools.vcf.gz

bcftools index ${SAMPLE}_bcftools.vcf.gz
```

### 4.2 FreeBayes

```bash
freebayes -f ${REF} -F 0.01 -C 5 ${BAM} > ${SAMPLE}_freebayes.vcf

# Compress and index
bgzip ${SAMPLE}_freebayes.vcf
bcftools index ${SAMPLE}_freebayes.vcf.gz
```

### 4.3 GATK HaplotypeCaller

```bash
gatk HaplotypeCaller     -R ${REF}     -I ${BAM}     -O ${SAMPLE}_gatk.vcf.gz
```

### 4.4 Repeat for all conditions

For each sample, repeat variant calling for:
- Baseline: `${SAMPLE}_bwa_GRCh38.sort.bam`
- Base: `${SAMPLE}_novo_GRCh38.sort.bam`
- Consensus: `${SAMPLE}_novo_gnomAD_AF10.sort.bam`, `${SAMPLE}_novo_gnomAD_AF30.sort.bam`, etc.

> **Note:** We used QUAL ≥ 1 without post-call filtering. Many clinical pipelines apply calibrated filters or machine-learning models to improve precision, but such steps are orthogonal to the alignment strategy and were not applied here. The precision losses reported are therefore likely conservative upper bounds.

---

## Step 5: Benchmarking with hap.py

Benchmark variant call accuracy against GIAB truth sets using hap.py (vcfeval engine).

### 5.1 Prepare inputs

```bash
SAMPLE="HG001"
CALLER="bcftools"
CONDITION="gnomAD_AF10"  # or "bwa_GRCh38", "novo_GRCh38", etc.

QUERY_VCF="${SAMPLE}_${CALLER}_${CONDITION}.vcf.gz"
TRUTH_VCF="${SAMPLE}_GIAB_v4.2.1_truth.vcf.gz"
CONFIDENT_BED="HG001_GRCh38_confident_regions.bed"
REF="GRCh38.p13.fa"
STRAT_BED="GRCh38_stratifications.bed"
```

### 5.2 Run hap.py

```bash
singularity run hap.py_0.3.15.sif /opt/hap.py/bin/hap.py     ${TRUTH_VCF}     ${QUERY_VCF}     -f ${CONFIDENT_BED}     -r ${REF}     -X -V     -o ${SAMPLE}_${CALLER}_${CONDITION}_happy     --engine vcfeval     --stratification ${STRAT_BED}
```

### 5.3 hap.py parameters explained

| Parameter | Purpose |
|-----------|---------|
| `-f` | Confident regions BED (restricts evaluation to high-confidence loci) |
| `-r` | Reference FASTA |
| `-X` | Enable extended output (additional metrics) |
| `-V` | Verbose logging |
| `-o` | Output prefix |
| `--engine vcfeval` | Use vcfeval comparison engine (recommended for INDELs) |
| `--stratification` | Stratification BED for region-specific metrics |

### 5.4 Output files

hap.py produces several output files; the **summary CSV** is consumed by `calculate_fp_burden.py`:

| File | Description |
|------|-------------|
| `*.summary.csv` | **Main input for burden framework** — per-stratification metrics |
| `*.vcf.gz` | Annotated VCF with TP/FP/FN labels |
| `*.runinfo.json` | Run metadata and parameters |
| `*.metrics.json` | JSON-formatted metrics |

### 5.5 Required columns in summary CSV

The burden framework requires these columns from the hap.py summary CSV:

```
Type, Subtype, Subset, Filter, Genotype,
METRIC.Recall, METRIC.Precision,
QUERY.FP, QUERY.TP, TRUTH.TOTAL, Subset.IS_CONF.Size
```

---

## Complete Workflow Script

Below is a bash script that automates the full workflow for one sample across all conditions. Adapt paths and sample names as needed.

```bash
#!/bin/bash
# run_workflow.sh — Complete IUPAC consensus workflow for one sample

set -euo pipefail

# Configuration
SAMPLE="HG001"
R1="${SAMPLE}_R1.fastq.gz"
R2="${SAMPLE}_R2.fastq.gz"
REF="GRCh38.p13.fa"
TRUTH_VCF="${SAMPLE}_GIAB_v4.2.1_truth.vcf.gz"
CONFIDENT_BED="HG001_GRCh38_confident_regions.bed"
STRAT_BED="GRCh38_stratifications.bed"

# Step 1: Build consensus references (run once)
# (See Step 2 above — assumes consensus FASTAs and indices already built)

# Step 2: Align reads
echo "=== Alignment ==="

# Baseline: BWA-MEM / GRCh38
if [ ! -f ${SAMPLE}_bwa_GRCh38.sort.bam ]; then
    bwa mem -t 32 ${REF} ${R1} ${R2} | samtools view -bS - > ${SAMPLE}_bwa_GRCh38.bam
    samtools sort -@ 32 ${SAMPLE}_bwa_GRCh38.bam -o ${SAMPLE}_bwa_GRCh38.sort.bam
    samtools index ${SAMPLE}_bwa_GRCh38.sort.bam
    rm ${SAMPLE}_bwa_GRCh38.bam
fi

# Base: novoAlign / GRCh38
if [ ! -f ${SAMPLE}_novo_GRCh38.sort.bam ]; then
    novoalign -a --tune NOVASEQ -r Random -o BAM -c 32         -d GRCh38.p13.ndx -f ${R1} ${R2} > ${SAMPLE}_novo_GRCh38.bam
    novosort --pcrFree --md --kt -i -o ${SAMPLE}_novo_GRCh38.sort.bam ${SAMPLE}_novo_GRCh38.bam
    samtools index ${SAMPLE}_novo_GRCh38.sort.bam
    rm ${SAMPLE}_novo_GRCh38.bam
fi

# Consensus: novoAlign / IUPAC (example: gnomAD AF10)
for DB in gnomAD_AF10 gnomAD_AF30 1KG_AF10 1KG_AF30 dbSNP_CAF70 dbSNP_CAF90; do
    if [ ! -f ${SAMPLE}_novo_${DB}.sort.bam ]; then
        novoalign -a --tune NOVASEQ -r Random -o BAM -c 32             -d GRCh38_${DB}.ndx -f ${R1} ${R2} > ${SAMPLE}_novo_${DB}.bam
        novosort --pcrFree --md --kt -i -o ${SAMPLE}_novo_${DB}.sort.bam ${SAMPLE}_novo_${DB}.bam
        samtools index ${SAMPLE}_novo_${DB}.sort.bam
        rm ${SAMPLE}_novo_${DB}.bam
    fi
done

# Step 3: Variant calling
echo "=== Variant Calling ==="

for BAM in ${SAMPLE}_bwa_GRCh38.sort.bam ${SAMPLE}_novo_GRCh38.sort.bam ${SAMPLE}_novo_*.sort.bam; do
    CONDITION=$(echo ${BAM} | sed 's/.*_novo_//' | sed 's/.sort.bam//' | sed 's/.*_bwa_/bwa_/')

    for CALLER in bcftools freebayes gatk; do
        OUT_VCF="${SAMPLE}_${CALLER}_${CONDITION}.vcf.gz"

        if [ -f ${OUT_VCF} ]; then
            echo "Skipping ${OUT_VCF} (exists)"
            continue
        fi

        case ${CALLER} in
            bcftools)
                bcftools mpileup -Ou -f ${REF} ${BAM} |                     bcftools call -f GQ -mv -Oz -o ${OUT_VCF}
                bcftools index ${OUT_VCF}
                ;;
            freebayes)
                freebayes -f ${REF} -F 0.01 -C 5 ${BAM} |                     bgzip > ${OUT_VCF}
                bcftools index ${OUT_VCF}
                ;;
            gatk)
                gatk HaplotypeCaller -R ${REF} -I ${BAM} -O ${OUT_VCF}
                ;;
        esac

        echo "Called: ${OUT_VCF}"
    done
done

# Step 4: Benchmarking
echo "=== Benchmarking ==="

for VCF in ${SAMPLE}_*_*.vcf.gz; do
    # Parse caller and condition from filename
    PARTS=($(echo ${VCF} | tr '_' ' ' | sed 's/.vcf.gz//'))
    CALLER=${PARTS[1]}
    CONDITION=${PARTS[2]}

    OUT_PREFIX="${SAMPLE}_${CALLER}_${CONDITION}_happy"

    if [ -f ${OUT_PREFIX}.summary.csv ]; then
        echo "Skipping ${OUT_PREFIX} (exists)"
        continue
    fi

    singularity run hap.py_0.3.15.sif /opt/hap.py/bin/hap.py         ${TRUTH_VCF} ${VCF}         -f ${CONFIDENT_BED}         -r ${REF}         -X -V         -o ${OUT_PREFIX}         --engine vcfeval         --stratification ${STRAT_BED}

    echo "Benchmarked: ${OUT_PREFIX}"
done

echo "=== Workflow complete for ${SAMPLE} ==="
```

---

## Output File Tree

After running the complete workflow, the directory structure should look like:

```
project/
├── references/
│   ├── GRCh38.p13.fa                    # Baseline reference
│   ├── GRCh38.p13.ndx                   # novoAlign index
│   ├── GRCh38_gnomAD_AF10.fa          # IUPAC consensus
│   ├── GRCh38_gnomAD_AF10.ndx         # Consensus index
│   ├── GRCh38_gnomAD_AF30.fa
│   ├── GRCh38_1KG_AF10.fa
│   ├── GRCh38_1KG_AF30.fa
│   ├── GRCh38_dbSNP_CAF70.fa
│   └── GRCh38_dbSNP_CAF90.fa
├── alignments/
│   ├── HG001_bwa_GRCh38.sort.bam      # Baseline alignment
│   ├── HG001_novo_GRCh38.sort.bam     # Base alignment
│   ├── HG001_novo_gnomAD_AF10.sort.bam
│   └── ...
├── variants/
│   ├── HG001_bcftools_bwa_GRCh38.vcf.gz
│   ├── HG001_bcftools_novo_gnomAD_AF10.vcf.gz
│   ├── HG001_freebayes_bwa_GRCh38.vcf.gz
│   └── ...
├── happy/
│   ├── HG001_bcftools_bwa_GRCh38_happy.summary.csv
│   ├── HG001_bcftools_novo_gnomAD_AF10_happy.summary.csv
│   └── ...
└── burden/
    └── (output from calculate_fp_burden.py)
```

---

## Troubleshooting

### novoAlign: "ambiguous base in reference"

- **Cause:** BWA-MEM index files mixed with novoAlign.
- **Fix:** Ensure `.ndx` files were built with `novoindex`, not `bwa index`.

### hap.py: "No confident regions overlap"

- **Cause:** Chromosome naming mismatch (e.g., `chr1` vs `1`).
- **Fix:** Ensure reference FASTA, BED files, and VCFs use consistent chromosome names (`chr1`–`chr22`, `chrX`, `chrY`, `chrM`).

### BCFtools mpileup: slow performance

- **Cause:** No region restriction; processing entire genome.
- **Fix:** Add `-R` with target regions if analysing subset, or ensure sufficient RAM (>32 GB recommended).

### FreeBayes: out of memory

- **Cause:** Large BAM files with high coverage.
- **Fix:** Use `--region` to parallelise by chromosome, or increase system RAM.

### GATK: "contig not in reference dictionary"

- **Cause:** BAM header does not match reference FASTA.
- **Fix:** Regenerate BAM with correct reference, or use `samtools reheader`.

### IUPAC consensus: coordinate drift

- **Cause:** novoUtil substituted indels or complex variants.
- **Fix:** Ensure `--exclude-types indels` was used in Step 1. Verify with the length check in Step 2.3.

---

## Next Steps

After generating hap.py summary CSVs, proceed to the burden calculation framework:

```bash
python scripts/calculate_fp_burden.py \
    --manifest manifests/samples_manifest.tsv \
    --output-prefix burden_results/ \
    --regions AllAutosomes lowmappabilityall segdups CMRG MHC \
    --filters ALL
```

See [`docs/usage_manual.md`](usage_manual.md) for full documentation.
