# 🧬 AistSeq

**An in-house easy-to-purify Tn5-based plasmid sequencing platform using a compact benchtop sequencer**

---

## 📖 Overview

**AistSeq** is a streamlined, laboratory-scale plasmid sequencing pipeline designed for **rapid and cost-effective sequence verification** using Tn5-based library preparation and short-read sequencing (e.g., Illumina iSeq100).

This repository provides a **fully automated computational workflow** for:

* Reference-guided plasmid validation
* De novo plasmid assembly
* Coverage and quality assessment
* Mutation detection
* Functional annotation

The pipeline is optimized for **in-house workflows**, minimizing complexity while maintaining analytical rigor.

---

## 🎯 Key Features

* 🧪 **Dual-mode analysis**

  * Reference-guided assembly
  * De novo assembly (Unicycler)

* ⚙️ **End-to-end automation**

  * FASTQ → Assembly → QC → Annotation → Report

* 📊 **Comprehensive QC metrics**

  * Coverage depth & uniformity
  * GC bias analysis
  * Assembly error rate

* 📈 **Interactive visualization**

  * Coverage plots
  * Mutation plots
  * Summary statistics

* 📄 **Single HTML report per plasmid**

  * Integrates all results for easy interpretation

---

## 🧬 Pipeline Architecture

### 1. Reference-guided workflow

Implemented in:
📄 

**Steps:**

1. Quality control

   * `fastp` trimming and filtering

2. Alignment

   * `Bowtie2` mapping to reference

3. Coverage analysis

   * `samtools depth`
   * Low-coverage region detection

4. Variant calling

   * `bcftools mpileup + call + filter`

5. Consensus generation

   * `bcftools consensus`

6. Annotation

   * `pLannotate`

7. Report generation

   * HTML integration of all outputs

---

### 2. De novo workflow

Also implemented in:
📄 

**Steps:**

1. Quality control (`fastp`)
2. Assembly:

   * Conservative → Normal → Bold (fallback strategy)
3. Read alignment (`BWA`)
4. Coverage analysis
5. Annotation (`pLannotate`)
6. Report generation

---

## 📂 Repository Structure

```
AistSeq/
│
├── Plasmid_assembly-2-Options-comand.sh   # Main pipeline
├── generate_plots.py                      # Coverage & mutation visualization
├── compute_coverage_stats.py              # QC metrics + GC bias
├── generate_summary_plot-V6.py            # Summary bar plots
├── calculate_fasta.py                     # Sequence length calculation
│
├── Adapter_and_reference_database/        # Adapter + reference DB
└── output/                                # Results
```

---

## ⚙️ Requirements

### Software

* Python ≥ 3.8
* Conda
* fastp
* Bowtie2
* BWA
* samtools
* bcftools
* bedtools
* Unicycler
* pLannotate

---

## 📥 Input Format

### 1. ZIP file (reads)

Contains paired-end reads:

```
R1.fastq.gz
R2.fastq.gz
```

---

### 2. Metadata file (TSV)

Example:

```
Primary key	read1	read2	expected_reads	plasmid_name	output_directory	Type_of_run	Reference	size
97	R1.fastq.gz	R2.fastq.gz	10000	Mapping01	Sample1	y	pUC18.fa	3000
160	R1.fastq.gz	R2.fastq.gz	10000	Denovo01	Sample2	n		15000
```

**Key fields:**

| Column      | Description                           |
| ----------- | ------------------------------------- |
| Type_of_run | `y` = reference, `n` = de novo        |
| Reference   | Required for reference mode           |
| size        | Optional (auto-calculated if missing) |

---

## 🚀 Usage

### Basic command

```bash
bash Plasmid_assembly-2-Options-comand.sh \
    -i input.zip \
    -t metadata.txt \
    -o output_dir
```

---

## 📊 Output

For each plasmid:

### 1. Final report

```
report_of_<plasmid>.html
```

Includes:

* QC summary
* Coverage plot
* Mutation plot
* Annotation
* Assembly statistics

---

### 2. Key files

| File                  | Description              |
| --------------------- | ------------------------ |
| `*_consensus.fa`      | Final assembled sequence |
| `*.gbk`               | Annotated plasmid        |
| `coverage.txt`        | Depth per position       |
| `mutations.vcf`       | Variant calls            |
| `coverage_stats.html` | QC metrics               |

---

## 📈 Visualization Scripts

### 1. Coverage & Mutation Plot

📄 

```bash
python generate_plots.py coverage.txt coverage.html
python generate_plots.py mutations.vcf mutations.html
```

---

### 2. Coverage QC + GC Bias

📄 

```bash
python compute_coverage_stats.py \
    coverage.txt \
    plasmid.fa \
    coverage_stats.html \
    GC_bias.png \
    stats.tsv
```

---

### 3. Summary Bar Plots

📄 

```bash
python generate_summary_plot-V6.py output_dir \
    <12 count values>
```

Generates:

* `<5 kb`, `5–10 kb`, `>10 kb`
* Overall success distribution

---

### 4. Sequence Size Calculation

📄 

```bash
python calculate_fasta.py plasmid.fa
```

---

## 📊 Key Metrics Interpretation

### Coverage

* ≥20×: reliable sequencing
* ≥95% breadth: near-complete coverage

### Error rate

* <1e-3 → high-quality assembly
* 1e-3–1e-2 → acceptable
* > 1e-2 → potential issues

### GC bias

* Pearson r ≈ 0 → minimal bias
* High |r| → sequencing bias

---

## 🧪 Experimental Context

This pipeline supports the AistSeq framework:

* Tn5-based plasmid library preparation
* Minimal purification design (His10 + Protein A)
* Direct sequencing from *E. coli* lysate
* Compact benchtop sequencing (iSeq100)

---

## ⚠️ Notes

* Designed for **plasmid-scale genomes (~1–20 kb)**
* Best performance with:

  * High coverage
  * Single-contig assembly
* Multi-contig output indicates:

  * Assembly ambiguity
  * Low coverage
  * Structural complexity

---

## 📌 Summary

AistSeq provides:

* A **compact, reproducible plasmid sequencing workflow**
* Integrated **assembly + QC + annotation**
* Optimized for **in-house synthetic biology pipelines**

---
