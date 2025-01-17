Mouse RNAseq Workflow
================
Jessie Tignor
January 16th, 2024

### Introduction

This repository provides a comprehensive step-by-step guide for analyzing mouse RNA-seq data to identify up and downregulated genes. The workflow assumes no prior coding experience and includes all commands, package installations, and reasoning behind chosen input options.

### 1. Environment Setup

#### 1a. Install Conda Package Manager

Conda simplifies package installation and environment management.

``` bash
# Download Miniconda installer
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

# Run installer
bash Miniconda3-latest-Linux-x86_64.sh

# Follow the prompts and restart the terminal
source ~/.bashrc
```

#### 1b. Create and Activate Conda Environment

``` bash
# Create an environment named 'rna_seq'
conda create -n rna_seq python=3.8 -y

# Activate the environment
conda activate rna_seq
```

### 2. Package Installation

#### 2a. Install Required Tools and R Packages

``` bash
# Install fastp for trimming and filtering
conda install -c bioconda fastp -y

# Install HISAT2 for alignment
conda install -c bioconda hisat2 -y

# Install SAMtools for file conversion and indexing
conda install -c bioconda samtools -y

# Install featureCounts for gene count matrix generation
conda install -c bioconda subread -y

# Install R and RStudio (if not pre-installed)
conda install -c r r-base r-essentials -y

# Install DESeq2 from Bioconductor in R
R
> if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
> BiocManager::install("DESeq2")
> install.packages("ggplot2")
> install.packages("pheatmap")
> BiocManager::install("clusterProfiler")
> BiocManager::install("org.Mm.eg.db")  # For mouse annotation
q()
```
### 3. Quality Control
Trim and filter raw reads use fastp to remove low-quality reads and trim adapters.
``` bash
# Run fastp for paired-end data
fastp -i AI1_R1.fastq.gz -I AI1_R2.fastq.gz \
      -o AI1_R1_trimmed.fastq.gz -O AI1_R2_trimmed.fastq.gz \
      --html AI1_fastp_report.html --json AI1_fastp_report.json

# Repeat for all samples (AI4, AI7, AI9, parental)
```
### 4. Read Alignment
Align Reads
``` bash
# Download genome index
wget -O grcm38_index.tar.gz https://example.com/path_to_hisat2_index.tar.gz
tar -xvzf grcm38_index.tar.gz

# Align reads using HISAT2
hisat2 -x grcm38/genome \
       -1 AI1_R1_trimmed.fastq.gz -2 AI1_R2_trimmed.fastq.gz \
       -S AI1_aligned.sam --summary-file AI1_alignment_summary.txt

# Repeat for all samples
```















