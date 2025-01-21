Mouse RNAseq Workflow
================
Jessie Tignor
January 21st, 2024

### Introduction

This repository documents the steps to process and analyze mouse RNA sequencing data for identifying upregulated and downregulated genes. The workflow assumes no prior coding experience and includes all commands and package installations.

The analysis was performed using R and RStudio on a Linux Rocky 9 system with 256GB of RAM and 16 CPUs. Raw sequencing data was provided by the Illumina sequencing core. The samples include 4 mutants (AI1, AI4, AI7, AI9) and one control (parental). To replicate this workflow, replace the sample file names with your own.

### 1. Environment Setup

#### 1a. Install Conda Package Manager

Conda simplifies package installation and environment management. These commands download the latest Linux installer, rename it to a shorter file name, perform a silent install, and then delete the installer:

``` bash
mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm ~/miniconda3/miniconda.sh
# Refresh your terminal
source ~/miniconda3/bin/activate
```

#### 1b. Create and Activate Conda Environment

Creating a conda environment ensures that your analysis is isolated from other projects, helps manage dependencies more efficiently, and guarantees reproducibility. To create an environment:

``` bash
conda create -n rna_seq_analysis r-base=4.2
conda activate rna_seq_analysis
```

### 2. Package Installation

#### 2a. Install Required Tools and R Packages

``` bash
conda install -c bioconda fastp hisat2 samtools subread -y
```

### 3. Quality Control
To ensure the quality of your raw sequencing data, use fastp for trimming and filtering reads. For paired end data (gzip compressed):
``` bash
fastp -i AI1_R1.fastq.gz -I AI1_R2.fastq.gz -o AI1_R1_trimmed.fastq.gz -O AI1_R2_trimmed.fastq.gz -h AI1.fastp.html -j AI1.fastp.html

# Repeat for all samples (AI4, AI7, AI9, parental).
```
HTML and JSON reports are generated as well. Both provide a summary of filtering results.

### 4. Read Alignment

#### 4a. Download Genome Index
Download the pre-built genome_tran index for GRCm38 from the HISAT2 website. This index includes both the genome and transcript annotations.
``` bash
wget -c https://genome-idx.s3.amazonaws.com/hisat/genome_tran.tar.gz
tar -xvzf genome_tran.tar.gz
mv genome_tran GRCm38_index
```
#### 4b. Align Reads
HISAT2 maps RNA-seq reads to a reference genome for accurate alignment.
``` bash
hisat2 -p 16 -x GRCm38_index/genome_tran -1 AI1_R1_trimmed.fastq.gz -2 AI1_R2_trimmed.fastq.gz -S AI1_aligned.sam

# Repeat for all samples
```






``` bash
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("DESeq2")
```








