Mouse RNAseq Workflow
================
Jessie Tignor
January 16th, 2024

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
To ensure the quality of your raw sequencing data, use fastp for trimming and filtering reads. For paired end data (gzip compressed):
``` bash
fastp -i AI1_R1.fastq.gz -I AI1_R2.fastq.gz -o AI1_R1_trimmed.fastq.gz -O AI1_R2_trimmed.fastq.gz -h AI1.fastp.html -j AI1.fastp.html

# Repeat for all samples (AI4, AI7, AI9, parental).
```
HTML and JSON reports are generated as well. Both provide a summary of filtering results.

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















