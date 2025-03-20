Mouse RNAseq Workflow
================
Jessie Tignor
January 21st, 2024

### Introduction
This repository documents the steps to process and analyze mouse RNA sequencing data for identifying upregulated and downregulated genes. The workflow assumes no prior coding experience and includes all commands and package installations.

The analysis was performed using R and RStudio on a Linux Rocky 9 system with 256GB of RAM and 16 CPUs. Raw sequencing data was provided by the Illumina sequencing core. The samples include 4 mutants (AI1, AI4, AI7, AI9) and one control (parental). To replicate this workflow, replace the sample file names with your own.

### 1. Environment Setup
#### 1a. Install Conda Package Manager
`Anaconda` simplifies package installation and environment management. I installed `Miniconda` instead of the full distribution to keep the environment lightweight and avoid unnecessary pre-installed packages. These commands download the latest Linux installer, rename it to a shorter file name, perform a silent install, and then delete the installer:
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
**HTML** and **JSON** reports are generated as well. Both provide a summary of filtering results.

### 4. Read Alignment
#### 4a. Download Genome Index
Download the pre-built genome_tran index for **GRCm38** from the `HISAT2` website. This index includes both the genome and transcript annotations.
``` bash
wget -c https://genome-idx.s3.amazonaws.com/hisat/genome_tran.tar.gz
tar -xvzf genome_tran.tar.gz
mv genome_tran GRCm38_index
```

#### 4b. Align Reads
HISAT2 maps RNA-seq reads to a reference genome for accurate alignment.
``` bash
hisat2 -p 15 -x GRCm38_index/genome_tran -1 AI1_R1_trimmed.fastq.gz -2 AI1_R2_trimmed.fastq.gz -S AI1_aligned.sam

# Repeat for all samples
```
`-p` specifies the number of CPUS threads to use. I used 15 CPU threads because I have 16 CPUS and left one available for system processes. The alignment process is much faster when distributed across multiple cores.

### 5. Convert and Sort BAM Files  
After alignment, convert SAM files to BAM format, sort them, and create index files using `SAMtools`:  
```bash
samtools view -bS AI1_aligned.sam > AI1.bam
samtools sort AI1.bam -o AI1_sorted.bam
samtools index AI1_sorted.bam

# Repeat for all samples
```
Sorting optimizes read lookup speeds, and indexing allows for effician downstream analyses.

### 6. Gene Quantification with featureCounts  
#### 6a. Download Annotation File  
The GTF annotation file for **GRCm38** was downloaded from the **Ensembl database** (version 101) to ensure accurate gene mapping.  
```bash
wget ftp://ftp.ensembl.org/pub/release-101/gtf/mus_musculus/Mus_musculus.GRCm38.101.gtf.gz
gunzip Mus_musculus.GRCm38.101.gtf.gz
mv Mus_musculus.GRCm38.101.gtf GRCm38_annotation.gtf
```
This file containts gene annotations needed for quantifying mapped reads.

#### 6b. Generate Raw Gene Counts
Use `featureCounts` to create a count matrix by mapping aligned reads to annotated genes:
``` bash
featureCounts -T 15 -a GRCm38_annotation.gtf -o gene_counts.txt AI1_sorted.bam AI4_sorted.bam AI7_sorted.bam AI9_sorted.bam Parental_sorted.bam
```
`-T` specifies CPU threads.

The output, `gene_counts.txt` contains raw gene counts for all samples

### Citations:
1. Shifu Chen. 2023. Ultrafast one-pass FASTQ data preprocessing, quality control, and deduplication using fastp. iMeta 2: e107. Available here: <https://github.com/OpenGene/fastp> <https://doi.org/10.1002/imt2.107>.
