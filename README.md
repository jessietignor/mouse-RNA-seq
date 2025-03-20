Mouse RNAseq Workflow
================
Jessie Tignor
March 20th, 2025

### Introduction
This repository documents the steps to process and analyze mouse RNA sequencing data for identifying upregulated and downregulated genes. The workflow assumes no prior coding experience and includes all commands and package installations.

The analysis was performed using R and RStudio on a **Linux Rocky 9** system with 256GB of RAM and 16 CPUs. Raw sequencing data was provided by the **Illumina** sequencing core. The samples include 4 mutants (AI1, AI4, AI7, AI9) and one control (parental). To replicate this workflow, replace the sample file names with your own.

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
`HISAT2` maps RNA-seq reads to a reference genome for accurate alignment.
``` bash
hisat2 -p 15 -x GRCm38_index/genome_tran -1 AI1_R1_trimmed.fastq.gz -2 AI1_R2_trimmed.fastq.gz -S AI1_aligned.sam

# Repeat for all samples
```
`-p` specifies the number of CPUS threads to use. I used 15 CPU threads because I have 16 CPUS and left one available for system processes. The alignment process is much faster when distributed across multiple cores.

### 5. Convert and Sort BAM Files  
After alignment, convert **SAM** files to **BAM** format, sort them, and create index files using `SAMtools`:  
```bash
samtools view -bS AI1_aligned.sam > AI1.bam
samtools sort AI1.bam -o AI1_sorted.bam
samtools index AI1_sorted.bam

# Repeat for all samples
```
Sorting optimizes read lookup speeds, and indexing allows for effician downstream analyses.

### 6. Gene Quantification with featureCounts  
#### 6a. Download Annotation File  
The **GTF** annotation file for **GRCm38** was downloaded from the **Ensembl database** (Ensembl Release 102) to ensure accurate gene mapping.  
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
1. conda contributors. conda: A system-level, binary package and environment manager running on all major operating systems and platforms. Available here: <https://github.com/conda/conda> & <https://www.anaconda.com/docs/getting-started/miniconda/install#quickstart-install-instructions>.
2. Shifu Chen. 2023. Ultrafast one-pass FASTQ data preprocessing, quality control, and deduplication using fastp. iMeta 2: e107. Available here: <https://github.com/OpenGene/fastp> & <https://doi.org/10.1002/imt2.107>.
3. Kim, D., Langmead, B. & Salzberg, S. HISAT: a fast spliced aligner with low memory requirements. *Nat Methods* 12, 357–360 (2015). Available here: <https://doi.org/10.1038/nmeth.3317> & <https://daehwankimlab.github.io/hisat2/main/>.
4. Petr Danecek, James K Bonfield, Jennifer Liddle, John Marshall, Valeriu Ohan, Martin O Pollard, Andrew Whitwham, Thomas Keane, Shane A McCarthy, Robert M Davies, Heng Li, Twelve years of SAMtools and BCFtools, *GigaScience*, Volume 10, Issue 2, February 2021, giab008. Available here: <https://doi.org/10.1093/gigascience/giab008> & <https://www.htslib.org/>.
5. Peter W Harrison, M Ridwan Amode, Olanrewaju Austine-Orimoloye, Andrey G Azov, Matthieu Barba, If Barnes, Arne Becker, Ruth Bennett, Andrew Berry, Jyothish Bhai, Simarpreet Kaur Bhurji, Sanjay Boddu, Paulo R Branco Lins, Lucy Brooks, Shashank Budhanuru Ramaraju, Lahcen I Campbell, Manuel Carbajo Martinez, Mehrnaz Charkhchi, Kapeel Chougule, Alexander Cockburn, Claire Davidson, Nishadi H De Silva, Kamalkumar Dodiya, Sarah Donaldson, Bilal El Houdaigui, Tamara El Naboulsi, Reham Fatima, Carlos Garcia Giron, Thiago Genez, Dionysios Grigoriadis, Gurpreet S Ghattaoraya, Jose Gonzalez Martinez, Tatiana A Gurbich, Matthew Hardy, Zoe Hollis, Thibaut Hourlier, Toby Hunt, Mike Kay, Vinay Kaykala, Tuan Le, Diana Lemos, Disha Lodha, Diego Marques-Coelho, Gareth Maslen, Gabriela Alejandra Merino, Louisse Paola Mirabueno, Aleena Mushtaq, Syed Nakib Hossain, Denye N Ogeh, Manoj Pandian Sakthivel, Anne Parker, Malcolm Perry, Ivana Piližota, Daniel Poppleton, Irina Prosovetskaia, Shriya Raj, José G Pérez-Silva, Ahamed Imran Abdul Salam, Shradha Saraf, Nuno Saraiva-Agostinho, Dan Sheppard, Swati Sinha, Botond Sipos, Vasily Sitnik, William Stark, Emily Steed, Marie-Marthe Suner, Likhitha Surapaneni, Kyösti Sutinen, Francesca Floriana Tricomi, David Urbina-Gómez, Andres Veidenberg, Thomas A Walsh, Doreen Ware, Elizabeth Wass, Natalie L Willhoft, Jamie Allen, Jorge Alvarez-Jarreta, Marc Chakiachvili, Bethany Flint, Stefano Giorgetti, Leanne Haggerty, Garth R Ilsley, Jon Keatley, Jane E Loveland, Benjamin Moore, Jonathan M Mudge, Guy Naamati, John Tate, Stephen J Trevanion, Andrea Winterbottom, Adam Frankish, Sarah E Hunt, Fiona Cunningham, Sarah Dyer, Robert D Finn, Fergal J Martin, Andrew D Yates, Ensembl 2024, *Nucleic Acids Research*, Volume 52, Issue D1, 5 January 2024, Pages D891–D899. Available here: <https://doi.org/10.1093/nar/gkad1049> & <https://nov2020.archive.ensembl.org/Mus_musculus/Info/Index>.
6. Liao Y, Smyth GK and Shi W (2014). featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. *Bioinformatics*, 30(7):923-30. Available here: <https://pubmed.ncbi.nlm.nih.gov/24227677/> & <https://subread.sourceforge.net/featureCounts.html>.
7. Love, M.I., Huber, W. & Anders, S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biol 15, 550 (2014). Available here: <https://doi.org/10.1186/s13059-014-0550-8> & <https://bioconductor.org/packages/release/bioc/html/DESeq2.html>.
8. 
