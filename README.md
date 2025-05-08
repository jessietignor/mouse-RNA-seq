Mouse RNAseq Workflow
================
Jessie Tignor

May 8th, 2025

### Introduction
This repository documents the steps to process and analyze mouse RNA sequencing data for identifying upregulated and downregulated genes. The workflow assumes **no prior coding experience** and includes all commands and package installations.

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
conda install -c bioconda fastp hisat2 samtools subread gffread -y
```

### 3. Quality Control
To ensure the quality of your raw sequencing data, use `fastp` for trimming and filtering reads. For paired end data (gzip compressed):
``` bash
fastp -i AI1_R1.fastq.gz -I AI1_R2.fastq.gz -o AI1_R1_trimmed.fastq.gz -O AI1_R2_trimmed.fastq.gz -h AI1.fastp.html -j AI1.fastp.html

# Repeat for all samples (AI4, AI7, AI9, parental).
```
Flags:
- `-i/-I`: Input R1 and R2 FASTQ files
- `-o/-O`: Output trimmed R1 and R2 FASTQ files
- `-h/-j`: Generate HTML and JSON reports for quality metrics

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
Flags:
- `-p`: Number of threads
- `-x`: Index path
- `-1/-2`: Paired-end reads
- `-S`: Output SAM file

I used 15 CPU threads because I have 16 CPUS and left one available for system processes. The alignment process is much faster when distributed across multiple cores.

### 5. Convert and Sort BAM Files  
After alignment, convert **SAM** files to **BAM** format, sort them, and create index files using `SAMtools`:  
```bash
samtools view -bS AI1_aligned.sam > AI1.bam
samtools sort AI1.bam -o AI1_sorted.bam
samtools index AI1_sorted.bam

# Repeat for all samples
```
Sorting optimizes read lookup speeds, and indexing allows for efficient downstream analyses.

### 6. Check Strandedness
You typically only need to check strandedness on **one representative sample**, such as AI1, assuming all samples were prepared with the same library construction protocol.
#### 6a. Convert GTF to BED
We'll use `gffread` and `awk` to prepare a BED file from the GTF annotation. This is required for the `infer_experiment.py` tool:
```bash
gffread Mus_musculus.GRCm38.101.gtf -T -o gffread_annotation.gtf
awk '$3 == "exon" {print $1"\t"$4-1"\t"$5"\t"$10"\t0\t"$7}' Mus_musculus.GRCm38.101.gtf | tr -d '";' > annotation.bed
```
Flags:
-  `tr -d '";'`: Removes quotes and semicolons to meet BED format requirements
#### 6b. Run Strandedness Check
```bash
infer_experiment.py -i AI1_sorted.bam -r annotation.bed
```
Interpret the output to determine the library type. Look at the **fraction of reads** explained by each orientation:
- If the majority matches `1++,1--,2+-,2-+` → `-s 1` (forward stranded)
- If the majority matches `1+-,1-+,2++,2--` → `-s 2` (reverse stranded)
- If reads are evenly split or ambiguous → `-s 0` (unstranded)

For `featureCounts`, there are values that correspond to each strandedness:
- Unstranded → use `-s 0`
- Stranded (forward) → use `-s 1`
- Stranded (reverse) → use `-s 2`

In our case, the data was **reversely stranded**, so we’ll set `-s 2` in the next step.
> **Note:** This value is important when using `featureCounts`. Incorrect values will lead to undercounting or misassignment of reads.

#### 7. Generate Raw Gene Counts
Use `featureCounts` to create a count matrix by mapping aligned reads to annotated genes:
```bash
featureCounts -T 15 -s 2 -p -a Mus_musculus.GRCm38.101.gtf -o gene_counts.txt \
  AI1_sorted.bam AI4_sorted.bam AI7_sorted.bam AI9_sorted.bam Parental_sorted.bam
```
Flags:
- `-T`: Number of threads
- `-s 2`: Strandedness
- `-p`: Paired-end reads
- `-a`: GTF annotation file
- `-o`: Output count file
The output file `gene_counts.txt` contains raw gene counts for each sample and will be used for differential expression analysis. You can open this file in any text editor or spreadsheet program to inspect the counts.

### 8. Differential Expression Analysis in R (DESeq2)
Switch to RStudio and load the necessary libraries:
```r
library(DESeq2)
library(biomaRt)
```

#### 8a. Load Counts and Metadata
```r
# Read in count matrix
counts <- read.delim("gene_counts.txt", comment.char="#", check.names=FALSE)
counts <- counts[, -c(2:6)]  # Remove extra columns (Chr, Start, End, Strand, Length)
rownames(counts) <- counts$Geneid
counts <- counts[, -1]

# Rename columns if needed to match your metadata
colnames(counts) <- c("AI1", "AI4", "AI7", "AI9", "Parental")

# Create metadata
sample_info <- data.frame(
  row.names = colnames(counts),
  condition = c("AI", "AI", "AI", "AI", "Parental")
)
```

#### 8b. Run DESeq2
```r
dds <- DESeqDataSetFromMatrix(countData = counts, colData = sample_info, design = ~ condition)
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "AI", "Parental"))
res_df <- as.data.frame(res)
res_df$ensembl_id <- rownames(res_df)
```

#### 8c. Annotate Results with Gene Name and Symbol
```r
mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
gene_info <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id", "external_gene_name"),
                   filters = "ensembl_gene_id", values = res_df$ensembl_id, mart = mart)

# Merge annotations
res_annotated <- merge(res_df, gene_info, by.x="ensembl_id", by.y="ensembl_gene_id", all.x=TRUE)
res_annotated <- res_annotated[, c("ensembl_id", "entrezgene_id", "external_gene_name", colnames(res_df)[1:6])]
colnames(res_annotated)[1:3] <- c("ENSEMBL", "ENTREZID", "GENENAME")
```

#### 8d. Export to CSV
```r
# All genes
write.csv(res_annotated, "AllAI_vs_Parental_all_genes.csv", row.names=FALSE)

# Upregulated (log2FoldChange > 0)
write.csv(res_annotated[res_annotated$log2FoldChange > 0, ], "AllAI_vs_Parental_upregulated.csv", row.names=FALSE)

# Downregulated (log2FoldChange < 0)
write.csv(res_annotated[res_annotated$log2FoldChange < 0, ], "AllAI_vs_Parental_downregulated.csv", row.names=FALSE)
```

Each file contains all relevant statistics, including raw counts (from input), log2 fold change, p-values, ENSEMBL ID, gene name, and Entrez ID.

### 9. Create Curated Gene Set Heatmaps in R
Switch to RStudio and load the necessary libraries:
```r
library(ggplot2)
library(reshape2)
library(RColorBrewer)
```

#### 9a. Define Gene Sets
Replace or append additional gene lists as needed.
```r
androgen_metabolic <- c("Hsd3b1", "Cyp11a1", "Cyp17a1", "Srd5a1", "Srd5a2")
ar_target <- c("Olfm1", "Slc38a3", "Filip1", "Hapln1", "Prdm1", "Plekha6", "Prg4", "Pygm", "Abtb2", 
               "Tnni1", "Speg", "Rarb", "Nsd1", "Fam180a", "Aknad1", "Nat8l", "Smox", "Car11", "Hepacam",
               "Gck", "Slc7a8", "Ankrd2", "Dlx3", "Sox10", "Lgi3", "Wnt3")
senescence <- c("Ezh2", "Lmnb1", "Cdk2", "Cdk6")
androgen_response <- c("Nkx3-1", "Fkbp5", "Tgm4", "Klk3")
wnt_genes <- c("Wnt2", "Wnt3", "Wnt4", "Wnt6", "Yap1", "Wwtr1")
```

#### 9b. Generate Fold Change Table
Assuming `full_annotated2` contains gene counts per sample:
```r
# Combine all genes into one data frame
curated_genes <- unique(c(androgen_metabolic, ar_target, senescence, androgen_response, wnt_genes))
subset_df <- full_annotated2[full_annotated2$SYMBOL %in% curated_genes, ]

# Calculate log2 fold change (mutant vs parental)
log2_fc <- data.frame(
  gene = subset_df$SYMBOL,
  AI1 = log2((subset_df$AI1 + 1) / (subset_df$Parental + 1)),
  AI4 = log2((subset_df$AI4 + 1) / (subset_df$Parental + 1)),
  AI7 = log2((subset_df$AI7 + 1) / (subset_df$Parental + 1)),
  AI9 = log2((subset_df$AI9 + 1) / (subset_df$Parental + 1))
)
```

#### 9c. Plot Heatmaps by Pathway
Replace `genes` with any of the defined gene sets above.
```r
plot_pathway_heatmap <- function(gene_set, title) {
  df <- full_annotated2[full_annotated2$SYMBOL %in% gene_set, ]
  log2_fc <- data.frame(
    gene = df$SYMBOL,
    AI1 = log2((df$AI1 + 1) / (df$Parental + 1)),
    AI4 = log2((df$AI4 + 1) / (df$Parental + 1)),
    AI7 = log2((df$AI7 + 1) / (df$Parental + 1)),
    AI9 = log2((df$AI9 + 1) / (df$Parental + 1))
  )
  melted <- melt(log2_fc, id.vars = "gene", variable.name = "mutant", value.name = "log2 fold change")
  melted$gene <- factor(melted$gene, levels = rev(unique(log2_fc$gene)))

  ggplot(melted, aes(x = mutant, y = gene, fill = `log2 fold change`)) +
    geom_tile(color = "white") +
    scale_fill_distiller(palette = "YlOrRd", direction = 1) +
    theme_minimal() +
    labs(title = title, x = "Mutant", y = "Gene") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 8),
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
}

# Example heatmap
plot_pathway_heatmap(ar_target, "AR Target Gene Expression")
```
Repeat the `plot_pathway_heatmap()` call for other pathways by passing the appropriate gene list and title.
```r
plot_pathway_heatmap(senescence, "Senescence Gene Expression")
plot_pathway_heatmap(androgen_response, "Androgen Response Gene Expression")
plot_pathway_heatmap(androgen_metabolic, "Androgen Metabolic Gene Expression")
plot_pathway_heatmap(wnt_genes, "WNT Signaling Gene Expression")
```

### Citations:
1. conda contributors. conda: A system-level, binary package and environment manager running on all major operating systems and platforms. Available here: <https://github.com/conda/conda> & <https://www.anaconda.com/docs/getting-started/miniconda/install#quickstart-install-instructions>.

2. Chen, S. (2023). Ultrafast one-pass FASTQ data preprocessing, quality control, and deduplication using fastp. *iMeta*, 2: e107. Available here: <https://github.com/OpenGene/fastp> & <https://doi.org/10.1002/imt2.107>.

3. Kim, D., Paggi, J.M., Park, C., Bennett, C., & Salzberg, S.L. (2019). Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype. *Nature Biotechnology*, 37(8), 907–915. Available here: <https://doi.org/10.1038/s41587-019-0201-4> & <https://daehwankimlab.github.io/hisat2/main/>.

4. Danecek, P., Bonfield, J.K., Liddle, J., Marshall, J., Ohan, V., Pollard, M.O., et al. (2021). Twelve years of SAMtools and BCFtools. *GigaScience*, 10(2), giab008. Available here: <https://doi.org/10.1093/gigascience/giab008> & <https://www.htslib.org/>.

5. Harrison, P.W., Amode, M.R., Austine-Orimoloye, O.G., et al. (2024). Ensembl 2024. *Nucleic Acids Research*, 52(D1), D891–D899. Available here: <https://doi.org/10.1093/nar/gkad1049> & <https://nov2020.archive.ensembl.org/Mus_musculus/Info/Index>.

6. Pertea, G., & Pertea, M. (2020). GFF Utilities: GffRead and GffCompare. *F1000Research*, 9:304. Available here: <https://doi.org/10.12688/f1000research.23297.2> & <https://github.com/gpertea/gffread>.

7. Liao, Y., Smyth, G.K., & Shi, W. (2014). featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. *Bioinformatics*, 30(7), 923–930. Available here: <https://pubmed.ncbi.nlm.nih.gov/24227677/> & <https://subread.sourceforge.net/featureCounts.html>.

8. Love, M.I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*, 15(12), 550. Available here: <https://doi.org/10.1186/s13059-014-0550-8> & <https://bioconductor.org/packages/release/bioc/html/DESeq2.html>.
