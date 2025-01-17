Mouse RNAseq Workflow
================
Jessie Tignor
January 16th, 2024

### Introduction

This repository provides a comprehensive step-by-step guide for analyzing mouse RNA-seq data to identify up and downregulated genes. The workflow assumes no prior coding experience and includes all commands, package installations, and reasoning behind chosen input options.

### 1. Environment Setup

#### Install Conda Package Manager

Conda simplifies package installation and environment management.

The easiest way to set up the RNAseq environment is by using a bio-cookie template. This template will download all the necessary tools into a new virtual environment, so that there are no disruptions with your current computer packaging versions. It is high reccomended to run on all of these commands on a high performance cluster.

The template will also ask if you would like to download a human or mouse genome for transcriptome alignment. If you want to align to a genome which is neither human or murine, then you must select "None" in the option and download the genome and annotation file separately. Below are a list of the tools which are downloaded during the installation steps. See here for a list of availabel genomes (<https://genome.ucsc.edu/cgi-bin/hgGateway> / <http://useast.ensembl.org/info/data/ftp/index.html> )

detailed pipeline tailored to mouse RNA sequencing analysis to identify pathways associated with up and downregulated genes
Mouse RNA-Seq Analysis Workflow
# Introduction
This repository provides a comprehensive step-by-step guide for analyzing RNA-seq data to identify up- and downregulated genes and their associated pathways. The workflow assumes no prior coding experience and includes all commands, package installations, and reasoning behind chosen input options.

System and Software Setup

System Specifications:

OS: Linux Rocky 9

Memory: 256GB RAM

CPUs: 16 cores

Data Source: Raw sequencing reads and quality metrics were provided by the Illumina sequencing core (DRAGEN fastqc reports, read metrics, and read instrument analytics).

Samples: 4 mutants (AI1, AI4, AI7, AI9) and 1 parental control.
