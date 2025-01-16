---
title: "Mouse RNAseq"
author: "Jessie Tignor"
date: "January 16th, 2025"
output: 
  github_document:
    fig_width: 5
    fig_height: 5
---
```{r setup, include=FALSE}
knitr::opts_chunk$set(base.dir = TRUE)
```


### Introduction 
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
