# R for Bioinformatics

## Project 1: Phylogenetic Analysis of Coronavirus Genomes and Spike Proteins (phylogenetic_analysis_coronavirus.R)

### Overview:
This R project performs a phylogenetic analysis of SARS-CoV-2 and MERS-CoV genomes and spike protein sequences. 
It includes multiple sequence alignment (MSA), tree construction using Neighbor Joining (NJ) and Maximum Likelihood (ML), and model selection using JC69, GTR, and JTT. 
The trees are visualized with bootstrap analysis and compared using tanglegrams.

### Requirements

Install the required R packages:
install.packages(c("ape", "seqinr", "Biostrings", "tidyverse", "msa", "phangorn", "phytools", "stringr", "taxize"))

### Data

Genome Sequences: SARS-CoV-2 and MERS-CoV genomes (DNA).
Protein Sequences: Spike protein sequences (Amino Acids).
Data is loaded from URLs or local FASTA files.

### Analysis Steps

Data Loading: Read genome and protein sequences from FASTA files.
Multiple Sequence Alignment: Perform MSA using ClustalW.
Phylogenetic Trees: Build NJ trees and compute distance matrices (JC69 for genomes, JTT for proteins).
Bootstrap Analysis: Assess tree reliability with bootstrap resampling.
Maximum Likelihood Trees: Generate and optimize ML trees.
Model Comparison: Compare JC69, GTR, and JTT models using AIC.
Visualization: Export trees and tanglegrams as PDF files.
Outputs

Phylogenetic Trees: NJ and ML trees with bootstrap values.
Tanglegrams: Comparison of NJ vs ML trees for genomes and proteins.

### Significance

This project helps us understand how SARS-CoV-2 (the virus that causes COVID-19) and MERS-CoV (another coronavirus) are related to each other and how they’ve evolved over time. By looking at both the virus genomes and their spike proteins, the analysis shows how different virus strains are connected and how they’ve changed. This kind of information is important because it can help scientists develop better treatments and vaccines.

By using different methods to create evolutionary trees, this project gives a clearer picture of the viruses’ genetic history. It also ensures that the results are reliable by comparing different models and running tests to confirm the findings. In simple terms, this project helps scientists track how the viruses have evolved, which is crucial for fighting current and future outbreaks.
