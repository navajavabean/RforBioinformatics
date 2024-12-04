# Load necessary libraries for sequence analysis, phylogenetics, and plotting
lapply(c("ape", "seqinr", "Biostrings", "tidyverse", "msa", "phangorn", "phytools", "stringr", "taxize"),
       library, character.only = TRUE)

# ########## Question 1 ##########

# Read genome and spike protein sequences from online fasta files
genomeSeq <- readDNAStringSet("https://raw.githubusercontent.com/idohatam/Biol-3315-files/main/SARS_MERS_coronavirus.raw_sequence.fasta")
proteinSeq <- readAAStringSet("https://raw.githubusercontent.com/idohatam/Biol-3315-files/main/spike.fa")

# Perform multiple sequence alignment (MSA) using the ClustalW method
genomeSeq_msa <- msa(genomeSeq, method = 'ClustalW')  # MSA for genome sequences
proteinSeq_msa <- msa(proteinSeq, method = 'ClustalW')  # MSA for protein sequences

# Convert the aligned sequences into 'phyDat' objects, which are compatible with phylogenetic analysis
genomeSeq_phy <- as.phyDat(msaConvert(genomeSeq_msa, "seqinr::alignment"), type = "DNA")
proteinSeq_phy <- as.phyDat(msaConvert(proteinSeq_msa, "seqinr::alignment"), type = "AA")

# Create distance matrices for DNA (genomes) and amino acid (proteins) sequences
genomeSeq_dml <- dist.ml(genomeSeq_phy, model = "JC69")  # JC69 model for DNA
proteinSeq_dml <- dist.ml(proteinSeq_phy, model = "JTT")  # JTT model for proteins

# Generate Neighbor Joining (NJ) trees from the distance matrices
genomeSeq_nj <- phangorn::NJ(genomeSeq_dml)  # NJ tree for genome sequences
proteinSeq_nj <- phangorn::NJ(proteinSeq_dml)  # NJ tree for protein sequences

# Root the trees using midpoint rooting (useful for unrooted trees)
genomeSeq_nj_mid <- midpoint(genomeSeq_nj)
proteinSeq_nj_mid <- midpoint(proteinSeq_nj)

# Set seed for reproducibility in the bootstrap analysis
set.seed(100)

# Perform bootstrap analysis for the NJ trees to assess tree stability
genomeSeq_nj_bs <- bootstrap.phyDat(genomeSeq_phy, FUN = function(x) NJ(dist.ml(x, model = "JC69")), bs=1000)
proteinSeq_nj_bs <- bootstrap.phyDat(proteinSeq_phy, FUN = function(x) NJ(dist.ml(x, model = "JTT")), bs=1000)

# Generate Maximum Likelihood (ML) trees using a more complex model (GTR for genomes, WAG for proteins)
genomeSeq_nj_pml <- pml(genomeSeq_nj, genomeSeq_phy, model="GTR", k=4, inv=.2)
proteinSeq_nj_pml <- pml(proteinSeq_nj, proteinSeq_phy, model="WAG", k=4)

# Optimize ML trees by adjusting various parameters for best fitting
genomeSeq_nj_pml_opt <- optim.pml(genomeSeq_nj_pml, optNni=TRUE, optBf=TRUE, optQ=TRUE, optInv=TRUE, optGamma=TRUE, optEdge=TRUE)
proteinSeq_nj_pml_opt <- optim.pml(proteinSeq_nj_pml, optNni=TRUE, optBf=TRUE, optQ=TRUE, optInv=TRUE, optGamma=TRUE, optEdge=TRUE)

# Perform bootstrap analysis for the ML trees to assess their stability
genomeSeq_pml_bs <- bootstrap.pml(genomeSeq_nj_pml_opt, bs=100, trees=TRUE, optNni=TRUE)
proteinSeq_pml_bs <- bootstrap.pml(proteinSeq_nj_pml_opt, bs=100, trees=TRUE, optNni=TRUE)

# Root the optimized ML trees using midpoint rooting
genomeSeq_pml_mid <- midpoint.root(genomeSeq_nj_pml_opt$tree)
proteinSeq_pml_mid <- midpoint.root(proteinSeq_nj_pml_opt$tree)

# Create tanglegrams to compare NJ and ML trees for genome and protein sequences
obj_genome <- cophylo(genomeSeq_nj_mid, genomeSeq_pml_mid)  # NJ vs ML for genomes
obj_protein <- cophylo(proteinSeq_nj_mid, proteinSeq_pml_mid)  # NJ vs ML for proteins
obj_comparison <- cophylo(genomeSeq_pml_mid, proteinSeq_pml_mid)  # ML genomes vs ML proteins

# ########## Question 2 ##########

# Save the NJ trees with bootstrap values for genome sequences as PDFs
pdf("genomeSeq_nj_tree.pdf")
plotBS(genomeSeq_nj_mid, genomeSeq_nj_bs, type = "phylogram", main = "NJ Tree of Coronavirus Genomes")
dev.off()

# Save the NJ trees with bootstrap values for protein sequences as PDFs
pdf("proteinSeq_nj_tree.pdf")
plotBS(proteinSeq_nj_mid, proteinSeq_nj_bs, type = "phylogram", main = "NJ Tree of Coronavirus Spike Proteins")
dev.off()

# Export tanglegrams comparing NJ and ML trees as PDFs
# NJ vs ML for genomes
pdf("tanglegram_genomes_NJ_ML.pdf", width = 10, height = 12)
par(mar = c(5, 5, 6, 5))  # Adjust plot margins
plot(obj_genome)
title("NJ vs ML Comparison for Coronavirus Genomes", line=-1)
dev.off()

# NJ vs ML for proteins
pdf("tanglegram_proteins_NJ_ML.pdf", width=10, height=12)
par(mar=c(5, 5, 6, 5))  # Adjust plot margins
plot(obj_protein)
title("NJ vs ML Comparison for Coronavirus Spike Proteins", line=-1)
dev.off()

# ML genomes vs ML proteins
pdf("tanglegram_genomes_ML_ML.pdf", width=10, height=12)
par(mar=c(5, 5, 6, 5))  
plot(obj_comparison)  
title("ML vs ML Comparison for Coronavirus Protein and Genome", line =-1)
dev.off()


# ########## Question 3 ##########

### A ###

# Observations on gene trees vs genome trees: 
# No noticeable differences between the gene trees and the genome trees tanglegrams. 
# This indicates that the dataset is reliable and the sequences accurately reflect the evolutionary relationships.

### B ###

# Spike protein tree based on nucleic acids:
# If the spike protein tree was made using nucleic acid sequences instead of amino acids, 
# we might see some differences. This is because DNA mutations and protein mutations occur differently,
# and DNA-based trees might show different evolutionary relationships compared to protein-based trees.

### C ###

# Model comparison for genome trees:
# The GTR model was used for analyzing the genome data, which is more complex and assumes that 
# DNA changes happen at different rates. This model provides a more realistic representation of DNA evolution 
# than the simpler JC69 model, which assumes all changes happen at the same rate.

### D ### 

# Model testing for the genome and protein sequences:
genome_model_test <- modelTest(genomeSeq_dml, model = "all")
spike_model_test <- modelTest(proteinSeq_dml, model = "all")

# Display the best model based on AIC (Akaike Information Criterion)
genome_model_test$Model[genome_model_test$AIC == min(genome_model_test$AIC)]
spike_model_test$Model[spike_model_test$AIC == min(spike_model_test$AIC)]

# JC-69 model did not perform as well as the GTR model. 
# The GTR model is more detailed and accounts for different types of DNA changes, 
# whereas the JC-69 model treats all changes equally, which is less realistic.
