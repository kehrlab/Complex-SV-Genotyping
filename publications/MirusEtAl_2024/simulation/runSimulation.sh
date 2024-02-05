#!/bin/bash

# Tim Mirus
# Simulate 7 different data sets containing 5 samples each.
# Simulate 902 variants per data set.
# Genotype all variants within each data set with different algorithms:
# GGTyper, BayesTyper, Paragraph

###################### Variables #######################################

genome=/misc/rci-rg/ag_kehr/data/reference/hg38/BWA_index/genome.fa
nThreads=20
coverage=15 # haplotype coverage; -> total coverage 30x

###################### Data simulation #################################

echo "Get reference sequences..."
./scripts/getReferenceChromosomes.sh $genome 20 21

echo "Simulate variant locations and structures..."
# requires installation of sv_simulator (https://git.uni-regensburg.de/mit16436/sv_simulator)
./scripts/simulateVariants.sh 5 7

echo "Build samples..."
./scripts/simulateAndAlignReads.sh $genome $coverage $nThreads

######################## Genotyping #####################################

echo "Genotype wth ggtyper"
./scripts/ggtypeSimData.sh $nThreads

echo "Genotype with paragraph (VCF)"
./scripts/paragraphSimData_vcf.sh $genome $nThreads

echo "Genotype with paragraph (graph input)"
./scripts/paragraphSimData_graph.sh $genome $nThreads

echo "Genotype with BayesTyper"
./scripts/bayestypeSimData.sh $nThreads

########################### Evaluation ##################################

./scripts/gather_paragraph_results_vcf.sh
./scripts/gather_bayestyper_results.sh
python scripts/extractGenotypeInfo.py data/simulatedSamples

Rscript scripts/gather_genotype_results.R
Rscript scripts/gather_paragraph_results_vcf.R
Rscript scripts/gather_paragraph_results_graph.R
Rscript scripts/gather_bayestyper_results.R


Rscript scripts/create_simulation_plots.R