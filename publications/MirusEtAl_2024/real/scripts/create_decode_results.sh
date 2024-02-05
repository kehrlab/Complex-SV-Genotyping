#!/bin/bash

# Tim Mirus
# Process the raw results obtained by genotyping 20 variants in Icelandic genomes.
# Raw results are not provided here due to privacy concerns.

if [ ! $# -eq 1 ]; then
	echo "Correct call: ./create_decode_results.sh <results_file.txt>"
	exit 1
fi

results=$1

if [ -f decode_results/decode_genotype_results.tsv ]; then
	rm decode_results/decode_genotype_results.tsv
fi

if [ ! -d decode_results/plots ]; then
	mkdir decode_results/plots
fi

if [ ! -d decode_results/files ]; then 
	mkdir decode_results/files
fi

echo "Split results file into trio and genotype information"
echo -e "Variant\tSampleID\tGenotype\tMean_Quality\tLower_Bound\tUpper_Bound\tCertainty\tTotalReads\tOutliers\tAvgMapQ\tMinMapQ\tMaxMapQ\tFilterPass" > decode_results/decode_genotype_results.tsv

cat $results | sed '1,3440d' | sed 's/ /\t/g' >> decode_results/decode_genotype_results.tsv
head -n 3440 $results | sed 's/ /\t/g' > decode_results/decode_trios.tsv


echo "Merge decode calls and trios"
Rscript decode_results/scripts/gather_trio_calls.R

echo "Calculate transmission stats"
Rscript decode_results/scripts/create_transmission_stats.R

echo "Calculate HWE stats"
Rscript decode_results/scripts/create_hwe_stats.R
