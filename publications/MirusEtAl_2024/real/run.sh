#!/bin/bash

# Tim Mirus
# execute the script required to reproduce the results in order

########################## set variables ###################################

genomePath=/misc/rci-rg/ag_kehr/data/reference/hg38/BWA_index/genome.fa
nThreads=12

# genotyping results for 20 variants using GGTyper on Icelandic data
# no access to samples, results provided by Bjarni V. Halldorsson
decodeResultsFile=decode_results/rerun231121.txt

# bayestyper results on polaris data (20 variants)
# results obtained by Birte Kehr, following this script:
# The actual analysis was performed as multiple jobs
# ./scripts/bayestypePolaris.sh data/fullPolarisFiles.txt results/bayestyperPolaris $nThreads
bayesTyperResultsFile=polarisResults/bayestyper_polaris_unfiltered.vcf

###########################################################################

echo "Genotype polaris data..."
./scripts/genotypePolarisData.sh $genomePath $nThreads
./scripts/gather_bayestyper_results.sh ${bayesTyperResultsFile}

echo "Create decode results (statistics and plots)"
./scripts/create_decode_results.sh ${decodeResultsFile}

echo "Create polaris and comparison results and plots..."
Rscript scripts/gather_polaris_info.R
Rscript scripts/create_polaris_stats.R
Rscript scripts/create_plots.R
Rscript scripts/create_bootstrap_plot.R
Rscript scripts/create_hwe_status_plot.R

echo "Create example plot of posterior quality distribution..."
sampleFile=data/fullSamples.txt
if [ -f $sampleFile ]; then
	rm $sampleFile
fi
ls /misc/rci-rg/ag_kehr/projects/Polaris_Diversity_Cohort/bam_exit/ERR195539?/*.bam > ${sampleFile}
./scripts/genotypeFullSamples.sh ${nThreads}

echo "Benchmark the runtime of the different algorithms on 10 samples..."
sampleFile=data/fullPolarisFiles.txt
if [ -f $sampleFile ]; then
	rm $sampleFile
fi
ls /misc/rci-rg/ag_kehr/projects/Polaris_Diversity_Cohort/bam/*.bam | head -n 10 > ${sampleFile}
./genotypePolarisData_runtime.sh ${sampleFile} ${genomePath} ${nThreads}