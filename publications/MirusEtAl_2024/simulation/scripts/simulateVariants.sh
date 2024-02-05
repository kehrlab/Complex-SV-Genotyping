#!/bin/bash

# Tim Mirus
# Call svsim for <nRuns> times to randomly distribute 
# complex SVs on the two chromosomes (20, 21),
# write the variant descriptions to disk and create
# <nSamples> different samples with distinct haplotypes (sequences + genotype information)

if [ ! $# -eq 2 ]; then
    echo "Usage: ./simulateVariants.sh <nSamples> <nRuns>"
    exit 1
fi

nSamples=$1
nRun=$2

varConfFile=config/variantConfig.json
varTempDir=config/variantTemplates

pySim="python -m svsim"

if [ ! -d data/simulatedSamples ]; then
    mkdir data/simulatedSamples
fi

echo "Start variant simulation..."
for i in `seq $nRun`; do
    echo "Repetition: ${i}"
    outDir=data/simulatedSamples/${i}
    if [ ! -d $outDir ]; then
        mkdir $outDir
    fi

    $pySim --templateDir $varTempDir \
        --variantConfig $varConfFile \
        --nSamples $nSamples \
        --outputDir $outDir \
        --jsonOutput $outDir/simulated_variants.json \
        --statFile $outDir/variant_stats.txt \
        --chromosomeDir data/ \
        --chromosomes 20 21 \
        --nThreads 1 > logs/sv_sim.log \
        --exclude "chr21:10864561-12915808"
done
