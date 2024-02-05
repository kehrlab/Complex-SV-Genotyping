#!/bin/bash

# Tim Mirus
# Extract specified chromosome sequences from reference genome for generation of modified chromosome sequences.

if [ $# -lt 2 ]; then
    echo "Correct Call:  ./getReferenceChromosomes.sh <bwa-genome> <chr1 chr2 ...>"
    exit 1
fi

genome=$1
simRef=data/simRef.fa

if [ -f ${simRef} ]; then
	rm ${simRef}
fi
touch ${simRef}

for chr in "$@"; do
    if [ $chr == $genome ]; then
        continue
    fi
    samtools faidx $genome chr${chr} > data/chr${chr}.fa
    echo data/chr${chr}.fa >> ${simRef}
    gzip data/chr$chr.fa
done
