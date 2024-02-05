#!/bin/bash

# Tim Mirus
# Genotype variant Spink14 in some Polaris samples with GGTyper
# and call R-script to produce plot of posterior distribution of genotype quality

if [ ! $# -eq 2 ]; then
	echo "Usage: ./genotypeFullSamples.sh <sampleFile> <nThreads>"
	exit 1
fi

sampleFile=$1
nThreads=$2

ggtyperBin=~/Git/Complex-SV-Genotyping/ggtyper

variantFile=variant_descriptions/variants.json

##### ggtyper ############

echo "GGTyper"

if [ ! -d data/fullProfiles ]; then
	mkdir data/fullProfiles
fi

if [ ! -d variant_descriptions/variantProfiles ]; then
	mkdir variant_descriptions/variantProfiles
fi

if [ ! -d fullResults ]; then
	mkdir fullResults
fi

$ggtyperBin profile-samples $sampleFile data/fullProfiles.txt data/fullProfiles -f -T $nThreads
$ggtyperBin profile-variants $variantFile variant_descriptions/variantProfiles.txt variant_descriptions/variantProfiles -S data/fullProfiles.txt -m 500 -T $nThreads
$ggtyperBin genotype variant_descriptions/variantProfiles.txt data/fullProfiles.txt fullResults/ggtyper -e -d -T $nThreads

Rscript scripts/create_bootstrap_plot.R