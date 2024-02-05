#!/bin/bash

# Tim Mirus
# Genotype polaris samples (file locations supplied as first argument) with BayesTyper. 
# Note that the script does not consider BayesTyper's limit of 30 samples and
# therefore must be called repeatedly with different sample lists.

if [ ! $# -eq 3 ]; then
	echo "Usage: ./bayestypeSimData.sh <sampleList> <outputDir> <nThreads>"
	exit 1
fi

sampleFile=$1
outputDir=$2
nThreads=$3

# KMC binary
kmcPath=~/Git/KMC3.2.2.linux.x64/bin/kmc

# Directory containing BayesTyper and BayesTyperTools
bayesDir=~/Git/bayesTyper_v1.5_linux_x86_64/bin

# Bayestyper reference and decoy sequences; bundle needs to be downloaded separately (https://drive.google.com/file/d/1ioTjLFkfmvOMsXubJS5_rwpfajPv5G1Q/view?pli=1)
bayestyperDecoy=/misc/rci-rg/ag_kehr/projects/2023_12_01_Genotyping_Sim/cxsv-paper-code/bayestyper_GRCh38_bundle_v1.3/GRCh38_decoy.fa
bayestyperCanon=/misc/rci-rg/ag_kehr/projects/2023_12_01_Genotyping_Sim/cxsv-paper-code/bayestyper_GRCh38_bundle_v1.3/GRCh38_canon.fa

# python script for converting variant json file to VCF
convertScript=/misc/rci-rg/ag_kehr/projects/2023_12_01_Genotyping_Sim/cxsv-paper-code/scripts/VariantConverter/convert_to_vcf.py

# reference genome used for alignment
genome=/misc/rci-rg/ag_kehr/data/reference/hg38/BWA_index/genome.fa

if [ ! -d ${outputDir} ]; then
	mkdir -p ${outputDir}
fi
if [ ! -d ${outputDir}/tmp ]; then
	mkdir ${outputDir}/tmp
fi

#############################################
# create input variant file (VCF from json) #
############################################

variantPrefix=variant_descriptions/bayestyper_variants

echo "Convert json file to VCF..."
python ${convertScript} variant_descriptions/variants.json ${genome} /dev/null 0 1 ${variantPrefix}_temp.vcf	
bcftools sort -o ${variantPrefix}_sorted.vcf ${variantPrefix}_temp.vcf
bcftools norm -f ${genome} ${variantPrefix}_sorted.vcf > ${variantPrefix}.vcf
rm ${variantPrefix}_sorted.vcf
rm ${variantPrefix}_temp.vcf

##############################################################
# create bloom filters and bayestyper-compatible sample list #
##############################################################
bayestyperFile=${outputDir}/bayestyperSamples.txt
if [ -f ${bayestyperFile} ]; then
	rm ${bayestyperFile}
fi
touch ${bayestyperFile}

for s in `cat ${sampleFile}`; do
	sampleName=${s##*/}    # delete longest match of pattern from the beginning
	sampleName=${sampleName%.GRCh38.bam}     # delete shortest match of pattern from the end

	bloomDir=${outputDir}/bloom
	if [ ! -d ${bloomDir} ]; then
		mkdir ${bloomDir}
	fi

	# sex information is required by bayestyper; extract from polaris metadata by sample name
	sex=`cat data/polaris_sample_info.tsv | grep ${sampleName} | cut -f4`
	
	kmcPre=${bloomDir}/${sampleName}_kmc
	echo -e "${sampleName}\t${sex}\t${kmcPre}" >> $bayestyperFile

	${kmcPath} -k55 -ci1 -fbam ${s} ${kmcPre} ${outputDir}/tmp
	${bayesDir}/bayesTyperTools makeBloom -k ${kmcPre} -p $nThreads
done

#########################
# the actual genotyping #
#########################
if [ -d bayestyper_cluster_data ]; then
	rm -r bayestyper_cluster_data
fi
if [ -d bayestyper_unit_1 ]; then
	rm -r bayestyper_unit_1
fi

echo "Cluster variants..." # creates directories bayestyper_cluster_data and bayestyper_unit_1
${bayesDir}/bayesTyper cluster -v ${variantPrefix}.vcf -s ${bayestyperFile} -g ${bayestyperCanon} -d ${bayestyperDecoy} -p $nThreads

echo "Genotype variants..." # needss to be called separately for every bayestyper_unit - directory (may be more than one for many variants)
${bayesDir}/bayesTyper genotype -v bayestyper_unit_1/variant_clusters.bin -c bayestyper_cluster_data/ -s ${bayestyperFile} -g ${bayestyperCanon} -d ${bayestyperDecoy} -o bayestyper_unit_1/bayestyper_genotypes_polaris -z -p $nThreads --min-genotype-posterior 0.1

${bayesDir}/bayesTyperTools filter -v bayestyper_unit_1/bayestyper_genotypes_polaris.vcf -o results/bayestyper_unfiltered --min-homozygote-genotypes 0 --min-genotype-posterior 0 --min-number-of-kmers 0 --kmer-coverage-file ""