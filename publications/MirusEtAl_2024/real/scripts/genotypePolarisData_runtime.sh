#!/bin/bash

source /misc/rci-rg/ag_kehr/user/mit16436/anaconda3/etc/profile.d/conda.sh
conda activate cxsvPaper

if [ ! $# -eq 3 ]; then
	echo "Usage: ./genotypePolarisData <sampleFile> <genome> <nThreads>"
	exit 1
fi

sampleFile=$1
tempSampleFile=${sampleFile}_head.txt
cat ${sampleFile} | head -n 1 > ${tempSampleFile}
referenceFile=$2
nThreads=$3

outputDir=`pwd`/polarisResults/runtime

ggtyperBin=~/Git/Complex-SV-Genotyping/ggtyper
paragraphPy=~/paragraph/bin/multigrmpy.py
bayesDir=~/Git/bayesTyper_v1.5_linux_x86_64/bin

kmcPath=~/Git/KMC3.2.2.linux.x64/bin/kmc

bayestyperDecoy=/misc/rci-rg/ag_kehr/projects/2023_12_01_Genotyping_Sim/cxsv-paper-code/bayestyper_GRCh38_bundle_v1.3/GRCh38_decoy.fa
bayestyperCanon=/misc/rci-rg/ag_kehr/projects/2023_12_01_Genotyping_Sim/cxsv-paper-code/bayestyper_GRCh38_bundle_v1.3/GRCh38_canon.fa

variantFile=variant_descriptions/variants.json

vcfConverter="python ../VariantConverter/convert_to_vcf.py"

resourceTracker=../scripts/track_resources.sh

##########################

if [ ! -d data/polarisProfiles ]; then
	mkdir data/polarisProfiles
fi

if [ ! -d variant_descriptions/variantProfiles ]; then
	mkdir variant_descriptions/variantProfiles
fi

if [ ! -d ${outputDir} ]; then
	mkdir -p ${outputDir}
fi

##### GGTyper ############

echo "GGTyper"

if [ -f polarisResults/runtime/ggtyper_resources.tsv ]; then
	rm polarisResults/runtime/ggtyper_resources.tsv
fi

$resourceTracker ggtyper 1 polarisResults/runtime/ggtyper_resources.tsv &

# GGTyper v0.1.0 (changed profile-variants arguments)
$ggtyperBin profile-samples ${tempSampleFile} data/polarisProfiles.txt data/polarisProfiles -f -T $nThreads
$ggtyperBin profile-variants $variantFile variant_descriptions/variantProfiles.txt variant_descriptions/variantProfiles data/polarisProfiles.txt -T $nThreads
$ggtyperBin genotype variant_descriptions/variantProfiles.txt data/polarisProfiles.txt polarisResults/runtime/ggtyper -T $nThreads

killall track_resources.sh

## paragraph (vcf) ###

echo "Paragraph (VCF input)"
resDir=polarisResults/runtime/paragraph_results_vcf
if [ -f data/polarisSamples_full.txt ]; then
	rm data/polarisSamples_full.txt
fi

echo -e "id\tpath\tdepth\tread length" > data/polarisSamples_full.txt
for f in `cat ${tempSampleFile}`; do
	echo -e "`samtools samples $f | awk 'BEGIN{}{print $1}END{}'`\t$f\t30\t150" >> data/polarisSamples_full.txt # need to get idxdepth for the samples!
done

echo "Create VCF file..."
if [ ! -f variant_descriptions/real_variants.vcf ]; then
	conda activate py
	$vcfConverter $variantFile $referenceFile /dev/null 0 1 variant_descriptions/variants.vcf
	conda activate cxsvPaper
 	bcftools sort -o variant_descriptions/real_variants.vcf variant_descriptions/variants.vcf
fi

conda activate py

if [ -f polarisResults/runtime/paragraph_resources_vcf.tsv ]; then
	rm polarisResults/runtime/paragraph_resources_vcf.tsv
fi

$resourceTracker grmpy 1 polarisResults/runtime/paragraph_resources_vcf.tsv &
python $paragraphPy -i variant_descriptions/real_variants.vcf -m data/polarisSamples_full.txt -r $referenceFile -o $resDir -t $nThreads
killall track_resources.sh


### bayestyper ###
conda activate cxsvPaper

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
bcftools norm -f ${referenceFile} variant_descriptions/real_variants.vcf > ${variantPrefix}.vcf

##############################################################
# create bloom filters and bayestyper-compatible sample list #
##############################################################
echo "Bayestyper"
bayestyperFile=${outputDir}/bayestyperSamples.txt
if [ -f ${bayestyperFile} ]; then
	rm ${bayestyperFile}
fi
touch ${bayestyperFile}

SECONDS=0
for s in `cat ${tempSampleFile}`; do
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

	$resourceTracker kmc 1 polarisResults/runtime/kmc_resources_vcf.tsv &
	${kmcPath} -k55 -ci1 -t${nThreads} -v -fbam ${s} ${kmcPre} ${outputDir}/tmp
	killall track_resources.sh
	$resourceTracker bayesTyperTools 1 polarisResults/runtime/makeBloom_resources_vcf.tsv &
	${bayesDir}/bayesTyperTools makeBloom -k ${kmcPre} -p $nThreads
	killall track_resources.sh
done
echo -e "kmer_bloom\t${SECONDS}" > polarisResults/runtime/bayestyper_runtime.txt

if [ -d bayestyper_cluster_data ]; then
	rm -r bayestyper_cluster_data
fi
if [ -d bayestyper_unit_1 ]; then
	rm -r bayestyper_unit_1
fi

$resourceTracker bayesTyper 1 polarisResults/runtime/bayestyper_resources_vcf.tsv &
echo "Cluster variants..." # creates directories bayestyper_cluster_data and bayestyper_unit_1
${bayesDir}/bayesTyper cluster -v ${variantPrefix}.vcf -s ${bayestyperFile} -g ${bayestyperCanon} -d ${bayestyperDecoy} -p $nThreads
echo -e "clustering\t${SECONDS}" >> polarisResults/runtime/bayestyper_runtime.txt

echo "Genotype variants..." # needss to be called separately for every bayestyper_unit - directory (may be more than one for many variants)
${bayesDir}/bayesTyper genotype -v bayestyper_unit_1/variant_clusters.bin -c bayestyper_cluster_data/ -s ${bayestyperFile} -g ${bayestyperCanon} -d ${bayestyperDecoy} -o bayestyper_unit_1/bayestyper_genotypes_polaris -z -p $nThreads --min-genotype-posterior 0.1
killall track_resources.sh

${bayesDir}/bayesTyperTools filter -v bayestyper_unit_1/bayestyper_genotypes_polaris.vcf.gz -o polarisResults/runtime/bayestyper_unfiltered --min-homozygote-genotypes 0 --min-genotype-posterior 0 --min-number-of-kmers 0 --kmer-coverage-file ""
echo -e "genotyping\t${SECONDS}" >> polarisResults/runtime/bayestyper_runtime.txt
