#!/bin/bash

# Tim Mirus
# Genotype 20 variants in 199 Polaris genomes with algorithms that accept
# BAM files that have been reduced to regions of interest (all but BayesTyper).
# These BAM files must be located in 'real/data/polaris'

if [ ! $# -eq 2 ]; then
	echo "Usage: ./genotypePolarisData <genome> <nThreads>"
	exit 1
fi

referenceFile=$1
nThreads=$2

ggtyperBin=~/Git/Complex-SV-Genotyping/ggtyper
paragraphPy=~/paragraph/bin/multigrmpy.py
paragraphBin=~/paragraph/bin/grmpy

variantFile=variant_descriptions/variants.json

vcfConverter="python ../scripts/VariantConverter/convert_to_vcf.py"
graphConverter="python ../scripts/VariantConverter/convert_to_graphs.py"

resourceTracker=../scripts/track_resources.sh

##### GGTyper ############

echo "GGTyper"

if [ ! -d data/polarisProfiles ]; then
	mkdir data/polarisProfiles
fi

if [ ! -d variant_descriptions/variantProfiles ]; then
	mkdir variant_descriptions/variantProfiles
fi

if [ ! -d polarisResults ]; then
	mkdir polarisResults
fi

if [ -f polarisResults/ggtyper_resources.tsv ]; then
	rm polarisResults/ggtyper_resources.tsv
fi


sampleFile=data/polarisFiles.txt
if [ -f $sampleFile ]; then
	rm $sampleFile
fi
touch $sampleFile
for f in `find data/polaris/ -name "*.bam"`; do
	echo $f >> $sampleFile
done

$resourceTracker ggtyper 1 polarisResults/ggtyper_resources.tsv &

$ggtyperBin profile-samples $sampleFile data/polarisProfiles.txt data/polarisProfiles -w -f -T $nThreads
$ggtyperBin profile-variants $variantFile variant_descriptions/variantProfiles.txt variant_descriptions/variantProfiles -S data/polarisProfiles.txt -m 500 -T $nThreads
$ggtyperBin genotype variant_descriptions/variantProfiles.txt data/polarisProfiles.txt polarisResults/ggtyper -e -d -T $nThreads

killall track_resources.sh

### Paragraph (vcf) ######

echo "Paragraph (VCF input)"
resDir=polarisResults/paragraph_results_vcf
if [ -f data/polarisSamples.txt ]; then
	rm data/polarisSamples.txt
fi

echo -e "id\tpath\tdepth\tread length" > data/polarisSamples.txt
for f in `find data/polaris -name "*.bam"`; do 
	echo -e "`samtools samples $f | awk 'BEGIN{}{print $1}END{}'`\t$f\t30\t150" >> data/polarisSamples.txt # need to get idxdepth for the samples!
done

echo "Create VCF file..."
if [ ! -f variant_descriptions/real_variants.vcf ]; then
	$vcfConverter $variantFile $referenceFile /dev/null 150 1 variant_descriptions/variants.vcf
 	bcftools sort -o variant_descriptions/real_variants.vcf variant_descriptions/variants.vcf
fi

$resourceTracker grmpy 1 paragraph_resources_vcf.tsv &
python $paragraphPy -i variant_descriptions/real_variants.vcf -m data/polarisSamples.txt -r $referenceFile -o $resDir -t $nThreads
killall track_resources.sh
mv paragraph_resources_vcf.tsv $resDir/

./scripts/gather_paragraph_results_vcf.sh

### Paragraph (graph) ####

echo "Paragraph (Graph input)"
resDir=polarisResults/paragraph_results_graph
if [ ! -d $resDir/logs ]; then
	mkdir -p $resDir/logs
fi

echo "Create graphs..."
$graphConverter $variantFile $referenceFile variant_descriptions/graphs/

$resourceTracker grmpy 1 paragraph_resources_graph.tsv &
for g in `ls variant_descriptions/graphs/`; do
	echo $g
	path=variant_descriptions/graphs/${g}
	vName=${g%_graph.json}
	$paragraphBin -g $path -m data/polarisSamples.txt -r $referenceFile -O $resDir -t $nThreads -M 1000 2> $resDir/logs/paragraph_log_${g}.txt
done
killall track_resources.sh
mv paragraph_resources_graph.tsv ${resDir}/
python scripts/extractGenotypeInfo.py ${resDir} ${resDir}/paragraph_graph_genotypes.tsv
