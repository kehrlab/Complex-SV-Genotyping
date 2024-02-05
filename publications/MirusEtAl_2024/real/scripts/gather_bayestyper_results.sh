#!/bin/bash

# Tim Mirus
# Extract genotype information for all variants and samples from VCF file
# produced by BayesTyper

if [ ! $# -eq 1 ]; then
	echo "Usage: ./gather_bayestyper_results.sh <bayesTyperOutput.vcf>"
	exit 1
fi
resFile=$1

bayesToolsBin=~/Git/bayesTyper_v1.5_linux_x86_64/bin/bayesTyperTools

outFile=polarisResults/bayestyper_gtInfo.tsv
if [ -f $outFile ]; then
	rm $outFile
fi
touch $outFile
echo -e "Sample\tVariant\tGenotype" > $outFile

tmpFile=polaris_sampleNames.tmp

bcftools view ${resFile} | grep \#CHROM | awk '{for(i=10;i<=NF;i++) printf $i"\n"; print ""}' > ${tmpFile}

for s in `cat ${tmpFile}`; do
	bcftools view ${resFile} -H -s $s | awk -v sample=$s '{printf sample"\t"$3"\t"$10"\n"}' >> $outFile
done
rm ${tmpFile}

Rscript scripts/gather_bayestyper_results.R