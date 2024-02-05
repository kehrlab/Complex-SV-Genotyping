#!/bin/bash

# Tim Mirus
# Extract genotype information for all variants and samples from VCF files
# produced by Paragraph for all 7 data sets

outFile=results/paragraph_vcf_gtInfo.tsv
if [ -f $outFile ]; then
	rm $outFile
fi
touch $outFile
echo -e "Sample\tVariant\tRun\tGenotype" > $outFile

for i in `ls data/simulatedSamples`; do
	resFile=data/simulatedSamples/${i}/paragraph_vcf_results/genotypes.vcf.gz
	
	bcftools view ${resFile} | grep \#CHROM | awk '{for(i=10;i<=NF;i++) printf $i"\n"; print ""}' > paragraph_sampleNames.tmp

	for s in `cat paragraph_sampleNames.tmp`; do
		bcftools view ${resFile} -H -s $s | awk -v sample=$s -v run=$i '{printf sample"\t"$3"\t"run"\t"$10"\n"}' >> $outFile
	done
	rm paragraph_sampleNames.tmp
done
Rscript scripts/gather_paragraph_results_vcf.R
