#!/bin/bash

# Tim Mirus
# Extract genotype information for all variants and samples from VCF file
# produced by Paragraph (VCF)

resDir=polarisResults/paragraph_results_vcf
bcftools view ${resDir}/genotypes.vcf.gz | grep \#CHROM | awk '{for(i=10;i<NF;i++) printf $i"\n"; print ""}' > ${resDir}/paragraph_sampleNames.txt

outFile=${resDir}/paragraph_gtInfo.tsv
if [ -f $outFile ]; then
	rm $outFile
fi
touch $outFile
echo -e "Sample\tVariant\tGenotype" > $outFile

for s in `cat ${resDir}/paragraph_sampleNames.txt`; do
	bcftools view ${resDir}/genotypes.vcf.gz -H -s $s | awk -v sample=$s '{printf sample"\t"$3"\t"$10"\n"}' >> $outFile
done

Rscript scripts/gather_paragraph_results_vcf.R
