#!/bin/bash

# Tim Mirus
# Extract genotype information for all variants and samples from VCF file
# produced by BayesTyper for all 7 data sets

bayesToolsBin=~/Git/bayesTyper_v1.5_linux_x86_64/bin/bayesTyperTools

outFile=results/bayestyper_gtInfo.tsv
if [ -f $outFile ]; then
	rm $outFile
fi
touch $outFile
echo -e "Sample\tVariant\tRun\tGenotype" > $outFile

wd=`pwd`
for i in `ls data/simulatedSamples`; do
	d=`pwd`/data/simulatedSamples/${i}

	resFile=${d}/bayestyper_unit_1/bayestyper_genotypes.vcf.gz
	tmpFile=bayestyper_sampleNames.tmp

	bcftools view ${resFile} | grep \#CHROM | awk '{for(i=10;i<=NF;i++) printf $i"\n"; print ""}' > ${tmpFile}

	for s in `cat ${tmpFile}`; do
		bcftools view ${resFile} -H -s $s | awk -v sample=$s -v run=$i '{printf sample"\t"$3"\t"run"\t"$10"\n"}' >> $outFile
	done
	rm ${tmpFile}
done
Rscript scripts/gather_bayestyper_results.R
