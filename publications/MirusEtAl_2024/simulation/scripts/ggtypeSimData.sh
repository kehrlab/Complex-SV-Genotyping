#!/bin/bash

# Tim Mirus
# Genotype all simulated variants in all 7 data sets at all coverages using GGTyper

if [ ! $# -eq 1 ]; then
    echo "Usage: ./ggtypeSimData.sh <nThreads>"
    exit 1
fi

ggtypeBin=/misc/rci-rg/ag_kehr/user/mit16436/Git/Complex-SV-Genotyping/ggtyper
regionFile=regions.txt
nThreads=$1

for i in `ls data/simulatedSamples`; do
	echo "Repetition ${i}"
	d=data/simulatedSamples/${i}
	
	for cov in 30 25 20 10; do
		covDir=${d}/reads/alignment/${cov}
		echo "Coverage ${cov}"
		mkdir ${d}/variantProfiles
    		
		if [ -f ${covDir}/samples.txt ]; then
			rm ${covDir}/samples.txt
		fi
		touch ${covDir}/samples.txt
		for f in `find ${covDir} -name "*.bam"`; do
			echo $f >> ${covDir}/samples.txt
		done

		mkdir ${covDir}/sampleProfiles

		$ggtypeBin profile-samples ${covDir}/samples.txt ${covDir}/sampleProfiles.txt ${covDir}/sampleProfiles -r $regionFile -T $nThreads
		$ggtypeBin profile-variants ${d}/simulated_variants.json ${d}/variantProfiles.txt ${d}/variantProfiles -S ${covDir}/sampleProfiles.txt -m 500 -T $nThreads
		$ggtypeBin genotype ${d}/variantProfiles.txt ${covDir}/sampleProfiles.txt ${covDir}/ggtyper -T $nThreads -e -d
		rm -r ${d}/variantProfiles
	done
done
