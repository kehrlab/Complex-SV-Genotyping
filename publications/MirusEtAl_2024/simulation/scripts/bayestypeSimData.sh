#!/bin/bash

# Tim Mirus
# Genotype all simulated variants in all 7 data sets at all coverages using BayesTyper

if [ ! $# -eq 1 ]; then
	echo "Usage: ./bayestypeSimData.sh <nThreads>"
	exit 1
fi

nThreads=$1

kmcPath=~/Git/KMC3.2.2.linux.x64/bin/kmc
bayesDir=~/Git/bayesTyper_v1.5_linux_x86_64/bin
bayestyperDecoy=/misc/rci-rg/ag_kehr/projects/2023_12_01_Genotyping_Sim/cxsv-paper-code/bayestyper_GRCh38_bundle_v1.3/GRCh38_decoy.fa
convertScript=/misc/rci-rg/ag_kehr/projects/2023_12_01_Genotyping_Sim/cxsv-paper-code/scripts/VariantConverter/convert_to_vcf.py
simRef=/misc/rci-rg/ag_kehr/projects/2023_12_01_Genotyping_Sim/cxsv-paper-code/simulation/data/simRef.fa
genome=/misc/rci-rg/ag_kehr/data/reference/hg38/BWA_index/genome.fa

# for each simulated data set (1-7)
wd=`pwd`
for i in `ls data/simulatedSamples/`; do
	echo "Repetition ${i}"
	
	d=${wd}/data/simulatedSamples/${i}
	sampleFile=${d}/bayestyperSamples.txt

	if [ -f ${sampleFile} ]; then
		rm ${sampleFile}
	fi
	touch ${sampleFile}

	# Calculate k-mer bloom filters for samples at all coverages
	# and create list of samples
	for cov in 10 20 25 30; do
		echo "Coverage: ${cov}"
		covDir=${d}/reads/alignment/${cov}
		cd ${covDir}

		echo "Create k-mer bloom filters"
		for s in 0 1 2 3 4; do 
			echo -e "\tSample ${s}..."
			echo -e "Sample_${s}_${cov}\tF\t${covDir}/Sample_${s}_kmc" >> ${sampleFile}
			if [ -f Sample_${s}_kmc.bloomData ]; then
				continue
			fi
			${kmcPath} -k55 -ci1 -fbam Sample_${s}_aln_sorted.bam Sample_${s}_kmc ./
			${bayesDir}/bayesTyperTools makeBloom -k Sample_${s}_kmc -p $nThreads
		done
	done

	# convert variant file (json) to VCF
	variantPrefix=${d}/bayestyper_variants
	if [ ! -f ${variantPrefix}.vcf ]; then
		echo "Convert json file to VCF..."
	
		python ${convertScript} ${d}/simulated_variants.json ${genome} /dev/null 0 1 ${variantPrefix}_temp.vcf

		bcftools sort -o ${variantPrefix}_sorted.vcf ${variantPrefix}_temp.vcf
		bcftools norm -f ${genome} ${variantPrefix}_sorted.vcf > ${variantPrefix}.vcf
		rm ${variantPrefix}_sorted.vcf
		rm ${variantPrefix}_temp.vcf
	fi

	cd ${d}
	if [ -d bayestyper_cluster_data ]; then
		rm -r bayestyper_cluster_data
	fi
	if [ -d bayestyper_unit_1 ]; then
		rm -r bayestyper_unit_1
	fi

	# genotype all samples belonging to the data set
	echo "Cluster variants..."
	${bayesDir}/bayesTyper cluster -v ${variantPrefix}.vcf -s ${sampleFile} -g ${simRef} -d ${bayestyperDecoy} -p $nThreads

	echo "Genotype variants..."
	${bayesDir}/bayesTyper genotype -v bayestyper_unit_1/variant_clusters.bin -c bayestyper_cluster_data/ -s ${sampleFile} -g ${simRef} -d ${bayestyperDecoy} -o bayestyper_unit_1/bayestyper_genotypes -z -p $nThreads
	
	cd $wd
done