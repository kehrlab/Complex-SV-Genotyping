#!/bin/bash

# Tim Mirus
# Genotype all simulated variants in all 7 data sets at all coverages using Paragraph (VCF input)

if [ ! $# -eq 2 ]; then
    echo "Usage: ./paragraphSimData.sh <bwaGenome> <nThreads>"
    exit 1
fi
paragraph=~/paragraph/bin/multigrmpy.py
converter=../scripts/VariantConverter/convert_to_vcf.py

genome=$1
nThreads=$2


for i in `ls data/simulatedSamples`; do
    if [ $i == "1" ]; then
	    continue
    fi

    echo "Repetition ${i}"
    d=data/simulatedSamples/${i}
    
    # create vcf file for paragraph
    echo "Create VCF file"
    varFile=${d}/paragraph_variants.vcf
    if [ ! -f $varFile ]; then
    	python3 $converter ${d}/simulated_variants.json $genome /dev/null 0 1 $varFile
    fi

    echo "Gather samples"
    # create list of samples for paragraph
    sampleFile=${d}/paragraphSamples.txt
    echo -e "id\tpath\tdepth\tread length" > $sampleFile
    for cov in `ls data/simulatedSamples/${i}/reads/alignment/`; do
        for f in `ls data/simulatedSamples/${i}/reads/alignment/${cov} | grep .bam | grep -v ".bai" | grep -v distributions`; do
            sampleName=${f%_aln_sorted.bam}_${i}_${cov}
            path="`pwd`/data/simulatedSamples/${i}/reads/alignment/${cov}/${f}"
            cov20=`samtools coverage ${path} -r chr20 -H`
            cov21=`samtools coverage ${path} -r chr21 -H`
            d20=`echo ${cov20} | cut -d " " -f7`
            d21=`echo ${cov21} | cut -d " " -f7`
            l20=`echo ${cov20} | cut -d " " -f3`
            l21=`echo ${cov21} | cut -d " " -f3`
            coverage=`python scripts/combineDepths.py ${d20} ${l20} ${d21} ${l21}`
            echo -e "${sampleName}\t${path}\t${coverage}\t150" >> ${sampleFile}
        done
    done

    echo "Genotype"
    # run paragraph
    python3 $paragraph -i $varFile -m $sampleFile -r $genome -o ${d}/paragraph_vcf_results -t $nThreads 2> logs/paragraph_vcf.log
done
