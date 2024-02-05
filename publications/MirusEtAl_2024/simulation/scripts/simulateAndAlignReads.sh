#!/bin/bash

# Tim Mirus
# Given fasta files containing 2 chromosome haplotypes each,
# simulate reads with ART (30x coverage), align them with BWA and sort the BAM files
# before subsampling them to different lower coverages (10x, 20x, 25x)

if [ ! $# -eq 3 ]; then
    echo "Usage: ./simulateAndAlignReads.sh <bwaGenome> <coverage> <nThreads>"
    exit 1
fi

genome=$1
coverage=$2
nThreads=$3

logDir=logs

cov=$((coverage + coverage))

for run in `ls data/simulatedSamples`; do
    echo "Repetition ${run}" 
    d=data/simulatedSamples/${run}
    mkdir -p ${d}/reads/alignment/${cov}
    for s in `ls ${d} | grep .fa`; do
        SECONDS=0
        sampleName=${s%.fa}
        echo -e "\t${s}"
        prefix=${d}/reads/${sampleName}_
        echo -e "\t\tSimulate reads..."
        art_illumina \
           -f $coverage \
	       -i ${d}/${s} \
	       -l 150 -m 420 -s 120 \
	       -ss HS25 \
	       -p -q -rs 1234 \
	       --noALN \
	       -o ${prefix} 2> /dev/null 1> /dev/null
        
	    echo -e "\t\tAlign reads..."
        bwa mem \
		-M \
		-t $nThreads \
        	-R "@RG\tID:${sampleName}_${run}.1\tSM:${sampleName}_${run}" \
		$genome \
		${prefix}1.fq ${prefix}2.fq \
		2> logs/bwa-${s}-c${coverage}-hap12.err \
        | samtools view -Sb \
        | samtools sort -@ $nThreads -o ${d}/reads/alignment/${cov}/${sampleName}_aln_sorted.bam

        echo -e "\t\tIndex bam file..."
        samtools index ${d}/reads/alignment/${cov}/${sampleName}_aln_sorted.bam
        echo -e "\tFinished after ${SECONDS}s"

        echo -e "\t\tSubsampling..."
        mkdir ${d}/reads/alignment/25
        mkdir ${d}/reads/alignment/20
        mkdir ${d}/reads/alignment/10
        echo -e "\t\t\t25x"
        samtools view -hb -s 1234.83333 -@ $nThreads -o ${d}/reads/alignment/25/${sampleName}_aln_sorted.bam ${d}/reads/alignment/${cov}/${sampleName}_aln_sorted.bam
        samtools index ${d}/reads/alignment/25/${sampleName}_aln_sorted.bam
        echo -e "\t\t\t20x"
        samtools view -hb -s 1234.66667 -@ $nThreads -o ${d}/reads/alignment/20/${sampleName}_aln_sorted.bam ${d}/reads/alignment/${cov}/${sampleName}_aln_sorted.bam
        samtools index ${d}/reads/alignment/20/${sampleName}_aln_sorted.bam
        echo -e "\t\t\t10x"
        samtools view -hb -s 1234.33333 -@ $nThreads -o ${d}/reads/alignment/10/${sampleName}_aln_sorted.bam ${d}/reads/alignment/${cov}/${sampleName}_aln_sorted.bam
        samtools index ${d}/reads/alignment/10/${sampleName}_aln_sorted.bam
    done
done
