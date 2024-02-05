#!/bin/bash

# Tim Mirus
# Genotype all simulated variants in all 7 data sets at all coverages using Paragraph (Graph input)

if [ ! $# -eq 2 ]; then
    echo "Usage: ./paragraphSimData.sh <bwaGenome> <nThreads>"
    exit 1
fi

paragraph=~/paragraph/bin/grmpy
converter=../scripts/VariantConverter/convert_to_graphs.py

genome=$1
nThreads=$2

for i in `ls data/simulatedSamples`; do
    echo "Repetition ${i}"
    d=data/simulatedSamples/${i}
    
    # create graphs for paragraph
    echo "Create graph representations"

    graphDir=${d}/variantGraphs
    mkdir ${graphDir}
    python3 $converter ${d}/simulated_variants.json ${genome} ${graphDir}

    # create list of samples for paragraph
    sampleFile=${d}/paragraphSamples.txt
    
    echo "Genotype"
    # run paragraph
    for g in `ls ${graphDir}`; do
	echo -e "\t${g}"
	path=${graphDir}/${g}
	vName=${g%_graph.json}
    	$paragraph -g $path -m $sampleFile -r $genome -O ${d}/paragraph_graph_results -t $nThreads --log-level error --log-file logs/paragraph_graph_${g}_${i}.log
    done
done
