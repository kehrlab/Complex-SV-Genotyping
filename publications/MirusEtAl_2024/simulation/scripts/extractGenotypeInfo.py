#!/usr/bin/python3

# Tim Mirus
# Extract relevant genotype information from JSON files produced by Paragraph (graph input)
# for all 7 data sets

import json
import sys
import os

if len(sys.argv) != 2:
    print("Usage: python extractGenotypeInfo.py <inputDir>")
    exit(1)

inputDir = sys.argv[1]
outPath = "results/paragraph_graph_gtInfo.tsv"

runs = [r for r in os.listdir(inputDir) if os.path.isdir(inputDir + '/' + r)]

with open(outPath, 'w') as outFile:
    header = 'Sample\tVariant\tRun\tGenotype\tQuality\tReads\n'
    for r in runs:
        print('Repetition: ', r)
        files = [f for f in os.listdir(inputDir + '/' + r + '/paragraph_graph_results') if os.path.isfile(inputDir + '/' + r + '/paragraph_graph_results/' + f) and '.json' in f]
        print('Extracting information from ', len(files), ' files')

        outFile.write(header)
        for f in files:
            path = inputDir + '/' + r + '/paragraph_graph_results/' + f
            variant = f.replace('_graph.json', '')
            with open(path) as inFile:
                variantData = json.load(inFile)
            for sample in variantData['samples']:
                sampleData = variantData['samples'][sample]
                gtInfo = sampleData['gt']
                genotype = gtInfo['GT']
                if genotype == '.':
                    quality = 0
                    nR = 0
                else:
                    quality = gtInfo['GQ']
                    nR = gtInfo['num_reads']
                line = sample + '\t' + variant + '\t' + r + '\t' + genotype + '\t' + str(quality) + '\t' + str(nR) + '\n'
                outFile.write(line)
