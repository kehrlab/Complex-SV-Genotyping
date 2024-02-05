#!/usr/bin/python3

# Tim Mirus
# Extract relevant genotype information from JSON file produced by Paragraph (graph input)

import json
import sys
import os

if len(sys.argv) != 3:
    print("Usage: python extractGenotypeInfo.py <inputDir> <output.tsv>")
    exit(1)

inputDir = sys.argv[1]
outPath = sys.argv[2]
files = [f for f in os.listdir(inputDir) if os.path.isfile(inputDir + '/' + f) and '.json' in f]
print('Extracting information from ', len(files), ' files')

with open(outPath, 'w') as outFile:
    header = 'Sample\tVariant\tGenotype\tQuality\tReads\n'
    outFile.write(header)
    for f in files:
        path = inputDir + '/' + f
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
            line = sample + '\t' + variant + '\t' + genotype + '\t' + str(quality) + '\t' + str(nR) + '\n'
            outFile.write(line)
