#!/usr/bin/python

# written by Tim Mirus
# Takes a JSON file containing GGTyper variant descriptions and a BED file as input.
# Output: tsv file containing for each variant information, whether at least on breakpoint is contained in a region given in the BED file.

import json
from intervaltree import Interval, IntervalTree
import sys


def createIntervalTrees(bedFile, margin = 0, excludeRegions = None, verbose = False):
    chromosomeTrees = {}

    print("Build chromosome trees")
    n = 0
    with open(bedFile, 'r') as f:
        for line in f:
            n += 1
            content = line.strip().split()
            chrom = content[0]
            start = int(content[1])
            end = int(content[2])

            if not chrom in chromosomeTrees:
                chromosomeTrees[chrom] = IntervalTree()

            chromosomeTrees[chrom][start:end] = 'R'
    print("Added " + str(n) + " regions to Interval Trees")
    return chromosomeTrees


if __name__ == "__main__":

    if len(sys.argv) != 4:
        print("Usage: python find_rep_vars.py <variants.json> <repeats.bed> <outFile>")
        exit(1)

    variantFile = sys.argv[1]
    bedFile = sys.argv[2]
    outFile = sys.argv[3]

    repeatVariants = 0

    with open(variantFile) as jsonFile:
        variantData = json.load(jsonFile)

    with open(outFile, 'w') as f:
        chromosomeTrees = createIntervalTrees(bedFile)

        for v in variantData:
            inRepeat = False
            variant = variantData[v]

            for chrom in variant["VAR"]:
                chromosome = variant["VAR"][chrom]

                for j in chromosome:
                    junction = chromosome[j]

                    chr1 = junction["rNameLeft"]
                    chr2 = junction["rNameRight"]
                    xLeft = junction["xLeft"]
                    xRight = junction["xRight"]

                    if chr1 in chromosomeTrees:
                        if len(chromosomeTrees[chr1][xLeft]) != 0:
                            inRepeat = True
                            break
                    if chr2 in chromosomeTrees:
                        if len(chromosomeTrees[chr2][xRight]) != 0:
                            inRepeat = True
                            break
                
                if inRepeat:
                    break
            if inRepeat:
                repeatVariants += 1
            f.write(v + "\t" + str(inRepeat) + "\n")
    print("Variants in repeat regions: ", repeatVariants)
