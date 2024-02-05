from Bio import SeqIO
from Bio.Seq import Seq
import vcfpy
import json
import sys
import gc
import copy

##########################################################
## Function Definitions                                  #
##########################################################

def createChromosomeVarSequence(structure, chromosome, referenceIndex, margin, zeroBased = False):
    sequence = ""
    shift = 1
    if zeroBased:
        shift = 0

    keys = [int(k) for k in structure.keys()]
    keys.sort()
    # left 
    chrLeft = structure[str(keys[0])]["rNameLeft"]
    xLeft = structure[str(keys[0])]["xLeft"]
    reverse = False
    if structure[str(keys[0])]["directionLeft"] == "left":
        begin = xLeft - margin - shift
        end = xLeft - shift
    else:
        begin = xLeft - shift
        end = xLeft + margin - shift
        reverse = True
    tempSeq = referenceIndex[chrLeft]
    tempSeq = tempSeq[begin:(end+1)]
    if reverse:
        tempSeq = tempSeq.reverse_complement()
    sequence += tempSeq.seq

    # middle
    for i in range(0, len(keys) - 1):
        junctionLeft = structure[str(keys[i])]
        junctionRight = structure[str(keys[i+1])]
        begin = min(junctionLeft["xRight"], junctionRight["xLeft"]) - shift
        end = max(junctionLeft["xRight"], junctionRight["xLeft"]) - shift
        chrom = junctionLeft["rNameRight"]
        tempSeq = referenceIndex[chrom]
        tempSeq = tempSeq[begin:(end+1)]
        if junctionLeft["directionRight"] == "left":
            tempSeq = tempSeq.reverse_complement()
        sequence += tempSeq.seq

    # right
    chrRight = structure[str(keys[-1])]["rNameRight"]
    xRight = structure[str(keys[-1])]["xRight"]
    reverse = False
    if structure[str(keys[-1])]["directionRight"] == "right":
        begin = xRight - shift
        end = xRight + margin - shift
    else:
        begin = xRight - margin - shift
        end = xRight - shift
        reverse = True
    tempSeq = referenceIndex[chrLeft]
    tempSeq = tempSeq[begin:(end+1)]
    if reverse:
        tempSeq = tempSeq.reverse_complement()
    sequence += tempSeq.seq

    return sequence

def createChromosomeRefSequence(structure, chromosome, referenceIndex, margin, zeroBased = False):
    keys = [int(k) for k in structure.keys()]
    keys.sort()

    shift = 1
    if zeroBased:
        shift = 0

    begin = structure[str(keys[0])]["xLeft"] - margin - shift
    end = structure[str(keys[-1])]["xRight"] + margin - shift
    sequence = referenceIndex[chromosome][begin:(end+1)]
    return (begin + 1, sequence.seq)

def adjustJunctions(structure, chromosome):
    if len(structure) != 1:
        return structure
    k = list(structure.keys())[0]

    if (structure[k]["rNameLeft"] != chromosome or structure[k]["rNameRight"] != chromosome):
        return None
    if (structure[k]["xLeft"] <= structure[k]["xRight"]):
        return structure
    
    structure[str(int(k) + 1)] = copy.deepcopy(structure[k])
    structure[str(int(k) + 2)] = copy.deepcopy(structure[k])
    structure[k]["xLeft"] = structure[k]["xRight"] - 1
    structure[str(int(k) + 2)]["xRight"] = structure[str(int(k) + 1)]["xLeft"] + 1
    return structure

def writeSequences(sequenceDict, filename, id_ext):
    recs = [SeqIO.SeqRecord(Seq(sequenceDict[s]), id = s + '|' + id_ext) for s in sequenceDict]
    SeqIO.write(recs, filename, "fasta")

def createVcfHeader():
    lines = [
        vcfpy.HeaderLine("fileformat", "VCFv4.2"),
        vcfpy.HeaderLine("filedate", "20231411"),
        vcfpy.HeaderLine("source", "Hand"),
        vcfpy.HeaderLine("reference", "/home/tim/Git/Data/genome.fa")
    ]
    header = vcfpy.Header(lines, vcfpy.SamplesInfos([]))
    header.add_filter_line(vcfpy.OrderedDict([('ID', 'PASS'),('Description', 'All filters passed')]))
    return header

def writeVcf(header, records, filename):
    writer = vcfpy.Writer.from_path(filename, header)

##########################################################
## End Function Definitions
##########################################################

# parse parameters
if len(sys.argv) < 6 or len(sys.argv) > 7:
    print("Usage: python3 convert_variants.py <variantFile> <referenceFile> <outputDir> <margin> <zeroBased> [optional: <vcfFile>]")
    exit(1)


variantFile = sys.argv[1]
referenceFile = sys.argv[2]
outputDir = sys.argv[3]
margin = int(sys.argv[4])
zeroBased = (sys.argv[5] == "0")

vcfFile = ""
if len(sys.argv) == 7:
    vcfFile = sys.argv[6]

## create reference index
print("Use 0-based indices: ", zeroBased)
print("Create reference index...")
referenceIndex = SeqIO.index(referenceFile, "fasta")

# create VCF header
header = createVcfHeader()

# load the variant descriptions
with open(variantFile) as jsonFile:
    variantData = json.load(jsonFile)

# for each variant
records = []
contigs = set()
counter = 1
for v in variantData:
    percent = counter / len(variantData) * 100
    # get the chromosome descriptions
    variantDescription = variantData[v]["VAR"]
    variantSequences = dict()
    referenceSequences = dict()
    gc.collect()
    
    # create chromosome sequences (region) from junctions for each chromosome
    suitable = (len(variantDescription) == 1)
    for chrom in variantDescription:
        if not suitable:
            print("Skip due to more than one involved chromosome \t(", v, ")", sep = "", file=sys.stderr)
            break
        junctions = variantDescription[chrom]
        junctions = adjustJunctions(junctions, chrom)
        if junctions == None:
            print("Skip due to None \t(", v, ")", sep = "", file=sys.stderr)
            continue
        variantSequences[chrom] = createChromosomeVarSequence(junctions, chrom, referenceIndex, margin, zeroBased)
        refSeq = createChromosomeRefSequence(junctions, chrom, referenceIndex, margin, zeroBased)
        referenceSequences[chrom] = refSeq[1]
        if variantSequences[chrom][0] != referenceSequences[chrom][0]:
            print("ERROR: Padding bases do not match (beginning)\t(", v, ")", sep = "", file=sys.stderr)
        if variantSequences[chrom][len(variantSequences[chrom])-1] != referenceSequences[chrom][len(referenceSequences[chrom])-1]:
            print("ERROR: Padding bases do not match (end)\t(", v, ")", sep = "", file=sys.stderr)
        # add VCF record information
        if suitable and vcfFile != "": 
            if not chrom in contigs:
                contigs.add(chrom)
                header.add_contig_line(
                    vcfpy.OrderedDict([('ID', chrom), ('length', len(referenceIndex[chrom])), ('assembly', 'hg38'), ('species', 'Homo sapiens')])
                )
            record = dict()
            record["CHR"] = chrom 
            record["POS"] = refSeq[0]
            record["ID"] = v 
            record["REF"] = referenceSequences[chrom]
            record["ALT"] = variantSequences[chrom]
            records.append(record)
    
    if outputDir != "/dev/null":
        # write the new sequences to fasta files
        varFileName = outputDir + '/' + v + '_VAR.fa'
        refFileName = outputDir + '/' + v + '_REF.fa'
        writeSequences(variantSequences, varFileName, 'VAR')
        writeSequences(referenceSequences, refFileName, 'REF')
    
    print("Progress: ", int(percent), "%\t(",v,")               ", end = "\r", sep = "")
    counter += 1

# write to VCF
print("")
if vcfFile != "":
    print("Writing VCF...")
    writeVcf(header, records, vcfFile)
    with open(vcfFile, "a") as vcf:
        for record in records:
            rString = record["CHR"] + "\t" + str(record["POS"]) + "\t" + record["ID"] + "\t"
            rString += record["REF"] + "\t" + record["ALT"] + "\t.\tPASS\t.\n"
            vcf.write(str(rString))
