# Complex-SV-Genotyping

## Description

This repository holds an implementation of a genotyping algorithm for complex structural variants using paired-end short-read data.
The aim is to enable genotyping of complex SVs of any structure.

## Installation

Download the source code and navigate to the corresponding directory:
```
git clone https://github.com/kehrlab/Complex-SV-Genotyping
cd Complex-SV-Genotyping
```

### Requirements
- g++
- Seqan v2.4
- htslib v1.17
- zlib
- openmp
- eigen
- nlohmann/json (included)

All requirements can be easily installed using conda:
```
conda create -c bioconda -c conda-forge -n genotyping cxx-compiler seqan=2.4 zlib intel-openmp htslib=1.17 eigen=3.4.0
conda activate genotyping
```

If the required libraries have not been installed system-wide, the following two lines in the Makefile need to be uncommented and adjusted to the correct paths:
```
# INCLUDE_PATH=/home/tim/.conda/envs/genotyping/include
# LIB_PATH=/home/tim/.conda/envs/genotyping/lib
```  
before calling `make`.  
Alternatively, the parameters can be supplied when calling make:
```
make CXXFLAGS="-I${INCLUDE_PATH} -L${LIB_PATH}"
```  

## Usage

The program is split in 3 commands:  
`ggtyper profile-samples`, `ggtyper profile-variants` and `ggtyper genotype`, which should be called in that order.  
To show the help, use `ggtyper -h/--help`:
```
GGTyper - Genotyping of complex structural variants
===================================================

SYNOPSIS
    ggtyper COMMAND [ARGUMENTS] [OPTIONS]

COMMANDS
    profile-samples       Create profiles of given bam files and write them to disk.
    profile-variants      Create profiles of given variants and write them to disk.
    genotype              Genotype given variants (specified by profiles) in all samples (specified by profiles).

VERSION
    GGTyper version: 0.0.1-9836abc
    Last update on 2024-02-26 15:48:21
    Contact: Tim Mirus (Tim.Mirus[at]ukr.de)

Try 'ggtyper COMMAND --help' for more information on each command.
```

## profile-samples

The `profile-samples` command requires as first argument the `.bam` file or a list of bam files as a `.txt` file.
The second argument is the location where the list of created sample profiles should be written. This must be a `.txt` file.
The third argument specifies the folder in which the created profiles are to be stored.  

A simple call would be
```
./ggtyper profile-samples [bamFiles.txt] sampleProfiles.txt ./ -T 10
```
where 
```
[bamFiles.txt]
filePath1.bam
filePath2.bam
...
```

*Note:* profile-samples expects one read group per BAM file. If more (or less) are present in the header, the command will fail. To force the use of such BAM files, use flag `-f`. All read groups within one file will be merged in that case.

## profile-variants

The `profile-variants` command is similar to sample-variants, but takes as first input either a `json` file of variant descriptions or a `txt` file containing one such json file per line.
Details on that format may be found below. The second and third argument once again denote the target location of the profile list and the target directory of the profiles, respectively.  
Variant profiles need to be created with library parameters (insert size range and read size) that match the samples to be genotyped, wich are determined from a list of sample profiles (fourth argument).

### Variant Description
(Complex) Structural Variants can be described by their novel junctions or break-ends. We use the JSON format to specify the novel junctions of the alternate alleles for each affected chromosome. It is imperative that these junctions are supplied in the order they are found in on the variant allele.
  
```json
{
    "[Variant1]" : {
        "[AlleleName]" : {
            "[chromosomeA]" : {
                "1" : {
                    "xLeft" : [position of left breakpoint],
                    "xRight" : [position of right breakpoint],
                    "rNameLeft" : [chromosome of left breakpoint],
                    "rNameRight" : [chromosome of right breakpoint],
                    "directionLeft" : ["left"/"right"; is the sequence to the left or right of the left breakpoint joined?],
                    "directionRight" : ["left"/"right"; is the sequence to the left or right of the right breakpoint joined?]
                },
                "2" : {
                    ...
                }
            }, 
            "[chromosomeB]" : ...
        }
    },
    "[Variant2]" : ...
}
```

Some examples of valid variant descriptions can be found [here](examples/). 

A simple call would be
```
./ggtyper profile-variants [variantDescriptions.json] variantProfiles.txt ./ sampleProfiles.txt
```

## genotype

The `genotype` command requires as arguments a list of variant profiles as returned by `profile-variants`, a list
of sample profiles as returned by `profile-samples` and the desired prefix/path of the output file.
The results will be written to `${PREFIX}_genotype_results.tsv`.  

We can genotype the variants and samples specified in variantProfiles.txt and sampleProfiles.txt for with
```
./ggtyper genotype variantProfiles.txt sampleProfiles.txt ./out -T 10 -e
```

### Other parameters
- The number of threads (OpenMP) can be specified with `-T`.
- Option `-d` writes created distributions and results into a directory created in the same directory as the bam file, .i.e., "[bamFile.bam]\_distributions"
- minQ specifies the minimum mapping quality required for reads to be considered. Default 0. 
- Estimation of Variant difficulty is enabled with `-e`, and returns for each genotype of a variant the number of reads required to genotype the variant with expected quality Q.
- Option `-P` allows specification of genotype priors. This currently works only for bi-allelic variants. There must be one variant per line. Each line has the format: [variantName] \t [p(G_0)] \t [p(G_1)] \t [p(g_2)]
- VCF output is optional and can be enable with `-V`. This split all variants into breakends. The feature is experimental and not all quality metrics are reported in the file.

## Output 
The output of `genotype` is a tab-separated file with the following columns:

- Variant: Name of the variant as given in the variant description.
- Sample: Name of the sample as read from the BAM header.
- File: Name of the BAM file.
- Genotype: Called Genotype, in the form [AlleleName1]/[AlleleName2]. 
Reference allele name is `REF`, variant allele name is given in variant description.
- Mean_Quality: Average genotype quality (determined by bootstrapping).
- Lower_Bound: Lower boundary of 95-% CI (determined by bootstrapping).
- Upper_Bound: Upper boundary of 95-% CI (determined by bootstrapping).
- Certainty: Genotype certainty (determined by bootstrapping).
- TotalReads: Number of read pairs the genotype call is based on.
- Outliers: Number of observed read pairs that could not be assigned to a group present in the profile.
- AvgMapQ: Average mapping quality of the used reads.
- MinMapQ: Lowest mapping quality among the used reads.
- MaxMapQ: Largest mapping quality among the used reads
- QualityPass: Flag indicating whether a certain quality (Certainty + AvgMapQ) was reached.
- Difficulty: Difficulty of this variant in this sample (derived from Kullback-Leibler divergence).

## Citation
Tim Mirus, Robert Lohmayer, Clementine Döhring, Bjarni V Halldórsson, Birte Kehr, GGTyper: genotyping complex structural variants using short-read sequencing data, Bioinformatics, Volume 40, Issue Supplement_2, September 2024, Pages ii11–ii19, https://doi.org/10.1093/bioinformatics/btae391

## Version and License
```
Last update: 21.03.2024
GGTyper version: 0.1.0
SeqAn version: 2.4
HTSLib version: 1.17
Author: Tim Mirus (Tim.Mirus[at]ukr.de)
```

GGTyper is distributed under the GPL-3.0. Consult the accompanying LICENSE file for more details.
