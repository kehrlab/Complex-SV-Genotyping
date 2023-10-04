# Complex-SV-Genotyping

## Description

This project is an implementation of a genotyping algorithm for complex structural variants using paired-end short-read data.
The aim is to enable genoyping of complex SVs of any structure.

## Installation

### Requirements
- g++
- Seqan v2.4
- htslib v1.17
- zlib
- boost
- openmp
- eigen

All requirements can be easily installed using conda:
```
conda create -c bioconda -c conda-forge -n genotyping cxx-compiler seqan=2.4 zlib boost intel-openmp htslib=1.17 eigen=3.4.0
conda activate genotyping
```

If the required libraries have not been installed system-wide, the following two lines in the Makefile need to be uncommented and adjusted to the correct paths:
```
# INCLUDE_PATH=/home/tim/.conda/envs/genotyping/include
# LIB_PATH=/home/tim/.conda/envs/genotyping/lib
```  
or call make with addition parameters:
```
make CXXFLAGS="-I${INCLUDE_PATH} -L${LIB_PATH}"
```  
  
Afterwards, installation is straightforward:
```
git clone https://github.com/kehrlab/Complex-SV-Genotyping
cd Complex-SV-Genotyping
make
sudo make install   # install in /usr/bin, optional
```

## Usage

The program can be run with the command `genotype`. To show the help, use the command line argument `-h`:
```
genotype -h
```
This will show the following help page:

```text
SYNOPSIS

DESCRIPTION

OPTIONS
    -h, --help
          Display the help message.
    -i, --input-file INPUT_FILE
          BAM input file Valid filetype is: .bam.
    -F, --file-list INPUT_FILE
          text fiel containing list of bam files Valid filetype is: .txt.
    -o, --output-file OUTPUT_FILE
          output file for genotype information
    -C, --config-file INPUT_FILE
          json file containing novel junctions Valid filetypes are: .json and .JSON.
    -V, --output-vcf OUTPUT_FILE
          VCF output for variant breakends Valid filetype is: .vcf.
    -G, --reference-genome INPUT_FILE
          fasta file containing reference sequences
    -s, --sequence-directory OUTPUT_DIRECTORY
          directory that created allele sequences are written to
    -R, --sampling-regions INPUT_FILE
          text file containing all regions used for sampling default insert sizes,
          one region per line
    -T, --num-threads INTEGER
          max number of threads to use
    -v, --verbose
          print additional information
    -w, --whole-genome
          use all records instead of sampling to determine insert size distribution
    -d, --output-distributions
          write insert size distributions for all genotypes to text files
    -p, --profile
          profile performance
    -e, --estimate-difficulty INTEGER
          get an estimate of how many reads are needed to genotype variant with a
          certain quality
    -k, --no-split
          do not use split reads
    -l, --no-spanning
          do not use split reads
    -m, --no-standard
          do not use split reads
    -n, --no-insert-sizes
          combine each insert size distribution into one probability for the
          corresponding read pair class
    -minQ, --minMapQ INTEGER
          minimum mapping quality used to filter records (inclusive)
    -mode, --distribution-mode INTEGER
          internal storing of insert size distributions (0: one distribution for
          split and spanning, each; 1: split and spanning parted by orientation; 2:
          like 0, but location is considered; 3: insert size is only used for RF,
          otherwise like 0)
    -gc, --gc-correct
          perform gc bias correction based on observed (sampled) records ind the bam
          file
    -M, --load-to-memory
          load reference sequences of sampled chromosomes into memory
    -z, --coverage INTEGER
          average sampling coverage for determination of theoretical distribution;
          if not given no sampling is performed
    -x, --stats
          record stats per read group
    -L, --legacy-mode
          do not calculate profile matrices
VERSION
    Last update: 
    genotype version: 
    SeqAn version: 2.4.
```


Basic usage requires specification of a target bam file, the reference genome and a description of the variants to genotype:  
```
./genotype -i [sample.bam] -C [variantDescription.json] -G [reference.fa]
```

### Bam File
Using the option `-i` the user can specify a single bam file to be genotyped. However, in order to genotype multiple samples in parallel, the option `-F` allows specification of a text file containing a list of files to be processed, with one filepath per line, e.g.

```
./genotype -F bamFiles.txt ...
```

where bamFiles.txt is of the form
```
bamFile1.bam
bamFile2.bam
...
```
If not specified otherwise, it is assumed that the bam file contains alignment data for the whole genome. If that is not the case, the option `-R` should be used to specify the regions available (format chr:start-end). If only a small region is contained in the bam file, one should use the flag `-w` to specify that all reads in the file need to be used for calculation of the library distribution.

### Reference genome
Option `-G` is used to specify the filepath of the reference genome. This has to be a fasta file containing the sequence of the reference genome used during alignment of the read data. It is used among other things to create the sequence of the variant allele from the description and to calculate sample-specific correction factors for biases.

The current version of the human genome (hg38) can be found [here](https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz).

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

### Other parameters
- The number of threads (OpenMP) can be specified with `-T`.  
- GC Bias correction is enabled with `-gc`.   
- Prefix of the output files is specified with `-o`.  
- Option `-d` writes created distributions and results into a directory created in the same directory as the bam file, .i.e., "[bamFile.bam]\_distributions".  
- minQ specifies the minimum mapping quality required for reads to be considered. Default 0.  
- Importance sampling is enabled with `-z <sampling coverage>`. Can't be used in conjunction with GC bias correction.  
- Estimation of Variant difficulty is enabled with `-e <Q>`, and returns for each genotype of a variant the number of reads required to genotype the variant with expected quality Q.
- VCF output is enabled with the parameter `-V <filename>`.
