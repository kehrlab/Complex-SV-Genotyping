#!/bin/bash
source /misc/rci-rg/ag_kehr/user/mit16436/anaconda3/etc/profile.d/conda.sh

# script to download VCF file from gnomAD-SV and create json file with all complex SVs
cd data
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/genome_sv/gnomad.v4.1.sv.sites.vcf.gz
gunzip gnomad.v4.1.sv.sites.vcf.gz

cat gnomad.v4.1.sv.sites.vcf | grep "^#" > gnomad.v4.1.sv.sites.filtered.vcf
cat gnomad.v4.1.sv.sites.vcf | grep SVTYPE=CPX | grep -v "<BND>" | grep PASS | grep -vE "^chr[X,Y]" >> gnomad.v4.1.sv.sites.filtered.vcf
cat gnomad.v4.1.sv.sites.vcf | grep "^#" > gnomad.v4.1.sv.sites.unfiltered.vcf
cat gnomad.v4.1.sv.sites.vcf | grep SVTYPE=CPX | grep -v "<BND>" | grep -vE "^chr[X,Y]" >> gnomad.v4.1.sv.sites.unfiltered.vcf


cd ..
mkdir data/gnomad_variants

conda activate py
python3 scripts/convertGnomadVCF.py data/gnomad.v4.1.sv.sites.filtered.vcf data/gnomad_variants/gnomad_cpx_filtered.json
python3 scripts/convertGnomadVCF.py data/gnomad.v4.1.sv.sites.unfiltered.vcf data/gnomad_variants/gnomad_cpx_unfiltered.json

#python3 scripts/find_variants_in_region.py data/gnomad_variants/gnomad_cpx_filtered.json data/hg38_centromeres.bed data/gnomad_variants/gnomad_centromere_info.tsv
#python3 scripts/find_variants_in_region.py data/gnomad_variants/gnomad_cpx_filtered.json data/hg38_problems.bed data/gnomad_variants/gnomad_problem_info.tsv
#python3 scripts/find_variants_in_region.py data/gnomad_variants/gnomad_cpx_filtered.json data/hg38_seg_dups.bed data/gnomad_variants/gnomad_seg_dups_info.tsv
#python3 scripts/find_variants_in_region.py data/gnomad_variants/gnomad_cpx_filtered.json data/hg38_repeats.bed data/gnomad_variants/gnomad_repeats_info.tsv

if [ -d polarisResults/ggtyper_gnomad_bih ]; then
	conda activate cxsvPaper
	Rscript scripts/eval_gnomad.R > polarisResults/gnomad_cpx_evaluation.txt
fi

