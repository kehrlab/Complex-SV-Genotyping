# Tim Mirus
# Reformat Paragraph (VCF input) results for easier joining and analysis.

library(tidyverse)

convert_gt <- function(gt) {
	if (gt == "0/0") {
		return("REF/REF")
	} else if (gt == "0/1" || gt == "1/0") {
		return("REF/VAR")
	} else if (gt == "1/1") {
		return ("VAR/VAR")
	}
	return("")
}

extract_quality <- function(qString) {
	qualities <- as.numeric(strsplit(qString, split = ",", fixed = TRUE)[[1]])
	if (any(is.na(qualities))) {
		return (0)
	}
	idx <- which.max(qualities)
	for (i in 1:length(qualities)) {
		if (qualities[i] < qualities[idx] && qualities[i] > 0) {
			idx <- i
		}
	}
	return(qualities[idx])
}

gtData <- read.table("results/paragraph_vcf_gtInfo.tsv", sep = "\t", header = TRUE)

gtInfo <- strsplit(gtData$Genotype, split = ":", fixed = TRUE)

genotypes <- sapply(gtInfo, function(x) convert_gt(x[1]))
qualities <- sapply(gtInfo, function(x) extract_quality(x[length(x)]))

sampleInfo <- strsplit(gtData$Sample, split = "_")
samples <- sapply(sampleInfo, function(x) paste0(x[1], "_", x[2]))
coverage <- sapply(sampleInfo, function(x) x[4])

gtData$Genotype <- genotypes
gtData$Quality <- qualities
gtData$Coverage <- coverage
gtData$Sample <- samples

gtData <- as_tibble(gtData)

gtData %>% write_tsv("results/paragraph_vcf_genotypes.tsv")
