# Tim Mirus
# Reformat Paragraph (VCF) results for easier joining and analysis.

library(tidyverse)

resDir <- "polarisResults/paragraph_results_vcf"
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

gtData <- read.table(paste0(resDir, "/paragraph_gtInfo.tsv"), header = TRUE, sep = "\t")

gtInfo <- strsplit(gtData$Genotype, split = ":")

genotypes <- sapply(gtInfo, function(x) convert_gt(x[1]))
qualities <- sapply(gtInfo, function(x) extract_quality(x[length(x)]))

gtData$Genotype <- genotypes
gtData$Quality <- qualities
gtData <- as_tibble(gtData)

gtData %>% write_tsv(paste0(resDir, "/paragraph_genotypes.tsv"))
