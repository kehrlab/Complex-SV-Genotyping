# Tim Mirus
# Reformat BayesTyper results for easier joining and analysis.

library(tidyverse)

convert_gt <- function(gt) {
	if (gt == "0/0") {
		return("REF/REF")
	} else if (gt == "0/1" || gt == "1/0") {
		return("REF/VAR")
	} else if (gt == "1/1") {
		return ("VAR/VAR")
	}
	return("Undetermined")
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

gtData <- read.table("polarisResults/bayestyper_gtInfo.tsv", sep = "\t", header = TRUE)

gtInfo <- strsplit(gtData$Genotype, split = ":", fixed = TRUE)

genotypes <- sapply(gtInfo, function(x) convert_gt(x[1]))
qualities <- sapply(gtInfo, function(x) as.numeric(x[5]))

gtData$Genotype <- genotypes
gtData$Quality <- qualities

gtData <- as_tibble(gtData)

gtData %>% write_tsv("polarisResults/bayestyper_genotypes.tsv")
