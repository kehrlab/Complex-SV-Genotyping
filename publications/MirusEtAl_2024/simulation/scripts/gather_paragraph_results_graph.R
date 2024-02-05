# Tim Mirus
# Reformat Paragraph (Graph input) results for easier joining and analysis.

library(tidyverse)

convert_gt <- function(gt) {
	if (gt == "REF/REF") {
		return("REF/REF")
	} else if (gt == "REF/ALT" || gt == "ALT/REF") {
		return("REF/VAR")
	} else if (gt == "ALT/ALT") {
		return ("VAR/VAR")
	}
	return(NA)
}

gtData <- read.table("results/paragraph_graph_gtInfo.tsv", header = TRUE, sep = "\t")

gtData$Genotype <- sapply(gtData$Genotype, function(x) {convert_gt(x)})
sampleInfo <- strsplit(gtData$Sample, split = "_")
coverage <- sapply(sampleInfo, function(x) x[4])
sampleName <- sapply(sampleInfo, function(x) {paste0(x[1], "_", x[2])})
gtData$Coverage <- coverage
gtData$Sample <- sampleName
gtData <- as_tibble(gtData)

gtData %>% write_tsv("results/paragraph_graph_genotypes.tsv")
