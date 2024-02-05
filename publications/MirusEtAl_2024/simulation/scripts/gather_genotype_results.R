# Tim Mirus
# Gather genotype calls for GGTyper across all 7 data sets and all coverages
# and store them in a new data frame.
# Print some preliminary statistics (precision / recall).

library(tidyverse)

translateGenotypes <- function(gts)
{
	gts[gts == "REF/REF"] <- "HomRef"
	gts[gts == "REF/VAR"] <- "HetVar"
	gts[gts == "VAR/VAR"] <- "HomVar"
	return(gts)
}

translateActualGenotypes <- function(gts) 
{
	gts[gts == "0/0"] <- "HomRef"
	gts[gts == "0/1"] <- "HetVar"
	gts[gts == "1/1"] <- "HomVar"
	return (gts)
}

gtInfo <- data.frame()
sampleInfo <- data.frame()
for (i in 1:7) 
{
	for (cov in c(30, 25, 20, 10))
	{
		d <- paste0("data/simulatedSamples/", i, "/reads/alignment/", cov, "/")
		temp <- read.table(paste0(d, "ggtyper_genotype_results.tsv"), header = TRUE, sep = "\t")
		gtInfo <- rbind(
				gtInfo,
				cbind(temp, Coverage = cov, Run = i)
				)
	}
	for (s in 0:4) 
	{
		temp <- read.table(paste0("data/simulatedSamples/", i, "/Sample_", s, "_genotypes.tsv"), sep = "\t", header = TRUE)
		colnames(temp) <- c("Variant", "RealGenotype")
		temp <- cbind(temp, Sample = paste0("Sample_", s), Run = i)
		sampleInfo <- rbind(sampleInfo, temp)
	}
}
gtInfo <- as_tibble(gtInfo) %>%
	mutate(Sample = sapply(Sample, function(s){paste0(strsplit(s, split = "_")[[1]][1], "_", strsplit(s, split = "_")[[1]][2])}))
sampleInfo <- as_tibble(sampleInfo) %>%
	mutate(RealGenotype = translateActualGenotypes(RealGenotype))

gtInfo %>% write_tsv("results/ggtyper_genotypes.tsv")
sampleInfo %>% write_tsv("results/sample_info.tsv")

# Print some performance statistics

gtInfo <- gtInfo %>%
	inner_join(sampleInfo, c("Variant", "Sample", "Run")) %>%
	mutate(PredictedGenotype = translateGenotypes(Genotype)) %>%
	mutate(CorrectPrediction = (PredictedGenotype == RealGenotype))

options(pillar.sigfig = 5)
cat("\n\nFraction of variants with AvgMapQ >= 40:\n\n")
gtInfo %>%
	summarize(GoodFraction = mean(AvgMapQ >= 40)) %>%
	print()

cat("\n\nPrecision and Recall after filtering by AvgMapQ:\n")
gtInfo %>%
	filter(grepl(Variant, pattern = "csvT")) %>%
	group_by(Run) %>%
	summarize(PrecisionRun = mean(CorrectPrediction[AvgMapQ >= 40]), RecallRun = mean(AvgMapQ >= 40), .groups = "keep") %>%
	print()

cat("\n\nPrecision and Recall after filtering by Certainty:\n")
gtInfo %>%
	group_by(Run) %>%
	summarize(PrecisionRun = mean(CorrectPrediction[Certainty > 0.9]), RecallRun = mean(Certainty > 0.9)) %>%
	ungroup() %>%
	summarize(Precision = mean(PrecisionRun), P_SE = sd(PrecisionRun) / sqrt(length(unique(Run))), Recall = mean(RecallRun), R_SE = sd(RecallRun) / sqrt(length(unique(Run)))) %>%
	print()


gtInfo <-  gtInfo %>%
	mutate(VariantType = sapply(Variant, function(s) strsplit(s, split = "_")[[1]][1])) %>%
	mutate(PassFilter = Certainty >= 0.9 & AvgMapQ >= 40) %>%
	mutate(Complex = ifelse(grepl(x = VariantType, pattern = "csv"), "Complex", "Canonical"))

cat("\n\nComplex and canonic error rate and recall without mapQ filter:\n\n")
gtInfo %>% 
	filter(PassFilter) %>%
	filter(Coverage == 30) %>%
	group_by(Run, Complex) %>%
	summarize(Correct = mean(CorrectPrediction), .groups = "keep") %>%
	ungroup() %>%
	group_by(Complex) %>%
	summarize(MeanAccuracy = mean(Correct), seAccuracy = sd(Correct) / sqrt(7)) %>%
	print()

gtInfo %>%	
	filter(Coverage == 30) %>%
	group_by(Complex) %>%
	summarize(Recall = mean(PassFilter)) %>%	
	print()

cat("\n\nIndividual error rate after filtering:\n\n")
gtInfo %>% 
	filter(PassFilter) %>%
	filter(Coverage == 30) %>%
	group_by(Run, VariantType) %>%
	summarize(Correct = mean(CorrectPrediction), .groups = "keep") %>%
	ungroup() %>%
	group_by(VariantType) %>%
	summarize(MeanAccuracy = mean(Correct), seAccuracy = sd(Correct) / sqrt(7)) %>%
	print()

cat("\n\nIndividual error rate without filtering:\n\n")
gtInfo %>%	
	filter(Coverage == 30) %>%
	group_by(Run, VariantType) %>%
	summarize(Correct = mean(CorrectPrediction), .groups = "keep") %>%
	ungroup() %>%
	group_by(VariantType) %>%
	summarize(MeanAccuracy = mean(Correct), seAccuracy = sd(Correct) / sqrt(7)) %>%
	print()

cat("\n\nOverall error rate after filtering:\n\n")
gtInfo %>%	
	filter(PassFilter) %>%
	filter(Coverage == 30) %>%
	group_by(Run) %>%
	summarize(Correct = mean(CorrectPrediction), .groups = "keep") %>%	
	ungroup() %>%
	summarize(MeanAccuracy = mean(Correct), seAccuracy = sd(Correct) / sqrt(7)) %>%
	print()

cat("\n\nOverall error rate without filtering:\n\n")
gtInfo %>%	
	group_by(Run) %>%
	filter(Coverage == 30) %>%
	summarize(Correct = mean(CorrectPrediction), .groups = "keep") %>%	
	ungroup() %>%
	summarize(MeanAccuracy = mean(Correct), seAccuracy = sd(Correct) / sqrt(7)) %>%
	print()

cat("\n\nOverall recall:\n\n")
gtInfo %>%	
	filter(Coverage == 30) %>%
	group_by(Run) %>%
	summarize(RecallRun = mean(PassFilter)) %>%	
	ungroup() %>%
	summarize(Recall = mean(RecallRun), RecallSE = sd(RecallRun) / sqrt(length(unique(Run)))) %>%
	print()

cat("\n\nIndividual Recall:\n\n")
gtInfo %>%
	filter(Coverage == 30) %>%
	group_by(Run, VariantType) %>%
	summarize(Recall = mean(PassFilter), .groups = "keep") %>%
	group_by(VariantType) %>%
	summarize(MeanRecall = mean(Recall), SeRecall = sd(Recall) / sqrt(length(unique(Run))), .groups = "keep") %>%
	print()
