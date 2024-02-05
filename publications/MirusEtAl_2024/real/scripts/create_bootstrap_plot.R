# Tim Mirus
#
# Select genotype call with certainty < 1 (Spink14 in some Polaris genomes) and plot posterior distribution
# of genotype quality 

library(tidyverse)
library(Cairo)
source("../scripts/R-scripts/theme.R")

genotypes <- read_tsv("fullResults/ggtyper_genotype_results.tsv") %>%
	filter(Variant == "SPINK14")


certainties <- genotypes %>% pull(Certainty)
minIdx <- which.min(certainties)
certainty <- min(certainties)

filepaths <- genotypes %>% pull(File)
filepath <- filepaths[minIdx]

samples <- genotypes %>% pull(Sample)
sample <- samples[minIdx]


cat("Sample: ", sample, "\tCertainty: ", certainty, "\n")

bootstrap_path <- paste0(filepath, "_distributions/SPINK14/SPINK14_bootstrap_data.txt")
bootstrap_data <- read_tsv(bootstrap_path)

# determine best and second best genotype
likelihoods <- bootstrap_data %>%
	select("REF/REF", "REF/VAR", "VAR/VAR") %>%
	rename("G0" = "REF/REF", "G1" = "REF/VAR", "G2" = "VAR/VAR") %>%
	pivot_longer(cols = c("G0", "G1", "G2"), names_to = "Genotype", values_to = "Likelihood") %>%
	group_by(Genotype) %>%
	summarize(MeanLikelihood = mean(Likelihood))

minIdx <- which.min(likelihoods %>% pull(MeanLikelihood))
maxIdx <- which.max(likelihoods %>% pull(MeanLikelihood))
secondIdx <- c(1,2,3)[!c(1, 2, 3) %in% c(minIdx, maxIdx)]

gts <- likelihoods %>% pull(Genotype)
bestGt <- gts[minIdx]
secondGt <- gts[secondIdx]

density_data <- density(bootstrap_data %>% pull(Quality))
density_df <- data.frame(x = density_data$x, y = density_data$y)

negFraction <- sum(bootstrap_data %>% pull(Quality) < 0) / nrow(bootstrap_data)
posFraction <- sum(bootstrap_data %>% pull(Quality) >= 0) / nrow(bootstrap_data)

# plot distribution of qualities with markings
p <- bootstrap_data %>%
	ggplot() +
	geom_density(aes(x = Quality)) +
	geom_vline(aes(xintercept = 0), linetype = "dashed") +
	geom_bar(data = density_df[which(density_df$x < 0),], aes(x = x, y = y), fill = "#f58080", alpha = 0.8, stat = "identity") +
	custom_theme_large +
	ylab("Density") +
	xlab(expression(Q["k,l"](S))) +
	annotate(geom = "text", label = paste0(bestGt, "\n", 100*round(posFraction, 3), "%"), x = 45, y = 0.0005, size = 10) +
	annotate(geom = "text", label = paste0(secondGt, "\n", 100*round(negFraction, 3), "%"), x = -45, y = 0.0005, size = 10) +
	theme(
		axis.title = element_text(size = 26),
		axis.text = element_text(size = 22)
	)
CairoPDF("plots/bootstrap_distribution.pdf", width = 16, height = 7, version = "1.5")
plot(p)
dev.off()
