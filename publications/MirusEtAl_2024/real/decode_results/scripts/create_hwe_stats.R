# load original decode data frame, 
# combine with trios, remove trios where child is a parent of another trio
# then calculate HWE statistic on remaining parents

library(tidyverse)
library(gridExtra)

hwe <- function(genotypes)
{
	g0obs <- sum(genotypes == "REF/REF")
	g1obs <- sum(genotypes == "REF/VAR")
	g2obs <- sum(genotypes == "VAR/VAR")

	N <- length(genotypes)
	a <- (2*g0obs + g1obs) / (2*N)
	b <- (2*g2obs + g1obs) / (2*N)

	g0exp <- a*a*N
	g2exp <- b*b*N
	g1exp <- 2*a*b*N

	chi <- (g0exp - g0obs)^2 / (g0exp) + (g1exp - g1obs)^2 / (g1exp) + (g2exp - g2obs)^2 / (g2exp)

	df <- data.frame(
		Genotype = c("G0", "G1", "G2", "G0", "G1", "G2"),
		Count = c(g0obs, g1obs, g2obs, g0exp, g1exp, g2exp), 
		Type = c(rep("Observed", 3), rep("Expected", 3))
		)
	p <- as_tibble(df) %>%
		ggplot(aes(x = Genotype, y = Count, group = Type, fill = Type)) +
		geom_bar(stat = "identity", position = "dodge") +
		theme_bw()
	return (list(p = p, data = df, chi = chi))
}

gts <- read_tsv("decode_results/decode_genotype_results.tsv") 

trios <- read.table("decode_results/decode_trios.tsv", sep = "\t", header = FALSE)
colnames(trios) <- c("child", "mother", "father")

children <- as_tibble(trios) %>%
	select(child)

gts <- gts %>%
	filter(!SampleID %in% children$child)

variants <- unique(gts$Variant)

variant_plots <- list()
chi_values <- c()
all_data <- data.frame()
for (v in variants) {
	temp <- gts %>%
		filter(Variant == v)
	hwe_result <- hwe(temp$Genotype)
	variant_plots <- append(variant_plots, list(hwe_result$p + ggtitle(v)))
	chi_values <- c(chi_values, hwe_result$chi)
	all_data <- rbind(
		all_data,
		cbind(hwe_result$data, Variant = v)
	)
}
names(chi_values) <- variants

variant_plots <- variant_plots[order(chi_values)]
chi_values <- sort(chi_values)
names(variant_plots) <- names(chi_values)

chi_squares <- data.frame(Chi = chi_values, Variant = names(chi_values), Significant = chi_values > 10.83)
as_tibble(chi_squares) %>% write_tsv("decode_results/files/chi_values.tsv")
as_tibble(all_data) %>% write_tsv("decode_results/files/hwe_counts.tsv")
