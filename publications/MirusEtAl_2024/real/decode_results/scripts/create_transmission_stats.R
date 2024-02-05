# Load output of 'gather_trio_calls.R' and calculate transmission rate 
# per variant and overall, with and without filtering

library(tidyverse)

transmissions <- read_tsv("decode_results/files/transmission_info.tsv")

transmissions %>%
	filter(Transmission >= 0) %>%
	group_by(Variant) %>%
	summarize(
		N = n(), 
		TR_Unfiltered = mean(Transmission), 
		TR_Filtered = mean(Transmission[MinCertainty >= 0.9 & MapQ >= 40]), 
		.groups = "keep"
		) %>%
	pivot_longer(
		cols = c("TR_Unfiltered", "TR_Filtered"), 
		names_to = "Status", 
		values_to = "TR"
		) %>%
	ungroup() %>% 
	write_tsv("decode_results/files/transmission_per_variant.tsv")

transmissions <- transmissions %>%
	filter(Transmission >= 0)

nMax <- nrow(transmissions)
transmission_filtered <- data.frame()
for (q in seq(0, max(transmissions$MinQ), by = 2))
{
	temp <- transmissions %>%
		filter(MinQ >= q) %>%
		summarize(N = n(), Percentage = n() / nMax, TR = mean(Transmission))
	transmission_filtered <- rbind(
				       transmission_filtered,
				       data.frame(N = temp$N, Percentage = temp$Percentage, TR = temp$TR, Threshold = q, Filter = "Quality")
				       )
}
for (ce in sort(unique(transmissions$MinCertainty), decreasing = FALSE))
{
	temp <- transmissions %>%
		filter(MinCertainty >= ce) %>%
		summarize(N = n(), Percentage = n() / nMax, TR = mean(Transmission))

		transmission_filtered <- rbind(
				       transmission_filtered,
				       data.frame(N = temp$N, Percentage = temp$Percentage, TR = temp$TR, Threshold = ce, Filter = "Certainty")
				       )
}

as_tibble(transmission_filtered) %>% write_tsv("decode_results/files/transmission_filtered.tsv")