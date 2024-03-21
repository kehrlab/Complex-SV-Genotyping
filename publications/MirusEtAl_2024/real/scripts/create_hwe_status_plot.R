library(tidyverse)
library(Cairo)
source("../scripts/R-scripts/theme.R")
source("scripts/id_assignment.R")

get_status <- function(chi_value)
{
	if (is.na(chi_value)) {
		return("NO_ALT")
	} else if (as.numeric(chi_value) < 10.83) {
		return("HWE")
	}
	return("NO_HWE")
}

variant_order <- read.table("data/variant_order.tsv", sep = "\n")[[1]]

hwe_counts <- read_tsv("polarisResults/hwe_counts.tsv") %>%
	group_by(Algorithm, Variant) %>%
	summarize(Chi = unique(Chi), .groups = "keep") %>%
	mutate(HWE_STATUS = get_status(Chi)) %>%
	mutate(Size = 1 / Chi + 1) %>%
	mutate(CALL_STATUS = "OK") %>%
	bind_rows(tibble(Algorithm = "BayesTyper", Variant = "NA_variant_7", Chi = NA, HWE_STATUS = "NOT_FOUND", Size = NA, CALL_STATUS = "NO_CALLS"))

hwe_counts %>% write_tsv("polarisResults/hwe_counts_chi.tsv")
hwe_counts %>% write_tsv("../figure_data/Fig_8.tsv")

p <- hwe_counts %>%
	ggplot(aes(x = Variant, y = Algorithm, fill = CALL_STATUS)) +
	geom_tile(alpha = 0.8) +
	custom_theme_large +
	theme(
		axis.title = element_blank(),
		axis.text.x = element_text(angle = 30, hjust = 0.9, vjust = 0.9),
		legend.position = "bottom",
		#legend.title = element_blank(),
		panel.grid = element_blank(),
		#legend.background = element_rect(colour = "light grey"),
		legend.box.just = "left",
		legend.key.size = unit(2, "line"),
	) +
	scale_colour_manual(limits = c("HWE", "NO_HWE", "NO_ALT"), values = c("#12becc", "#f58080", "black"), name = "", na.value = "grey") +
	scale_fill_manual(limits = c("NO_CALLS"), values = c("grey"), na.value = "white", name = "") + 
	scale_x_discrete(limits = variant_order, labels = ids[variant_order]) +
	geom_point(aes(size = Chi + 1, col = HWE_STATUS, shape = HWE_STATUS)) +
	scale_shape_manual(limits = c("HWE", "NO_HWE", "NO_ALT"), values = c(20, 20, 4), name = "") +
	scale_size_area(max_size = 20, trans = "log", na.value = 10, breaks = log(c(6, 11, 51)), labels = c(5, 10, 50), name = expression(chi^2)) +
	guides(
	       size = guide_legend(title.vjust = 1, order = 1),
	       color = guide_legend(override.aes = list(size = 10), title = NULL, order = 2),
	       shape = guide_legend(title = NULL, order = 2), 
	       fill = guide_legend(title = NULL, order = 3)
	)
CairoPDF("plots/hwe_status.pdf", width = 22, height = 5, version = "1.5")
plot(p)
dev.off()
