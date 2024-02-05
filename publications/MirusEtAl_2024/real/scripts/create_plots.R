# Tim Mirus

library(tidyverse)
library(gridExtra)
library(Cairo)
source("../scripts/R-scripts/theme.R")

########################################################################
# Decode plot with Transmission Rate, MIER and Chi^2-Value per Variant #
########################################################################

transmission_decode <- read_tsv("decode_results/files/transmission_per_variant.tsv") %>%
	mutate(Filter = gsub(Status, pattern = "TR_", replacement = "")) %>%
	rename(Value = TR) %>%
  	mutate(Value = Value * 100) %>%
	mutate(Type = "TR") %>%
	select(-N, -Status)
decode_data <- read_tsv("decode_results/files/transmission_info.tsv")
chi_decode <- read_tsv("decode_results/files/chi_values.tsv") %>%
	mutate(Value = log10(Chi)) %>%
	mutate(Type = "log(Chi^2)") %>%
	select(-Significant, -Chi) %>%
	mutate(Filter = "Unfiltered")

variants <- chi_decode$Variant
variant_order <- variants

write.table(variant_order, file = "data/variant_order.tsv", sep = "\n", row.names = FALSE, col.names = FALSE)

mier_decode <- decode_data %>%
  group_by(Variant) %>%
  summarize(
    MIER_Unfiltered = mean(MendelianError) * 100, 
    MIER_Filtered = mean(MendelianError[MinCertainty >= 0.9 & MapQ >= 40]) * 100
    ) %>%
  pivot_longer(
    cols = c(MIER_Unfiltered, MIER_Filtered),
    names_to = "Filter",
    names_prefix = "MIER_",
    values_to = "Value"
    ) %>%
  mutate(Type = "MIER")

decode_stats <- bind_rows(
	transmission_decode,
	chi_decode, 
	mier_decode
)


type.labs <- c("TR", "MIER", expression(log[10](chi^2)))
line_data <- data.frame(Type = factor(c("log(Chi^2)", "TR"), levels = c("TR", "MIER", "log(Chi^2)"), labels = type.labs), y = c(log10(10.83), 50))
p <- decode_stats %>%
  mutate(Type = factor(Type, levels = c("TR", "MIER", "log(Chi^2)"), labels = type.labs)) %>%
	ggplot(aes(x = Variant, y = Value, group = Filter, fill = Filter)) +
	facet_wrap(
	  ~Type,
	  ncol = 1,
	  scales = "free_y", 
	  strip.position = "left",
	  labeller = label_parsed
	  ) +
	geom_bar(stat = "identity", position = "dodge2", width = 0.5) +
	scale_x_discrete(limits = variants) +
	custom_theme +
	theme(
	      axis.title.x = element_blank(), 
	      axis.title.y = element_blank(),
	      legend.title = element_blank(),
	      legend.position = "top",
	      legend.spacing.x = unit(1, "lines"),
	      strip.background = element_rect(fill = "white", colour = "grey"),
	      axis.text.x = element_blank(),
	      panel.spacing = unit(1, "lines")
	) +
	geom_hline(data = line_data, aes(yintercept = y), linetype = "dashed") +
	scale_fill_custom_d("lit_colours")
CairoPDF("plots/decode_stats.pdf", width = 21, height = 6, version = "1.5")
plot(p)
dev.off()

########################################################################
# Combined Transmission plot for all algorithms, data sets and filters #
########################################################################

decode_transmission <- read_tsv("decode_results/files/transmission_filtered.tsv") %>%
  mutate(Algorithm = "GGTyper") %>%
  mutate(Data = "Icelandic")
polaris_transmission <- read_tsv("polarisResults/transmission_rates.tsv") %>%
  mutate(Data = "Polaris")
transmissions <- bind_rows(
  decode_transmission,
  polaris_transmission
)
p_trans <- transmissions %>%
  ggplot(aes(
    x = 100 * Percentage, 
    y = TR * 100, 
    group = interaction(Algorithm, Filter, Data), 
    col = interaction(Algorithm, Filter, drop = TRUE, sep = " - "), 
    linetype = Data)
    ) +
  geom_line(linewidth = 3) +
  geom_hline(yintercept = 50, linetype = "dashed") +
  ylim(25, 75) +
  xlim(1, 100) +
  xlab("% of Trios") +
  ylab("Transmission Rate [%]") +
  custom_theme_large +
  scale_linetype_manual(breaks = c("Polaris", "Icelandic"), values = c("solid", "dashed")) +
  guides(linetype = guide_legend(order = 1), col = guide_legend("Algorithm & Filter", order = 2)) +
  theme(	
  	legend.position = c(0.98, 0.98),
  	legend.background = element_rect(colour = "light grey"),
	legend.key.width=unit(1.1,"inches"),
	legend.justification = c("right", "top"), 
	legend.box = "horizontal"
  ) +
  scale_colour_custom_d("lit_colours")

CairoPDF("plots/transmission_rate.pdf", width = 22, height = 9, version = "1.5")
plot(p_trans)
dev.off()

######################################################################################
# Observed and expected genotype counts for Polaris (GGTyper, Paragraph, BayesTyper) #
######################################################################################

hwe_polaris <- read_tsv("polarisResults/hwe_counts.tsv")
plots <- list()
for (v in variants) {
  p_hwe <- hwe_polaris %>%
    filter(Variant == v) %>%
    #filter(Algorithm == "GGTyper" | Algorithm == "BayesTyper" | Algorithm == "Paragraph(Graph)") %>%
    ggplot(aes(x = Genotype, y = Count, group = Type, fill = Type)) +
    facet_grid(cols = vars(Algorithm)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.8) +
    custom_theme +
    theme(
	  legend.position = "None", 
	  strip.text.x = element_text(size = 12), 
	  plot.title = element_text(size = 18),
	  axis.title = element_text(size = 14),
	  axis.text = element_text(size = 12)
	  ) +
    ggtitle(v) +
    scale_fill_custom_d("lit_colours")
  plots[[v]] <- p_hwe
}


CairoPDF("plots/polaris_hwe.pdf", width = 21, height = 30, version = "1.5")
grid.arrange(grobs = plots, nrow = 7, ncol = 3)
#plot(p_hwe_1)
dev.off()


#############################################################################
# Observed and expected genotype counts for GGTyper (Decode and Polaris)    #
#############################################################################

hwe_polaris <- read_tsv("polarisResults/hwe_counts.tsv") %>%
  filter(Algorithm == "GGTyper") %>%
  mutate(Data = "Polaris Diversity") %>%
  group_by(Variant, Type) %>%
  mutate(Fraction = Count / sum(Count))
hwe_decode <- read_tsv("decode_results/files/hwe_counts.tsv") %>%
  mutate(Data = "Icelandic") %>%
  group_by(Variant, Type) %>%
  mutate(Fraction = Count / sum(Count))
hwe <- bind_rows(
  hwe_polaris,
  hwe_decode
)
plots <- list()
for (v in variants) {
  p_hwe_gg <- hwe %>%
    filter(Variant == v) %>%
    ggplot(aes(x = Genotype, y = Fraction, group = Type, fill = Type)) +
    facet_grid(cols = vars(Data)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.5) +
    custom_theme +
    theme(legend.position = "None", plot.title = element_text(size = 18), strip.text.x = element_text(size = 16)) +
    ggtitle(v) +
    scale_fill_custom_d("lit_colours")
  plots[[v]] <- p_hwe_gg
}
CairoPDF("plots/ggtyper_hwe.pdf", width = 20, height = 30, version = "1.5")
grid.arrange(grobs = plots, nrow = 7, ncol = 3)
dev.off()

################################################################################
# Compare GGtyper's allele frequencies for decode against diversity populations#
################################################################################

hwe_polaris_eur <- read_tsv("polarisResults/hwe_counts_eur.tsv") %>%
  filter(Algorithm == "GGTyper") %>%
  group_by(Variant, Type) %>%
  mutate(Fraction = Count / sum(Count)) %>%
  rename(EurFraction = Fraction)
hwe_polaris_eas <- read_tsv("polarisResults/hwe_counts_eas.tsv") %>%
  filter(Algorithm == "GGTyper") %>%
  group_by(Variant, Type) %>%
  mutate(Fraction = Count / sum(Count)) %>%
  rename(EasFraction = Fraction)
hwe_polaris_afr <- read_tsv("polarisResults/hwe_counts_afr.tsv") %>%
  filter(Algorithm == "GGTyper") %>%
  group_by(Variant, Type) %>%
  mutate(Fraction = Count / sum(Count)) %>%
  rename(AfrFraction = Fraction)

hwe_fraction <- hwe_decode %>%
  inner_join(hwe_polaris_eur, by = c("Variant", "Genotype", "Type")) %>%
  inner_join(hwe_polaris_afr, by = c("Variant", "Genotype", "Type")) %>%
  inner_join(hwe_polaris_eas, by = c("Variant", "Genotype", "Type")) %>%
  filter(Type == "Observed") %>%
  mutate(eur_sq = (Fraction - EurFraction)^2, afr_sq = (Fraction - AfrFraction)^2, eas_sq = (Fraction - EasFraction)^2) %>%
  group_by(Variant) %>%
  summarize(AFR = sum(afr_sq), EAS = sum(eas_sq), EUR = sum(eur_sq)) %>%
  pivot_longer(cols = c(AFR, EAS, EUR), names_to = "Population", values_to = "SSQ")
  
overall_fraction <- hwe_fraction %>%
  group_by(Population) %>%
  summarize(SSQ = sum(SSQ)) %>%
  ggplot(aes(x = Population, y = SSQ)) +
  geom_col(position = "dodge")

p_ssq1 <- hwe_fraction %>%
  rename(Value = SSQ) %>%
  mutate(Type = "SSQ") %>%
  ggplot(aes(x = Variant, y = Value, group = Population, fill = Population)) +
  facet_wrap(~Type, strip.position = "left", ncol = 1, labeller = label_parsed) +
  geom_bar(stat = "identity", position = "dodge2", width = 0.5) +
  scale_x_discrete(limits = variant_order) +
  custom_theme +
  theme(
    axis.title.x = element_blank(), 
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 0.95, vjust = 0.95),
    legend.title = element_blank(),
    legend.position = "top",
    legend.spacing.x = unit(1, "lines"),
    strip.background = element_rect(fill = "white", colour = "grey"),
    strip.text.y = element_text(margin = margin(0.2, 0.28, 0.2, 0.28, "cm"))
    ) +
    scale_fill_manual(values = c("#0059ec", "#fcab79", "#6a00e5"))

CairoPDF("plots/ssq_decode_pop.pdf", width = 21, height = 3.5, version = "1.5")
plot(p_ssq1)
dev.off()
