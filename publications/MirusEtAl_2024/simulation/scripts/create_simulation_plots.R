# Tim Mirus
# Evaluate genotype calls on simulated data and create plots.

library(tidyverse)
library(Cairo)
source("../scripts/R-scripts/theme.R")

dir.create("plots")

translateGenotypes <- function(gts)
{
        gts[gts == "REF/REF"] <- "HomRef"
        gts[gts == "REF/VAR"] <- "HetVar"
        gts[gts == "VAR/VAR"] <- "HomVar"
        return(gts)
}

translateGenotypes2 <- function(gts)
{
  gts[gts == "HomRef"] <- "G_0"
  gts[gts == "HetVar"] <- "G_1"
  gts[gts == "HomVar"] <- "G_2"
  return(gts)
}

###################### read and called genotypes ###############################
gq <- 20
ggtyper_gts <- read_tsv("results/ggtyper_genotypes.tsv") %>%
	mutate(Algorithm = "GGTyper") %>%
	mutate(PassFilter = !is.na(Genotype) & Certainty >= 0.9 & AvgMapQ >= 40) %>%
	rename(Quality = Mean_Quality)
paragraph_vcf_gts <- read_tsv("results/paragraph_vcf_genotypes.tsv") %>%
	mutate(Algorithm = "Paragraph(VCF)") %>%
	mutate(PassFilter = !is.na(Genotype) & Quality >= gq)
paragraph_graph_gts <- read_tsv("results/paragraph_graph_genotypes.tsv") %>%
	mutate(Algorithm = "Paragraph(Graph)") %>%
	mutate(PassFilter = !is.na(Genotype) & Quality >= gq)
bayestyper_gts <- read_tsv("results/bayestyper_genotypes.tsv") %>%
	mutate(Algorithm = "BayesTyper") %>%
	mutate(PassFilter = (Genotype != "Undetermined") & !is.na(Genotype) & Quality >= gq)

gts <- ggtyper_gts %>%
	bind_rows(paragraph_vcf_gts) %>%
	bind_rows(paragraph_graph_gts) %>%
	bind_rows(bayestyper_gts) %>%
	mutate(
		VariantType = sapply(Variant, function(v){strsplit(v, split = "_")[[1]][1]})
	)
sampleInfo <- read_tsv("results/sample_info.tsv")

variantTypes <- gts$VariantType
variantTypes[variantTypes == "TandemDuplication"] <- "Tandem Duplication"
gts$VariantType <- variantTypes

gts <- gts %>% 
	inner_join(sampleInfo, by = c("Variant", "Run", "Sample")) %>%
	mutate(PredictedGenotype = translateGenotypes(Genotype)) %>%
	mutate(CorrectPrediction = (PredictedGenotype == RealGenotype))
if (any(is.na(gts$CorrectPrediction))) {
	isCorrect <- gts$CorrectPrediction
	isCorrect[is.na(isCorrect)] <- FALSE
	gts$CorrectPrediction <- isCorrect
}

### create plot for overall error rate at different filter strengths ##########
overall_unfiltered <- gts %>%
	filter(Coverage == 30) %>%
	group_by(Run, Algorithm) %>%
	summarize(ErrorRate = 1 - mean(CorrectPrediction), .groups = "keep") %>%
	group_by(Algorithm) %>%
	summarize(
	  Error = mean(ErrorRate), 
	  SE = sd(ErrorRate) / sqrt(length(unique(Run))), 
	  .groups = "keep"
	  ) %>%
	mutate(Variant = "Overall") %>%
	mutate(Filter = "None")
overall_filtered_ggtyper <- gts %>%
	filter(Algorithm == "GGTyper") %>%
	filter(Coverage == 30) %>%
	filter(Certainty >= 0.9) %>%
	filter(AvgMapQ >= 40) %>%
	group_by(Run, Algorithm) %>%
	summarize(
	  ErrorRate = 1 - mean(CorrectPrediction), 
	  .groups = "keep"
	  ) %>%
	group_by(Algorithm) %>%
	summarize(
	  Error = mean(ErrorRate), 
	  SE = sd(ErrorRate) / sqrt(length(unique(Run))), 
	  .groups = "keep"
	  ) %>%
	mutate(Variant = "Overall") %>%
	mutate(Filter = "AvgMapQ + Certainty")
overall_filtered_other <- gts %>%
	filter(Coverage == 30) %>%
	filter(Quality >= 20) %>% # look into the paper, do they filter their results?
	group_by(Run, Algorithm) %>%
	summarize(ErrorRate = 1 - mean(CorrectPrediction), .groups = "keep") %>%
	group_by(Algorithm) %>%
	summarize(
	  Error = mean(ErrorRate), 
	  SE = sd(ErrorRate) / sqrt(length(unique(Run))), 
	  .groups = "keep"
	  ) %>%
	mutate(Variant = "Overall") %>%
	mutate(Filter = "Quality")

overall_performance <- bind_rows(
	overall_unfiltered,
	overall_filtered_ggtyper,
	overall_filtered_other
)
p1 <- overall_performance %>%
	ggplot(aes(x = Algorithm, y = Error * 100, group = factor(Filter, levels = c("None", "Quality", "AvgMapQ + Certainty")), fill = Filter)) +
	geom_bar(stat = "identity", position = "dodge") + 
	geom_errorbar(aes(ymin = 100 * (Error - SE), ymax = 100 * (Error + SE)), position = "dodge") +
	ylab("Error Rate [%]") +
	custom_theme +
	theme(axis.text.x = element_text(angle = 20, vjust = 0.9, hjust = 0.9)) +
	scale_fill_custom_d("lit_colours")

CairoPDF("plots/overall_performance.pdf", width = 16, height = 9, version = "1.5")
plot(p1)
dev.off()

##### create plot of performance per variant type for all algorithms ##########
individual_performance <- gts %>%
	filter(Coverage == 30) %>%
	group_by(Run, VariantType, Algorithm) %>%
	summarize(ErrorRate = 1 - mean(CorrectPrediction), .groups = "keep") %>%
	group_by(VariantType, Algorithm) %>%
	summarize(Error = mean(ErrorRate), SE = sd(ErrorRate) / sqrt(length(unique(Run))), .groups = "keep") %>%
	mutate(Filter = "None")
individual_performance_filtered_ggtyper <- gts %>%
	filter(Algorithm == "GGTyper") %>%
	filter(AvgMapQ >= 40 & Certainty >= 0.9) %>%
	filter(Coverage == 30) %>%
	group_by(Run, VariantType, Algorithm) %>%
	summarize(ErrorRate = 1 - mean(CorrectPrediction), .groups = "keep") %>%
	group_by(VariantType, Algorithm) %>%
	summarize(Error = mean(ErrorRate), SE = sd(ErrorRate) / sqrt(length(unique(Run))), .groups = "keep") %>%
	mutate(Filter = "AvgMapQ + Certainty")
individual_performance_filtered_other <- gts %>%
	filter(Quality >= 20) %>%
	filter(Coverage == 30) %>%
	group_by(Run, VariantType, Algorithm) %>%
	summarize(ErrorRate = 1 - mean(CorrectPrediction), .groups = "keep") %>%
	group_by(VariantType, Algorithm) %>%
	summarize(Error = mean(ErrorRate), SE = sd(ErrorRate) / sqrt(length(unique(Run))), .groups = "keep") %>%
	mutate(Filter = "Quality")

individual_performance <- individual_performance %>%
	bind_rows(
		  individual_performance_filtered_ggtyper, 
		  individual_performance_filtered_other
	) %>%
	bind_rows(overall_performance %>% rename(VariantType = Variant))
individual_performance %>% write_tsv("results/variant_performances.tsv")

individual_performance <- individual_performance %>%
  bind_rows(
    data.frame(VariantType = "Translocation", Algorithm = "Paragraph(VCF)", Error = 1, SE = NA, Filter = "None"),
    data.frame(VariantType = "Translocation", Algorithm = "Paragraph(Graph)", Error = 1, SE = NA, Filter = "None"),
    data.frame(VariantType = "csvT", Algorithm = "Paragraph(VCF)", Error = 1, SE = NA, Filter = "None"),
    data.frame(VariantType = "csvT", Algorithm = "Paragraph(Graph)", Error = 1, SE = NA, Filter = "None"),
    data.frame(VariantType = "Translocation", Algorithm = "BayesTyper", Error = 1, SE = NA, Filter = "None"),
    data.frame(VariantType = "csvT", Algorithm = "BayesTyper", Error = 1, SE = NA, Filter = "None")
  )

p1 <- individual_performance %>%
	ggplot(aes(x = factor(Filter, levels = c("None", "Quality", "AvgMapQ + Certainty")), y = 100 * Error, fill = Filter)) +
	facet_grid(rows = vars(VariantType), cols = vars(Algorithm), scales = "free_y") +
	geom_bar(stat = "identity", position = "dodge") +
	geom_errorbar(aes(ymin = 100 * (Error - SE), ymax = 100 * (Error + SE)), position = "dodge") +
	ylab("Error Rate [%]") +
	xlab("Filter") +
	custom_theme +
	theme(
	      axis.text.x = element_text(angle = 30, hjust = 0.9, vjust = 0.9), 
	      strip.text.x = element_text(size = 18)
	) +
	scale_fill_custom_d("lit_colours")

p2 <- individual_performance %>%
	filter(Algorithm == "GGTyper") %>%
	ggplot(aes(x = VariantType, y = 100 * Error, group = factor(Filter, levels = c("None", "Quality", "AvgMapQ + Certainty")), fill = factor(Filter, levels = c("None", "Quality", "AvgMapQ + Certainty")))) +
	geom_bar(stat = "identity", position = "dodge") +
	geom_errorbar(aes(ymin = 100 * (Error - SE), ymax = 100 * (Error + SE)), position = "dodge") +
	ylab("Error Rate [%]") +
	custom_theme +
	theme(
	      axis.title.x = element_blank(),
	      axis.text.x = element_text(angle = 30, hjust = 0.9, vjust = 0.9)
	      ) +
	scale_x_discrete(limits = c("Overall", "csvA", "csvB", "csvC", "csvD", "csvI", "csvT", "Deletion", "Inversion", "Duplication", "Tandem Duplication", "Translocation")) +
	geom_vline(xintercept = c(1.5, 7.5), linewidth = 2) +
	guides(fill = guide_legend(title = "Filter")) +
	scale_fill_custom_d("lit_colours")

p3 <- individual_performance %>%
  filter(Filter == "None") %>%
  ggplot(aes(x = VariantType, y = 100 * (1 - Error), group = Algorithm, fill = Algorithm)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_errorbar(aes(ymin = 100 * (1 - Error - SE), ymax = 100 * (1 - Error + SE)), position = "dodge", width = 0.7) +
  ylab("Precision [%]") +
  custom_theme_large +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 0.95, vjust = 0.95),
    legend.position = "top",
    legend.title = element_blank(),
    legend.spacing.x = unit(1, "lines")
  ) +
  scale_x_discrete(limits = c("Overall", "csvA", "csvB", "csvC", "csvD", "csvI", "csvT", "Deletion", "Inversion", "Duplication", "Tandem Duplication", "Translocation")) +
  geom_vline(xintercept = c(1.5, 7.5), linewidth = 2) +
  scale_fill_custom_d("lit_colours")

CairoPDF("plots/individual_performance.pdf", width = 20, height = 34, version = "1.5")
plot(p1)
dev.off()


CairoPDF("plots/individual_performance_ggtyper.pdf", width = 20, height = 10, version = "1.5")
plot(p2)
dev.off()

CairoPDF("plots/individual_performance_unfiltered.pdf", width = 16, height = 8, version = "1.5")
plot(p3)
dev.off()

############ quality / certainty distributions per variant type ###############
certainties <- gts %>%
	filter(Coverage == 30) %>%
	filter(Algorithm == "GGTyper") %>%
	ggplot(aes(x = VariantType, y = Certainty, col = VariantType)) +
	geom_violin()

qualities <- gts %>%
	filter(Coverage == 30) %>%
	filter(Algorithm != "GGTyper") %>%
	ggplot(aes(x = VariantType, y = Quality, col = VariantType)) +
	facet_grid(rows = vars(Algorithm)) +
	geom_violin()

CairoPDF("plots/certainty_per_type.pdf", width = 30, height = 10, version = "1.5")
plot(certainties)
plot(qualities)
dev.off()

###################### precision-recall curves at 30x ##########################

precision_recall <- tibble()
for (ce in sort(unique(gts$Certainty), decreasing = TRUE)) {
  pr <- gts %>%
  	filter(Algorithm == "GGTyper") %>%
  	filter(Coverage == 30) %>%
  	mutate(PassFilter = Certainty >= ce) %>%
  	group_by(Run, Algorithm) %>%
  	summarize(Recall = sum(PassFilter) / n(), Precision = mean(CorrectPrediction[PassFilter]), .groups = "keep") %>%
  	group_by(Algorithm) %>%
  	summarize(MeanRecall = mean(Recall), SeRecall = sd(Recall) / sqrt(length(unique(Run))), MeanPrecision = mean(Precision), SePrecision = sd(Precision) / sqrt(length(unique(Run)))) %>%
  	mutate(Filter = "Certainty") %>%
  	mutate(Threshold = ce) %>%
  	rename(Recall = MeanRecall, Precision = MeanPrecision)
  precision_recall <- bind_rows(precision_recall, pr)
  pr <- gts %>%
  	filter(Algorithm == "GGTyper") %>%
  	filter(Coverage == 30) %>%
  	mutate(PassFilter = Certainty >= ce & AvgMapQ >= 40) %>%
  	group_by(Run, Algorithm) %>%
  	summarize(Recall = sum(PassFilter) / n(), Precision = mean(CorrectPrediction[PassFilter]), .groups = "keep") %>%
  	group_by(Algorithm) %>%
  	summarize(MeanRecall = mean(Recall), SeRecall = sd(Recall) / sqrt(length(unique(Run))), MeanPrecision = mean(Precision), SePrecision = sd(Precision) / sqrt(length(unique(Run)))) %>%
  	mutate(Filter = "AvgMapQ + Certainty") %>%
  	mutate(Threshold = ce) %>%
  	rename(Recall = MeanRecall, Precision = MeanPrecision)
  	precision_recall <- bind_rows(precision_recall, pr)
}
for (q in seq(0, max(floor(max(gts %>% pull(Quality)))), by = 2)) {
  pr <- gts %>%
  	filter(Coverage == 30) %>%
  	mutate(PassFilter = Quality >= q) %>%
  	group_by(Run, Algorithm) %>%
  	summarize(Recall = sum(PassFilter) / n(), Precision = mean(CorrectPrediction[PassFilter]), .groups = "keep") %>%
  	group_by(Algorithm) %>%
  	summarize(MeanRecall = mean(Recall), SeRecall = sd(Recall) / sqrt(length(unique(Run))), MeanPrecision = mean(Precision), SePrecision = sd(Precision) / sqrt(length(unique(Run)))) %>%
  	mutate(Filter = "Quality") %>%
  	mutate(Threshold = q) %>%
  	rename(Recall = MeanRecall, Precision = MeanPrecision)
  precision_recall <- bind_rows(precision_recall, pr)
}

p4 <- precision_recall %>%
	ggplot(aes(x = Recall, y = Precision, group = interaction(Algorithm, Filter), col = Algorithm, linetype = Filter)) +
	geom_path() +
	scale_x_reverse(limits = c(1, 0.5)) +
	custom_theme +
	scale_colour_custom_d("lit_colours")
p5 <- precision_recall %>%
	filter(Algorithm == "GGTyper") %>%
	ggplot(aes(x = Recall * 100, y = Precision * 100, group = Filter, col = Filter)) +
	geom_path(linewidth = 2) +
	xlim(70, 100) +
	custom_theme_large +
	scale_colour_custom_d("lit_colours") +
	theme(
		legend.position = c(0.01, 0.01), 
    		legend.background = element_rect(colour = "light grey"),
		legend.justification = c("left", "bottom")
	) + 
	ylab("Precision [%]") +
	xlab("Recall [%]")

CairoPDF("plots/pr_filter.pdf", width = 13, height = 5, version = "1.5")
plot(p4)
plot(p5)
dev.off()

############# precision and recall vs coverage for all algorithms ##############

gts <- gts %>%
	mutate(
	  Complex = ifelse(
	    grepl(x = VariantType, pattern = "csv"), 
	    "Complex", 
	    "Canonical"
	    )
	  )

pr <- gts %>%
	group_by(Coverage, Run, Complex, Algorithm) %>%
	summarize(
	  Recall = sum(PassFilter) / n(), 
	  Precision = mean(CorrectPrediction[PassFilter]), 
	  .groups = "keep"
	  ) %>%
	ungroup() %>%
	group_by(Coverage, Complex, Algorithm) %>%
	summarize(
	  MeanRecall = mean(Recall), 
	  MeanPrecision = mean(Precision), 
	  SeRecall = sd(Recall) / sqrt(length(unique(Run))), 
	  SePrecision = sd(Precision) / sqrt(length(unique(Run))), 
	  .groups = "keep"
	  ) %>%
	rename(Recall = MeanRecall, Precision = MeanPrecision) %>%
	pivot_longer(
	  cols = c("Recall", "Precision"), 
	  values_to = "Value", 
	  names_to = "Type"
	  )
errors <- pr %>% pull(SeRecall)
for (i in 1:nrow(pr))
{
	if (pr$Type[i] == "Precision")
	{
		errors[i] = pr$SePrecision[i]
	}
}
pr$Error <- errors


p6 <- pr %>%
	rename(Metric = Type) %>%
	rename(Type = Complex) %>%
	ggplot(aes(x = Coverage, y = Value, group = interaction(Metric, Type), col = Metric, linetype = Type)) +
	facet_grid(rows = vars(Algorithm)) +
	geom_line() +
	geom_errorbar(aes(ymin = Value - Error, ymax = Value + Error)) +
	custom_theme +
	theme(axis.title.y = element_blank()) +
	scale_colour_custom_d("lit_colours") 

p7 <- pr %>%
	rename(Metric = Type) %>%
	rename(Type = Complex) %>%
	filter(Algorithm == "GGTyper") %>%
	ggplot(aes(x = Coverage, y = Value, group = interaction(Metric, Type), col = Metric, linetype = Type)) +
	geom_line(linewidth = 2) +
	custom_theme +
	theme(axis.title.y = element_blank()) +
	scale_colour_custom_d("lit_colours")

p8_complete <- pr %>%
  rename(Metric = Type) %>%
  rename(Type = Complex) %>%
  ggplot(aes(x = Coverage, y = Value, group = interaction(Algorithm, Type), col = Algorithm)) +
  facet_grid(cols = vars(Metric), scales = "free") +
  geom_line(aes(linetype = Type), linewidth = 3) +
  geom_errorbar(aes(ymin = Value - Error, ymax = Value + Error), width = 0.25) +
  custom_theme_large +
  scale_linetype_manual(
    breaks = c("Complex", "Canonical"), 
    values = c("solid", "dashed")
    ) +
  theme(
    axis.title.y = element_blank(), 
    legend.position = c(0.99, 0.01),
    legend.justification = c("right", "bottom"), 
    legend.background = element_rect(colour = "light grey"),
    legend.box = "horizontal"
    ) +
  scale_colour_custom_d("lit_colours")

p8_select <- pr %>%
  rename(Metric = Type) %>%
  rename(Type = Complex) %>%
  filter(Algorithm == "GGTyper" | Algorithm == "BayesTyper") %>%
  ggplot(aes(x = Coverage, y = Value, group = interaction(Algorithm, Type), col = Algorithm)) +
  facet_grid(cols = vars(Metric), scales = "free") +
  geom_line(aes(linetype = Type), linewidth = 2.5) +
  geom_errorbar(aes(ymin = Value - Error, ymax = Value + Error), width = 0.25) +
  custom_theme_large +
  scale_linetype_manual(
    breaks = c("Complex", "Canonical"), 
    values = c("solid", "dashed")
    ) +
  theme(
    axis.title.y = element_blank(), 
    legend.position = c(0.01, 0.01),
    legend.justification = c("left", "bottom"), 
    legend.background = element_rect(colour = "light grey"),
    legend.box = "horizontal"
    ) +
  scale_colour_custom_d("lit_colours")

p8_precision <- pr %>%
	rename(Metric = Type) %>%
	rename(Type = Complex) %>%
	filter(Metric == "Precision") %>%
	filter(Algorithm == "GGTyper" | Algorithm == "BayesTyper") %>%
	ggplot(aes(x = Coverage, y = Value, group = interaction(Algorithm, Type), col = Algorithm)) +
	geom_line(aes(linetype = Type), linewidth = 2.5) +
	geom_errorbar(aes(ymin = Value - Error, ymax = Value + Error), width = 0.25) +
	custom_theme_large +
	scale_linetype_manual(breaks = c("Complex", "Canonical"), values = c("solid", "dashed")) +
	ylab("Precision") +
	scale_colour_custom_d("lit_colours")

CairoPDF("plots/pr_complex.pdf", width = 16, height = 21, version = "1.5")
plot(p6)
dev.off()

CairoPDF("plots/pr_complex_ggtyper.pdf", width = 16, height = 9, version = "1.5")
plot(p7)
dev.off()

CairoPDF("plots/pr_complex_wide.pdf", width = 21, height = 7, version = "1.5")
plot(p8_complete)
plot(p8_select)
dev.off()

CairoPDF("plots/pr_precision.pdf", width = 16, height = 9, version = "1.5")
plot(p8_precision)
dev.off()

pr %>%
	rename(Metric = Type) %>%
	rename(Type = Complex) %>%
	filter(Coverage == 30) %>%
	write_tsv("results/precision_recall.tsv")

############### plot certainty and precision vs difficulty #####################
cov_labs <- c("10x", "20x", "25x", "30x")
covs <- c(10, 20, 25, 30)
names(cov_labs) <- covs

difficulties <- gts %>%
	filter(Algorithm == "GGTyper") %>%
	group_by(Coverage, VariantType) %>%
	summarize(Certainty = mean(Certainty), Difficulty = mean(Difficulty), Precision = mean(CorrectPrediction), .groups = "keep") %>%
	ungroup() %>%
	pivot_longer(cols = c("Precision", "Certainty"), names_to = "Type", values_to = "Value")
p9 <- difficulties %>%
	ggplot(aes(x = Difficulty, y = Value, col = str_wrap(VariantType, width = 12))) +
	facet_grid(cols = vars(Coverage), rows = vars(Type), switch = "y", labeller = labeller(Coverage = cov_labs)) +
	geom_point(size = 12, pch = "x") +
	geom_smooth(method = "lm", se = FALSE, col = "black") +
	custom_theme_large +
	theme(
	  axis.text.x = element_text(angle = 30, vjust = 0.9, hjust = 0.9), 
	  axis.title.y = element_blank(), 
	  legend.title = element_blank(), 
	  legend.spacing.y = unit(1, "lines") 
	) +
	scale_x_continuous(breaks = c(5000, 10000, 15000)) +
	guides(col = guide_legend(byrow = TRUE))

CairoPDF("plots/difficulties.pdf", width = 20, height = 10, version = "1.5")
plot(p9)
dev.off()

##################### confusion plots per algorithm ############################

n_per_genotype <- gts %>%
  mutate(RealGenotype = translateGenotypes2(RealGenotype)) %>%
  filter(Coverage == 30) %>%
  group_by(VariantType, Algorithm, RealGenotype) %>%
  summarize(N_gt = n(), .groups = "keep")

gts %>%
	filter(Coverage == 30) %>%
	group_by(Algorithm, VariantType) %>%
	summarize(N = n(), .groups = "keep") %>%
	write_tsv("results/n_per_alg.tsv")

p <- gts %>%
  mutate(RealGenotype = translateGenotypes2(RealGenotype)) %>%
  mutate(PredictedGenotype = translateGenotypes2(PredictedGenotype)) %>%
  filter(Coverage == 30 & !is.na(Genotype)) %>%
  group_by(VariantType, Algorithm, RealGenotype, PredictedGenotype) %>%
  summarize(
    N = n(),
    .groups = "keep"
    ) %>%
  inner_join(n_per_genotype, by = c("VariantType", "Algorithm", "RealGenotype")) %>%
  mutate(Fraction = N / N_gt) %>%
  ggplot(aes(x = RealGenotype, y = PredictedGenotype)) +
  facet_grid(cols = vars(Algorithm), rows = vars(VariantType)) +
  geom_tile(aes(fill = Fraction)) +
  geom_label(aes(label = N), alpha = 0.5) +
  scale_x_discrete(limits = c("G_0", "G_1", "G_2"), labels = c(expression("G"["0"]), expression("G"["1"]), expression("G"["2"]))) +
  scale_y_discrete(limits = c("G_0", "G_1", "G_2"), labels = c(expression("G"["0"]), expression("G"["1"]), expression("G"["2"]))) +
  custom_theme +
  theme(
    legend.position = "None",
    strip.text.y = element_text(size = 16)
  ) +
  ylab("Called") +
  xlab("Real")

CairoPDF("plots/n_per_gt.pdf", width = 16, height = 24, version = "1.5")
plot(p)
dev.off()
