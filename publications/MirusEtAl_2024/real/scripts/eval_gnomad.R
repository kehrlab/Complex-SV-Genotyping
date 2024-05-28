# evaluate gnomadSV genotyping results

# Tim Mirus
# Join information from Polaris Diversity and Kids cohort and store relevant attributes.

library(tidyverse)
library(Cairo)

source("../scripts/R-scripts/theme.R")

ceMin <- 0.9
qMin <- 40
minFreq <- 0.005


kids_data <- read.table(file = "data/kids_cohort.csv", sep = ",", header = TRUE)
diversity_data <- read.table(file = "data/diversity_cohort.csv", sep = ",", header = TRUE)

kids_data <- as_tibble(kids_data) %>%
  rename(SampleID = Sample.ID, FamilyID = Family.ID) %>%
  select(SampleID, FamilyID, Sex, Superpopulation) %>%
  mutate(Role = "Child")
polaris_data <- as_tibble(diversity_data) %>%
  rename(SampleID = Sample.ID, FamilyID = Family.ID) %>%
  select(SampleID, FamilyID, Sex, Superpopulation) %>%
  mutate(Role = ifelse(Sex == "F", "Mother", "Father")) %>%
  bind_rows(kids_data)

roles <- polaris_data$Role
roles[is.na(polaris_data$FamilyID)] <- NA
polaris_data$Role <- roles

polaris_data %>%
  filter(is.na(Role) | Role != "Child") %>%
  group_by(Superpopulation) %>% 
  summarize(N = n(), .groups = "keep")

#polaris_data %>%
#  write_tsv("~/polaris_info.tsv")


#########################
## function definitions #
#########################

adjust_genotypes <- function(gt_vector)
{
  gt_vector <- sapply(gt_vector, function(g) {
    if (g == "REF/REF")
    {
      return (g)
    } else if (g == "ALT/REF" || g == "REF/ALT" || g == "VAR/REF")
    {
      return ("REF/VAR")
    } else if (g == "ALT/ALT")
    {
      return ("VAR/VAR")
    } else {
      return (g)
    }
  })
}

mendelian_error <- function(role, genotype)
{
  if (length(role) != 3 || length(genotype) != 3)
  {
    return(NA)
  }
  
  child <- genotype[which(role == "Child")]
  mother <- genotype[which(role == "Mother")]
  father <- genotype[which(role == "Father")]
  child_alleles <- strsplit(child, split = "/")[[1]]
  mother_alleles <- strsplit(mother, split = "/")[[1]]
  father_alleles <- strsplit(father, split = "/")[[1]]
  
  parent_alleles <- c(mother_alleles, father_alleles)
  
  a1_mother <- sum(mother_alleles == child_alleles[1])
  a2_mother <- sum(mother_alleles == child_alleles[2])
  a1_father <- sum(father_alleles == child_alleles[1])
  a2_father <- sum(father_alleles == child_alleles[2])
  
  if (a1_mother == 0 && a1_father == 0)
  {
    return (TRUE)
  }
  
  if (a2_mother == 0 && a2_father == 0)
  {
    return (TRUE)
  }
  
  if (a1_mother == 0 && a2_mother == 0)
  {
    return (TRUE)
  }
  
  if (a1_father == 0 && a2_father == 0)
  {
    return (TRUE)
  }
  
  return (FALSE)
}

transmission <- function(role, genotype)
{
  if (length(role) != 3 || length(genotype) != 3)
  {
    return(NA)
  }
  
  child <- genotype[which(role == "Child")]
  mother <- genotype[which(role == "Mother")]
  father <- genotype[which(role == "Father")]
  
  child_alleles <- strsplit(child, split = "/")[[1]]
  mother_alleles <- strsplit(mother, split = "/")[[1]]
  father_alleles <- strsplit(father, split = "/")[[1]]
  
  parent_alleles <- c(mother_alleles, father_alleles)
  pass <- -1
  for (a in unique(parent_alleles))
  {
    if (sum(parent_alleles == a) == 1)
    {
      pass <- 0
      if (sum(child_alleles == a) == 1)
      {
        pass <- 1
      }
      break
    }
  }
  return(pass)
}

calculate_allele_freqs <- function(data, pop = "all")
{
  if (pop != "all") {
    data <- data %>% filter(Population == pop)
  }
  n <- length(unique(data %>% pull(SampleID)))
  f <- data %>%
    group_by(Variant, Algorithm, VariantType) %>%
    summarize(
      AlleleFreq = (sum(Genotype == "REF/VAR") + 2*sum(Genotype == "VAR/VAR")) / (2*n),
      N_Het = sum(Genotype == "REF/VAR"),
      .groups = "keep"
    ) %>%
    arrange(desc(AlleleFreq))
  return (f)
}

calculate_pop_hwe <- function(data, population = "all")
{
  hw_genotypes <- genotype_results %>%
    filter(is.na(Role) | Role != "Child") %>%
    filter(Chromosome != "chrX" & Chromosome != "chrY")
  if (population != "all")
  {
    hw_genotypes <- hw_genotypes %>%
      filter(Population == population)
  }
  
  gts <- hw_genotypes$Genotype
  gts[gts == "REF/REF"] <- "G0"
  gts[gts == "REF/VAR"] <- "G1"
  gts[gts == "VAR/VAR"] <- "G2"
  hw_genotypes$Genotype <- gts
  
  observed_counts <- hw_genotypes %>%
    group_by(Algorithm, Variant, Genotype, VariantType) %>%
    summarize(Observed = n(), .groups = "keep") %>%
    pivot_wider(
      names_from = Genotype, 
      names_prefix = "Observed_", 
      values_from = Observed, 
      values_fill = 0
    ) 
  
  expected <- t(apply(
    observed_counts[, c("Observed_G0", "Observed_G1", "Observed_G2")], 
    1, 
    function(x) {
      N <- sum(x)
      aObs <- (2*x[1] + x[2]) / (2*N)
      bObs <- (2*x[3] + x[2]) / (2*N)
      
      g0Exp <- aObs^2 * N
      g2Exp <- bObs^2 * N
      g1Exp <- 2*aObs*bObs*N
      
      chi <- (x[1] - g0Exp)^2/g0Exp + (x[2] - g1Exp)^2/g1Exp + (x[3] - g2Exp)^2/g2Exp
      values <- c(g0Exp, g1Exp, g2Exp, chi)
      names(values) <- c("Expected_G0", "Expected_G1", "Expected_G2", "Chi")
      return(values)
    }
  ))
  
  hwe_counts <- observed_counts %>% 
    bind_cols(expected) %>%
    pivot_longer(
      cols = c(
        Observed_G0, Observed_G1, Observed_G2, 
        Expected_G0, Expected_G1, Expected_G2
      ), 
      names_sep = "_", 
      names_to = c("Type", "Genotype"), 
      values_to = "Count"
    )
  hwe_counts <- hwe_counts %>% filter(!is.na(Chi))
  
  return(hwe_counts)
}

determine_multiallelic <- function(VariantNames)
{
  varNames <- gsub(x = unique(VariantNames), pattern = "_[0-9]+$", replacement = "")
  dupNames <- unique(varNames[which(duplicated(varNames))])
  return(gsub(x = VariantNames, pattern = "_[0-9]+$", replacement = "") %in% dupNames)
}

load_ggtyper_results <- function(d, population_data)
{
  results <- tibble()
  for (f in list.files(path = d, pattern = "*_genotype_results.tsv", full.names = TRUE))
  {
    results <- results %>%
      bind_rows(read_tsv(f))
  }
  
  results <- results %>%
    filter(TotalReads > 0) %>%
    rename(SampleID = Sample, Quality = Mean_Quality) %>%
    select(SampleID, Variant, Genotype, Quality, Certainty, AvgMapQ) %>%
    mutate(Algorithm = "GGTyper") %>%
    mutate(SampleID = gsub(pattern = "-N1-DNA1-WGS1", replacement = "", x = SampleID)) %>%
    mutate(Chromosome = gsub(Variant, pattern = "_.*", replacement = "")) %>%
    mutate(VariantType = gsub(x = Variant, pattern = "chr[0-9,X,Y]+_", replacement = ""))%>%
    mutate(VariantType = gsub(pattern = "[0-9,_]+_", replacement = "", x = VariantType)) %>%
    mutate(VariantType = gsub(x = VariantType, pattern = "_[0-9]+", replacement = "")) %>%
    mutate(Multi = determine_multiallelic(Variant)) %>%
    filter(!is.na(AvgMapQ))
    

  return(results)
}

load_gnomad_info <- function(d)
{
  gnomad <- tibble()
  for (f in list.files(path = d, pattern = "gnomad_structural_variants_[0-9]+-", full.names = TRUE))
  {
      temp <- read_csv(f, progress = FALSE)
      gnomad <- gnomad %>%
        bind_rows(temp)
  }
  return(gnomad)
}

# load GGTyper results
ggtyper_res <- load_ggtyper_results("polarisResults/ggtyper_gnomad_bih/", polaris_data)
genotype_results <- ggtyper_res
rm(ggtyper_res)

gnomad_info <- load_gnomad_info("data/gnomad_structural_variants_download/") %>%
  select("Allele Frequency", Size, "Variant ID", "Homozygote Count") %>%
  rename(AlleleFrequency = "Allele Frequency", VariantID = "Variant ID", HomCount = "Homozygote Count") %>%
  filter(!duplicated(VariantID))

name_assignments <- read_tsv("data/gnomad_variants/gnomad_cpx_name_assignments.tsv", col_names = c("VariantID", "Variant")) %>%
  mutate(VariantID = gsub(VariantID, pattern = "gnomAD-SV_v3_", replacement = "")) %>%
  mutate(VariantID = toupper(VariantID))

gnomad_info <- gnomad_info %>%
  inner_join(name_assignments, by = "VariantID")

# additional files containing information about repeat / problematic region status for each SV (based on UCSC table browser bed files)
# created with scripts/find_variants_in_region and data/hg38*.bed
#repeat_info <- read_tsv("data/gnomad_variants/gnomad_repeats_info.tsv", col_names = c("Variant", "Repeat"))
#seg_dup_info <- read_tsv("data/gnomad_variants/gnomad_seg_dups_info.tsv", col_names = c("Variant", "SegDup"))
#problem_info <- read_tsv("data/gnomad_variants/gnomad_problem_info.tsv", col_names = c("Variant", "ProblemRegion"))

genotype_results <- genotype_results %>%
#  inner_join(repeat_info, by = "Variant") %>%
#  inner_join(seg_dup_info, by = "Variant") %>%
#  inner_join(problem_info, by = "Variant") %>%
  inner_join(gnomad_info, by = "Variant")

genotype_results %>% write_tsv("polarisResults/gnomadCPX_all_genotypes.tsv")

genotype_results <- genotype_results %>% filter(AlleleFrequency >= minFreq) 



# load Possible results from other algorithms here...


# Join with Polaris population info
genotype_results <- genotype_results %>%
  inner_join(polaris_data, by = "SampleID") %>%
  rename(Population = Superpopulation) %>%
  select(-Sex) %>%
  mutate(Genotype = adjust_genotypes(Genotype))

# genotype_results <- genotype_results %>%
#   filter(!Multi & !ProblemRegion)

# determine for each variant the number of samples it can be found in
variantCounts <- genotype_results %>%
  group_by(Algorithm, Variant) %>%
  summarize(Poly = sum(Genotype != "REF/REF"), .groups = "keep")


no_alleles <- sum(genotype_results$Variant %in% (variantCounts %>% filter(Poly == 0) %>% pull(Variant))) / 199


# filter out variants not contained in the polaris data
genotype_results <- genotype_results %>%
  filter(Variant %in% (variantCounts %>% filter(Poly > 0) %>% pull(Variant)))


# determine for each variant type the number of variants
genotype_results <- genotype_results %>% 
  mutate(VariantType = gsub(x = Variant, pattern = "chr[0-9,X,Y]+_", replacement = ""))%>%
  mutate(VariantType = gsub(pattern = "[0-9,_]+_", replacement = "", x = VariantType)) %>%
  mutate(VariantType = gsub(x = VariantType, pattern = "_[0-9]+", replacement = ""))

# determine for each variant type the number of remaining variants
typeCounts <- genotype_results %>%
  group_by(Algorithm, VariantType, Genotype) %>%
  summarize(Unfiltered = n(), Filtered = sum(Certainty >= ceMin & AvgMapQ >= qMin), .groups = "keep") %>%
  pivot_longer(cols = c("Unfiltered", "Filtered"), names_to = "Type", values_to = "N")

########################################################################
# Calculate allele frequencies and compare them to gnomAD frequencies  #
########################################################################

allele_freqs <- calculate_allele_freqs(genotype_results %>% filter(Certainty >= ceMin & AvgMapQ >= qMin), "all")
allele_freqs_gnomad <- allele_freqs %>% inner_join(gnomad_info %>% select(Variant, AlleleFrequency), by = "Variant")
allele_correlation <- cor(allele_freqs_gnomad$AlleleFreq, allele_freqs_gnomad$AlleleFrequency)^2


allele_freqs_unfiltered <- calculate_allele_freqs(genotype_results, "all")
allele_freqs_filtered <- calculate_allele_freqs(genotype_results %>% filter(Certainty >= ceMin & AvgMapQ >= qMin), "all")
allele_freqs <-bind_rows(
     allele_freqs_unfiltered %>% mutate(Filtered = "Unfiltered"),
     allele_freqs_filtered %>% mutate(Filtered = "Filtered")
   )

#######################################################################################################
# Calculate Observed and Expected allele frequencies and Chi-Squared values per Variant and Algorithm #
#######################################################################################################

hwe_counts_all <- calculate_pop_hwe(genotype_results, "all")
hwe_counts_eur <- calculate_pop_hwe(genotype_results, "EUR")
hwe_counts_afr <- calculate_pop_hwe(genotype_results, "AFR")
hwe_counts_eas <- calculate_pop_hwe(genotyoe_results, "EAS")

hwe_counts <-
  bind_rows(
    hwe_counts_all %>% mutate(Population = "All"),
    hwe_counts_afr %>% mutate(Population = "AFR"),
    hwe_counts_eas %>% mutate(Population = "EAS"),
    hwe_counts_eur %>% mutate(Population = "EUR")
  )

hwe_counts_all <- calculate_pop_hwe(genotype_results %>% filter(Certainty > minCe), "all")
hwe_counts_eur <- calculate_pop_hwe(genotype_results %>% filter(Certainty > minCe), "EUR")
hwe_counts_afr <- calculate_pop_hwe(genotype_results %>% filter(Certainty > minCe), "AFR")
hwe_counts_eas <- calculate_pop_hwe(genotyoe_results %>% filter(Certainty > minCe), "EAS")

hwe_counts_filtered <-
  bind_rows(
    hwe_counts_all %>% mutate(Population = "All"),
    hwe_counts_afr %>% mutate(Population = "AFR"),
    hwe_counts_eas %>% mutate(Population = "EAS"),
    hwe_counts_eur %>% mutate(Population = "EUR")
  )

hwe_counts_type_eur <- hwe_counts %>%
  filter(Population == "EUR") %>%
  mutate(VariantType = gsub(x = Variant, pattern = "chr[0-9,X,Y]+_", replacement = ""))%>%
  mutate(VariantType = gsub(pattern = "[0-9,_]+_", replacement = "", x = VariantType)) %>%
  mutate(VariantType = gsub(x = VariantType, pattern = "_[0-9]+", replacement = "")) %>%
  group_by(Algorithm, VariantType) %>%
  summarize(HWE_count = sum(Chi <= 10.83) / 6, N = n() / 6, Fraction = HWE_count / N, .groups = "keep")


hwe_values <- hwe_counts %>%
  group_by(Algorithm, Variant, Population) %>%
  summarize(Chi = unique(Chi), .groups = "keep")


#############################
# Calculate Trio statistics #
#############################


# create transmission and inheritance information
trio_info <- genotype_results %>%
  filter(Chromosome != "chrY" & Chromosome != "chrX") %>%
  filter(!is.na(FamilyID)) %>%
  filter(Genotype != "" & !is.na(Genotype)) %>%
  group_by(Variant, Algorithm, FamilyID, Population, VariantType) %>%
  summarize(
    TR = transmission(Role, Genotype), 
    MendelianError = mendelian_error(Role, Genotype), 
    MinQ = min(Quality), 
    MinCertainty = min(Certainty),
    MinMapQ = min(AvgMapQ),
    .groups = "keep"
  ) %>%
  filter(!is.na(TR))
#trio_info %>% write_tsv("gnomad_trio_info.tsv")

# TR curve
transmission_info <- trio_info %>% 
  filter(TR >= 0) %>%
  filter(MinMapQ >= qMin)
N_per_algorithm <- transmission_info %>%
  group_by(Algorithm) %>%
  summarize(nMax = n(), .groups = "keep")


transmission_rates <- tibble()
for (q in seq(0, max(trio_info$MinQ), by = 2))
{
  tr <- transmission_info %>%
    filter(MinQ >= q) %>%
    group_by(Algorithm) %>%
    summarize(N = n(), Threshold = q, TR = mean(TR), Filter = "Quality", .groups = "keep")
  tr <- tr %>%
    inner_join(N_per_algorithm, by = "Algorithm") %>%
    mutate(Percentage = N / nMax) %>%
    select(-nMax)
  transmission_rates <- transmission_rates %>%
    bind_rows(tr)
}
for (ce in sort(unique(transmission_info %>% filter(!is.na(MinCertainty)) %>% pull(MinCertainty))))
{
  tr <- transmission_info %>%
    filter(!is.na(MinCertainty)) %>%
    filter(MinCertainty >= ce) %>%
    group_by(Algorithm) %>%
    summarize(N = n(), Threshold = ce, TR = mean(TR), Filter = "Certainty", .groups = "keep")
  tr <- tr %>%
    inner_join(N_per_algorithm, by = "Algorithm") %>%
    mutate(Percentage = N / nMax) %>%
    select(-nMax)
  transmission_rates <- transmission_rates %>%
    bind_rows(tr)
}

type_transmission_rates <- tibble()
N_per_algorithm_vtype <- transmission_info %>%
  group_by(Algorithm, VariantType) %>%
  summarize(nMax = n(), .groups = "keep")
for (q in seq(0, max(trio_info$MinQ), by = 2))
{
  tr <- transmission_info %>%
    filter(MinQ >= q) %>%
    group_by(Algorithm, VariantType) %>%
    summarize(N = n(), Threshold = q, TR = mean(TR), Filter = "Quality", .groups = "keep")
  tr <- tr %>%
    inner_join(N_per_algorithm_vtype, by = c("Algorithm", "VariantType")) %>%
    mutate(Percentage = N / nMax) %>%
    select(-nMax)
  type_transmission_rates <- type_transmission_rates %>%
    bind_rows(tr)
}
for (ce in sort(unique(transmission_info %>% filter(!is.na(MinCertainty)) %>% pull(MinCertainty))))
{
  tr <- transmission_info %>%
    filter(!is.na(MinCertainty)) %>%
    filter(MinCertainty >= ce) %>%
    group_by(Algorithm, VariantType) %>%
    summarize(N = n(), Threshold = ce, TR = mean(TR), Filter = "Certainty", .groups = "keep")
  tr <- tr %>%
    inner_join(N_per_algorithm_vtype, by = c("Algorithm", "VariantType")) %>%
    mutate(Percentage = N / nMax) %>%
    select(-nMax)
  type_transmission_rates <- type_transmission_rates %>%
    bind_rows(tr)
}
# transmission_rates %>% write_tsv("polarisResults/transmission_rates.tsv")

# trio_info <- trio_info %>%
#   filter(MinCertainty > 0.99)

# Calculate MIER per Variant and Algorithm
error_rates_var <- trio_info %>%
  group_by(Algorithm, Variant) %>%
  summarize(MIER = mean(MendelianError), .groups = "keep")

error_rates_var_filtered <- trio_info %>%
  filter(MinCertainty > ceMin & MinMapQ >= qMin) %>%
  group_by(Algorithm, Variant) %>%
  summarize(MIER = mean(MendelianError), .groups = "keep")

error_rates_var <- bind_rows(
    error_rates_var %>% mutate(Filter = "Unfiltered"),
    error_rates_var_filtered %>% mutate(Filter = "Filtered")
  )

error_rates_var_pop <- trio_info %>%
  group_by(Algorithm, Variant, Population) %>%
  summarize(MIER = mean(MendelianError), .groups = "keep")

error_rates_type <- trio_info %>%
  mutate(VariantType = gsub(x = Variant, pattern = "chr[0-9,X,Y]+_", replacement = ""))%>%
  mutate(VariantType = gsub(pattern = "[0-9,_]+_", replacement = "", x = VariantType)) %>%
  mutate(VariantType = gsub(x = VariantType, pattern = "_[0-9]+", replacement = "")) %>%
  group_by(Algorithm, VariantType) %>%
  summarize(MIER = mean(MendelianError), N = n(), .groups = "keep")

error_rates_type_filtered <- trio_info %>%
  filter(MinCertainty > ceMin & MinMapQ >= qMin) %>%
  mutate(VariantType = gsub(x = Variant, pattern = "chr[0-9,X,Y]+_", replacement = ""))%>%
  mutate(VariantType = gsub(pattern = "[0-9,_]+_", replacement = "", x = VariantType)) %>%
  mutate(VariantType = gsub(x = VariantType, pattern = "_[0-9]+", replacement = "")) %>%
  group_by(Algorithm, VariantType) %>%
  summarize(MIER = mean(MendelianError), N = n(), .groups = "keep")

error_rates_type <- bind_rows(
  error_rates_type %>% mutate(Filter = "Unfiltered"),
  error_rates_type_filtered %>% mutate(Filter = "Filtered")
)

error_rates_type_pop <- trio_info %>%
  mutate(VariantType = gsub(x = Variant, pattern = "chr[0-9,X,Y]+_", replacement = ""))%>%
  mutate(VariantType = gsub(pattern = "[0-9,_]+_", replacement = "", x = VariantType)) %>%
  mutate(VariantType = gsub(x = VariantType, pattern = "_[0-9]+", replacement = "")) %>%
  group_by(Algorithm, VariantType, Population) %>%
  summarize(MIER = mean(MendelianError), .groups = "keep")
  

error_rates_overall <- trio_info %>%
  group_by(Algorithm) %>%
  summarize(MIER = mean(MendelianError), .groups = "keep")

error_rates_overall_filtered <- trio_info %>%
  filter(MinCertainty >= ceMin & MinMapQ >= qMin) %>%
  group_by(Algorithm) %>%
  summarize(MIER = mean(MendelianError), .groups = "keep")

################## create plots #######################################

p_varFreq <- variantCounts %>%
  filter(Poly > 0) %>%
  ggplot(aes(x = Poly)) + 
  geom_histogram(binwidth = 1, col = "black", fill = "light pink") +
  custom_theme +
  xlab("# samples carrying at least one variant allele") +
  ylab("Number of variants")

p_typeFreq <- typeCounts %>%
  ggplot(aes(x = VariantType, y = N, fill = factor(Genotype, levels = c("VAR/VAR", "REF/VAR", "REF/REF")))) +
  geom_bar(stat = "identity", position = "stack") + # position = fill
  facet_grid(cols = vars(Type)) +
  custom_theme +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 0.9, vjust = 0.9),
    legend.title = element_blank()
  )

p_typeFreq2 <- typeCounts %>%
  ggplot(aes(x = VariantType, y = N, fill = factor(Genotype, levels = c("VAR/VAR", "REF/VAR", "REF/REF")))) +
  geom_bar(stat = "identity", position = "fill") + # position = fill
  facet_grid(cols = vars(Type)) +
  custom_theme +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 0.9, vjust = 0.9),
    legend.title = element_blank()
  )

p_chi <- hwe_counts %>% 
  ggplot(aes(x = Chi)) +
  geom_histogram(binwidth = 1, colour = "black", fill = "light pink") +
  geom_vline(xintercept = 10.83, col = "red", lty = "dashed") +
  facet_grid(cols = vars(factor(Population, levels = c("All", "AFR", "EAS", "EUR")))) +
  xlab(expression(chi^2)) +
  ylab("Count") +
  custom_theme +
  ggtitle("Distribution of Chi^2 values (Unfiltered)")

p_chi_filtered <- hwe_counts_filtered %>% 
  ggplot(aes(x = Chi)) +
  geom_histogram(binwidth = 1, colour = "black", fill = "light pink") +
  geom_vline(xintercept = 10.83, col = "red", lty = "dashed") +
  facet_grid(cols = vars(factor(Population, levels = c("All", "AFR", "EAS", "EUR")))) +
  xlab(expression(chi^2)) +
  ylab("Count") +
  custom_theme +
  ggtitle("Distribution of Chi^2 values (Filtered)")

p_error_type <- error_rates_type %>%
  ggplot(aes(x = VariantType, y = MIER * 100)) +
  geom_bar(stat = "identity", col = "black", fill = "#12becc") +
  ylab("MIER [%]") +
  xlab("") +
  facet_grid(cols = vars(Filter)) +
  custom_theme +
  theme(axis.text.x = element_text(angle = 30, vjust = 0.9, hjust = 0.9))

p_n_type_filtered <- error_rates_type %>%
  ggplot(aes(x = VariantType, y = N)) +
  geom_bar(stat = "identity", col = "black", fill = "#12becc") +
  ylab("Number of genotype calls") +
  facet_grid(cols = vars(Filter)) +
  custom_theme_large +
  theme(
    axis.text.x = element_text(angle = 30, vjust = 0.9, hjust = 0.9),
    axis.title.x = element_blank()
    )

p_n_type_filtered_vars <- error_rates_type %>%
  filter(Filter == "Unfiltered") %>%
  ggplot(aes(x = VariantType, y = N / 199)) +
  geom_bar(stat = "identity", col = "black", fill = "#12becc") +
  ylab("Number of distinct variants") +
  custom_theme_large +
  theme(
    axis.text.x = element_text(angle = 30, vjust = 0.9, hjust = 0.9),
    axis.title.x = element_blank()
    )

p_certainty_dist <- genotype_results %>%
  ggplot(aes(x = VariantType, y = Certainty, group = interaction(VariantType, Genotype), col = Genotype)) +
  geom_violin(col = "black") +
  geom_jitter(alpha = 0.1, col = "light pink") +
  custom_theme +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 0.9, vjust = 0.9)
  )

p_hwe_error <- error_rates_var %>%
  inner_join(hwe_values, by = "Variant") %>%
  mutate(Bin = floor(Chi / 10)) %>%
  ggplot(aes(x = Chi, y = MIER * 100, group = Bin)) +
  geom_point(col = "light pink", alpha = 0.3) +
  geom_boxplot() +
  facet_grid(cols = vars(Filter)) +
  geom_vline(xintercept = 10.83, lty = "dashed", col = "red") +
  ylab("MIER [%]") +
  custom_theme

p_allele <- allele_freqs %>%
  ggplot(aes(x = AlleleFreq, y = N_Het, col = VariantType)) +
  geom_point() +
  facet_grid(cols = vars(Filtered)) +
  custom_theme

p_af <- allele_freqs_gnomad %>%
  ggplot(aes(x = AlleleFrequency, y = AlleleFreq)) +
  geom_point(size = 3) +
  xlab("gnomAD-SV AF") +
  ylab("Polaris AF") +
  custom_theme +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 0.9, vjust = 0.9),
    panel.spacing = unit(2, "lines")
  ) +
  geom_smooth(method = "lm") + 
  geom_text(x = 0.75, y = 0.55, label = paste0("r=", round(cor(allele_freqs_gnomad$AlleleFrequency, allele_freqs_gnomad$AlleleFreq), 2)), size.unit = "pt", size = 20)

point1 <- transmission_rates %>%
  filter(Filter == "Certainty") %>%
  filter(Threshold >= 0.9) %>%
  arrange(Threshold)
point1 <- point1[1,]

p_tr <- transmission_rates %>%
  group_by(Filter, Algorithm) %>%
  ggplot(aes(x = Percentage * 100, y = TR * 100, col = Filter)) +
  geom_line(linewidth = 2) +
  geom_hline(yintercept = 50, lty = "dashed") +
  ylim(0, 100) +
  ylab("Transmission Rate [%]") +
  xlab("# of Trios [%]") +
  custom_theme +
  geom_point(aes(x = point1$Percentage * 100, y = point1$TR * 100), col = "black", size = 4, shape = 8)


#p_ttr <- type_transmission_rates %>%
#  group_by(Filter, Algorithm) %>%
#  ggplot(aes(x = Percentage, y = TR, col = Filter)) +
#  geom_line() +
#  geom_hline(yintercept = 0.5, lty = "dashed") +
#  ylim(0.25, 0.75) +
#  custom_theme +
#  facet_grid(cols = vars(VariantType))


####### write plots

CairoPDF("plots/gnomadSV_plots.pdf", width = 16, height = 9, version = "1.5")
plot(p_varFreq)
plot(p_typeFreq)
plot(p_typeFreq2)
plot(p_error_type)
plot(p_n_type_filtered)
plot(p_hwe_error)
plot(p_allele)
plot(p_chi)
plot(p_chi_filtered)
plot(p_certainty_dist)
plot(p_tr)
plot(p_af)
plot(p_n_type_filtered_vars)
dev.off()

#pdf("plots/gnomadSV_transmission_rates.pdf", width = 100, height = 9)
#plot(p_ttr)
#dev.off()

#pdf("plots/gnomadSV_VariantType_chis.pdf", width = 16, height = 9)
#for (vt in hwe_counts %>% pull(VariantType) %>% unique()) {
#  p_chi <- hwe_counts %>% 
#    filter(VariantType == vt) %>%
#    ggplot(aes(x = Chi)) +
#    geom_histogram(binwidth = 1, colour = "black", fill = "light pink") +
#    geom_vline(xintercept = 10.83, col = "red", lty = "dashed") +
#    facet_grid(cols = vars(factor(Population, levels = c("All", "AFR", "EAS", "EUR")))) +
#    xlab(expression(chi^2)) +
#    ylab("Count") +
#    custom_theme +
#    ggtitle(vt)
#  plot(p_chi)
#}
#dev.off()

cat("Remove due to no alleles: ", no_alleles, "\n")
cat("Without filtering: \n")
cat("Calls: ", nrow(genotype_results))
cat("#: ",length(genotype_results %>% pull(Variant) %>% unique), "\n")
cat("# per variant type: \n")
genotype_results %>%
  group_by(VariantType) %>%
  summarize(N = length(unique(Variant))) %>%
  print
cat("----------------------------------------\n")
cat("With filtering: \n")
cat("Calls: ", nrow(genotype_results %>% filter(Quality >= qMin & Certainty >= ceMin)), "\n")
cat("#: ",length(genotype_results %>% filter(Quality >= qMin & Certainty >= ceMin) %>% pull(Variant) %>% unique), "\n")
cat("# per variant type: \n")
genotype_results %>%
  filter(Quality >= qMin & Certainty >= ceMin) %>%
  group_by(VariantType) %>%
  summarize(N = length(unique(Variant))) %>%
  print
cat("----------------------------------------\n")
cat("r^2: ", allele_correlation, "\n")
cat("HWE:\n")
hwe_counts %>%
  group_by(Algorithm, Population) %>%
  summarize(HWE = sum(Chi[!is.na(Chi)] <= 10.83) / 6, Total = n() / 6, Percentage = HWE / Total * 100, .groups = "keep") %>%
  print
cat("Overall filtered MIER:\n")
print(error_rates_overall_filtered)
cat("Recall:\n")
genotype_results %>%
  summarize(Recall = mean(Certainty >= ceMin & AvgMapQ >= qMin)) %>%
  print
cat("TR for certainty >= 0.9:\n")
print(point1)


#################
vcf_af <- read_csv("data/gnomad_variants/gnomad_cpx_AF.csv") # these AFs are based on the VCF file!
p_af_vcf <- vcf_af %>% 
  mutate(eur_AN = fin_AN + nfe_AN, eur_AF = (fin_AF * fin_AN + nfe_AN * nfe_AF) / eur_AN) %>% 
  pivot_longer(cols = c("eas_AF", "afr_AF", "eur_AF"), values_to = "AF", names_to = "Population") %>%
  mutate(Population = sapply(Population, function(x){if (x == "eur_AF") return("EUR") else if (x == "eas_AF") return("EAS") else return("AFR")})) %>%
  select(-c("fin_AF", "fin_AN", "nfe_AN", "nfe_AF", "eas_AN", "eur_AN", "afr_AN")) %>%
  filter(AF > 0) %>%
  ggplot(aes(x = AF)) +
  facet_grid(cols = vars(Population)) +
  geom_histogram(col = "black", fill = "#12becc") +
  custom_theme +
  ylab("Count") +
  xlab("gnomAD-SV allele frequency (VCF)")


p_af_browser_1 <- gnomad_info %>%
  filter(AlleleFrequency > 0) %>%
  ggplot(aes(x = AlleleFrequency)) +
  geom_histogram(binwidth = 0.01, col = "black", fill = "#12becc") +
  custom_theme_large +
  ylab("Count") +
  xlab("gnomAD-SV allele frequency (Browser)")
p_af_browser_2 <- gnomad_info %>%
  filter(AlleleFrequency > 0.005) %>%
  ggplot(aes(x = AlleleFrequency)) +
  geom_histogram(binwidth = 0.01, col = "black", fill = "#12becc") +
  custom_theme_large +
  ylab("Count") +
  xlab("gnomAD-SV allele frequency (Browser)")

CairoPDF("plots/afs.pdf", width = 16, height = 9, version = "1.5")
plot(p_af_vcf)
plot(p_af_browser_1)
plot(p_af_browser_2)
dev.off()
