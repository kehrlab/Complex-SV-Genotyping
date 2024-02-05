# Tim Mirus
# Calculate population and inheritance statistics for the Polaris data.

library(tidyverse)

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


################################################################
# gather genotype calls and family information into one tibble #
################################################################

polarisInfo <- read_tsv("data/polaris_info.tsv") %>%
  mutate(SampleID = paste0(SampleID, "-N1-DNA1-WGS1"))
ggtyper_res <- read_tsv("polarisResults/ggtyper_genotype_results.tsv") %>%
  rename(SampleID = Sample, Quality = Mean_Quality) %>%
  select(SampleID, Variant, Genotype, Quality, Certainty, AvgMapQ) %>%
  mutate(Algorithm = "GGTyper")

paragraph_vcf <- read_tsv("polarisResults/paragraph_results_vcf/paragraph_genotypes.tsv") %>%
  rename(SampleID = Sample) %>%
  select(SampleID, Variant, Genotype, Quality) %>%
  mutate(Algorithm = "Paragraph(Graph)")

paragraph_graph <- read_tsv("polarisResults/paragraph_results_graph/paragraph_graph_genotypes.tsv") %>%
  rename(SampleID = Sample) %>%
  mutate(Algorithm = "Paragraph(VCF)")

bayestyper_res <- read_tsv("polarisResults/bayestyper_genotypes.tsv") %>%
	mutate(SampleID = paste0(Sample, "-N1-DNA1-WGS1")) %>%
	select(SampleID, Variant, Genotype, Quality) %>%
	mutate(Algorithm = "BayesTyper")

genotype_results <- bind_rows(ggtyper_res, paragraph_vcf, paragraph_graph, bayestyper_res)


rm(ggtyper_res)
rm(paragraph_graph)
rm(paragraph_vcf)
rm(bayestyper_res)
gc()

genotype_results <- genotype_results %>%
  inner_join(polarisInfo, by = "SampleID") %>%
  rename(Population = Superpopulation) %>%
  select(-Sex) %>%
  mutate(Genotype = adjust_genotypes(Genotype))

rm(polarisInfo)
gc()
genotype_results %>% write_tsv("polarisResults/all_results.tsv")

#############################
# Calculate Trio statistics #
#############################


# create transmission and inheritance information
trio_info <- genotype_results %>%
  filter(!is.na(FamilyID)) %>%
  filter(Genotype != "" & !is.na(Genotype)) %>%
  group_by(Variant, Algorithm, FamilyID) %>%
  summarize(
    TR = transmission(Role, Genotype), 
    MendelianError = mendelian_error(Role, Genotype), 
    MinQ = min(Quality), 
    MinCertainty = min(Certainty), 
    MinMapQ = min(AvgMapQ),
    .groups = "keep"
    ) %>%
  filter(!is.na(TR))
trio_info %>% write_tsv("polarisResults/trio_info.tsv")

# TR curve
transmission_info <- trio_info %>% 
  filter(TR >= 0)
N_per_algorithm <- transmission_info %>%
  group_by(Algorithm) %>%
  summarize(nMax = n())

transmission_rates <- tibble()
for (q in seq(0, max(trio_info$MinQ), by = 2))
{
  tr <- transmission_info %>%
    filter(MinQ >= q) %>%
    group_by(Algorithm) %>%
    summarize(N = n(), Threshold = q, TR = mean(TR), Filter = "Quality")
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
    summarize(N = n(), Threshold = ce, TR = mean(TR), Filter = "Certainty")
  tr <- tr %>%
    inner_join(N_per_algorithm, by = "Algorithm") %>%
    mutate(Percentage = N / nMax) %>%
    select(-nMax)
  transmission_rates <- transmission_rates %>%
    bind_rows(tr)
}
transmission_rates %>% write_tsv("polarisResults/transmission_rates.tsv")

# Calculate MIER per Variant and Algorithm
error_rates <- trio_info %>%
  group_by(Algorithm, Variant) %>%
  summarize(MIER = mean(MendelianError), .groups = "keep")
error_rates %>% write_tsv("polarisResults/mendelian_errors.tsv")

#######################################################################################################
# Calculate Observed and Expected allele frequencies and Chi-Squared values per Variant and Algorithm #
#######################################################################################################

hw_genotypes <- genotype_results %>%
  filter(is.na(Role) | Role != "Child")
gts <- hw_genotypes$Genotype
gts[gts == "REF/REF"] <- "G0"
gts[gts == "REF/VAR"] <- "G1"
gts[gts == "VAR/VAR"] <- "G2"
hw_genotypes$Genotype <- gts

observed_counts <- hw_genotypes %>%
  group_by(Algorithm, Variant, Genotype) %>%
  summarize(Observed = n()) %>%
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

hwe_counts %>% write_tsv("polarisResults/hwe_counts.tsv")

hwe_counts %>%
  group_by(Algorithm) %>%
  summarize(HWE = sum(Chi <= 10.83) / 6) %>%
  write_tsv("polarisResults/hwe_numbers.tsv")

######### HWE only on EUR population ##########################################

hw_genotypes <- genotype_results %>%
  filter(is.na(Role) | Role != "Child") %>%
  filter(Population == "EUR")
gts <- hw_genotypes$Genotype
gts[gts == "REF/REF"] <- "G0"
gts[gts == "REF/VAR"] <- "G1"
gts[gts == "VAR/VAR"] <- "G2"
hw_genotypes$Genotype <- gts

observed_counts <- hw_genotypes %>%
  group_by(Algorithm, Variant, Genotype) %>%
  summarize(Observed = n()) %>%
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

hwe_counts %>% write_tsv("polarisResults/hwe_counts_eur.tsv")

hwe_counts %>%
  group_by(Algorithm) %>%
  summarize(HWE = sum(Chi <= 10.83) / 6) %>%
  write_tsv("polarisResults/hwe_numbers_eur.tsv")

######### HWE only on AFR population ##########################################

hw_genotypes <- genotype_results %>%
  filter(is.na(Role) | Role != "Child") %>%
  filter(Population == "AFR")
gts <- hw_genotypes$Genotype
gts[gts == "REF/REF"] <- "G0"
gts[gts == "REF/VAR"] <- "G1"
gts[gts == "VAR/VAR"] <- "G2"
hw_genotypes$Genotype <- gts

observed_counts <- hw_genotypes %>%
  group_by(Algorithm, Variant, Genotype) %>%
  summarize(Observed = n()) %>%
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

hwe_counts %>% write_tsv("polarisResults/hwe_counts_afr.tsv")

hwe_counts %>%
  group_by(Algorithm) %>%
  summarize(HWE = sum(Chi <= 10.83) / 6) %>%
  write_tsv("polarisResults/hwe_numbers_afr.tsv")

######### HWE only on EAS population ##########################################

hw_genotypes <- genotype_results %>%
  filter(is.na(Role) | Role != "Child") %>%
  filter(Population == "EAS")
gts <- hw_genotypes$Genotype
gts[gts == "REF/REF"] <- "G0"
gts[gts == "REF/VAR"] <- "G1"
gts[gts == "VAR/VAR"] <- "G2"
hw_genotypes$Genotype <- gts

observed_counts <- hw_genotypes %>%
  group_by(Algorithm, Variant, Genotype) %>%
  summarize(Observed = n()) %>%
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

hwe_counts %>% write_tsv("polarisResults/hwe_counts_eas.tsv")

hwe_counts %>%
  group_by(Algorithm) %>%
  summarize(HWE = sum(Chi <= 10.83) / 6) %>%
  write_tsv("polarisResults/hwe_numbers_eas.tsv")
