# Tim Mirus
# Join information from Polaris Diversity and Kids cohort and store relevant attributes.

library(tidyverse)

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
  summarize(N = n())

polaris_data %>%
	write_tsv("data/polaris_info.tsv")
