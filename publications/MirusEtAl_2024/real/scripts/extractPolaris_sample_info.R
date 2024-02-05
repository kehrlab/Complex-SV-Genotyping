# Tim Mirus
# Join information from Polaris Diversity and Kids cohort and store relevant attributes.

library(tidyverse)

kids <- read_csv("data/kids_cohort.csv", name_repair = "universal")
diversity <- read_csv("data/diversity_cohort.csv", name_repair = "universal")

all <- kids %>%
  bind_rows(diversity)

sampleInfo <- all %>% 
  mutate(ENA = sapply(strsplit(ENA.link, split = "/"), function(x) {x[length(x)]})) %>%
  select(Sample.ID, ENA, Family.ID, Sex)

sampleInfo %>% write_tsv("data/polaris_sample_info.tsv")
