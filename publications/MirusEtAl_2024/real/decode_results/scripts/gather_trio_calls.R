# combine genotype calls on Icelandic genomes with Family info,
# calculate transmission and inheritance information and
# store resulting data frame for further analysis

library(tidyverse)

trios <- read.table("decode_results/decode_trios.tsv", sep = "\t", header = FALSE)
colnames(trios) <- c("child", "father", "mother")

genotype_data <- as_tibble(read.table("decode_results/decode_genotype_results.tsv", sep = "\t", header = TRUE))

trio_pass <- function(child, mother, father)
{
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


mendelian_error <- function(child, mother, father)
{
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

trio_calls <- data.frame()
for (i in 1:nrow(trios))
{
	for (variant in unique(genotype_data$Variant))
	{
		temp <- genotype_data %>%
			filter(Variant == variant) %>%
			filter(SampleID %in% trios[i,])
		cert <- min(temp$Certainty)
		q <- min(temp$Mean_Quality)
		mapQ <- min(temp$AvgMapQ)
		r <- mean(temp$TotalReads)
		
		child <- temp %>%
			filter(SampleID == trios[i,1]) %>%
			pull(Genotype)
		father <- temp %>%
			filter(SampleID == trios[i, 2]) %>%
			pull(Genotype)
		mother <- temp %>%
			filter(SampleID == trios[i, 3]) %>%
			pull(Genotype)
		
		trio_calls <- rbind(
				    trio_calls,
				    data.frame(
  						Variant = variant, 
  						TrioID = i, 
  						MinCertainty = cert, 
  						MinQ = q, 
  						MapQ = mapQ, 
  						AvgReads = r, 
  						Child = child, 
  						Mother = mother, 
  						Father = father, 
  						Transmission = trio_pass(child, mother, father), 
  						MendelianError = mendelian_error(child, mother, father)
  						)
				    )
	}
}

as_tibble(trio_calls) %>% write_tsv("decode_results/files/transmission_info.tsv")
